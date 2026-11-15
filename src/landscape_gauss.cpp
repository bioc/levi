#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// -----------------------------------------------------------------------------
// landscape_gauss(): builds the landscape by normalised convolution, that is,
// with the Nadaraya-Watson estimator using a Gaussian kernel.
//
//            sum_i w_i(x) * s_i                          ( -||x - p_i||^2 )
//   S(x) = ----------------------  ,   w_i(x) = exp      ( -------------- )
//               sum_i w_i(x)                             (    2 sigma^2   )
//
// Dividing by sum_i w_i(x) is the key step: it makes S(x) a convex average of
// the signals, hence always contained in [min s_i, max s_i]. Data with no
// variation yields the neutral value for any sigma and any network density,
// which was not the case when background points with signal 0 entered the
// average.
//
// The denominator D(x) = sum_i w_i(x) measures network occupancy and defines
// the silhouette: where D is low there is no network, and the cell gets NA
// instead of 0. This keeps "no network here" from being confused with "network
// present, with minimal expression".
//
// Cost: the Gaussian is separable, so the 2-D convolution becomes two 1-D
// passes, O(res^2 * r) with r ~ 3*sigma, instead of O(res^2 * n). The number
// of points enters only in the initial deposit, which is O(n).
// -----------------------------------------------------------------------------

// Separable blur. The kernel is deliberately not renormalised by the local
// weight: numerator and denominator undergo the same truncation near the
// border, so their ratio remains a legitimate weighted average.
//
// Both passes add the kernel terms in ascending offset t, starting from 0,
// so every cell receives exactly the same sequence of floating-point
// operations as the earlier cell-by-cell loop and the results are
// bit-identical; the loops are only reordered so that the inner one runs over
// contiguous memory.
static void gaussBlur(std::vector<double>& g, std::vector<double>& tmp,
                      int res, const std::vector<double>& kern, int r)
{
    const std::size_t R = (std::size_t)res;
    std::fill(tmp.begin(), tmp.end(), 0.0);
    for (int i = 0; i < res; ++i) {                       // pass along rows
        double* o = &tmp[(std::size_t)i * R];
        const int lo = std::max(-r, -i), hi = std::min(r, res - 1 - i);
        for (int t = lo; t <= hi; ++t) {
            const double k = kern[(std::size_t)(t + r)];
            const double* in = &g[(std::size_t)(i + t) * R];
            for (int j = 0; j < res; ++j) o[j] += k * in[j];
        }
    }
    for (int i = 0; i < res; ++i) {                       // pass along columns
        const double* in = &tmp[(std::size_t)i * R];
        double* o = &g[(std::size_t)i * R];
        for (int j = 0; j < res; ++j) {
            const int lo = std::max(-r, -j), hi = std::min(r, res - 1 - j);
            double acc = 0.0;
            for (int t = lo; t <= hi; ++t)
                acc += kern[(std::size_t)(t + r)] * in[j + t];
            o[j] = acc;
        }
    }
}

// Bilinear deposit: spreads the value over the 4 neighbouring cells, so the
// result does not jump when a point falls between grid cells.
static inline void splat(std::vector<double>& num, std::vector<double>* den,
                         int res, double fi, double fj, double s, double wt)
{
    const int i0 = (int)std::floor(fi), j0 = (int)std::floor(fj);
    const double di = fi - i0, dj = fj - j0;
    for (int a = 0; a <= 1; ++a) {
        for (int b = 0; b <= 1; ++b) {
            const int i = i0 + a, j = j0 + b;
            if (i < 0 || i >= res || j < 0 || j >= res) continue;
            const double w = (a ? di : 1.0 - di) * (b ? dj : 1.0 - dj) * wt;
            if (w <= 0.0) continue;
            const std::size_t p = (std::size_t)i * res + j;
            num[p] += w * s;
            if (den) (*den)[p] += w;
        }
    }
}

//' Landscape by normalised convolution with a Gaussian kernel
//'
//' @param coord n x 2 matrix with the normalised coordinates of the points
//'   (nodes followed by edge midpoints).
//' @param SignalOut n x 1 matrix with the combined signal (test vs control).
//' @param signalExp n x 1 matrix with the test signal.
//' @param signalCtrl n x 1 matrix with the control signal.
//' @param resolutionValue side of the square output grid.
//' @param zoomValue coordinate of the lower-left corner of the grid.
//' @param increase grid step, in coordinate units.
//' @param sigma kernel width, in grid cells.
//' @param occFrac fraction of the occupancy produced by an isolated point
//'   below which the cell is considered background and gets NA.
//' @param weights optional n-vector of non-negative support weights, one per
//'   point. A point of weight w counts as w copies of itself in both the
//'   numerator and the occupancy; weight 0 removes it. Empty means all ones.
//' @return List with m1 (combined signal), m2 (test), m3 (control) and occ
//'   (relative network occupancy), all of them resolutionValue x
//'   resolutionValue matrices.
//' @keywords internal
// [[Rcpp::export(rng = false)]]
List landscape_gauss(NumericMatrix coord,
                     NumericMatrix SignalOut,
                     NumericMatrix signalExp,
                     NumericMatrix signalCtrl,
                     int resolutionValue,
                     double zoomValue,
                     double increase,
                     double sigma,
                     double occFrac,
                     NumericVector weights = NumericVector(0))
{
    const int res = resolutionValue;
    if (res < 1) stop("resolutionValue must be >= 1");
    const int n = coord.nrow();
    if (n < 1) stop("coord must have at least one row");
    if (SignalOut.nrow() != n || signalExp.nrow() != n || signalCtrl.nrow() != n)
        stop("coord and the signal matrices must have the same number of rows");
    if (!(increase > 0.0)) stop("increase must be > 0");
    const bool weighted = weights.size() > 0;
    if (weighted && weights.size() != n)
        stop("weights must have one value per point");

    const std::size_t N = (std::size_t)res * res;

    if (!(sigma > 0.35)) sigma = 0.35;   // below this the kernel degenerates
    int r = (int)std::ceil(3.0 * sigma);
    if (r < 1) r = 1;
    if (r > res) r = res;

    // 1-D kernel, normalised to sum to 1
    std::vector<double> kern((std::size_t)(2 * r + 1));
    double ksum = 0.0;
    for (int t = -r; t <= r; ++t) {
        const double v = std::exp(-(double)t * t / (2.0 * sigma * sigma));
        kern[(std::size_t)(t + r)] = v;
        ksum += v;
    }
    for (std::size_t t = 0; t < kern.size(); ++t) kern[t] /= ksum;

    std::vector<double> numOut(N, 0.0), numExp(N, 0.0), numCtrl(N, 0.0),
                        den(N, 0.0), tmp(N, 0.0);

    // 1. deposit -- O(n)
    for (int m = 0; m < n; ++m) {
        const double fi = (coord(m, 0) - zoomValue) / increase;
        const double fj = (coord(m, 1) - zoomValue) / increase;
        if (!R_finite(fi) || !R_finite(fj)) continue;
        const double wt = weighted ? weights[m] : 1.0;
        if (!(wt > 0.0)) continue;
        splat(numOut,  &den, res, fi, fj, SignalOut(m, 0), wt);
        splat(numExp,  NULL, res, fi, fj, signalExp(m, 0), wt);
        splat(numCtrl, NULL, res, fi, fj, signalCtrl(m, 0), wt);
    }

    // 2. separable blur -- O(res^2 * r)
    gaussBlur(numOut,  tmp, res, kern, r);
    gaussBlur(numExp,  tmp, res, kern, r);
    gaussBlur(numCtrl, tmp, res, kern, r);
    gaussBlur(den,     tmp, res, kern, r);

    // 3. division and mask -- O(res^2)
    // An isolated point produces, in its own cell, occupancy ~ kern[r]^2.
    const double isolated = kern[(std::size_t)r] * kern[(std::size_t)r];
    const double occMin = occFrac * isolated;

    NumericMatrix mOut(res, res), mExp(res, res), mCtrl(res, res), mOcc(res, res);
    for (int i = 0; i < res; ++i) {
        for (int j = 0; j < res; ++j) {
            const std::size_t p = (std::size_t)i * res + j;
            const double d = den[p];
            mOcc(i, j) = (isolated > 0.0) ? d / isolated : 0.0;
            if (d > occMin && d > 0.0) {
                mOut(i, j)  = numOut[p]  / d;
                mExp(i, j)  = numExp[p]  / d;
                mCtrl(i, j) = numCtrl[p] / d;
            } else {
                mOut(i, j) = mExp(i, j) = mCtrl(i, j) = NA_REAL;
            }
        }
    }

    return List::create(_["m1"] = mOut, _["m2"] = mExp,
                        _["m3"] = mCtrl, _["occ"] = mOcc);
}

// Kernel and half-width used by landscape_gauss(), shared by the one-channel
// entry points below so that all of them blur with identical weights.
static std::vector<double> makeKernel(double& sigma, int res, int& r)
{
    if (!(sigma > 0.35)) sigma = 0.35;
    r = (int)std::ceil(3.0 * sigma);
    if (r < 1) r = 1;
    if (r > res) r = res;
    std::vector<double> kern((std::size_t)(2 * r + 1));
    double ksum = 0.0;
    for (int t = -r; t <= r; ++t) {
        const double v = std::exp(-(double)t * t / (2.0 * sigma * sigma));
        kern[(std::size_t)(t + r)] = v;
        ksum += v;
    }
    for (std::size_t t = 0; t < kern.size(); ++t) kern[t] /= ksum;
    return kern;
}

//' Smoothed occupancy of the network, for repeated one-channel landscapes
//'
//' The occupancy depends only on the coordinates and the support weights, not
//' on the signal, so the permutation loops compute it once and reuse it
//' through \code{landscape_m1()}.
//'
//' @inheritParams landscape_gauss
//' @return Resolution x resolution matrix with the smoothed occupancy where it
//'   exceeds the silhouette threshold and NA elsewhere.
//' @keywords internal
// [[Rcpp::export(rng = false)]]
NumericMatrix landscape_occupancy(NumericMatrix coord,
                                  int resolutionValue,
                                  double zoomValue,
                                  double increase,
                                  double sigma,
                                  double occFrac,
                                  NumericVector weights = NumericVector(0))
{
    const int res = resolutionValue;
    if (res < 1) stop("resolutionValue must be >= 1");
    const int n = coord.nrow();
    if (n < 1) stop("coord must have at least one row");
    if (!(increase > 0.0)) stop("increase must be > 0");
    const bool weighted = weights.size() > 0;
    if (weighted && weights.size() != n)
        stop("weights must have one value per point");
    const std::size_t N = (std::size_t)res * res;
    int r;
    const std::vector<double> kern = makeKernel(sigma, res, r);
    std::vector<double> den(N, 0.0), unused(N, 0.0), tmp(N, 0.0);
    for (int m = 0; m < n; ++m) {
        const double fi = (coord(m, 0) - zoomValue) / increase;
        const double fj = (coord(m, 1) - zoomValue) / increase;
        if (!R_finite(fi) || !R_finite(fj)) continue;
        const double wt = weighted ? weights[m] : 1.0;
        if (!(wt > 0.0)) continue;
        splat(unused, &den, res, fi, fj, 0.0, wt);
    }
    gaussBlur(den, tmp, res, kern, r);
    const double isolated = kern[(std::size_t)r] * kern[(std::size_t)r];
    const double occMin = occFrac * isolated;
    NumericMatrix out(res, res);
    for (int i = 0; i < res; ++i) {
        for (int j = 0; j < res; ++j) {
            const double d = den[(std::size_t)i * res + j];
            out(i, j) = (d > occMin && d > 0.0) ? d : NA_REAL;
        }
    }
    return out;
}

//' Combined-signal landscape with a precomputed occupancy
//'
//' Returns exactly the \code{m1} matrix of \code{landscape_gauss()} for the
//' same arguments, at the cost of one blur instead of four.
//'
//' @inheritParams landscape_gauss
//' @param signal numeric vector with the combined signal of every point.
//' @param occupancy matrix from \code{landscape_occupancy()} for the same
//'   coordinates, grid, sigma and weights.
//' @return Resolution x resolution matrix; NA outside the silhouette.
//' @keywords internal
// [[Rcpp::export(rng = false)]]
NumericMatrix landscape_m1(NumericMatrix coord,
                           NumericVector signal,
                           NumericMatrix occupancy,
                           int resolutionValue,
                           double zoomValue,
                           double increase,
                           double sigma,
                           NumericVector weights = NumericVector(0))
{
    const int res = resolutionValue;
    if (res < 1) stop("resolutionValue must be >= 1");
    const int n = coord.nrow();
    if (n < 1) stop("coord must have at least one row");
    if (signal.size() != n)
        stop("signal must have one value per point");
    if (occupancy.nrow() != res || occupancy.ncol() != res)
        stop("occupancy must be a resolution x resolution matrix");
    if (!(increase > 0.0)) stop("increase must be > 0");
    const bool weighted = weights.size() > 0;
    if (weighted && weights.size() != n)
        stop("weights must have one value per point");
    const std::size_t N = (std::size_t)res * res;
    int r;
    const std::vector<double> kern = makeKernel(sigma, res, r);
    std::vector<double> num(N, 0.0), tmp(N, 0.0);
    for (int m = 0; m < n; ++m) {
        const double fi = (coord(m, 0) - zoomValue) / increase;
        const double fj = (coord(m, 1) - zoomValue) / increase;
        if (!R_finite(fi) || !R_finite(fj)) continue;
        const double wt = weighted ? weights[m] : 1.0;
        if (!(wt > 0.0)) continue;
        splat(num, NULL, res, fi, fj, signal[m], wt);
    }
    gaussBlur(num, tmp, res, kern, r);
    NumericMatrix out(res, res);
    for (int i = 0; i < res; ++i) {
        for (int j = 0; j < res; ++j) {
            const double d = occupancy(i, j);
            out(i, j) = R_IsNA(d) ? NA_REAL : num[(std::size_t)i * res + j] / d;
        }
    }
    return out;
}

//' Index of the node nearest to each grid cell
//'
//' Used by the GUI to report which gene a region selected on the landscape
//' belongs to. Brute force in O(res^2 * n), exact and cheap: for a network of
//' 1200 nodes on a 240 x 240 grid this is about 70 million comparisons.
//'
//' @param coordNodes n x 2 matrix with the normalised coordinates of the nodes.
//' @param resolutionValue side of the square grid.
//' @param zoomValue coordinate of the lower-left corner of the grid.
//' @param increase grid step.
//' @return Integer matrix resolutionValue x resolutionValue with the
//'   (1-based) index of the node nearest to each cell.
//' @keywords internal
// [[Rcpp::export(rng = false)]]
IntegerMatrix nearest_node_grid(NumericMatrix coordNodes,
                                int resolutionValue,
                                double zoomValue,
                                double increase)
{
    const int res = resolutionValue;
    const int n = coordNodes.nrow();
    if (res < 1) stop("resolutionValue must be >= 1");
    if (n < 1) stop("coordNodes must have at least one row");

    std::vector<double> gx((std::size_t)res), gy((std::size_t)res);
    for (int i = 0; i < res; ++i) gx[(std::size_t)i] = zoomValue + i * increase;
    for (int j = 0; j < res; ++j) gy[(std::size_t)j] = zoomValue + j * increase;

    IntegerMatrix out(res, res);
    for (int i = 0; i < res; ++i) {
        const double px = gx[(std::size_t)i];
        for (int j = 0; j < res; ++j) {
            const double py = gy[(std::size_t)j];
            double best = R_PosInf;
            int arg = NA_INTEGER;
            for (int m = 0; m < n; ++m) {
                const double dx = px - coordNodes(m, 0);
                const double dy = py - coordNodes(m, 1);
                const double d2 = dx * dx + dy * dy;
                if (d2 < best) { best = d2; arg = m + 1; }
            }
            out(i, j) = arg;
        }
    }
    return out;
}
