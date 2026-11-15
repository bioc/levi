// Largest regional excess mass of a landscape, for the permutation nulls.
// Cells whose value exceeds 0.5 by `threshold` (over) or falls below it by
// the same amount (under) are joined by eight-connectivity; a component with
// at least `min_cells` cells has mass = sum(excess) * increase^2, and the
// maximum over components is returned for each direction. This is the C++
// counterpart of the R reference .maximumRegionMassR().
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

static int rm_find(std::vector<int>& parent, int a) {
    while (parent[a] != a) {
        parent[a] = parent[parent[a]];
        a = parent[a];
    }
    return a;
}

static double side_mass(const NumericMatrix& z, bool over, double threshold,
                        int min_cells, double increase)
{
    const int nr = z.nrow(), nc = z.ncol();
    const std::size_t N = (std::size_t)nr * nc;
    std::vector<char> on(N, 0);
    std::vector<double> excess(N, 0.0);
    bool any = false;
    for (int j = 0; j < nc; ++j) {
        for (int i = 0; i < nr; ++i) {
            const double v = z(i, j);
            if (!R_finite(v)) continue;
            const std::size_t p = (std::size_t)j * nr + i;
            if (over ? (v >= 0.5 + threshold) : (v <= 0.5 - threshold)) {
                on[p] = 1;
                const double e = over ? v - 0.5 - threshold : 0.5 - v - threshold;
                excess[p] = e > 0.0 ? e : 0.0;
                any = true;
            }
        }
    }
    if (!any) return 0.0;

    std::vector<int> parent(N);
    for (std::size_t p = 0; p < N; ++p) parent[p] = (int)p;
    // Forward neighbours (down, right, down-right, down-left) cover all eight.
    static const int di[4] = {1, 0, 1, 1};
    static const int dj[4] = {0, 1, 1, -1};
    for (int j = 0; j < nc; ++j) {
        for (int i = 0; i < nr; ++i) {
            const std::size_t p = (std::size_t)j * nr + i;
            if (!on[p]) continue;
            for (int k = 0; k < 4; ++k) {
                const int ii = i + di[k], jj = j + dj[k];
                if (ii < 0 || ii >= nr || jj < 0 || jj >= nc) continue;
                const std::size_t q = (std::size_t)jj * nr + ii;
                if (!on[q]) continue;
                const int a = rm_find(parent, (int)p), b = rm_find(parent, (int)q);
                if (a != b) parent[b] = a;
            }
        }
    }
    std::vector<double> mass(N, 0.0);
    std::vector<int> cells(N, 0);
    for (std::size_t p = 0; p < N; ++p) {
        if (!on[p]) continue;
        const int r = rm_find(parent, (int)p);
        mass[r] += excess[p];
        cells[r] += 1;
    }
    double best = 0.0;
    for (std::size_t p = 0; p < N; ++p)
        if (on[p] && parent[p] == (int)p && cells[p] >= min_cells && mass[p] > best)
            best = mass[p];
    return best * increase * increase;
}

//' Maximum regional excess mass in both directions
//'
//' @param z landscape matrix, NA outside the silhouette.
//' @param threshold distance from the neutral value 0.5 a cell must exceed.
//' @param min_cells minimum number of cells of a region.
//' @param increase grid step, in coordinate units.
//' @return Numeric vector c(over, under).
//' @keywords internal
// [[Rcpp::export(rng = false)]]
NumericVector max_region_mass(NumericMatrix z, double threshold,
                              int min_cells, double increase)
{
    NumericVector out = NumericVector::create(
        _["over"]  = side_mass(z, true,  threshold, min_cells, increase),
        _["under"] = side_mass(z, false, threshold, min_cells, increase));
    return out;
}
