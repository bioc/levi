// One side of threshold-free cluster enhancement on a graph, with connected
// components maintained by union-find while thresholds descend, so a call
// costs O(n_steps * n) instead of rebuilding an induced subgraph at every
// threshold. It reproduces the R reference (.tfceR in leviNetworkStatistics.R)
// to rounding error: the thresholds are seq(dh, top, length.out = n_steps),
// a node is active at h when v >= h, and each node accumulates
// size^E * h^H * dh over the thresholds in ascending order.
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

static int uf_find(std::vector<int>& parent, int a) {
    while (parent[a] != a) {
        parent[a] = parent[parent[a]];
        a = parent[a];
    }
    return a;
}

// v: non-negative statistic per node (NA and negative values are never
// active); edges: m x 2 matrix of 0-based node indices.
// [[Rcpp::export(rng = false)]]
NumericVector tfce_side(NumericVector v, IntegerMatrix edges, double E,
                        double H, int n_steps) {
    const int n = v.size();
    const int m = edges.nrow();
    NumericVector ans(n, 0.0);
    if (n_steps < 1) return ans;

    std::vector<double> val(n);
    double top = 0.0;
    for (int i = 0; i < n; ++i) {
        const double x = v[i];
        val[i] = (R_finite(x) && x > 0.0) ? x : -1.0;
        if (val[i] > top) top = val[i];
    }
    if (!(top > 0.0)) return ans;

    const double dh = top / n_steps;
    std::vector<double> hs(n_steps);
    // seq.default(dh, top, length.out = n_steps) is c(dh, dh + (1:(n - 2)) * by,
    // top): the last value is `top` itself, not dh + (n - 1) * by.
    const double by = (n_steps == 1) ? 0.0 : (top - dh) / (n_steps - 1);
    for (int k = 0; k < n_steps; ++k) hs[k] = dh + k * by;
    if (n_steps >= 2) hs[n_steps - 1] = top;

    // Adjacency in compressed form.
    std::vector<int> start(n + 1, 0);
    for (int e = 0; e < m; ++e) {
        start[edges(e, 0) + 1]++;
        start[edges(e, 1) + 1]++;
    }
    for (int i = 0; i < n; ++i) start[i + 1] += start[i];
    std::vector<int> adj(2 * static_cast<std::size_t>(m));
    std::vector<int> fill(start.begin(), start.end() - 1);
    for (int e = 0; e < m; ++e) {
        adj[fill[edges(e, 0)]++] = edges(e, 1);
        adj[fill[edges(e, 1)]++] = edges(e, 0);
    }

    std::vector<int> order(n);
    for (int i = 0; i < n; ++i) order[i] = i;
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return val[a] > val[b]; });

    std::vector<int> parent(n), size(n, 1);
    for (int i = 0; i < n; ++i) parent[i] = i;
    std::vector<char> on(n, 0);
    // Component size of every active node at every threshold, so that the
    // sums below run in ascending threshold order like the R code.
    std::vector<int> comp(static_cast<std::size_t>(n_steps) * n, 0);

    std::size_t next = 0;
    for (int k = n_steps - 1; k >= 0; --k) {
        const double h = hs[k];
        while (next < static_cast<std::size_t>(n) && val[order[next]] >= h) {
            const int a = order[next++];
            on[a] = 1;
            for (int q = start[a]; q < start[a + 1]; ++q) {
                const int b = adj[q];
                if (!on[b]) continue;
                int ra = uf_find(parent, a), rb = uf_find(parent, b);
                if (ra != rb) {
                    if (size[ra] < size[rb]) std::swap(ra, rb);
                    parent[rb] = ra;
                    size[ra] += size[rb];
                }
            }
        }
        for (std::size_t j = 0; j < next; ++j) {
            const int a = order[j];
            comp[static_cast<std::size_t>(k) * n + a] = size[uf_find(parent, a)];
        }
    }
    for (int k = 0; k < n_steps; ++k) {
        const double w = std::pow(hs[k], H) * dh;
        for (int a = 0; a < n; ++a) {
            const int s = comp[static_cast<std::size_t>(k) * n + a];
            if (s) ans[a] += std::pow(static_cast<double>(s), E) * w;
        }
    }
    return ans;
}
