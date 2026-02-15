#include "grouping_graph_solver.hpp"
#include "datatype.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct ArcPrimitive {
    std::array<double, 2> p0;
    std::array<double, 2> p1;
    std::array<double, 2> dir;
    std::array<double, 2> mid;
    double length;
    int polarity;
};

struct EdgeFeature {
    int u;
    int v;
    std::array<double, 5> feat;
    double weight;
};

struct DSU {
    std::vector<int> p;
    std::vector<int> r;
    explicit DSU(int n) : p(n), r(n, 0) {
        std::iota(p.begin(), p.end(), 0);
    }
    int find(int x) {
        if (p[x] == x) return x;
        p[x] = find(p[x]);
        return p[x];
    }
    void unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a == b) return;
        if (r[a] < r[b]) std::swap(a, b);
        p[b] = a;
        if (r[a] == r[b]) r[a]++;
    }
};

double clamp01(double x) {
    return std::max(0.0, std::min(1.0, x));
}

double norm2(const std::array<double, 2>& a) {
    return std::sqrt(a[0] * a[0] + a[1] * a[1]);
}

std::array<double, 2> sub(const std::array<double, 2>& a, const std::array<double, 2>& b) {
    return {a[0] - b[0], a[1] - b[1]};
}

double dot(const std::array<double, 2>& a, const std::array<double, 2>& b) {
    return a[0] * b[0] + a[1] * b[1];
}

double cross(const std::array<double, 2>& a, const std::array<double, 2>& b) {
    return a[0] * b[1] - a[1] * b[0];
}

std::array<double, 2> normalize(const std::array<double, 2>& v) {
    const double n = norm2(v);
    if (n < 1e-9) return {0.0, 0.0};
    return {v[0] / n, v[1] / n};
}

std::array<double, 5> parseWeights() {
    std::array<double, 5> w{{1.00, 0.90, 0.60, 0.55, 0.70}};
    const char* env = std::getenv("LU_GROUPING_EDGE_WEIGHTS");
    if (!env) return w;
    std::stringstream ss(env);
    std::string token;
    int i = 0;
    while (std::getline(ss, token, ',') && i < 5) {
        w[i++] = std::atof(token.c_str());
    }
    return w;
}

std::vector<ArcPrimitive> buildArcs(double* lines, int line_num) {
    std::vector<ArcPrimitive> arcs(line_num);
    for (int i = 0; i < line_num; ++i) {
        ArcPrimitive a;
        a.p0 = {lines[i * TUPLELENGTH + 0], lines[i * TUPLELENGTH + 1]};
        a.p1 = {lines[i * TUPLELENGTH + 2], lines[i * TUPLELENGTH + 3]};
        a.dir = normalize({lines[i * TUPLELENGTH + 4], lines[i * TUPLELENGTH + 5]});
        a.length = std::max(1e-6, lines[i * TUPLELENGTH + 6]);
        a.polarity = static_cast<int>(lines[i * TUPLELENGTH + 7]);
        a.mid = {(a.p0[0] + a.p1[0]) * 0.5, (a.p0[1] + a.p1[1]) * 0.5};
        arcs[i] = a;
    }
    return arcs;
}

bool intersectNormalLines(const ArcPrimitive& a, const ArcPrimitive& b, std::array<double,2>* center) {
    std::array<double,2> na = {-a.dir[1] * a.polarity, a.dir[0] * a.polarity};
    std::array<double,2> nb = {-b.dir[1] * b.polarity, b.dir[0] * b.polarity};
    const double det = cross(na, nb);
    if (std::abs(det) < 1e-6) return false;
    auto d = sub(b.mid, a.mid);
    const double t = cross(d, nb) / det;
    (*center) = {a.mid[0] + t * na[0], a.mid[1] + t * na[1]};
    return true;
}

EdgeFeature computeEdgeFeature(const ArcPrimitive& a, const ArcPrimitive& b, int u, int v, double spatial_radius) {
    EdgeFeature ef;
    ef.u = u;
    ef.v = v;

    const double d00 = norm2(sub(a.p0, b.p0));
    const double d01 = norm2(sub(a.p0, b.p1));
    const double d10 = norm2(sub(a.p1, b.p0));
    const double d11 = norm2(sub(a.p1, b.p1));
    const double min_end = std::min(std::min(d00, d01), std::min(d10, d11));
    const double endpoint_proximity = clamp01(1.0 - min_end / spatial_radius);

    const double tcont = clamp01((std::abs(dot(a.dir, b.dir)) - 0.2) / 0.8);

    auto dm = sub(b.mid, a.mid);
    const double dmid = std::max(1e-6, norm2(dm));
    auto dm_n = normalize(dm);
    const double ka = std::abs(cross(a.dir, dm_n)) / dmid;
    const double kb = std::abs(cross(b.dir, dm_n)) / dmid;
    const double curvature_compat = std::exp(-10.0 * std::abs(ka - kb));

    double shared_center = 0.0;
    double local_fit_gain = 0.0;
    std::array<double,2> c{};
    if (intersectNormalLines(a, b, &c)) {
        const double ra = norm2(sub(a.mid, c));
        const double rb = norm2(sub(b.mid, c));
        shared_center = std::exp(-std::abs(ra - rb) / (0.35 * spatial_radius + 1e-6));

        std::array<double, 4> rs{{norm2(sub(a.p0, c)), norm2(sub(a.p1, c)), norm2(sub(b.p0, c)), norm2(sub(b.p1, c))}};
        const double mean_r = (rs[0] + rs[1] + rs[2] + rs[3]) * 0.25;
        double var = 0.0;
        for (size_t i = 0; i < rs.size(); ++i) {
            const double dr = rs[i] - mean_r;
            var += dr * dr;
        }
        var /= 4.0;
        local_fit_gain = std::exp(-var / (0.12 * spatial_radius * spatial_radius + 1e-6));
    }

    ef.feat = {{endpoint_proximity, tcont, curvature_compat, shared_center, local_fit_gain}};
    ef.weight = 0.0;
    return ef;
}

std::vector<EdgeFeature> buildSparseEdges(const std::vector<ArcPrimitive>& arcs, int imgx, int imgy) {
    const int n = static_cast<int>(arcs.size());
    const double diag = std::sqrt(static_cast<double>(imgx * imgx + imgy * imgy));
    const double spatial_radius = std::max(16.0, 0.1 * diag);
    const int kmax = 8;

    std::vector<EdgeFeature> edges;
    for (int i = 0; i < n; ++i) {
        std::vector<std::pair<double, int>> neighbors;
        neighbors.reserve(n);
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            const double d = norm2(sub(arcs[i].mid, arcs[j].mid));
            if (d <= spatial_radius) neighbors.push_back({d, j});
        }
        std::sort(neighbors.begin(), neighbors.end(), [](const std::pair<double, int>& a, const std::pair<double, int>& b){ return a.first < b.first; });
        const int use = std::min(kmax, static_cast<int>(neighbors.size()));
        for (int t = 0; t < use; ++t) {
            int j = neighbors[t].second;
            if (i < j) {
                edges.push_back(computeEdgeFeature(arcs[i], arcs[j], i, j, spatial_radius));
            }
        }
    }
    return edges;
}

std::vector<int> labelsFromDSU(DSU& dsu, int n) {
    std::map<int, int> remap;
    std::vector<int> labels(n, -1);
    int next = 0;
    for (int i = 0; i < n; ++i) {
        int root = dsu.find(i);
        if (!remap.count(root)) remap[root] = next++;
        labels[i] = remap[root];
    }
    return labels;
}

double localNodeEnergy(int node, int candidate_label, const std::vector<int>& labels, const std::vector<EdgeFeature>& edges) {
    double e = 0.0;
    for (size_t i = 0; i < edges.size(); ++i) {
        const EdgeFeature& ed = edges[i];
        int other = -1;
        if (ed.u == node) other = ed.v;
        else if (ed.v == node) other = ed.u;
        if (other < 0) continue;
        const bool same = (candidate_label == labels[other]);
        if (ed.weight >= 0) {
            e += same ? 0.0 : ed.weight;
        } else {
            e += same ? -ed.weight : 0.0;
        }
    }
    return e;
}

std::vector<int> solveGreedy(int n, std::vector<EdgeFeature> edges) {
    std::sort(edges.begin(), edges.end(), [](const EdgeFeature& a, const EdgeFeature& b) { return a.weight > b.weight; });
    DSU dsu(n);
    for (size_t i = 0; i < edges.size(); ++i) {
        if (edges[i].weight <= 0) break;
        if (dsu.find(edges[i].u) != dsu.find(edges[i].v)) {
            dsu.unite(edges[i].u, edges[i].v);
        }
    }
    return labelsFromDSU(dsu, n);
}

std::vector<int> solveMulticutApprox(int n, const std::vector<EdgeFeature>& edges) {
    std::vector<int> labels = solveGreedy(n, edges);
    int next_label = *std::max_element(labels.begin(), labels.end()) + 1;

    for (int iter = 0; iter < 6; ++iter) {
        bool changed = false;
        for (int v = 0; v < n; ++v) {
            std::vector<int> candidates;
            candidates.push_back(labels[v]);
            for (size_t i = 0; i < edges.size(); ++i) {
                if (edges[i].u == v) candidates.push_back(labels[edges[i].v]);
                if (edges[i].v == v) candidates.push_back(labels[edges[i].u]);
            }
            candidates.push_back(next_label);
            std::sort(candidates.begin(), candidates.end());
            candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());

            int best_label = labels[v];
            double best_energy = localNodeEnergy(v, labels[v], labels, edges);
            for (size_t c = 0; c < candidates.size(); ++c) {
                double e = localNodeEnergy(v, candidates[c], labels, edges);
                if (e + 1e-9 < best_energy) {
                    best_energy = e;
                    best_label = candidates[c];
                }
            }
            if (best_label != labels[v]) {
                labels[v] = best_label;
                if (best_label == next_label) ++next_label;
                changed = true;
            }
        }
        if (!changed) break;
    }

    std::map<int, int> remap;
    int nxt = 0;
    for (int i = 0; i < n; ++i) {
        if (!remap.count(labels[i])) remap[labels[i]] = nxt++;
        labels[i] = remap[labels[i]];
    }
    return labels;
}

GroupingErrorBreakdown evaluateErrors(const std::vector<int>& labels, const std::vector<EdgeFeature>& edges) {
    GroupingErrorBreakdown eb;
    for (size_t i = 0; i < edges.size(); ++i) {
        const EdgeFeature& ed = edges[i];
        const bool same = labels[ed.u] == labels[ed.v];
        if (ed.weight >= 0) {
            if (!same) eb.attractive_cut_error += ed.weight;
        } else {
            if (same) eb.repulsive_join_error += -ed.weight;
        }
    }
    eb.total_energy = eb.attractive_cut_error + eb.repulsive_join_error;
    return eb;
}

} // namespace

void groupLSsGraphSolver(
    double *lines,
    int line_num,
    int imgx,
    int imgy,
    std::vector<std::vector<int>> *groups,
    const std::string& optimizer,
    GroupingRuntimeBreakdown* runtime,
    GroupingErrorBreakdown* error) {

    auto t0 = std::chrono::high_resolution_clock::now();
    if (line_num <= 0) {
        groups->clear();
        if (runtime) runtime->total_ms = 0.0;
        if (error) *error = GroupingErrorBreakdown();
        return;
    }

    std::vector<ArcPrimitive> arcs = buildArcs(lines, line_num);
    auto t1 = std::chrono::high_resolution_clock::now();

    std::vector<EdgeFeature> edges = buildSparseEdges(arcs, imgx, imgy);
    auto weights = parseWeights();
    for (size_t i = 0; i < edges.size(); ++i) {
        double score = 0.0;
        for (int k = 0; k < 5; ++k) score += weights[k] * edges[i].feat[k];
        edges[i].weight = score - 1.55;
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    std::vector<int> labels;
    if (optimizer == "greedy") labels = solveGreedy(line_num, edges);
    else labels = solveMulticutApprox(line_num, edges);
    auto t3 = std::chrono::high_resolution_clock::now();

    std::map<int, std::vector<int>> by_group;
    for (int i = 0; i < line_num; ++i) by_group[labels[i]].push_back(i);

    groups->clear();
    groups->reserve(by_group.size());
    for (std::map<int, std::vector<int>>::iterator it = by_group.begin(); it != by_group.end(); ++it) {
        groups->push_back(it->second);
    }

    if (runtime) {
        runtime->feature_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        runtime->graph_build_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
        runtime->optimization_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();
        runtime->total_ms = std::chrono::duration<double, std::milli>(t3 - t0).count();
    }
    if (error) {
        *error = evaluateErrors(labels, edges);
    }
}
