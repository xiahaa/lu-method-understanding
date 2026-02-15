#ifndef GROUPING_GRAPH_SOLVER_HPP
#define GROUPING_GRAPH_SOLVER_HPP

#include <vector>
#include <string>

struct GroupingRuntimeBreakdown {
    double feature_ms = 0.0;
    double graph_build_ms = 0.0;
    double optimization_ms = 0.0;
    double total_ms = 0.0;
};

struct GroupingErrorBreakdown {
    double attractive_cut_error = 0.0;
    double repulsive_join_error = 0.0;
    double total_energy = 0.0;
};

void groupLSsGraphSolver(
    double *lines,
    int line_num,
    int imgx,
    int imgy,
    std::vector<std::vector<int>> *groups,
    const std::string& optimizer,
    GroupingRuntimeBreakdown* runtime,
    GroupingErrorBreakdown* error);

#endif
