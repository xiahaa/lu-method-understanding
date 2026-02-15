#ifndef LOGGING_STATS_HPP
#define LOGGING_STATS_HPP

#include <string>
#include <vector>

struct TriggerFeatureStats {
    double avg_arc_length = 0.0;
    double avg_curvature_fluctuation = 0.0;
    double avg_gradient_direction_variance = 0.0;
    double avg_occlusion_ratio = 0.0;
};

struct CandidateStageStats {
    int arc_count = 0;
    int group_count = 0;
    int validated_count = 0;
};

struct ErrorStats {
    int count = 0;
    TriggerFeatureStats trigger_features;
    std::string rule;
};

struct PipelineLogStats {
    std::string dataset_label;
    std::string scenario_label;
    CandidateStageStats stage_curve;
    ErrorStats over_split;
    ErrorStats over_merge;
    int final_ellipse_count = 0;
};

void writePipelineStatsJson(const PipelineLogStats& stats, const std::string& output_path);

#endif
