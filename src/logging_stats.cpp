#include "logging_stats.hpp"
#include <fstream>
#include <iomanip>

void writePipelineStatsJson(const PipelineLogStats& stats, const std::string& output_path) {
    std::ofstream os(output_path.c_str());
    if (!os.is_open()) {
        return;
    }

    os << std::fixed << std::setprecision(6);
    os << "{\n";
    os << "  \"schema_version\": \"ellipse_pipeline_stats/v1\",\n";
    os << "  \"dataset_label\": \"" << stats.dataset_label << "\",\n";
    os << "  \"scenario_label\": \"" << stats.scenario_label << "\",\n";
    os << "  \"stage_curve\": {\n";
    os << "    \"arc_count\": " << stats.stage_curve.arc_count << ",\n";
    os << "    \"group_count\": " << stats.stage_curve.group_count << ",\n";
    os << "    \"validated_count\": " << stats.stage_curve.validated_count << "\n";
    os << "  },\n";
    os << "  \"errors\": {\n";
    os << "    \"over_split\": {\n";
    os << "      \"count\": " << stats.over_split.count << ",\n";
    os << "      \"trigger_rule\": \"" << stats.over_split.rule << "\",\n";
    os << "      \"trigger_features\": {\n";
    os << "        \"arc_length\": " << stats.over_split.trigger_features.avg_arc_length << ",\n";
    os << "        \"curvature_fluctuation\": " << stats.over_split.trigger_features.avg_curvature_fluctuation << ",\n";
    os << "        \"gradient_direction_variance\": " << stats.over_split.trigger_features.avg_gradient_direction_variance << ",\n";
    os << "        \"occlusion_ratio\": " << stats.over_split.trigger_features.avg_occlusion_ratio << "\n";
    os << "      }\n";
    os << "    },\n";
    os << "    \"over_merge\": {\n";
    os << "      \"count\": " << stats.over_merge.count << ",\n";
    os << "      \"trigger_rule\": \"" << stats.over_merge.rule << "\",\n";
    os << "      \"trigger_features\": {\n";
    os << "        \"arc_length\": " << stats.over_merge.trigger_features.avg_arc_length << ",\n";
    os << "        \"curvature_fluctuation\": " << stats.over_merge.trigger_features.avg_curvature_fluctuation << ",\n";
    os << "        \"gradient_direction_variance\": " << stats.over_merge.trigger_features.avg_gradient_direction_variance << ",\n";
    os << "        \"occlusion_ratio\": " << stats.over_merge.trigger_features.avg_occlusion_ratio << "\n";
    os << "      }\n";
    os << "    }\n";
    os << "  },\n";
    os << "  \"final_ellipse_count\": " << stats.final_ellipse_count << "\n";
    os << "}\n";
}
