#include "ellipse_detection.hpp"
#include <opencv2/core/core.hpp>
#include <algorithm>
#include <cmath>
#include <vector>
#include <numeric>

using namespace cv;

// Helper to compute dot product
static double dotProduct(const Point2d& a, const Point2d& b) {
    return a.x * b.x + a.y * b.y;
}

// Main ellipse detection
void ellipseDetection(
    const std::vector<point5d>& candidates_input,
    const std::vector<Point2d>& points,
    const std::vector<Point2d>& normals,
    double distance_tolerance,
    double normal_tolerance,
    double Tmin,
    double angleCoverage,
    const Mat& E,
    std::vector<int>& mylabels,
    std::vector<int>& labels,
    std::vector<point5d>& ellipses,
    std::vector<std::vector<Point2d>>& pointcloud,
    PipelineLogStats* pipeline_stats) {
    
    int n_points = points.size();
    labels.assign(n_points, 0);
    mylabels.assign(n_points, 0);
    ellipses.clear();
    pointcloud.clear();

    auto clampValue = [](double value, double lower, double upper) {
        return std::max(lower, std::min(upper, value));
    };

    auto candidateFeatures = [](const point5d& e) {
        TriggerFeatureStats f;
        const double major = std::max(e.a, e.b);
        const double minor = std::max(1e-6, std::min(e.a, e.b));
        const double perimeter = M_PI * (1.5 * (major + minor) - std::sqrt(major * minor));
        f.avg_arc_length = perimeter;
        f.avg_curvature_fluctuation = std::fabs(major - minor) / std::max(major, 1e-6);
        f.avg_gradient_direction_variance = std::pow(f.avg_curvature_fluctuation, 2);
        f.avg_occlusion_ratio = clampValue(1.0 - (minor / major), 0.0, 1.0);
        return f;
    };

    auto accumulateFeatures = [&](const std::vector<point5d>& source, const std::vector<bool>& mask) {
        TriggerFeatureStats agg;
        int cnt = 0;
        for (size_t i = 0; i < source.size() && i < mask.size(); ++i) {
            if (!mask[i]) continue;
            TriggerFeatureStats f = candidateFeatures(source[i]);
            agg.avg_arc_length += f.avg_arc_length;
            agg.avg_curvature_fluctuation += f.avg_curvature_fluctuation;
            agg.avg_gradient_direction_variance += f.avg_gradient_direction_variance;
            agg.avg_occlusion_ratio += f.avg_occlusion_ratio;
            cnt++;
        }
        if (cnt > 0) {
            agg.avg_arc_length /= cnt;
            agg.avg_curvature_fluctuation /= cnt;
            agg.avg_gradient_direction_variance /= cnt;
            agg.avg_occlusion_ratio /= cnt;
        }
        return agg;
    };

    std::vector<bool> merged_flags(candidates_input.size(), false);
    
    // Compute goodness for each candidate
    std::vector<double> goodness(candidates_input.size(), 0.0);
    
    for (size_t i = 0; i < candidates_input.size(); i++) {
        Point2d ellipseCenter(candidates_input[i].x, candidates_input[i].y);
        Point2d ellipseAxes(candidates_input[i].a, candidates_input[i].b);
        
        int tbins = std::min(180.0, floor(M_PI * (1.5 * (ellipseAxes.x + ellipseAxes.y) - 
                                                  sqrt(ellipseAxes.x * ellipseAxes.y)) * Tmin));
        
        // Find points within bounding box
        std::vector<int> s_dx;
        for (int j = 0; j < n_points; j++) {
            if (points[j].x >= (ellipseCenter.x - ellipseAxes.x - 1) &&
                points[j].x <= (ellipseCenter.x + ellipseAxes.x + 1) &&
                points[j].y >= (ellipseCenter.y - ellipseAxes.x - 1) &&
                points[j].y <= (ellipseCenter.y + ellipseAxes.x + 1)) {
                s_dx.push_back(j);
            }
        }
        
        // Filter using Rosin distance
        std::vector<Point2d> s_dx_points;
        for (int idx : s_dx) {
            s_dx_points.push_back(points[idx]);
        }
        
        std::vector<double> rosin_dists = dRosinSquare(candidates_input[i], s_dx_points);
        std::vector<int> inliers;
        for (size_t j = 0; j < s_dx.size(); j++) {
            if (rosin_dists[j] <= 1.0) {
                inliers.push_back(s_dx[j]);
            }
        }
        
        if (inliers.empty()) continue;
        
        // Compute ellipse normals
        std::vector<Point2d> inlier_points;
        for (int idx : inliers) {
            inlier_points.push_back(points[idx]);
        }
        
        std::vector<Point2d> ellipse_normals = computePointAngle(candidates_input[i], inlier_points);
        
        // Filter using normal consistency
        std::vector<double> p_dot_temp;
        for (size_t j = 0; j < inliers.size(); j++) {
            p_dot_temp.push_back(dotProduct(normals[inliers[j]], ellipse_normals[j]));
        }
        
        int p_cnt = 0;
        for (double dot : p_dot_temp) {
            if (dot > 0) p_cnt++;
        }
        
        std::vector<int> filtered_inliers;
        if (p_cnt > inliers.size() * 0.5) {
            for (size_t j = 0; j < inliers.size(); j++) {
                if (p_dot_temp[j] > 0 && p_dot_temp[j] >= 0.923879532511287) {
                    filtered_inliers.push_back(inliers[j]);
                }
            }
        } else {
            for (size_t j = 0; j < inliers.size(); j++) {
                if (p_dot_temp[j] < 0 && (-p_dot_temp[j]) >= 0.923879532511287) {
                    filtered_inliers.push_back(inliers[j]);
                }
            }
        }
        
        inliers = filtered_inliers;
        
        if (inliers.empty()) continue;
        
        // Take inliers
        std::vector<Point2d> inlier_pts;
        for (int idx : inliers) {
            inlier_pts.push_back(points[idx]);
        }
        
        std::vector<bool> take_mask = takeInliers(inlier_pts, ellipseCenter, tbins);
        std::vector<int> filtered_inliers2;
        for (size_t j = 0; j < inliers.size(); j++) {
            if (take_mask[j]) {
                filtered_inliers2.push_back(inliers[j]);
            }
        }
        inliers = filtered_inliers2;
        
        if (inliers.empty()) continue;
        
        // Calculate goodness
        double perimeter = M_PI * (1.5 * (ellipseAxes.x + ellipseAxes.y) - 
                                  sqrt(ellipseAxes.x * ellipseAxes.y));
        double support_inliers_ratio = inliers.size() / floor(perimeter);
        
        std::vector<Point2d> inlier_pts2;
        for (int idx : inliers) {
            inlier_pts2.push_back(points[idx]);
        }
        double completeness_ratio = calcuCompleteness(inlier_pts2, ellipseCenter, tbins) / 360.0;
        
        goodness[i] = sqrt(support_inliers_ratio * completeness_ratio);
    }
    
    // Sort candidates by goodness
    std::vector<size_t> goodness_index(candidates_input.size());
    std::iota(goodness_index.begin(), goodness_index.end(), 0);
    std::sort(goodness_index.begin(), goodness_index.end(),
              [&goodness](size_t i1, size_t i2) { return goodness[i1] > goodness[i2]; });
    
    std::vector<point5d> candidates;
    for (size_t idx : goodness_index) {
        if (goodness[idx] > 0) {
            candidates.push_back(candidates_input[idx]);
        }
    }
    
    // Process with different angle thresholds
    std::vector<double> angles = {300, 210, 150, 90};
    
    // Remove angles less than angleCoverage
    angles.erase(std::remove_if(angles.begin(), angles.end(),
                                [angleCoverage](double a) { return a < angleCoverage; }),
                 angles.end());
    
    if (angles.empty() || angles.back() != angleCoverage) {
        angles.push_back(angleCoverage);
    }
    
    for (size_t angleLoop = 0; angleLoop < angles.size(); angleLoop++) {
        // Find unlabeled points
        std::vector<int> idx;
        for (int i = 0; i < n_points; i++) {
            if (labels[i] == 0) {
                idx.push_back(i);
            }
        }
        
        if (idx.size() < 2 * M_PI * (6 * distance_tolerance) * Tmin) {
            break;
        }
        
        // Get subset of points and normals
        std::vector<Point2d> subset_points, subset_normals;
        for (int i : idx) {
            subset_points.push_back(points[i]);
            subset_normals.push_back(normals[i]);
        }
        
        // Run sub-detection
        std::vector<int> L2, L;
        std::vector<point5d> C;
        std::vector<bool> validCandidates;
        std::vector<std::vector<Point2d>> pcl;
        
        subEllipseDetection(candidates, subset_points, subset_normals,
                          distance_tolerance, normal_tolerance, Tmin,
                          angles[angleLoop], E, angleLoop,
                          L2, L, C, validCandidates, pcl);
        
        // Update candidates
        std::vector<point5d> new_candidates;
        for (size_t i = 0; i < candidates.size(); i++) {
            if (validCandidates[i]) {
                new_candidates.push_back(candidates[i]);
            }
        }
        candidates = new_candidates;
        
        // Process detected ellipses
        for (size_t i = 0; i < C.size(); i++) {
            bool flag = false;
            
            for (size_t j = 0; j < ellipses.size(); j++) {
                double center_dist = sqrt(pow(C[i].x - ellipses[j].x, 2) + 
                                         pow(C[i].y - ellipses[j].y, 2));
                double axes_dist = sqrt(pow(C[i].a - ellipses[j].a, 2) + 
                                       pow(C[i].b - ellipses[j].b, 2));
                double phi_diff = fabs(C[i].phi - ellipses[j].phi);
                
                if (center_dist <= distance_tolerance &&
                    axes_dist <= distance_tolerance &&
                    phi_diff <= 0.314159265358979) {
                    flag = true;
                    if (i < merged_flags.size()) {
                        merged_flags[i] = true;
                    }
                    
                    // Update labels
                    for (size_t k = 0; k < idx.size(); k++) {
                        if (L[k] == (int)(i + 1)) {
                            labels[idx[k]] = j + 1;
                        }
                        if (L2[k] == (int)(i + 1)) {
                            mylabels[idx[k]] = j + 1;
                        }
                    }
                    break;
                }
            }
            
            if (!flag) {
                // Add new ellipse
                for (size_t k = 0; k < idx.size(); k++) {
                    if (L[k] == (int)(i + 1)) {
                        labels[idx[k]] = ellipses.size() + 1;
                    }
                    if (L2[k] == (int)(i + 1)) {
                        mylabels[idx[k]] = ellipses.size() + 1;
                    }
                }
                
                ellipses.push_back(C[i]);
                if (i < pcl.size()) {
                    pointcloud.push_back(pcl[i]);
                }
            }
        }
    }

    if (pipeline_stats != nullptr) {
        std::vector<bool> split_flags(candidates_input.size(), false);
        for (size_t i = 0; i < candidates_input.size(); ++i) {
            if (goodness[i] > 0 && !merged_flags[i]) {
                split_flags[i] = true;
            }
        }

        int over_merge_count = 0;
        for (bool f : merged_flags) {
            if (f) over_merge_count++;
        }

        int over_split_count = 0;
        for (bool f : split_flags) {
            if (f) over_split_count++;
        }

        pipeline_stats->over_split.count = over_split_count;
        pipeline_stats->over_merge.count = over_merge_count;
        pipeline_stats->over_split.trigger_features = accumulateFeatures(candidates_input, split_flags);
        pipeline_stats->over_merge.trigger_features = accumulateFeatures(candidates_input, merged_flags);
        pipeline_stats->over_split.rule =
            "arc_length<th_perimeter || curvature_fluctuation>th_curvature || gradient_direction_variance>th_grad_var || occlusion_ratio>th_occlusion";
        pipeline_stats->over_merge.rule =
            "arc_length>th_perimeter && curvature_fluctuation<th_curvature && gradient_direction_variance<th_grad_var && occlusion_ratio<th_occlusion";
        pipeline_stats->final_ellipse_count = static_cast<int>(ellipses.size());
    }
}
