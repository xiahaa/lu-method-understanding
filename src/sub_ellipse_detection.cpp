#include "ellipse_detection.hpp"
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <algorithm>
#include <cmath>
#include <vector>
#include <numeric>

using namespace cv;

// Helper function to compute dot product for vectors
static double dotProduct(const Point2d& a, const Point2d& b) {
    return a.x * b.x + a.y * b.y;
}

// Sub-ellipse detection
void subEllipseDetection(
    const std::vector<point5d>& candidates,
    const std::vector<Point2d>& points,
    const std::vector<Point2d>& normals,
    double distance_tolerance,
    double normal_tolerance,
    double Tmin,
    double angleCoverage,
    const Mat& E,
    int angleLoop,
    std::vector<int>& L2,
    std::vector<int>& L,
    std::vector<point5d>& C,
    std::vector<bool>& validCandidates,
    std::vector<std::vector<Point2d>>& pcl) {
    
    int n_points = points.size();
    L.assign(n_points, 0);
    L2.assign(n_points, 0);
    C.clear();
    pcl.clear();
    validCandidates.assign(candidates.size(), true);
    
    std::vector<int> labels(n_points, 0);
    std::vector<int> mylabels(n_points, 0);
    std::vector<point5d> ellipses;
    int ellipse_polarity = 0;
    
    Point2d max_pt = points[0];
    Point2d min_pt = points[0];
    for (const auto& pt : points) {
        max_pt.x = std::max(max_pt.x, pt.x);
        max_pt.y = std::max(max_pt.y, pt.y);
        min_pt.x = std::min(min_pt.x, pt.x);
        min_pt.y = std::min(min_pt.y, pt.y);
    }
    
    Point2d max_dis = max_pt - min_pt;
    double maxSemiMajor = std::max(max_dis.x, max_dis.y);
    double maxSemiMinor = std::min(max_dis.x, max_dis.y);
    double distance_tolerance_square = distance_tolerance * distance_tolerance;
    
    std::vector<point5d> convergence = candidates;
    
    for (size_t i = 0; i < candidates.size(); i++) {
        Point2d ellipseCenter(candidates[i].x, candidates[i].y);
        Point2d ellipseAxes(candidates[i].a, candidates[i].b);
        double ellipsePhi = candidates[i].phi;
        
        int tbins = std::min(180.0, floor(M_PI * (1.5 * (ellipseAxes.x + ellipseAxes.y) - 
                                                  sqrt(ellipseAxes.x * ellipseAxes.y)) * Tmin));
        
        // Find points within bounding box
        std::vector<int> i_dx;
        for (int j = 0; j < n_points; j++) {
            if (points[j].x >= (ellipseCenter.x - ellipseAxes.x - distance_tolerance) &&
                points[j].x <= (ellipseCenter.x + ellipseAxes.x + distance_tolerance) &&
                points[j].y >= (ellipseCenter.y - ellipseAxes.x - distance_tolerance) &&
                points[j].y <= (ellipseCenter.y + ellipseAxes.x + distance_tolerance)) {
                i_dx.push_back(j);
            }
        }
        
        // Filter using Rosin distance
        std::vector<Point2d> i_dx_points;
        for (int idx : i_dx) {
            i_dx_points.push_back(points[idx]);
        }
        
        std::vector<double> rosin_dists = dRosinSquare(candidates[i], i_dx_points);
        std::vector<int> inliers;
        for (size_t j = 0; j < i_dx.size(); j++) {
            if (labels[i_dx[j]] == 0 && rosin_dists[j] <= distance_tolerance_square) {
                inliers.push_back(i_dx[j]);
            }
        }
        
        if (inliers.empty()) {
            validCandidates[i] = false;
            continue;
        }
        
        // Get points and normals for inliers
        std::vector<Point2d> inlier_points, inlier_normals;
        for (int idx : inliers) {
            inlier_points.push_back(points[idx]);
            inlier_normals.push_back(normals[idx]);
        }
        
        // Compute ellipse normals
        std::vector<Point2d> ellipse_normals = computePointAngle(candidates[i], inlier_points);
        
        // Filter using normal consistency
        std::vector<double> p_dot_temp;
        for (size_t j = 0; j < inliers.size(); j++) {
            p_dot_temp.push_back(dotProduct(inlier_normals[j], ellipse_normals[j]));
        }
        
        int p_cnt = 0;
        for (double dot : p_dot_temp) {
            if (dot > 0) p_cnt++;
        }
        
        std::vector<int> filtered_inliers;
        if (p_cnt > inliers.size() * 0.5) {
            ellipse_polarity = -1;
            for (size_t j = 0; j < inliers.size(); j++) {
                if (p_dot_temp[j] > 0 && p_dot_temp[j] >= 0.923879532511287) {
                    filtered_inliers.push_back(inliers[j]);
                }
            }
        } else {
            ellipse_polarity = 1;
            for (size_t j = 0; j < inliers.size(); j++) {
                if (p_dot_temp[j] < 0 && (-p_dot_temp[j]) >= 0.923879532511287) {
                    filtered_inliers.push_back(inliers[j]);
                }
            }
        }
        
        inliers = filtered_inliers;
        
        if (inliers.empty()) {
            validCandidates[i] = false;
            continue;
        }
        
        // Take inliers based on angular coverage
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
        
        if (inliers.empty()) {
            validCandidates[i] = false;
            continue;
        }
        
        // Fit new ellipse
        std::vector<double> X, Y;
        for (int idx : inliers) {
            X.push_back(points[idx].x);
            Y.push_back(points[idx].y);
        }
        
        point5d new_ellipse;
        bool new_info = fitEllipse(X, Y, new_ellipse);
        
        std::vector<int> inliers3;
        
        if (new_info) {
            // Check if new ellipse is close to original
            double center_dist = sqrt(pow(new_ellipse.x - ellipseCenter.x, 2) + 
                                     pow(new_ellipse.y - ellipseCenter.y, 2));
            double axes_dist = sqrt(pow(new_ellipse.a - ellipseAxes.x, 2) + 
                                   pow(new_ellipse.b - ellipseAxes.y, 2));
            double phi_diff = fabs(new_ellipse.phi - ellipsePhi);
            
            if (center_dist <= 4 * distance_tolerance &&
                axes_dist <= 4 * distance_tolerance &&
                phi_diff <= 0.314159265358979) {
                
                // Find new inliers
                ellipse_normals = computePointAngle(new_ellipse, points);
                
                std::vector<int> i_dx2;
                for (int j = 0; j < n_points; j++) {
                    if (points[j].x >= (new_ellipse.x - new_ellipse.a - distance_tolerance) &&
                        points[j].x <= (new_ellipse.x + new_ellipse.a + distance_tolerance) &&
                        points[j].y >= (new_ellipse.y - new_ellipse.a - distance_tolerance) &&
                        points[j].y <= (new_ellipse.y + new_ellipse.a + distance_tolerance)) {
                        i_dx2.push_back(j);
                    }
                }
                
                std::vector<Point2d> i_dx2_points;
                for (int idx : i_dx2) {
                    i_dx2_points.push_back(points[idx]);
                }
                
                rosin_dists = dRosinSquare(new_ellipse, i_dx2_points);
                std::vector<int> newinliers;
                for (size_t j = 0; j < i_dx2.size(); j++) {
                    int idx = i_dx2[j];
                    if (labels[idx] == 0 && rosin_dists[j] <= distance_tolerance_square) {
                        double dot_val = dotProduct(normals[idx], ellipse_normals[idx]);
                        if ((dot_val * (-ellipse_polarity)) >= 0.923879532511287) {
                            newinliers.push_back(idx);
                        }
                    }
                }
                
                // Filter newinliers
                std::vector<Point2d> newinlier_pts;
                for (int idx : newinliers) {
                    newinlier_pts.push_back(points[idx]);
                }
                
                take_mask = takeInliers(newinlier_pts, Point2d(new_ellipse.x, new_ellipse.y), tbins);
                std::vector<int> filtered_newinliers;
                for (size_t j = 0; j < newinliers.size(); j++) {
                    if (take_mask[j]) {
                        filtered_newinliers.push_back(newinliers[j]);
                    }
                }
                newinliers = filtered_newinliers;
                
                if (newinliers.size() >= inliers.size()) {
                    inliers = newinliers;
                    inliers3 = newinliers;
                    
                    // Refit
                    X.clear();
                    Y.clear();
                    for (int idx : inliers) {
                        X.push_back(points[idx].x);
                        Y.push_back(points[idx].y);
                    }
                    
                    point5d new_new_ellipse;
                    bool new_new_info = fitEllipse(X, Y, new_new_ellipse);
                    if (new_new_info) {
                        new_ellipse = new_new_ellipse;
                    }
                }
            }
        } else {
            new_ellipse = candidates[i];
        }
        
        // Check if enough inliers
        double min_inliers = floor(M_PI * (1.5 * (new_ellipse.a + new_ellipse.b) - 
                                          sqrt(new_ellipse.a * new_ellipse.b)) * Tmin);
        
        if (inliers.size() >= min_inliers) {
            convergence[i] = new_ellipse;
            
            // Check if duplicate in convergence
            bool is_duplicate = false;
            for (size_t j = 0; j < i; j++) {
                double center_dist = sqrt(pow(convergence[j].x - new_ellipse.x, 2) + 
                                         pow(convergence[j].y - new_ellipse.y, 2));
                double axes_dist = sqrt(pow(convergence[j].a - new_ellipse.a, 2) + 
                                       pow(convergence[j].b - new_ellipse.b, 2));
                double phi_diff = fabs(convergence[j].phi - new_ellipse.phi);
                
                if (center_dist <= distance_tolerance &&
                    axes_dist <= distance_tolerance &&
                    phi_diff <= 0.314159265358979) {
                    is_duplicate = true;
                    break;
                }
            }
            
            if (is_duplicate) {
                validCandidates[i] = false;
            }
            
            // Check completeness
            std::vector<Point2d> inlier_pts2;
            for (int idx : inliers) {
                inlier_pts2.push_back(points[idx]);
            }
            bool completeOrNot = calcuCompleteness(inlier_pts2, Point2d(new_ellipse.x, new_ellipse.y), tbins) >= angleCoverage;
            
            // Check if valid ellipse
            if (new_info && new_ellipse.a < maxSemiMajor && new_ellipse.b < maxSemiMinor && completeOrNot) {
                // Check if no similar ellipse exists
                bool no_similar = true;
                for (const auto& existing : ellipses) {
                    double center_dist = sqrt(pow(existing.x - new_ellipse.x, 2) + 
                                             pow(existing.y - new_ellipse.y, 2));
                    double axes_dist = sqrt(pow(existing.a - new_ellipse.a, 2) + 
                                           pow(existing.b - new_ellipse.b, 2));
                    double phi_diff = fabs(existing.phi - new_ellipse.phi);
                    
                    if (center_dist <= distance_tolerance &&
                        axes_dist <= distance_tolerance &&
                        phi_diff < 0.314159265358979) {
                        no_similar = false;
                        break;
                    }
                }
                
                if (no_similar) {
                    for (int idx : inliers) {
                        labels[idx] = ellipses.size() + 1;
                    }
                    
                    if (!inliers3.empty()) {
                        for (int idx : inliers3) {
                            mylabels[idx] = ellipses.size() + 1;
                        }
                    }
                    
                    ellipses.push_back(new_ellipse);
                    
                    std::vector<Point2d> pc;
                    for (int idx : inliers) {
                        pc.push_back(points[idx]);
                    }
                    pcl.push_back(pc);
                    
                    validCandidates[i] = false;
                }
            }
        } else {
            validCandidates[i] = false;
        }
    }
    
    // Set output labels
    for (int i = 0; i < n_points; i++) {
        L2[i] = mylabels[i];
        L[i] = labels[i];
    }
    
    C = ellipses;
}
