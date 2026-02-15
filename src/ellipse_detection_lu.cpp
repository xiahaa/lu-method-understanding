#include "ellipse_detection.hpp"
#include "generate_candidates.hpp"
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <iostream>
#include <chrono>

using namespace cv;

// Main ellipse detection LU method
void ellipseDetectionLU(
    const Mat& I,
    double Tac,
    double Tr,
    int specified_polarity,
    bool verbose,
    std::vector<point5d>& ellipses,
    std::vector<point5d>& posi,
    Mat& L,
    std::vector<std::vector<Point2d>>& points,
    PipelineLogStats* pipeline_stats) {
    
    double angleCoverage = Tac;  // default 165
    double Tmin = Tr;            // default 0.6
    double unit_dis_tolerance = 2;  // max([2, 0.005 * min([size(I, 1), size(I, 2)])])
    double normal_tolerance = M_PI / 9;  // 20 degrees = pi/9

    if (pipeline_stats != nullptr) {
        if (pipeline_stats->dataset_label.empty()) pipeline_stats->dataset_label = "unspecified";
        if (pipeline_stats->scenario_label.empty()) pipeline_stats->scenario_label = "unspecified";
    }
    
    auto t0 = std::chrono::high_resolution_clock::now();
    
    Mat gray_image;
    if (I.channels() > 1) {
        cvtColor(I, gray_image, COLOR_BGR2GRAY);
    } else {
        gray_image = I.clone();
    }
    
    // Generate ellipse candidates
    std::vector<point5d> candidates;
    Mat edge;
    std::vector<Point2d> normals;
    Mat lsimg;
    
    generateEllipseCandidatesStandalone(gray_image, 2, specified_polarity, 
                                       candidates, edge, normals, lsimg, pipeline_stats);
    
    auto t1 = std::chrono::high_resolution_clock::now();
    
    if (verbose) {
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
        std::cout << "The time of generating ellipse candidates: " 
                  << duration.count() / 1000.0 << "s" << std::endl;
        std::cout << "Init ellipse candidates: " << candidates.size() << std::endl;
    }
    
    if (candidates.empty()) {
        candidates.clear();
        ellipses.clear();
        posi.clear();
        L = Mat::zeros(I.size(), CV_32S);
        points.clear();
        if (pipeline_stats != nullptr) {
            pipeline_stats->final_ellipse_count = 0;
        }
        return;
    }
    
    posi = candidates;
    
    // Find edge points
    std::vector<Point2d> edge_points;
    for (int y = 0; y < edge.rows; y++) {
        for (int x = 0; x < edge.cols; x++) {
            if (edge.at<uchar>(y, x) > 0) {
                edge_points.push_back(Point2d(x, y));
            }
        }
    }
    
    // Run main detection
    std::vector<int> mylabels, labels;
    std::vector<std::vector<Point2d>> pcl;
    
    ellipseDetection(candidates, edge_points, normals, 
                    unit_dis_tolerance, normal_tolerance, Tmin, 
                    angleCoverage, gray_image,
                    mylabels, labels, ellipses, pcl, pipeline_stats);
    
    if (verbose) {
        auto t2 = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t0);
        std::cout << "-----------------------------------------------------------" << std::endl;
        std::cout << "Running time: " << total_duration.count() / 1000.0 << "s" << std::endl;
    }
    
    // Create label image
    L = Mat::zeros(I.rows, I.cols, CV_32S);
    for (size_t i = 0; i < edge_points.size(); i++) {
        int x = static_cast<int>(edge_points[i].x);
        int y = static_cast<int>(edge_points[i].y);
        if (x >= 0 && x < L.cols && y >= 0 && y < L.rows) {
            L.at<int>(y, x) = mylabels[i];
        }
    }
    
    points = pcl;
}
