#ifndef GENERATE_CANDIDATES_HPP
#define GENERATE_CANDIDATES_HPP

#include <vector>
#include "opencv2/core/core.hpp"
#include "datatype.h"
#include "logging_stats.hpp"

// Generate ellipse candidates from image
// Returns: candidates (vector of point5d), edge image, normals, and line segment image
void generateEllipseCandidatesStandalone(
    const cv::Mat& input_image,
    int edge_process_select,  // 1: sobel, 2: canny
    int specified_polarity,   // 1: positive, -1: negative, 0: all
    std::vector<point5d>& candidates,
    cv::Mat& edge_image,
    std::vector<cv::Point2d>& normals,
    cv::Mat& lsimg,
    PipelineLogStats* pipeline_stats = nullptr);

#endif // GENERATE_CANDIDATES_HPP
