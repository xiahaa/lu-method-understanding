#ifndef ELLIPSE_DETECTION_HPP
#define ELLIPSE_DETECTION_HPP

#include <vector>
#include "opencv2/core/core.hpp"
#include "datatype.h"

// Compute point normals for an ellipse
std::vector<cv::Point2d> computePointAngle(const point5d& ellipse, const std::vector<cv::Point2d>& points);

// Compute Rosin distance squared
std::vector<double> dRosinSquare(const point5d& param, const std::vector<cv::Point2d>& points);

// Take inliers based on angular coverage
std::vector<bool> takeInliers(const std::vector<cv::Point2d>& x, const cv::Point2d& center, int tbins);

// Calculate completeness ratio
double calcuCompleteness(const std::vector<cv::Point2d>& x, const cv::Point2d& center, int tbins);

// Fit ellipse using least squares
bool fitEllipse(const std::vector<double>& X, const std::vector<double>& Y, point5d& ellipse);

// Fit ellipse using Ahn's method
point5d fitAhn(const std::vector<double>& Xi, const std::vector<double>& Yi, const point5d& ellipara);

// Sub-ellipse detection
void subEllipseDetection(
    const std::vector<point5d>& candidates,
    const std::vector<cv::Point2d>& points,
    const std::vector<cv::Point2d>& normals,
    double distance_tolerance,
    double normal_tolerance,
    double Tmin,
    double angleCoverage,
    const cv::Mat& E,
    int angleLoop,
    std::vector<int>& L2,
    std::vector<int>& L,
    std::vector<point5d>& C,
    std::vector<bool>& validCandidates,
    std::vector<std::vector<cv::Point2d>>& pcl);

// Main ellipse detection
void ellipseDetection(
    const std::vector<point5d>& candidates,
    const std::vector<cv::Point2d>& points,
    const std::vector<cv::Point2d>& normals,
    double distance_tolerance,
    double normal_tolerance,
    double Tmin,
    double angleCoverage,
    const cv::Mat& E,
    std::vector<int>& mylabels,
    std::vector<int>& labels,
    std::vector<point5d>& ellipses,
    std::vector<std::vector<cv::Point2d>>& pointcloud);

// Main ellipse detection LU
void ellipseDetectionLU(
    const cv::Mat& I,
    double Tac,
    double Tr,
    int specified_polarity,
    bool verbose,
    std::vector<point5d>& ellipses,
    std::vector<point5d>& posi,
    cv::Mat& L,
    std::vector<std::vector<cv::Point2d>>& points);

// Draw ellipses on image
void drawEllipses(const std::vector<point5d>& ellipses, cv::Mat& image);

#endif // ELLIPSE_DETECTION_HPP
