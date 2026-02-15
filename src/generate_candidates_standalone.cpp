#include "generate_candidates.hpp"
#include "datatype.h"
#include "utils.hpp"
#include "image.hpp"
#include "elsd.hpp"
#include "group_forming.hpp"
#include "gen_init_set.hpp"
#include "clustering.hpp"
#include "gradient.hpp"
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <vector>
#include <cstring>

using namespace cv;

// Forward declarations from elsd.hpp and other headers
extern double * LineSegmentDetection(int * n_out, double * img, int X, int Y, double scale, double sigma_scale,
                                     double quant, double ang_th, double log_eps, double density_th, int n_bins,
                                     int ** reg_img, int * reg_x, int * reg_y);
extern void groupLSs(double *lines, int lines_num, int *reg, int reg_x, int reg_y, std::vector<std::vector<int>> *groups);
extern void calcuGroupCoverage(double *lines, int dim, std::vector<std::vector<int>> groups, double *&coverages);
extern void calculateGradient2(double *img, unsigned int imgx, unsigned int imgy, image_double *angles);
extern void calculateGradient3(double *img, unsigned int imgx, unsigned int imgy, image_double *angles);
extern PairGroupList * getValidInitialEllipseSet(double *lines, int linenum, std::vector<std::vector<int>> *groups,
                                                 double *coverages, image_double angles, double distance_tolerance,
                                                 int specified_polarity);
extern void generateEllipseCandidates(PairGroupList * pairGroupList, double distance_tolerance,
                                     double *&candidates, int *candidates_num);
extern void freePairGroupList(PairGroupList * pairGroupList);

// Standalone version of generateEllipseCandidates without MEX
void generateEllipseCandidatesStandalone(
    const Mat& input_image,
    int edge_process_select,
    int specified_polarity,
    std::vector<point5d>& candidates,
    Mat& edge_image,
    std::vector<Point2d>& normals,
    Mat& lsimg,
    PipelineLogStats* pipeline_stats) {
    
    candidates.clear();
    normals.clear();
    
    // Convert Mat to double array
    Mat gray;
    if (input_image.channels() > 1) {
        cvtColor(input_image, gray, COLOR_BGR2GRAY);
    } else {
        gray = input_image.clone();
    }
    
    if (gray.type() != CV_8U) {
        gray.convertTo(gray, CV_8U);
    }
    
    int imgy = gray.rows;
    int imgx = gray.cols;
    
    double *data = (double*)malloc(imgy * imgx * sizeof(double));
    for (int c = 0; c < imgx; c++) {
        for (int r = 0; r < imgy; r++) {
            data[c + r * imgx] = gray.at<uchar>(r, c);
        }
    }
    
    int n;
    std::vector<std::vector<int>> groups;
    double * coverages;
    int * reg;
    int reg_x;
    int reg_y;
    
    // LSD parameters
    double scale = 0.8;
    double sigma_scale = 0.6;
    double quant = 2.0;
    double ang_th = 22.5;
    double log_eps = 0.0;
    double density_th = 0.7;
    int n_bins = 1024;
    
    double* out = LineSegmentDetection(&n, data, imgx, imgy, scale, sigma_scale, quant,
                                       ang_th, log_eps, density_th, n_bins, &reg, &reg_x, &reg_y);
    groupLSs(out, n, reg, reg_x, reg_y, &groups);
    if (pipeline_stats != nullptr) {
        pipeline_stats->stage_curve.arc_count = n;
        pipeline_stats->stage_curve.group_count = static_cast<int>(groups.size());
    }
    free(reg);
    calcuGroupCoverage(out, n, groups, coverages);
    
    image_double angles;
    if (edge_process_select == 1)
        calculateGradient2(data, imgx, imgy, &angles);
    else
        calculateGradient3(data, imgx, imgy, &angles);
    
    PairGroupList * pairGroupList;
    double distance_tolerance = 2;
    double * candidates_data;
    int candidates_num = 0;
    
    pairGroupList = getValidInitialEllipseSet(out, n, &groups, coverages, angles, 
                                             distance_tolerance, specified_polarity);

    if (pipeline_stats != nullptr) {
        pipeline_stats->stage_curve.validated_count =
            (pairGroupList != NULL) ? pairGroupList->length : 0;
    }
    
    if (pairGroupList != NULL) {
        generateEllipseCandidates(pairGroupList, distance_tolerance, candidates_data, &candidates_num);
        
        // Convert to vector
        for (int i = 0; i < candidates_num; i++) {
            point5d ellipse;
            ellipse.x = candidates_data[i * 5 + 0];
            ellipse.y = candidates_data[i * 5 + 1];
            ellipse.a = candidates_data[i * 5 + 2];
            ellipse.b = candidates_data[i * 5 + 3];
            ellipse.phi = candidates_data[i * 5 + 4];
            candidates.push_back(ellipse);
        }
        
        freePairGroupList(pairGroupList);
        free(candidates_data);
    }
    
    // Create edge image
    edge_image = Mat::zeros(imgy, imgx, CV_8U);
    
    // Extract normals
    for (int c = 0; c < imgx; c++) {
        for (int r = 0; r < imgy; r++) {
            unsigned long addr = r * imgx + c;
            if (angles->data[addr] != NOTDEF) {
                edge_image.at<uchar>(r, c) = 255;
                double angle = angles->data[addr];
                normals.push_back(Point2d(cos(angle), sin(angle)));
            }
        }
    }
    
    // Create line segment image (optional, can be left empty for now)
    lsimg = Mat::zeros(imgy, imgx, CV_8U);
    
    // Cleanup
    free_image_double(angles);
    free(out);
    free(data);
    free(coverages);
}
