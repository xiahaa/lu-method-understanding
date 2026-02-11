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
extern double * mylsd(int * n_out, double * img, int X, int Y, int ** reg_img, int * reg_x, int * reg_y);
extern void groupLSs(double *lines, int lines_num, int *reg, int reg_x, int reg_y, std::vector<std::vector<int>> *groups);
extern void calcuGroupCoverage(double *lines, int dim, std::vector<std::vector<int>> &groups, double *&coverages);
extern void calculateGradient2(double *img, int imgx, int imgy, image_double *angles);
extern void calculateGradient3(double *img, int imgx, int imgy, image_double *angles);
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
    Mat& lsimg) {
    
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
    
    double* out = mylsd(&n, data, imgx, imgy, &reg, &reg_x, &reg_y);
    groupLSs(out, n, reg, reg_x, reg_y, &groups);
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
    } else {
        // No candidates found
        point5d zero_ellipse;
        zero_ellipse.x = 0;
        zero_ellipse.y = 0;
        zero_ellipse.a = 0;
        zero_ellipse.b = 0;
        zero_ellipse.phi = 0;
        candidates.push_back(zero_ellipse);
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
