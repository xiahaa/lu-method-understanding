#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "ellipse_detection.hpp"

using namespace cv;
using namespace std;

int main(int argc, char** argv) {
    // Parameters
    double Tac = 165;  // Elliptic angular coverage threshold
    double Tr = 0.6;   // Ratio of support inliers
    int specified_polarity = 0;  // 0: all, 1: positive, -1: negative
    
    // Parse command line arguments
    string filename;
    if (argc > 1) {
        filename = argv[1];
    } else {
        cout << "Usage: " << argv[0] << " <image_path> [Tac] [Tr] [polarity]" << endl;
        cout << "Example: " << argv[0] << " image.jpg 165 0.6 0" << endl;
        return -1;
    }
    
    if (argc > 2) Tac = atof(argv[2]);
    if (argc > 3) Tr = atof(argv[3]);
    if (argc > 4) specified_polarity = atoi(argv[4]);
    
    cout << "------read image------" << endl;
    Mat I = imread(filename);
    
    if (I.empty()) {
        cerr << "Error: Could not read image from " << filename << endl;
        return -1;
    }
    
    cout << "Image size: " << I.cols << " x " << I.rows << endl;
    cout << "Parameters: Tac=" << Tac << ", Tr=" << Tr << ", polarity=" << specified_polarity << endl;
    
    // Detecting ellipses
    cout << "------detecting ellipses------" << endl;
    std::vector<point5d> ellipses;
    std::vector<point5d> posi;
    Mat L;
    std::vector<std::vector<Point2d>> pcl;
    
    ellipseDetectionLU(I, Tac, Tr, specified_polarity, true, ellipses, posi, L, pcl);
    
    cout << "------drawing detected ellipses------" << endl;
    drawEllipses(ellipses, I);
    
    // Refine ellipses using Ahn's method
    cout << "------refining ellipses------" << endl;
    for (size_t i = 0; i < ellipses.size() && i < pcl.size(); i++) {
        if (pcl[i].empty()) continue;
        
        std::vector<double> Xi, Yi;
        for (const auto& pt : pcl[i]) {
            Xi.push_back(pt.x);
            Yi.push_back(pt.y);
        }
        
        ellipses[i] = fitAhn(Xi, Yi, ellipses[i]);
    }
    
    drawEllipses(ellipses, I);
    
    // Display results
    cout << "------results------" << endl;
    cout << "Ellipse parameters (x, y, a, b, phi_degrees):" << endl;
    for (size_t i = 0; i < ellipses.size(); i++) {
        double phi_deg = ellipses[i].phi * 180.0 / M_PI;
        cout << i + 1 << ": " 
             << ellipses[i].x << ", "
             << ellipses[i].y << ", "
             << ellipses[i].a << ", "
             << ellipses[i].b << ", "
             << phi_deg << endl;
    }
    cout << "The total number of detected ellipses: " << ellipses.size() << endl;
    
    // Wait for key press
    cout << "Press any key to exit..." << endl;
    waitKey(0);
    
    return 0;
}
