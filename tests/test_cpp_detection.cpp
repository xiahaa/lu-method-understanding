#include <iostream>
#include <fstream>
#include <iomanip>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "ellipse_detection.hpp"

using namespace cv;
using namespace std;

int main(int argc, char** argv) {
    // Test image path
    string test_image = "tests/test_data/test_ellipses.png";
    string output_file = "tests/test_data/cpp_results.csv";
    
    // Allow override from command line
    if (argc > 1) {
        test_image = argv[1];
    }
    if (argc > 2) {
        output_file = argv[2];
    }
    
    // Parameters (same as MATLAB test)
    double Tac = 165;  // Elliptic angular coverage threshold
    double Tr = 0.6;   // Ratio of support inliers
    int specified_polarity = 0;  // 0: all, 1: positive, -1: negative
    
    cout << "------C++ Ellipse Detection Test------" << endl;
    cout << "Reading image: " << test_image << endl;
    
    Mat I = imread(test_image);
    
    if (I.empty()) {
        cerr << "Error: Could not read image from " << test_image << endl;
        return -1;
    }
    
    cout << "Image size: " << I.cols << " x " << I.rows << endl;
    cout << "Detecting ellipses with parameters: Tac=" << Tac 
         << ", Tr=" << Tr << ", polarity=" << specified_polarity << endl;
    
    // Detecting ellipses
    std::vector<point5d> ellipses;
    std::vector<point5d> posi;
    Mat L;
    std::vector<std::vector<Point2d>> pcl;
    
    try {
        ellipseDetectionLU(I, Tac, Tr, specified_polarity, false, ellipses, posi, L, pcl);
        
        // Refine ellipses using Ahn's method
        cout << "Refining ellipses using Ahn's method" << endl;
        for (size_t i = 0; i < ellipses.size() && i < pcl.size(); i++) {
            if (pcl[i].empty()) continue;
            
            std::vector<double> Xi, Yi;
            for (const auto& pt : pcl[i]) {
                Xi.push_back(pt.x);
                Yi.push_back(pt.y);
            }
            
            ellipses[i] = fitAhn(Xi, Yi, ellipses[i]);
        }
        
        // Display results
        cout << "------C++ Results------" << endl;
        cout << "Detected " << ellipses.size() << " ellipses:" << endl;
        cout << "Format: (x, y, a, b, phi_degrees)" << endl;
        
        for (size_t i = 0; i < ellipses.size(); i++) {
            double phi_deg = ellipses[i].phi * 180.0 / M_PI;
            cout << i + 1 << ": " 
                 << fixed << setprecision(2)
                 << ellipses[i].x << ", "
                 << ellipses[i].y << ", "
                 << ellipses[i].a << ", "
                 << ellipses[i].b << ", "
                 << phi_deg << endl;
        }
        
        // Save results to CSV file
        cout << "Saving results to: " << output_file << endl;
        
        ofstream ofs(output_file);
        if (!ofs.is_open()) {
            cerr << "Error: Could not open output file " << output_file << endl;
            return -1;
        }
        
        // Write header
        ofs << "x,y,a,b,phi_degrees" << endl;
        
        // Write ellipse parameters
        ofs << fixed << setprecision(6);
        for (const auto& e : ellipses) {
            double phi_deg = e.phi * 180.0 / M_PI;
            ofs << e.x << "," << e.y << "," << e.a << "," 
                << e.b << "," << phi_deg << endl;
        }
        
        ofs.close();
        cout << "C++ test completed successfully!" << endl;
        
        // Optionally draw ellipses (skip display, just save image)
        if (!ellipses.empty()) {
            Mat I_copy = I.clone();
            for (const auto& e : ellipses) {
                cv::ellipse(I_copy, 
                           Point(e.x, e.y), 
                           Size(e.a, e.b), 
                           e.phi * 180.0 / M_PI, 
                           0, 360, 
                           Scalar(0, 0, 255), 
                           2);
            }
            imwrite("tests/test_data/cpp_detection_result.png", I_copy);
        }
        
    } catch (const exception& e) {
        cerr << "Error during C++ detection: " << e.what() << endl;
        return -1;
    }
    
    return 0;
}
