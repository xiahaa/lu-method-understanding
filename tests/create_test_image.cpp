#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <iostream>

using namespace cv;
using namespace std;

int main() {
    // Create a white image
    Mat img = Mat::ones(500, 500, CV_8UC3) * 255;
    
    // Convert to grayscale for edge-based detection
    Mat gray = Mat::ones(500, 500, CV_8UC1) * 255;
    
    // Draw thick black ellipses to create strong edges
    // Ellipse 1: center=(150, 150), axes=(60, 40), angle=30
    ellipse(gray, Point(150, 150), Size(60, 40), 30, 0, 360, Scalar(0), 5);
    
    // Ellipse 2: center=(350, 150), axes=(50, 30), angle=0  
    ellipse(gray, Point(350, 150), Size(50, 30), 0, 0, 360, Scalar(0), 5);
    
    // Ellipse 3: center=(250, 350), axes=(70, 50), angle=45
    ellipse(gray, Point(250, 350), Size(70, 50), 45, 0, 360, Scalar(0), 5);
    
    // Convert back to BGR for saving
    cvtColor(gray, img, cv::COLOR_GRAY2BGR);
    
    // Save the test image
    string output_path = "tests/test_data/test_ellipses.png";
    if (imwrite(output_path, img)) {
        cout << "Test image created successfully: " << output_path << endl;
        cout << "Ground truth ellipses:" << endl;
        cout << "1: center=(150, 150), a=60, b=40, angle=30" << endl;
        cout << "2: center=(350, 150), a=50, b=30, angle=0" << endl;
        cout << "3: center=(250, 350), a=70, b=50, angle=45" << endl;
    } else {
        cerr << "Failed to save test image" << endl;
        return -1;
    }
    
    return 0;
}
