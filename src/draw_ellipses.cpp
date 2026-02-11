#include "ellipse_detection.hpp"
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace cv;

// Draw ellipses on image
void drawEllipses(const std::vector<point5d>& ellipses, Mat& image) {
    // Create output image if needed
    Mat displayImage;
    if (image.channels() == 1) {
        cvtColor(image, displayImage, COLOR_GRAY2BGR);
    } else {
        displayImage = image.clone();
    }
    
    // Draw each ellipse
    for (const auto& ellipse : ellipses) {
        Point2d center(ellipse.x, ellipse.y);
        Size2d axes(ellipse.a, ellipse.b);
        double angle = ellipse.phi * 180.0 / M_PI; // Convert to degrees
        
        // Draw the ellipse with random color for better visibility
        Scalar color(rand() % 256, rand() % 256, rand() % 256);
        cv::ellipse(displayImage, center, axes, angle, 0, 360, color, 2, LINE_AA);
        
        // Draw the center
        circle(displayImage, center, 3, color, -1);
    }
    
    // Display the image
    imshow("Detected Ellipses", displayImage);
    waitKey(1);
}
