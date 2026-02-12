#include "ellipse_detection.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <opencv2/core/core.hpp>

using namespace cv;

// Compute the points' normals belonging to an ellipse, the normals have been already normalized
// param: [x0 y0 a b phi]
// points: [xi yi], n x 2
std::vector<Point2d> computePointAngle(const point5d& ellipse, const std::vector<Point2d>& points) {
    // Convert [x0 y0 a b phi] to Ax^2+Bxy+Cy^2+Dx+Ey+F = 0
    double a_square = ellipse.a * ellipse.a;
    double b_square = ellipse.b * ellipse.b;
    double sin_phi = sin(ellipse.phi);
    double cos_phi = cos(ellipse.phi);
    double sin_square = sin_phi * sin_phi;
    double cos_square = cos_phi * cos_phi;
    
    double A = b_square * cos_square + a_square * sin_square;
    double B = (b_square - a_square) * sin_phi * cos_phi * 2;
    double C = b_square * sin_square + a_square * cos_square;
    double D = -2 * A * ellipse.x - B * ellipse.y;
    double E = -2 * C * ellipse.y - B * ellipse.x;
    
    // Calculate points' normals to ellipse
    std::vector<Point2d> ellipse_normals;
    ellipse_normals.reserve(points.size());
    
    for (const auto& pt : points) {
        double angle = atan2(C * pt.y + B / 2 * pt.x + E / 2, 
                            A * pt.x + B / 2 * pt.y + D / 2);
        ellipse_normals.push_back(Point2d(cos(angle), sin(angle)));
    }
    
    return ellipse_normals;
}

// Compute Rosin distance squared
// param: [x0 y0 a b Phi], 1 x 5 or 5 x 1
// points: points to compute rosin distance, (xi, yi), size: n x 2
// dmin: minimum distance
std::vector<double> dRosinSquare(const point5d& param, const std::vector<Point2d>& points) {
    double ae2 = param.a * param.a;
    double be2 = param.b * param.b;
    double fe2 = ae2 - be2;
    
    std::vector<double> dmin;
    dmin.reserve(points.size());
    
    for (const auto& pt : points) {
        double x = pt.x - param.x;
        double y = pt.y - param.y;
        double xp = x * cos(-param.phi) - y * sin(-param.phi);
        double yp = x * sin(-param.phi) + y * cos(-param.phi);
        
        double X = xp * xp;
        double Y = yp * yp;
        double delta = (X + Y + fe2) * (X + Y + fe2) - 4 * fe2 * X;
        double A = (X + Y + fe2 - sqrt(delta)) / 2;
        double ah = sqrt(A);
        double bh2 = fe2 - A;
        double term = A * be2 + ae2 * bh2;
        double xi = ah * sqrt(ae2 * (be2 + bh2) / term);
        double yi = param.b * sqrt(bh2 * (ae2 - A) / term);
        
        double d[4];
        d[0] = (xp - xi) * (xp - xi) + (yp - yi) * (yp - yi);
        d[1] = (xp - xi) * (xp - xi) + (yp + yi) * (yp + yi);
        d[2] = (xp + xi) * (xp + xi) + (yp - yi) * (yp - yi);
        d[3] = (xp + xi) * (xp + xi) + (yp + yi) * (yp + yi);
        
        dmin.push_back(*std::min_element(d, d + 4));
    }
    
    return dmin;
}

// Calculate completeness
// x: all points of ellipse, 
// center: ellipse center (x, y)
// tbins: number of bins
double calcuCompleteness(const std::vector<Point2d>& x, const Point2d& center, int tbins) {
    std::vector<double> theta;
    theta.reserve(x.size());
    
    for (const auto& pt : x) {
        double angle = atan2(pt.y - center.y, pt.x - center.x);
        theta.push_back(angle);
    }
    
    double tmin = -M_PI;
    double tmax = M_PI;
    
    std::vector<int> tt;
    tt.reserve(theta.size());
    
    for (double angle : theta) {
        int bin = static_cast<int>(round((angle - tmin) / (tmax - tmin) * tbins + 0.5));
        if (bin < 1) bin = 1;
        if (bin > tbins) bin = tbins;
        tt.push_back(bin);
    }
    
    std::vector<int> h(tbins + 1, 0);
    for (int bin : tt) {
        h[bin]++;
    }
    
    int h_greatthanzero_num = 0;
    for (int i = 1; i <= tbins; i++) {
        if (h[i] > 0) h_greatthanzero_num++;
    }
    
    double completeness = h_greatthanzero_num * (360.0 / tbins);
    return completeness;
}

// Take inliers based on angular coverage
std::vector<bool> takeInliers(const std::vector<Point2d>& x, const Point2d& center, int tbins) {
    std::vector<double> theta;
    theta.reserve(x.size());
    
    for (const auto& pt : x) {
        double angle = atan2(pt.y - center.y, pt.x - center.x);
        theta.push_back(angle);
    }
    
    double tmin = -M_PI;
    double tmax = M_PI;
    
    std::vector<int> tt;
    tt.reserve(theta.size());
    
    for (double angle : theta) {
        int bin = static_cast<int>(round((angle - tmin) / (tmax - tmin) * tbins + 0.5));
        if (bin < 1) bin = 1;
        if (bin > tbins) bin = tbins;
        tt.push_back(bin);
    }
    
    std::vector<int> h(tbins + 1, 0);
    for (int bin : tt) {
        h[bin]++;
    }
    
    std::vector<int> mark(tbins + 1, 0);
    std::vector<int> compSize;
    int nComps = 0;
    std::vector<int> queue(tbins + 1);
    int du[2] = {-1, 1};
    
    for (int i = 1; i <= tbins; i++) {
        if (h[i] > 0 && mark[i] == 0) {
            nComps++;
            mark[i] = nComps;
            int front = 0, rear = 0;
            queue[front] = i;
            
            while (front <= rear) {
                int u = queue[front];
                front++;
                
                for (int j = 0; j < 2; j++) {
                    int v = u + du[j];
                    if (v == 0) v = tbins;
                    if (v > tbins) v = 1;
                    
                    if (mark[v] == 0 && h[v] > 0) {
                        rear++;
                        queue[rear] = v;
                        mark[v] = nComps;
                    }
                }
            }
            
            int size = 0;
            for (size_t k = 0; k < tt.size(); k++) {
                if (mark[tt[k]] == nComps) size++;
            }
            compSize.push_back(size);
        }
    }
    
    if (compSize.empty()) {
        return std::vector<bool>(x.size(), false);
    }
    
    int maxCompSize = *std::max_element(compSize.begin(), compSize.end());
    std::vector<int> validComps;
    for (size_t i = 0; i < compSize.size(); i++) {
        if (compSize[i] >= maxCompSize * 0.1 && compSize[i] > 10) {
            validComps.push_back(i + 1);
        }
    }
    
    std::vector<int> validBins;
    for (int i = 1; i <= tbins; i++) {
        if (std::find(validComps.begin(), validComps.end(), mark[i]) != validComps.end()) {
            validBins.push_back(i);
        }
    }
    
    std::vector<bool> idx(x.size());
    for (size_t i = 0; i < tt.size(); i++) {
        idx[i] = std::find(validBins.begin(), validBins.end(), tt[i]) != validBins.end();
    }
    
    return idx;
}
