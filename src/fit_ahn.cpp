#include "ellipse_detection.hpp"
#include <opencv2/core/core.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <vector>

using namespace cv;
using namespace Eigen;

// Helper function to get orthogonal point on ellipse
static MatrixXd getOrthoPoint(const VectorXd& xi, const VectorXd& yi, double a, double b) {
    int n = xi.size();
    MatrixXd result(2, n);
    
    // Orthogonal contacting point on ellipse
    for (int j = 0; j < n; j++) {
        double xk1_x = (xi(j) * a * b) / sqrt(b * b * xi(j) * xi(j) + a * a * yi(j) * yi(j));
        double xk1_y = (yi(j) * a * b) / sqrt(b * b * xi(j) * xi(j) + a * a * yi(j) * yi(j));
        
        double xk2_x, xk2_y;
        if (fabs(xi(j)) < a) {
            xk2_x = xi(j);
            xk2_y = (yi(j) >= 0 ? 1 : -1) * (b / a) * sqrt(a * a - xi(j) * xi(j));
        } else {
            xk2_x = (xi(j) >= 0 ? 1 : -1) * a;
            xk2_y = 0;
        }
        
        double x_e_x = 0.5 * (xk1_x + xk2_x);
        double x_e_y = 0.5 * (xk1_y + xk2_y);
        
        double x = x_e_x;
        double y = x_e_y;
        
        // Iterate 4 times
        for (int i = 0; i < 4; i++) {
            Matrix2d Qk;
            Qk << b * b * x, a * a * y,
                  (a * a - b * b) * y + b * b * yi(j), (a * a - b * b) * x - a * a * xi(j);
            
            Vector2d fk;
            fk << 0.5 * (a * a * y * y + b * b * x * x - a * a * b * b),
                  b * b * x * (yi(j) - y) - a * a * y * (xi(j) - x);
            
            Vector2d r = Vector2d(x_e_x, x_e_y) - Qk.completeOrthogonalDecomposition().solve(fk);
            x = r(0);
            y = r(1);
        }
        
        result(0, j) = x;
        result(1, j) = y;
    }
    
    return result;
}

// Helper function to calculate Jacobian matrix
static MatrixXd calcJacobianMatrix(double a, double b, double x, double y, double alpha, double xi, double yi) {
    double C = cos(alpha);
    double S = sin(alpha);
    
    Matrix2d R;
    R << C, S,
        -S, C;
    
    MatrixXd B(2, 5);
    
    // B1
    B(0, 0) = b * b * x * C - a * a * y * S;
    B(1, 0) = b * b * (yi - y) * C + a * a * (xi - x) * S;
    
    // B2
    B(0, 1) = b * b * x * S + a * a * y * C;
    B(1, 1) = b * b * (yi - y) * S - a * a * (xi - x) * C;
    
    // B3
    B(0, 2) = a * (b * b - y * y);
    B(1, 2) = 2 * a * y * (xi - x);
    
    // B4
    B(0, 3) = b * (a * a - x * x);
    B(1, 3) = -2 * b * x * (yi - y);
    
    // B5
    B(0, 4) = (a * a - b * b) * x * y;
    B(1, 4) = (a * a - b * b) * (x * x - y * y - x * xi + y * yi);
    
    Matrix2d Qk;
    Qk << b * b * x, a * a * y,
          (a * a - b * b) * y + b * b * yi, (a * a - b * b) * x - a * a * xi;
    
    MatrixXd r = R.transpose() * Qk.completeOrthogonalDecomposition().solve(B);
    
    return r;
}

// Fit ellipse using Ahn's method (orthogonal distance)
point5d fitAhn(const std::vector<double>& Xi, const std::vector<double>& Yi, const point5d& ellipara) {
    double Xc = ellipara.x;
    double Yc = ellipara.y;
    double a = ellipara.a;
    double b = ellipara.b;
    double alpha = ellipara.phi;
    
    double lambda = 0.1; // step size
    
    int n = Xi.size();
    MatrixXd XY(2, n);
    for (int i = 0; i < n; i++) {
        XY(0, i) = Xi[i];
        XY(1, i) = Yi[i];
    }
    
    for (int k = 0; k < 20; k++) {
        MatrixXd J(2 * n, 5);
        
        Matrix2d R;
        R << cos(alpha), sin(alpha),
            -sin(alpha), cos(alpha);
        
        MatrixXd r = R * (XY - Vector2d(Xc, Yc).replicate(1, n));
        MatrixXd x_new = getOrthoPoint(r.row(0), r.row(1), a, b);
        MatrixXd X_new = R.transpose() * x_new;
        X_new.row(0).array() += Xc;
        X_new.row(1).array() += Yc;
        MatrixXd X_new2 = XY - X_new;
        
        for (int i = 0; i < n; i++) {
            MatrixXd Ji = calcJacobianMatrix(a, b, x_new(0, i), x_new(1, i), alpha, r(0, i), r(1, i));
            J.block<2, 5>(i * 2, 0) = Ji;
        }
        
        VectorXd X_new2_vec(2 * n);
        for (int i = 0; i < n; i++) {
            X_new2_vec(i * 2) = X_new2(0, i);
            X_new2_vec(i * 2 + 1) = X_new2(1, i);
        }
        
        VectorXd r_vec = -J.completeOrthogonalDecomposition().solve(X_new2_vec);
        
        if ((lambda * r_vec).norm() < 1e-6) {
            break;
        }
        
        // Update
        Xc = Xc - lambda * r_vec(0);
        Yc = Yc - lambda * r_vec(1);
        a = a - lambda * r_vec(2);
        b = b - lambda * r_vec(3);
        alpha = alpha - lambda * r_vec(4);
    }
    
    point5d fit;
    fit.x = Xc;
    fit.y = Yc;
    fit.a = a;
    fit.b = b;
    fit.phi = alpha;
    
    return fit;
}
