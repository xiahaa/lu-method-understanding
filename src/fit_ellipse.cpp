#include "ellipse_detection.hpp"
#include <opencv2/core/core.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <vector>

using namespace cv;
using namespace Eigen;

// Fit ellipse using least squares method
// Returns true if successful, false otherwise
bool fitEllipse(const std::vector<double>& X, const std::vector<double>& Y, point5d& ellipse) {
    if (X.size() != Y.size() || X.size() < 6) {
        return false;
    }
    
    int n = X.size();
    
    // Build design matrix
    MatrixXd D(n, 6);
    for (int i = 0; i < n; i++) {
        double x = X[i];
        double y = Y[i];
        D(i, 0) = x * x;
        D(i, 1) = x * y;
        D(i, 2) = y * y;
        D(i, 3) = x;
        D(i, 4) = y;
        D(i, 5) = 1.0;
    }
    
    // Build scatter matrix
    MatrixXd S = D.transpose() * D;
    
    // Build 6x6 constraint matrix
    MatrixXd C = MatrixXd::Zero(6, 6);
    C(0, 2) = 2;
    C(1, 1) = -1;
    C(2, 0) = 2;
    
    // Solve generalized eigenvalue problem
    GeneralizedSelfAdjointEigenSolver<MatrixXd> solver(S, C);
    
    if (solver.info() != Success) {
        return false;
    }
    
    VectorXd eigenvalues = solver.eigenvalues();
    MatrixXd eigenvectors = solver.eigenvectors();
    
    // Find the positive (or near-zero) eigenvalue
    int I = -1;
    for (int i = 0; i < eigenvalues.size(); i++) {
        if (eigenvalues(i) > -1e-8 && !std::isinf(eigenvalues(i))) {
            I = i;
            break;
        }
    }
    
    if (I == -1) {
        return false;
    }
    
    // Extract eigenvector
    VectorXd A = eigenvectors.col(I).real();
    
    // Convert to geometric parameters
    std::vector<double> par(6);
    for (int i = 0; i < 6; i++) {
        par[i] = A(i);
    }
    
    double thetarad = 0.5 * atan2(par[1], par[0] - par[2]);
    
    double cost = cos(thetarad);
    double sint = sin(thetarad);
    double sin_squared = sint * sint;
    double cos_squared = cost * cost;
    double cos_sin = sint * cost;
    
    double Ao = par[5];
    double Au = par[3] * cost + par[4] * sint;
    double Av = -par[3] * sint + par[4] * cost;
    double Auu = par[0] * cos_squared + par[2] * sin_squared + par[1] * cos_sin;
    double Avv = par[0] * sin_squared + par[2] * cos_squared - par[1] * cos_sin;
    
    double tuCentre = -Au / (2.0 * Auu);
    double tvCentre = -Av / (2.0 * Avv);
    double wCentre = Ao - Auu * tuCentre * tuCentre - Avv * tvCentre * tvCentre;
    
    double uCentre = tuCentre * cost - tvCentre * sint;
    double vCentre = tuCentre * sint + tvCentre * cost;
    
    double Ru = -wCentre / Auu;
    double Rv = -wCentre / Avv;
    
    Ru = sqrt(fabs(Ru));
    Rv = sqrt(fabs(Rv));
    
    ellipse.x = uCentre;
    ellipse.y = vCentre;
    ellipse.a = Ru;
    ellipse.b = Rv;
    ellipse.phi = thetarad;
    
    // Ensure a >= b
    if (Ru < Rv) {
        ellipse.a = Rv;
        ellipse.b = Ru;
        if (thetarad < 0) {
            ellipse.phi = ellipse.phi + M_PI / 2;
        } else {
            ellipse.phi = ellipse.phi - M_PI / 2;
        }
    }
    
    return true;
}
