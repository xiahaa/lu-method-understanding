function fit = fitAhn( Xi, Yi, ellipara)
%FITWITHLSO Ellipse fitting with Least-Squares based on orthogonal distance
%
% [1] Sung Joon Ahn, W. Rauh, and M. Recknagel, "Ellipse fitting and
% parameter assessment of circular object targets for robot vision",
% Intelligent Robots and Systems, 1999.
%
% AUTHOR Sebastian Dingler <s.dingler@gmail.com>
%        Karlsruhe Institute of Technology (KIT), Germany
%
% DATE   22.12.2014

Xc = ellipara(1);
Yc = ellipara(2);
a = ellipara(3);
b = ellipara(4);
alpha = ellipara(5);

% step size
lambda=0.1;

for k=1:50
    J = [];
    X_new2 = [];
    for i=1:length(Xi)
        r = XY2xy([Xi(i);Yi(i)],alpha,[Xc;Yc]);
        x_new = getOrthoPoint(r(1),r(2),a,b);
        X_new = xy2XY2(x_new, alpha, [Xc;Yc]);
        xi = XY2xy([Xi(i);Yi(i)],alpha,[Xc;Yc]);
        J = [J;calcJacobianMatrix(a,b,x_new(1),x_new(2),alpha,xi(1),xi(2))];
        X_new2 = [X_new2; [Xi(i);Yi(i)]-X_new];
    end
    r=-pinv(J) * X_new2;
    % update 
    Xc=Xc-lambda*r(1);
    Yc=Yc-lambda*r(2);
    a=a-lambda*r(3);
    b=b-lambda*r(4);
    alpha=alpha-lambda*r(5);
end

fit = [Xc,Yc,a,b,alpha];
end

% Jacobian matrix at the orthogonal contacting point on ellipse
function r = calcJacobianMatrix(a, b, x, y, alpha, xi, yi)
C = cos(alpha);
S = sin(alpha);
R = [C S;-S C];
B1 = [b^2 * x * C-a^2 * y * S; b^2 * (yi-y) * C+a^2 * (xi-x) * S];
B2 = [b^2 * x * S+a^2 * y * C;b^2 * (yi-y) * S-a^2 * (xi-x) * C];
B3 = [a * (b^2-y^2); 2 * a * y * (xi-x)];
B4 = [b * (a^2-x^2);-2 * b * x * (yi-y)];
B5 = [(a^2-b^2) * x * y; (a^2-b^2) * (x^2 - y^2 - x * xi + y * yi)];
B = [B1 B2 B3 B4 B5];
Qk = [b^2 * x a^2*y; (a^2-b^2)*y+b^2*yi (a^2-b^2)*x-a^2*xi];
r = pinv(R)*pinv(Qk)*B;
end

function r = XY2xy(X, alpha, Xc)
%XY2XY Transforms ellipse that it is symmetric to x- and y-axis
%
% AUTHOR Sebastian Dingler <s.dingler@gmail.com>
%        Karlsruhe Institute of Technology (KIT), Germany
%
% DATE   22.12.2014
R = [cos(alpha) sin(alpha);-sin(alpha) cos(alpha)];
r = R*(X-Xc);
end

function r = getOrthoPoint(xi, yi, a, b)
%GETORTHOPOINT Compute orthogonal Point on ellipse
%
% [1] Sung Joon Ahn, W. Rauh, and M. Recknagel, "Ellipse fitting and
% parameter assessment of circular object targets for robot vision",
% Intelligent Robots and Systems, 1999.
%
% AUTHOR Sebastian Dingler <s.dingler@gmail.com>
%        Karlsruhe Institute of Technology (KIT), Germany
%
% DATE   22.12.2014

%Orthogonal contacting point on ellipse
xk1 = ([xi;yi]*a*b)/sqrt(b^2*xi^2+a^2*yi^2);
if abs(xi)<a
    xk2=[xi;sign(yi)*(b/a)*sqrt(a^2-xi^2)];
else
    xk2=[sign(xi)*a;0];
end

x_e = 0.5*(xk1+xk2); % x_0
for i=1:4
    x = x_e(1);
    y = x_e(2);
    Qk = [b^2*x a^2*y;(a^2-b^2)*y+b^2*yi (a^2-b^2)*x-a^2*xi];
    fk = [0.5*(a^2*y^2+b^2*x^2-a^2*b^2);b^2*x*(yi-y)-a^2*y*(xi-x)];
    x_e = x_e-pinv(Qk)*fk;
end
r = x_e;
end

function r = xy2XY2(x, alpha, Xc)
%XY2XY2 Rotates ellipse back (inverse function of XY2xy.m)
%
% AUTHOR Sebastian Dingler <s.dingler@gmail.com>
%        Karlsruhe Institute of Technology (KIT), Germany
%
% DATE   22.12.2014
R = [cos(alpha) sin(alpha);-sin(alpha) cos(alpha)];
r = pinv(R)*x+Xc;
end

