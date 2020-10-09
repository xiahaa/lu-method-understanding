
% done, take all points of ellipse, compute angle, make a histogram and count the number in each bin,
% the completeness is actually how many bins are nonzero, times a scalar
% 360, this is to convert to a completeness in terms of 360 coverage
function [completeness] = calcuCompleteness(x, center, tbins)
    [theta, ~] = cart2pol(x(:, 1) - center(1), x(:, 2) - center(2));%theta??(-pi,pi)????????num x 1
    tmin = -pi; tmax = pi;
    tt = round((theta - tmin) / (tmax - tmin) * tbins + 0.5);%theta????i????????????j??bin????tt??i????????j??????num x 1
    tt(tt < 1) = 1; tt(tt > tbins) = tbins;
    h = histc(tt, 1 : tbins);
    h_greatthanzero_num = sum(h>0);
    completeness = h_greatthanzero_num*(360 / tbins);
end