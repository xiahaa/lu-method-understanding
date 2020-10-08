function varargout = ellDetect(im)
    % parameters
    Tac = 165;
    Tr = 0.6;
    specified_polarity = 0;
    
    % detecting ellipses from real-world images
    % [ellipses, ~, posi] = ellipseDetectionByArcSupportLSs(I, Tac, Tr, specified_polarity);
    [ellipses] = ellipseDetectionLU(im, Tac, Tr, specified_polarity,false);
    
    % here I want to add a conversion function that both outputs canonical
    % parameters as well as conic parameters in matrix form.
    ellconics = [];
    for i = 1:size(ellipses,1)
        ellconics(i,:) = convertEllParam2Conic(ellipses(i,:));
    end
    varargout{1} = ellipses;
    varargout{2} = ellconics;
end

function ellconic = convertEllParam2Conic(ellParam)
    xc = ellParam(1);
    yc = ellParam(2);
    a = ellParam(3);
    b = ellParam(4);
    theta = ellParam(5);

	cos_theta = cos(theta); sin_theta = sin(theta); sin_2theta = sin(2 * theta);
	pow_cos_theta = cos_theta^2; pow_sin_theta = sin_theta^2;
	aa_inv = 1 / (a*a); bb_inv = 1 / (b*b);

	ellconic(1) = pow_cos_theta*aa_inv + pow_sin_theta*bb_inv;
	ellconic(2) = -0.5 * sin_2theta*(bb_inv - aa_inv);
	ellconic(3) = pow_cos_theta*bb_inv + pow_sin_theta*aa_inv;
	ellconic(4) = (-xc*pow_sin_theta + yc*sin_2theta / 2)*bb_inv - (xc*pow_cos_theta + yc*sin_2theta / 2)*aa_inv;
	ellconic(5) = (-yc*pow_cos_theta + xc*sin_2theta / 2)*bb_inv - (yc*pow_sin_theta + xc*sin_2theta / 2)*aa_inv;
	tmp1 = (xc*cos_theta + yc*sin_theta) / a; tmp2 = (yc*cos_theta - xc*sin_theta) / b;
	ellconic(6) = tmp1^2 + tmp2^2 - 1;

	k = ellconic(1) * ellconic(3) - ellconic(2) * ellconic(2);
	k = 1 / sqrt(abs(k));

    for i=1:6
        ellconic(i) = ellconic(i) * k;
    end
end