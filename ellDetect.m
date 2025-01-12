function varargout = ellDetect(im, varargin)
    % parameters
    Tac = 165;
    Tr = 0.6;
    specified_polarity = 0;
    
    % detecting ellipses from real-world images
    % [ellipses, ~, posi] = ellipseDetectionByArcSupportLSs(I, Tac, Tr, specified_polarity);
    [ellipses,~,~,pcl] = ellipseDetectionLU(im, Tac, Tr, specified_polarity,false);
    
    
    
    % finish traditional ellipse detection
    % try some refinement steps

    subpixel = 0;
    othfit = 0;
    if ~isempty(varargin)
        for i = 1:length(varargin)
            if ~isempty(strfind(varargin{i},'Subpixel'))
                subpixel = 1;
            end
            if ~isempty(strfind(varargin{i},'Ahn'))
                othfit = 1;
            end
        end
    end

    if ~isempty(ellipses)
        if size(ellipses,1) == 2
            if ellipses(1,3) > ellipses(2,3)
                id = 1;        
            else
                id = 2;
            end
        else
            id = 1;
        end
    
        % how about use subpixel edges, seems ok
        if subpixel == 1
            poly1 = param2ellipse(ellipses(id,:));
            addpath(genpath('C:\Users\xiahaa\Documents\DTU\publications\fore-end\ellipse_detection\src\3rdparty\Subpixel Matlab v2.11'));
            if size(im,3) == 3
                image = rgb2gray(im);
            end
            threshold = 10;
            edges = subpixelEdges(image, threshold); 
            rmpath(genpath('C:\Users\xiahaa\Documents\DTU\publications\fore-end\ellipse_detection\src\3rdparty\Subpixel Matlab v2.11'));
            x = edges.x; y = edges.y;
            xx = x.^2;xy = x.*y;yy = y.^2;
            f1 = poly1(1).*xx + poly1(2).*xy + poly1(3).*yy + poly1(4).*x + poly1(5).*y + poly1(6);
            valid = abs(f1) < 0.05;% | abs(f2) < 0.01;
        %     %% show edges
        %     visEdges(edges, 'showNormals', false);
            points = [edges.x(valid) edges.y(valid)];
        else
            points = pcl{id};
        end
    
        if othfit == 0
            if subpixel == 0
                % do nothing
            else
                new_ellipse = fitEllipse( points(:,1), points(:,2));
                ellipses(id,:) = new_ellipse;
            end
        else 
            new_ellipse = fitAhn( points(:,1), points(:,2), ellipses(id,:)');
            ellipses(id,:) = new_ellipse;
        end
    end
    % here I want to add a conversion function that both outputs canonical
    % parameters as well as conic parameters in matrix form.
    ellconics = [];
    for i = 1:size(ellipses,1)
        ellconics(i,:) = convertEllParam2Conic(ellipses(i,:));
    end
    varargout{1} = ellipses;
    varargout{2} = ellconics;
end


function poly = param2ellipse(param)
    if length(param)==6,
       error('For ellipse Coefficients to matrix: use coeffs2ellipse');
    end
    % param  = [ xc;yc;ax;bx;rho ]
    xc    = param(1);
    yc    = param(2);
    ax    = param(3);
    bx    = param(4);
    rho   = param(5);
    T     = [ cos(rho),-sin(rho),xc; sin(rho),cos(rho),yc; 0,0,1 ];
    C     = transpose(inv(T))*diag([1/ax^2;1/bx^2;-1])*inv(T);
    % C     = C/norm(C,'fro');
    % from xyabtheta to conic form
    poly = [C(1,1), C(1,2)*2, C(2,2) C(1,3)*2 C(2,3)*2 C(3,3)];
end

function []=simpledrawEllipse(ellipses_para,color)
    th=0:pi/180:2*pi;
    for i=1:size(ellipses_para,2)
        Semi_major= ellipses_para(3,i);
        Semi_minor= ellipses_para(4,i);
        x0= ellipses_para(1,i);
        y0= ellipses_para(2,i);
        Phi= ellipses_para(5,i);
        x=x0+Semi_major*cos(Phi)*cos(th)-Semi_minor*sin(Phi)*sin(th);
        y=y0+Semi_minor*cos(Phi)*sin(th)+Semi_major*sin(Phi)*cos(th);   

        plot(x,y,'Color', color, 'LineWidth',2);
    end
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