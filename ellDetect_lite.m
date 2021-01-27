function varargout = ellDetect_lite(im, varargin)
    % parameters
    Tac = 165;
    Tr = 0.6;
    specified_polarity = 0;
    
    % detecting ellipses from real-world images
%     [ellipses, ~, posi] = ellipseDetectionByArcSupportLSs(I, Tac, Tr, specified_polarity);
    [ellipses] = ellipseDetectionLU(im, Tac, Tr, specified_polarity,false);

    % finish traditional ellipse detection
    % try some refinement steps
    
    varargout{1} = ellipses;
end

