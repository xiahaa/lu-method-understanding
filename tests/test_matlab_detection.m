function test_matlab_detection()
% Test MATLAB ellipse detection and save results
% This script runs the MATLAB implementation and saves the detected ellipses
% to a CSV file for comparison with C++ implementation

clc; close all;

% Test image path
test_image = 'test_data/test_ellipses.png';

% Check if test image exists
if ~exist(test_image, 'file')
    error('Test image not found: %s', test_image);
end

% Parameters (same as default in C++ main)
Tac = 165;  % Elliptic angular coverage threshold
Tr = 0.6;   % Ratio of support inliers
specified_polarity = 0;  % 0: all, 1: positive, -1: negative

% Read image
fprintf('------MATLAB Ellipse Detection Test------\n');
fprintf('Reading image: %s\n', test_image);
I = imread(test_image);

% Detect ellipses using MATLAB implementation
fprintf('Detecting ellipses with parameters: Tac=%g, Tr=%g, polarity=%d\n', ...
    Tac, Tr, specified_polarity);

try
    [ellipses, ~, ~, pcl] = ellipseDetectionLU(I, Tac, Tr, specified_polarity, false);
    
    % Refine ellipses using Ahn's method
    fprintf('Refining ellipses using Ahn''s method\n');
    for i = 1:length(pcl)
        if ~isempty(pcl{i})
            points = pcl{i};
            ellipses(i,:) = fitAhn(points(:,1), points(:,2), ellipses(i,:));
        end
    end
    
    % Convert angle to degrees
    ellipses(:,5) = ellipses(:,5) ./ pi * 180;
    
    % Display results
    fprintf('------MATLAB Results------\n');
    fprintf('Detected %d ellipses:\n', size(ellipses, 1));
    fprintf('Format: (x, y, a, b, phi_degrees)\n');
    for i = 1:size(ellipses, 1)
        fprintf('%d: %.2f, %.2f, %.2f, %.2f, %.2f\n', ...
            i, ellipses(i,1), ellipses(i,2), ellipses(i,3), ...
            ellipses(i,4), ellipses(i,5));
    end
    
    % Save results to CSV file
    output_file = 'test_data/matlab_results.csv';
    fprintf('Saving results to: %s\n', output_file);
    
    % Create header
    fid = fopen(output_file, 'w');
    fprintf(fid, 'x,y,a,b,phi_degrees\n');
    for i = 1:size(ellipses, 1)
        fprintf(fid, '%.6f,%.6f,%.6f,%.6f,%.6f\n', ...
            ellipses(i,1), ellipses(i,2), ellipses(i,3), ...
            ellipses(i,4), ellipses(i,5));
    end
    fclose(fid);
    
    fprintf('MATLAB test completed successfully!\n');
    
    % Optionally draw ellipses
    drawEllipses(ellipses', I);
    saveas(gcf, 'test_data/matlab_detection_result.png');
    
catch ME
    fprintf('Error during MATLAB detection: %s\n', ME.message);
    fprintf('Stack trace:\n');
    disp(ME.stack);
    rethrow(ME);
end

end
