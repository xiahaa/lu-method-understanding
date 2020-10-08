clc;close all;clear all;
%image path
filename = './pics/35.jpg';

% read image 
disp('------read image------');
I = imread(filename);


%% detecting ellipses from real-world images
% [ellipses, ~, posi] = ellipseDetectionByArcSupportLSs(I, Tac, Tr, specified_polarity);
[ellipses,ellconics] = ellDetect(I);

disp('draw detected ellipses');
drawEllipses(ellipses',I);
% display
ellipses(:,5) = ellipses(:,5)./pi*180;
ellipses
disp(['The total number of detected ellipses??',num2str(size(ellipses,1))]);

% id = ellipses(:,5)<0;
% ellipses(id,5) = 360+ellipses(id,5);

%% i think the difference lays on the drawing of opencv functionality
for i = 1:size(ellipses)
    I = cv.ellipse(I, (ellipses(i,[1 2])),ellipses(i,[3 4]),'Angle',ellipses(5),'Color',[255,0,0]','Thickness',2);
    imshow(I);
end





