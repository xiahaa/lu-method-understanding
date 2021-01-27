
%% this script connects all part together: from object detection, tracking to classification
close all;clear;
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

detector=configEllipseDetector(1);
addpath(genpath(detector.path))


data_root_path = 'D:\dtu\aamed_ellipse_datasets\'; 

dataset_name = [{'Synthetic Images - Occluded Ellipses'},...
    {'Synthetic Images - Overlap Ellipses'},...
    {'Prasad Images - Dataset Prasad'},...
    {'Random Images - Dataset #1'},...
    {'Smartphone Images - Dataset #2'},...
    {'Concentric Ellipses - Dataset Synthetic'},...
    {'Concurrent Ellipses - Dataset Synthetic'},...
    {'Satellite Images - Dataset Meng #1'},...
    {'Satellite Images - Dataset Meng #2'}];

methods_name = 'Proposed';

gt_label = [{'occluded'},{'overlap'},{'prasad'},{'random'},{'smartphone'},...
    {'concentric'},{'concurrent'},{'satellite1'}, {'satellite2'}];

for dsi = [4] %,1,2,6,7,8,9 %,5,4,3 
    disp(['Evaluating dataset: ',dataset_name{dsi}]);
    imgname_path = [data_root_path,dataset_name{dsi},'\imagenames.txt'];
    fid = fopen(imgname_path,'r');
    imgnum = 0;
    imgname = [];
    while feof(fid) == 0
        imgnum = imgnum + 1;
        imgname{imgnum} = fgetl(fid);
    end
    fclose(fid);

    % folder to save result
%     if ~exist([data_root_path,dataset_name{dsi},'\',methods_name],'dir')
%         mkdir([data_root_path,dataset_name{dsi},'\',methods_name]);
%     end
    
    dirpath = [data_root_path,dataset_name{dsi},'\'];
    f = waitbar(0,['On ', dataset_name{dsi}]);
    gt_elps = Read_Ellipse_GT([dirpath,'gt\'], [dirpath,'images\'], imgname, gt_label{dsi});
    for i = 1:1:imgnum
        waitbar(i/imgnum,f,'Processing');
        
        imfilename = fullfile([data_root_path,dataset_name{dsi},'\images'],imgname{i});
        im = imread(imfilename);
        imshow(im);

        if ~exist(fullfile([data_root_path,dataset_name{dsi},'\Proposed'],[imgname{i},'.txt']), 'file')
            tic
            [ellipses] = detector.call(im,detector.options);
            time = toc;

            drawEllipses(ellipses,im,'r',false);
            gtelps = gt_elps{i};
            drawEllipses(gtelps(:,[1,2,3,4,5]),im,'g',false);
            hold off;
            fid = fopen(fullfile([data_root_path,dataset_name{dsi},'\Proposed'],[imgname{i},'.txt']),'w');
            if ~isempty(ellipses)
                fprintf(fid,'%d, %d, ', 1, size(ellipses,1));
                for j = 1:size(ellipses,1)-1
                    fprintf(fid,'%12.6f, %12.6f, %12.6f, %12.6f, %12.6f, ', ellipses(j,1),ellipses(j,2),ellipses(j,3),ellipses(j,4),ellipses(j,5));
                end
                j = size(ellipses,1);
                fprintf(fid,'%12.6f, %12.6f, %12.6f, %12.6f, %12.6f\n', ellipses(j,1),ellipses(j,2),ellipses(j,3),ellipses(j,4),ellipses(j,5));
            else
                fprintf(fid,'%d, %d, ', 0, -1);
            end
            pause(0.01);
            fclose(fid);
        else
%             continue;
            fid = fopen(fullfile([data_root_path,dataset_name{dsi},'\Proposed'],[imgname{i},'.txt']));
            data = str2num(fgetl(fid));
            if int32(data(1)) == 1
                numell = int32(data(2));
                ellipses = reshape(data(3:end),5,numell)';
                imshow(im);
                drawEllipses([ellipses(:,1:2) ellipses(:,[3,4]) ellipses(:,5)],im,'r',false);
                gtelps = gt_elps{i};
                drawEllipses(gtelps(:,[1,2,3,4,5]),im,'g',false);
                pause(0.01);
                hold off;
                ff=getframe;
                img = ff.cdata;
                imwrite(img,fullfile([data_root_path,dataset_name{dsi},'\Proposed\img'],[imgname{i},'-res.png']));
%                 clf;
            end
            fclose(fid);
        end
    end
    waitbar(1,f,'End...');
    close(f);


end


function [gt_elps, varargout] = Read_Ellipse_GT(gt_path, img_path, imgname, dataset)

imgnum = length(imgname);
gt_elps = cell(1, imgnum);
gt_size = zeros(imgnum,2);

if nargout == 1
    get_size = 0;
else
    get_size = 1;
end

if strcmp(dataset, 'prasad')
    for i = 1:imgnum
        
        if get_size  == 1
            img = imread([img_path,imgname{i}]);
            gt_size(i,1) = size(img,1); gt_size(i,2) = size(img,2);
        end
        
        
        fid_gt = fopen([gt_path,'gt_',imgname{i},'.txt'],'r');
        if fid_gt == -1
            error('gt file error');
        end
        str = fgetl(fid_gt);
        elp_num = str2double(str);
        elp_data = zeros(elp_num,5);
        for j = 1:elp_num
            str = fgetl(fid_gt);
            elp_data(j,:) = str2num(str);
        end
        gt_elps{i}=elp_data;
        fclose(fid_gt);
    end
    
    
    if get_size == 1
        varargout{1} = gt_size;
    end
    
    return;
end

if strcmp(dataset, 'random') || strcmp(dataset, 'smartphone') || strcmp(dataset, 'training')
    for i = 1:imgnum
        
        if get_size  == 1
            img = imread([img_path,imgname{i}]);
            gt_size(i,1) = size(img,1); gt_size(i,2) = size(img,2);
        end
        
        fid_gt = fopen([gt_path,'gt_',imgname{i},'.txt'],'r');
        if fid_gt == -1
            error('gt file error');
        end
        str = fgetl(fid_gt);
        elp_num = str2double(str);
        elp_data = zeros(elp_num,5);
        for j = 1:elp_num
            str = fgetl(fid_gt);
            elp_data(j,:) = str2num(str);
        end
        
        elp_data(:,1:2)=elp_data(:,1:2)+1;
        
        elp_data(:,5) = elp_data(:,5)/180*pi;
        gt_elps{i}=elp_data;
        fclose(fid_gt);
    end
    
    if get_size == 1
        varargout{1} = gt_size;
    end
    
    
    
    
    return;
end


if strcmp(dataset, 'overlap') || strcmp(dataset, 'occluded') || ...
        strcmp(dataset, 'occludedwithmultilines') || strcmp(dataset, 'overlapwithmultilines')
    
    for i = 1:imgnum
        if get_size  == 1
            img = imread([img_path,imgname{i}]);
            gt_size(i,1) = size(img,1); gt_size(i,2) = size(img,2);
        end
        
        gt = load([gt_path,imgname{i}(1:end-4), '.mat']);
        gt_elp = gt.ellipse_param';
        gt_elp(:,5) = gt_elp(:,5)/180*pi;
        
        temp = gt_elp(:,1); gt_elp(:,1) = gt_elp(:,2); gt_elp(:,2) = temp;
        temp = gt_elp(:,3); gt_elp(:,3) = gt_elp(:,4); gt_elp(:,4) = temp;
        temp = pi/2 - gt_elp(:,5);
        gt_elp(:,5) =  temp;
        
        gt_elps{i}=gt_elp;
    end
    
    
    if get_size == 1
        varargout{1} = gt_size;
    end
    
    return;
end


if strcmp(dataset, 'concentric') || strcmp(dataset, 'concurrent')
    
    for i = 1:imgnum
        
        if get_size  == 1
            img = imread([img_path,imgname{i}]);
            gt_size(i,1) = size(img,1); gt_size(i,2) = size(img,2);
        end
        
        fid_gt = fopen([gt_path,imgname{i},'.txt'],'r');
        if fid_gt == -1
            error('gt file error');
        end
        str = fgetl(fid_gt);
        elp_num = str2double(str);
        elp_data = zeros(elp_num,5);
        for j = 1:elp_num
            str = fgetl(fid_gt);
            tmp = str2num(str);
            elp_data(j,:) = [tmp(2), tmp(1), tmp(4), tmp(3), -tmp(5)/180*pi];
        end
        gt_elps{i}=elp_data;
        fclose(fid_gt);
    end
    
    if get_size == 1
        varargout{1} = gt_size;
    end
    return;
    
end



if strcmp(dataset, 'axisratiowithorientation') || strcmp(dataset, 'semimajorwithaxisratio')
    
    for i = 1:imgnum
        
        if get_size  == 1
            img = imread([img_path,imgname{i}]);
            gt_size(i,1) = size(img,1); gt_size(i,2) = size(img,2);
        end
        
        fid_gt = fopen([gt_path,imgname{i},'.txt'],'r');
        if fid_gt == -1
            error('gt file error');
        end
        str = fgetl(fid_gt);
        elp_num = str2double(str);
        elp_data = zeros(elp_num,5);
        for j = 1:elp_num
            str = fgetl(fid_gt);
            tmp = str2num(str);
            elp_data(j,:) = [tmp(2), tmp(1), tmp(4), tmp(3), -tmp(5)];
        end
        gt_elps{i}=elp_data;
        fclose(fid_gt);
    end
    
    if get_size == 1
        varargout{1} = gt_size;
    end
    return;
    
end

if strcmp(dataset, 'satellite1') || strcmp(dataset, 'satellite2') || strcmp(dataset, 'industrial') 
    A = load([gt_path(1:end-1),'\gt.mat']);
    gt_elps = A.gts;
    
    
    if get_size == 1
        varargout{1} = 0;
    end
    return;
end

error(['Current dataset does not exist: ', dataset]);




end


function [] = drawEllipses(ellipses_para,im,c,siwtch)
    hold on;
    size_im = size(im);
    if size(ellipses_para,1) ~= 5
        ellipses_para = ellipses_para';
    end

    th=0:pi/180:2*pi;
    for i=1:size(ellipses_para,2)
        Semi_major= ellipses_para(3,i);
        Semi_minor= ellipses_para(4,i);
        x0= ellipses_para(1,i);
        y0= ellipses_para(2,i);
        Phi= ellipses_para(5,i);
        x=x0+Semi_major*cos(Phi)*cos(th)-Semi_minor*sin(Phi)*sin(th);
        y=y0+Semi_minor*cos(Phi)*sin(th)+Semi_major*sin(Phi)*cos(th);
        if ~siwtch
            plot(x,y,'Color',c, 'LineWidth',1.5);
        else
            plot(y,x,'Color',c, 'LineWidth',1.5);
        end
    end
    if ~isempty(im)
        axis on; set(gca,'XTick',[],'YTick',[]);axis ij;axis equal;axis([0 size_im(2) 0 size_im(1)]);
    end
end







function detector=configEllipseDetector(id)
    detectors = { ...
        struct('path',['C:\Users\xiahaa\Documents\DTU\publications\fore-end\ellipse_detection\thirdparty\High-quality-ellipse-detection-master'], ...
        'name', 'lu', 'namepaper', 'LU', 'call', @ellDetect_lite,'options',''), ...
        struct('path',['C:\Users\xiahaa\Documents\DTU\publications\fore-end\ellipse_detection\thirdparty\High-quality-ellipse-detection-master'], ...
        'name', 'lu', 'namepaper', 'LU', 'call', @ellDetect,'options','Subpixel'), ...
        struct('path',['C:\Users\xiahaa\Documents\DTU\publications\fore-end\ellipse_detection\thirdparty\High-quality-ellipse-detection-master'], ...
        'name', 'lu', 'namepaper', 'LU', 'call', @ellDetect,'options','Ahn'), ...
        struct('path',['C:\Users\xiahaa\Documents\DTU\publications\fore-end\ellipse_detection\thirdparty\High-quality-ellipse-detection-master'], ...
        'name', 'lu', 'namepaper', 'LU', 'call', @ellDetect,'options','Subpixel, Ahn'), ...
        };
    detector = detectors{id};
end