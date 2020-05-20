#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <limits.h>
#include <float.h>
#include <iostream>
#ifdef __APPLE__
#else
#include "mex.h"
#endif
#include "opencv2/core/core.hpp" 
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/opencv.hpp"

#include "datatype.h"
#include "utils.hpp"
#include "image.hpp"
#include "elsd.hpp"
#include "group_forming.hpp"
#include "gen_init_set.hpp"
#include "clustering.hpp"
#include "gradient.hpp"

using namespace cv;

/*------------------------------------------------------------------------------------------------*/
/**
my code,Alan Lu
输入
img  : 输入图像的一维double型数组,大小为Y*X，按照行优先存储，传入前需要拥有内存
X    : 输入图像的columns
Y    ：输入图像的rows
输出
n_out: lsd算法检测得到的线段的数量n，return的返回值是n条线段，为一维double型数组，长度为8*n，每8个为一组，存着x1,y1,x2,y2,dx,dy,width,polarity
reg_img: 输出标记区域，是一维的int型数组，大小reg_y*reg_x,在相应的像素位置标记着它属于的线段(1,2,3,...n),如果值为0表示不属于任何线段.
         假如外部是int * region_img,则只需要 &region_img,就可以得到标记区域的返回，不需要时直接NULL传入
reg_x  : 输出标记区域的columns,不需要时直接NULL传入
reg_y  : 输出标记区域的rows,不需要时直接NULL传入
*/
double * mylsd(int * n_out, double * img, int X, int Y, int ** reg_img, int * reg_x, int * reg_y)
{
   /* LSD parameters */
  double scale = 0.8;       /* Scale the image by Gaussian filter to 'scale'. */
  double sigma_scale = 0.6; /* Sigma for Gaussian filter is computed as
                                sigma = sigma_scale/scale.                    */
  double quant = 2.0;       /* Bound to the quantization error on the
                                gradient norm.                                */
  double ang_th = 22.5;     /* Gradient angle tolerance in degrees.           */
  double log_eps = 0.0;     /* Detection threshold: -log10(NFA) > log_eps     */
  double density_th = 0.7;  /* Minimal density of region point2is in rectangle. */
  int n_bins = 1024;        /* Number of bins in pseudo-ordering of gradient
                               modulus.                                       */ 

  return LineSegmentDetection( n_out, img, X, Y, scale, sigma_scale, quant,
                               ang_th, log_eps, density_th, n_bins,
                               reg_img, reg_x, reg_y );
}







//==========================================END=======================================================================
/**
输入：
prhs[0]: 输入的灰度图像，单通道，大小是imgy x imgx
prhs[1]: 边缘提取选择，1 canny; 2 sobel
prhs[2]: 检测指定的椭圆极性
输出：
plhs[0]: 候选椭圆组合(xi,yi,ai,bi,phi_i)', 5 x m
plhs[1]: 边缘图，大小是imgy x imgx，设边缘点总数为 edgepix_n. 二值化，0 或者 255
plhs[2]: 边缘点的梯度向量矩阵，大小是 2 x edgepix_n, (cos(theta_rad),sin(theta_rad))'...
plhs[3]: 线段图，大小是imgy x imgx 
*/
/*
compile：
mex generateEllipseCandidates.cpp -IF:\OpenCV\opencv2.4.9\build\include -IF:\OpenCV\opencv2.4.9\build\include\opencv -IF:\OpenCV\opencv2.4.9\build\include\opencv2 -LF:\OpenCV\opencv2.4.9\build\x64\vc11\lib -IF:\Matlab\settlein\extern\include -LF:\Matlab\settlein\extern\lib\win64\microsoft -lopencv_core249 -lopencv_highgui249 -lopencv_imgproc249 -llibmwlapack.lib
*/
//======================================MEX function==================================================================
#ifdef __APPLE__
//
//  main.c
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

int main(int argc, const char * argv[]) {
    // insert code here...
    printf("Hello, World!\n");
    
    std::string images_folder = "/Volumes/document/fore-end/ellipse_detection/data/test_data";
    std::string out_folder = images_folder;
    std::vector<std::string> names;
    names.push_back(images_folder+"/"+"sample01.jpg");
    //glob(images_folder + "Lo3my4.*", names);
    
    for (const auto& image_name : names)
    {
        std::string name = image_name.substr(image_name.find_last_of("/") + 1);
        name = name.substr(0, name.find_last_of("."));
        
        Mat3b image = imread(image_name);
        Size sz = image.size();
        
        // Convert to grayscale
        Mat1b gray;//1b is one byte
        cvtColor(image, gray, CV_BGR2GRAY);
        
        int edge_process_select = 2;
        int specified_polarity = 0;
        
        int imgy = gray.rows;
        int imgx = gray.cols;
        
        uchar *ptr = (uchar*)gray.data;
        
        double *data=(double*)malloc(imgy*imgx*sizeof(double));//将输入矩阵中的图像数据转存到一维数组中
        for(int c=0;c<imgx;c++)
        {
            for(int r=0;r<imgy;r++)
            {
                data[c+r*imgx]=ptr[c+r*imgx];
            }
        }
        
        int n;//线段数量
        //int new_n;
        std::vector<std::vector<int>> groups;
        double * coverages;
        int * reg;
        int reg_x;
        int reg_y;
        double* out=mylsd(&n, data,imgx,imgy,&reg,&reg_x,&reg_y);
        groupLSs(out,n,reg,reg_x,reg_y,&groups);//分组， done.
        free(reg); //释放内存
        calcuGroupCoverage(out,n,groups,coverages);//计算每个组的覆盖角度, done
        
        printf("The number of output arc-support line segments: %i\n",n);
        printf("The number of arc-support groups:%i\n",groups.size());
        
        image_double angles;
        if(edge_process_select == 1)
            calculateGradient2(data,imgx,imgy,&angles); //version2, sobel; version 3 canny
        else
            calculateGradient3(data,imgx,imgy,&angles); //version2, sobel; version 3 canny
        PairGroupList * pairGroupList;
        double distance_tolerance = 2;//max( 2.0, 0.005*min(angles->xsize,angles->ysize) ); // 0.005%*min(xsize,ysize)
        double * candidates; //候选椭圆
        double * candidates_out;//输出候选椭圆指针
        int  candidates_num = 0;//候选椭圆数量
        //rejectShortLines(out,n,&new_n);
        pairGroupList = getValidInitialEllipseSet(out,n,&groups,coverages,angles,distance_tolerance,specified_polarity);
//        double *candidates_out = NULL;
        if(pairGroupList != NULL)
        {
            printf("The number of initial ellipses：%i \n",pairGroupList->length);
            generateEllipseCandidates(pairGroupList, distance_tolerance, candidates, &candidates_num);
            printf("The number of ellipse candidates: %i \n",candidates_num);
            
//            plhs[0] = mxCreateDoubleMatrix(5,candidates_num,mxREAL);
//            candidates_out = (double*)mxGetPr(plhs[0]);
            candidates_out = (double *)malloc(sizeof(double)*5*candidates_num);
            if (candidates_out == NULL)
            {
                std::cout << "initialize candidates_out failed!"<<std::endl;
                assert(1);
            }
            //候选圆组合(xi,yi,ai,bi,phi_i)', 5 x candidates_num, 复制到矩阵candidates_out中
            memcpy(candidates_out,candidates,sizeof(double)*5*candidates_num);
            
            freePairGroupList(pairGroupList);
            free(candidates);
        }
        else
        {
            printf("The number of initial ellipses：%i \n",0);
//            double *candidates_out;
//            plhs[0] = mxCreateDoubleMatrix(5,1,mxREAL);
//            candidates_out = (double*)mxGetPr(plhs[0]);
//            candidates_out[0] = candidates_out[1] = candidates_out[2] = candidates_out[3] = candidates_out[4] = 0;
        }
//        uchar *edgeimg_out;
        unsigned long edge_pixels_total_num = 0;//边缘总像素
        double *gradient_vec_out;
//        plhs[1] = mxCreateNumericMatrix(imgy,imgx,mxUINT8_CLASS,mxREAL);
//        edgeimg_out = (uchar*)mxGetData(plhs[1]);
        //将边缘图复制到矩阵edgeimg_out中
        //将梯度向量存到矩阵gradient_vec_out中
        cv::Mat edgeimg_out(imgy,imgx,CV_8UC1);
        unsigned long addr,g_cnt = 0;
        for ( int c = 0; c < imgx; c++ )
            for ( int r = 0; r < imgy; r++)
            {
                addr = r*imgx+c;
                if(angles->data[addr] == NOTDEF)
                    edgeimg_out.data[addr] = 0;
                else
                {
                    edgeimg_out.data[addr] = 255;//为边缘点，赋值为白色
                    //------------------------------------------------
                    edge_pixels_total_num++;
                }
            }
        cv::imshow("edge detection res", edgeimg_out);
        printf("edge pixel number: %i\n",edge_pixels_total_num);
        //申请edge_pixels_total_num x 2 来保存每一个边缘点的梯度向量，以列为优先，符合matlab的习惯
//        plhs[2] = mxCreateDoubleMatrix(2,edge_pixels_total_num,mxREAL);
//        gradient_vec_out = (double*)mxGetPr(plhs[2]);
//        for ( int c = 0; c < imgx; c++ )
//            for ( int r = 0; r < imgy; r++)
//            {
//                addr = r*imgx+c;
//                if(angles->data[addr] != NOTDEF)
//                {
//                    gradient_vec_out[g_cnt++] = cos(angles->data[addr]);
//                    gradient_vec_out[g_cnt++] = sin(angles->data[addr]);
//                }
//            }
        //---------------------------------------------------------------------
        //输出线段检测的图像
//        if(nlhs == 4)
        {
            Mat ls_mat = Mat::zeros(imgy,imgx,CV_8UC3);
            image.copyTo(ls_mat);
            for ( int i = 0; i<n ; i++)//draw lines
            {
                Point2d p1(out[8*i],out[8*i+1]),p2(out[8*i+2],out[8*i+3]);
                line(ls_mat,p1,p2,Scalar(255,0,0));
            }
            if(candidates_num > 0)//draw ellipses
            {
                for ( int i = 0; i<candidates_num; i++)
                    ellipse(ls_mat,cv::Point((int)candidates_out[i*5],(int)candidates_out[i*5+1]),cv::Size(candidates_out[i*5+2],candidates_out[i*5+3]),candidates_out[i*5+4]*180/M_PI,0,360,(Scalar(0,255,0)),1);
            }
            
            cv::imshow("seg detection res", ls_mat);
            cv::waitKey(0);
//            plhs[3] = mxCreateDoubleMatrix(imgy,imgx,mxREAL);
//            double * ls_img_out = (double*)mxGetPr(plhs[3]);
            //memcpy(ls_out_mat,ls_mat.data ,sizeof(unsigned char)*M*N);
//            for (int i = 0; i<imgx; i++)
//                for (int j = 0; j<imgy;j++)
//                    ls_img_out[i*imgy+j]=ls_mat.data[j*imgx+i];
        }
        //---------------------------------------------------------------------
        //这里的free是释放程序中用于产生候选圆所用到的一系列内存
        free(data);
        free(coverages);
        free(out);
        free_image_double(angles);
        
        //        Mat3b resultImage = image.clone();
        //        yaed.DrawDetectedEllipses(resultImage, ellsYaed);
        //
        //        imwrite(out_folder + name + ".png", resultImage);
        //
        //        imshow("Yaed", resultImage);
        //        waitKey();
        //
        //
        //        int yghds = 0;
    }
    
    return 0;
}



#else
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs!=3)
        mexErrMsgIdAndTxt( "MATLAB:revord:invalidNumInputs","One input required.");
    else if(nlhs > 4)
        mexErrMsgIdAndTxt( "MATLAB:revord:maxlhs","Too many output arguments.");
    uchar * inputimg = (uchar*)mxGetData(prhs[0]);
    int imgy,imgx;
    int edge_process_select = (int)mxGetScalar(prhs[1]);//边缘提取选择，1 canny; 2 sobel
    int specified_polarity  = (int)mxGetScalar(prhs[2]);//1,指定检测的椭圆极性要为正; -1指定极性为负; 0表示两种极性椭圆都检测
    imgy = (int)mxGetM(prhs[0]);
    imgx = (int)mxGetN(prhs[0]);
    double *data=(double*)malloc(imgy*imgx*sizeof(double));//将输入矩阵中的图像数据转存到一维数组中
    for(int c=0;c<imgx;c++)
    {
        for(int r=0;r<imgy;r++)
        {
            data[c+r*imgx]=inputimg[r+c*imgy];
        }
    }
    int n;//线段数量
    //int new_n;
    std::vector<std::vector<int>> groups;
    double * coverages;
    int * reg;
    int reg_x;
    int reg_y;
    double* out=mylsd(&n, data,imgx,imgy,&reg,&reg_x,&reg_y);
    groupLSs(out,n,reg,reg_x,reg_y,&groups);//分组， done.
    free(reg); //释放内存
    calcuGroupCoverage(out,n,groups,coverages);//计算每个组的覆盖角度, done

    printf("The number of output arc-support line segments: %i\n",n);
    printf("The number of arc-support groups:%i\n",groups.size());
  /*int groups_t = 0;
  for (int i = 0; i<groups.size(); i++)
  { 
    groups_t+= groups[i].size();
  }
  printf("Groups' total ls num:%i\n",groups_t);*/

    image_double angles;
    if(edge_process_select == 1)
        calculateGradient2(data,imgx,imgy,&angles); //version2, sobel; version 3 canny
    else
        calculateGradient3(data,imgx,imgy,&angles); //version2, sobel; version 3 canny
    PairGroupList * pairGroupList;
    double distance_tolerance = 2;//max( 2.0, 0.005*min(angles->xsize,angles->ysize) ); // 0.005%*min(xsize,ysize)
    double * candidates; //候选椭圆
    double * candidates_out;//输出候选椭圆指针
    int  candidates_num = 0;//候选椭圆数量
    //rejectShortLines(out,n,&new_n);
    pairGroupList = getValidInitialEllipseSet(out,n,&groups,coverages,angles,distance_tolerance,specified_polarity);
    if(pairGroupList != NULL)
    {
        printf("The number of initial ellipses：%i \n",pairGroupList->length);
        generateEllipseCandidates(pairGroupList, distance_tolerance, candidates, &candidates_num);
        printf("The number of ellipse candidates: %i \n",candidates_num);
    
        plhs[0] = mxCreateDoubleMatrix(5,candidates_num,mxREAL);
        candidates_out = (double*)mxGetPr(plhs[0]);
        //候选圆组合(xi,yi,ai,bi,phi_i)', 5 x candidates_num, 复制到矩阵candidates_out中
        memcpy(candidates_out,candidates,sizeof(double)*5*candidates_num);

        freePairGroupList(pairGroupList);
        free(candidates);
   }
   else
   {
        printf("The number of initial ellipses：%i \n",0);
        double *candidates_out;
        plhs[0] = mxCreateDoubleMatrix(5,1,mxREAL);
        candidates_out = (double*)mxGetPr(plhs[0]);
        candidates_out[0] = candidates_out[1] = candidates_out[2] = candidates_out[3] = candidates_out[4] = 0;
   }
   uchar *edgeimg_out;
   unsigned long edge_pixels_total_num = 0;//边缘总像素
   double *gradient_vec_out;
   plhs[1] = mxCreateNumericMatrix(imgy,imgx,mxUINT8_CLASS,mxREAL);
   edgeimg_out = (uchar*)mxGetData(plhs[1]);
   //将边缘图复制到矩阵edgeimg_out中
   //将梯度向量存到矩阵gradient_vec_out中
   unsigned long addr,g_cnt = 0;
   for ( int c = 0; c < imgx; c++ )
       for ( int r = 0; r < imgy; r++)
       {
           addr = r*imgx+c;
           if(angles->data[addr] == NOTDEF)
               edgeimg_out[c*imgy+r] = 0;
           else
           {
               edgeimg_out[c*imgy+r] = 255;//为边缘点，赋值为白色
               //------------------------------------------------
               edge_pixels_total_num++;
           }
       }
    printf("edge pixel number: %i\n",edge_pixels_total_num);
    //申请edge_pixels_total_num x 2 来保存每一个边缘点的梯度向量，以列为优先，符合matlab的习惯
    plhs[2] = mxCreateDoubleMatrix(2,edge_pixels_total_num,mxREAL);
    gradient_vec_out = (double*)mxGetPr(plhs[2]);
    for ( int c = 0; c < imgx; c++ )
        for ( int r = 0; r < imgy; r++)
        {
            addr = r*imgx+c;
            if(angles->data[addr] != NOTDEF)
            {
                gradient_vec_out[g_cnt++] = cos(angles->data[addr]);
                gradient_vec_out[g_cnt++] = sin(angles->data[addr]);
            }
        }
    //---------------------------------------------------------------------
    //输出线段检测的图像
    if(nlhs == 4)
    {
        Mat ls_mat = Mat::zeros(imgy,imgx,CV_8UC1);
        for ( int i = 0; i<n ; i++)//draw lines
        {
            Point2d p1(out[8*i],out[8*i+1]),p2(out[8*i+2],out[8*i+3]);
            line(ls_mat,p1,p2,Scalar(255,0,0));
        }
        if(candidates_num > 0)//draw ellipses
        {
            for ( int i = 0; i<candidates_num; i++)
                ellipse(ls_mat,cv::Point((int)candidates_out[i*5],(int)candidates_out[i*5+1]),cv::Size(candidates_out[i*5+2],candidates_out[i*5+3]),candidates_out[i*5+4]*180/M_PI,0,360,(Scalar(255,0,0)),1);
        }
        plhs[3] = mxCreateDoubleMatrix(imgy,imgx,mxREAL);
        double * ls_img_out = (double*)mxGetPr(plhs[3]);
        //memcpy(ls_out_mat,ls_mat.data ,sizeof(unsigned char)*M*N);
        for (int i = 0; i<imgx; i++)
            for (int j = 0; j<imgy;j++)
                ls_img_out[i*imgy+j]=ls_mat.data[j*imgx+i];
    }
    //---------------------------------------------------------------------
    //这里的free是释放程序中用于产生候选圆所用到的一系列内存
    free(data);
    free(coverages);
    free(out);
    free_image_double(angles);
}
#endif
