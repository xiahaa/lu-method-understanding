//
//  canny.cpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#include "canny.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <limits.h>
#include <float.h>
#include <iostream>
#include "opencv2/core/core.hpp"
//#include "opencv2/features2d/features2d.hpp"
//#include "opencv2/highgui/highgui.hpp"
#include "opencv2/opencv.hpp"

using namespace cv;

//=============================================================================
//需要包含如下头文件
//#include <opencv2\opencv.hpp>
//using namespace cv;
void cvCanny3(  const void* srcarr, void* dstarr,
              void* dxarr, void* dyarr,
              int aperture_size )
{
    //cv::Ptr<CvMat> dx, dy;
    cv::AutoBuffer<char> buffer;
    std::vector<uchar*> stack;
    uchar **stack_top = 0, **stack_bottom = 0;
    
    CvMat srcstub, *src = cvGetMat( srcarr, &srcstub );
    CvMat dststub, *dst = cvGetMat( dstarr, &dststub );
    
    CvMat dxstub, *dx = cvGetMat( dxarr, &dxstub );
    CvMat dystub, *dy = cvGetMat( dyarr, &dystub );
    
    CvSize size;
    int flags = aperture_size;
    int low, high;
    int* mag_buf[3];
    uchar* map;
    ptrdiff_t mapstep;
    int maxsize;
    int i, j;
    CvMat mag_row;
    
    if( CV_MAT_TYPE( src->type ) != CV_8UC1 ||
       CV_MAT_TYPE( dst->type ) != CV_8UC1 ||
       CV_MAT_TYPE( dx->type  ) != CV_16SC1 ||
       CV_MAT_TYPE( dy->type  ) != CV_16SC1 )
        CV_Error( CV_StsUnsupportedFormat, "" );
    
    if( !CV_ARE_SIZES_EQ( src, dst ))
        CV_Error( CV_StsUnmatchedSizes, "" );
    
    aperture_size &= INT_MAX;
    if( (aperture_size & 1) == 0 || aperture_size < 3 || aperture_size > 7 )
        CV_Error( CV_StsBadFlag, "" );
    
    
    size.width = src->cols;
    size.height = src->rows;
    
    //aperture_size = -1; //SCHARR
    cvSobel( src, dx, 1, 0, aperture_size );
    cvSobel( src, dy, 0, 1, aperture_size );
    
    Mat1f magGrad(size.height, size.width, 0.f);
    float maxGrad(0);
    float val(0);
    for(i=0; i<size.height; ++i)
    {
        float* _pmag = magGrad.ptr<float>(i);
        const short* _dx = (short*)(dx->data.ptr + dx->step*i);
        const short* _dy = (short*)(dy->data.ptr + dy->step*i);
        for(j=0; j<size.width; ++j)
        {
            val = float(abs(_dx[j]) + abs(_dy[j]));
            _pmag[j] = val;
            maxGrad = (val > maxGrad) ? val : maxGrad;
        }
    }
    
    //% Normalize for threshold selection
    //normalize(magGrad, magGrad, 0.0, 1.0, NORM_MINMAX);
    
    //% Determine Hysteresis Thresholds
    
    //set magic numbers
    const int NUM_BINS = 64;
    const double percent_of_pixels_not_edges = 0.9;
    const double threshold_ratio = 0.3;
    
    //compute histogram
    int bin_size = cvFloor(maxGrad / float(NUM_BINS) + 0.5f) + 1;
    if (bin_size < 1) bin_size = 1;
    int bins[NUM_BINS] = { 0 };
    for (i=0; i<size.height; ++i)
    {
        float *_pmag = magGrad.ptr<float>(i);
        for(j=0; j<size.width; ++j)
        {
            int hgf = int(_pmag[j]);
            bins[int(_pmag[j]) / bin_size]++;
        }
    }
    
    
    
    
    //% Select the thresholds
    float total(0.f);
    float target = float(size.height * size.width * percent_of_pixels_not_edges);
    int low_thresh, high_thresh(0);
    
    while(total < target)
    {
        total+= bins[high_thresh];
        high_thresh++;
    }
    high_thresh *= bin_size;
    low_thresh = cvFloor(threshold_ratio * float(high_thresh));
    
    if( flags & CV_CANNY_L2_GRADIENT )
    {
        Cv32suf ul, uh;
        ul.f = (float)low_thresh;
        uh.f = (float)high_thresh;
        
        low = ul.i;
        high = uh.i;
    }
    else
    {
        low = cvFloor( low_thresh );
        high = cvFloor( high_thresh );
    }
    
    
    buffer.allocate( (size.width+2)*(size.height+2) + (size.width+2)*3*sizeof(int) );
    mag_buf[0] = (int*)(char*)buffer;
    mag_buf[1] = mag_buf[0] + size.width + 2;
    mag_buf[2] = mag_buf[1] + size.width + 2;
    map = (uchar*)(mag_buf[2] + size.width + 2);
    mapstep = size.width + 2;
    
    maxsize = MAX( 1 << 10, size.width*size.height/10 );
    stack.resize( maxsize );
    stack_top = stack_bottom = &stack[0];
    
    memset( mag_buf[0], 0, (size.width+2)*sizeof(int) );
    memset( map, 1, mapstep );
    memset( map + mapstep*(size.height + 1), 1, mapstep );
    
    /* sector numbers
     (Top-Left Origin)
     
     1   2   3
     *  *  *
     * * *
     0*******0
     * * *
     *  *  *
     3   2   1
     */
    
#define CANNY_PUSH(d)    *(d) = (uchar)2, *stack_top++ = (d)
#define CANNY_POP(d)     (d) = *--stack_top
    
    mag_row = cvMat( 1, size.width, CV_32F );
    
    // calculate magnitude and angle of gradient, perform non-maxima supression.
    // fill the map with one of the following values:
    //   0 - the pixel might belong to an edge
    //   1 - the pixel can not belong to an edge
    //   2 - the pixel does belong to an edge
    for( i = 0; i <= size.height; i++ )
    {
        int* _mag = mag_buf[(i > 0) + 1] + 1;
        float* _magf = (float*)_mag;
        const short* _dx = (short*)(dx->data.ptr + dx->step*i);
        const short* _dy = (short*)(dy->data.ptr + dy->step*i);
        uchar* _map;
        int x, y;
        ptrdiff_t magstep1, magstep2;
        int prev_flag = 0;
        
        if( i < size.height )
        {
            _mag[-1] = _mag[size.width] = 0;
            
            if( !(flags & CV_CANNY_L2_GRADIENT) )
                for( j = 0; j < size.width; j++ )
                    _mag[j] = abs(_dx[j]) + abs(_dy[j]);
            
            else
            {
                for( j = 0; j < size.width; j++ )
                {
                    x = _dx[j]; y = _dy[j];
                    _magf[j] = (float)std::sqrt((double)x*x + (double)y*y);
                }
            }
        }
        else
            memset( _mag-1, 0, (size.width + 2)*sizeof(int) );
        
        // at the very beginning we do not have a complete ring
        // buffer of 3 magnitude rows for non-maxima suppression
        if( i == 0 )
            continue;
        
        _map = map + mapstep*i + 1;
        _map[-1] = _map[size.width] = 1;
        
        _mag = mag_buf[1] + 1; // take the central row
        _dx = (short*)(dx->data.ptr + dx->step*(i-1));
        _dy = (short*)(dy->data.ptr + dy->step*(i-1));
        
        magstep1 = mag_buf[2] - mag_buf[1];
        magstep2 = mag_buf[0] - mag_buf[1];
        
        if( (stack_top - stack_bottom) + size.width > maxsize )
        {
            int sz = (int)(stack_top - stack_bottom);
            maxsize = MAX( maxsize * 3/2, maxsize + 8 );
            stack.resize(maxsize);
            stack_bottom = &stack[0];
            stack_top = stack_bottom + sz;
        }
        
        for( j = 0; j < size.width; j++ )
        {
#define CANNY_SHIFT 15
#define TG22  (int)(0.4142135623730950488016887242097*(1<<CANNY_SHIFT) + 0.5)
            
            x = _dx[j];
            y = _dy[j];
            int s = x ^ y;
            int m = _mag[j];
            
            x = abs(x);
            y = abs(y);
            if( m > low )
            {
                int tg22x = x * TG22;
                int tg67x = tg22x + ((x + x) << CANNY_SHIFT);
                
                y <<= CANNY_SHIFT;
                
                if( y < tg22x )
                {
                    if( m > _mag[j-1] && m >= _mag[j+1] )
                    {
                        if( m > high && !prev_flag && _map[j-mapstep] != 2 )
                        {
                            CANNY_PUSH( _map + j );
                            prev_flag = 1;
                        }
                        else
                            _map[j] = (uchar)0;
                        continue;
                    }
                }
                else if( y > tg67x )
                {
                    if( m > _mag[j+magstep2] && m >= _mag[j+magstep1] )
                    {
                        if( m > high && !prev_flag && _map[j-mapstep] != 2 )
                        {
                            CANNY_PUSH( _map + j );
                            prev_flag = 1;
                        }
                        else
                            _map[j] = (uchar)0;
                        continue;
                    }
                }
                else
                {
                    s = s < 0 ? -1 : 1;
                    if( m > _mag[j+magstep2-s] && m > _mag[j+magstep1+s] )
                    {
                        if( m > high && !prev_flag && _map[j-mapstep] != 2 )
                        {
                            CANNY_PUSH( _map + j );
                            prev_flag = 1;
                        }
                        else
                            _map[j] = (uchar)0;
                        continue;
                    }
                }
            }
            prev_flag = 0;
            _map[j] = (uchar)1;
        }
        
        // scroll the ring buffer
        _mag = mag_buf[0];
        mag_buf[0] = mag_buf[1];
        mag_buf[1] = mag_buf[2];
        mag_buf[2] = _mag;
    }
    
    // now track the edges (hysteresis thresholding)
    while( stack_top > stack_bottom )
    {
        uchar* m;
        if( (stack_top - stack_bottom) + 8 > maxsize )
        {
            int sz = (int)(stack_top - stack_bottom);
            maxsize = MAX( maxsize * 3/2, maxsize + 8 );
            stack.resize(maxsize);
            stack_bottom = &stack[0];
            stack_top = stack_bottom + sz;
        }
        
        CANNY_POP(m);
        
        if( !m[-1] )
            CANNY_PUSH( m - 1 );
        if( !m[1] )
            CANNY_PUSH( m + 1 );
        if( !m[-mapstep-1] )
            CANNY_PUSH( m - mapstep - 1 );
        if( !m[-mapstep] )
            CANNY_PUSH( m - mapstep );
        if( !m[-mapstep+1] )
            CANNY_PUSH( m - mapstep + 1 );
        if( !m[mapstep-1] )
            CANNY_PUSH( m + mapstep - 1 );
        if( !m[mapstep] )
            CANNY_PUSH( m + mapstep );
        if( !m[mapstep+1] )
            CANNY_PUSH( m + mapstep + 1 );
    }
    
    // the final pass, form the final image
    for( i = 0; i < size.height; i++ )
    {
        const uchar* _map = map + mapstep*(i+1) + 1;
        uchar* _dst = dst->data.ptr + dst->step*i;
        
        for( j = 0; j < size.width; j++ )
        {
            _dst[j] = (uchar)-(_map[j] >> 1);
        }
    }
};

void Canny3(  InputArray image, OutputArray _edges,
            OutputArray _sobel_x, OutputArray _sobel_y,
            int apertureSize, bool L2gradient )
{
    Mat src = image.getMat();
    _edges.create(src.size(), CV_8U);
    _sobel_x.create(src.size(), CV_16S);
    _sobel_y.create(src.size(), CV_16S);
    
    
    CvMat c_src = src;
    CvMat c_dst = _edges.getMat();
    CvMat c_dx = _sobel_x.getMat();
    CvMat c_dy = _sobel_y.getMat();
    
    
    cvCanny3( &c_src, &c_dst,
             &c_dx, &c_dy,
             apertureSize + (L2gradient ? CV_CANNY_L2_GRADIENT : 0));
};
