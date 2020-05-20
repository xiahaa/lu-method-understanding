//
//  gradient.cpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
// done

#include "gradient.hpp"
#include <math.h>
#include <stdlib.h>
#include "opencv2/core/core.hpp"
#include "opencv2/opencv.hpp"
#include "canny.hpp"

using namespace cv;

void calculateGradient2( double * img_in, unsigned int imgx, unsigned int imgy, image_double * angles)
{
    if(img_in == NULL || imgx == 0 || imgy == 0)
        error("calculateGradient error!");
    image_double mod = new_image_double(imgx,imgy);
    (*angles) = new_image_double(imgx,imgy);
    unsigned int x,y,adr;
    double com1,com2;
    double gx,gy;
    double norm,norm_square;
    double threshold;
    double sum = 0;
    double value;
    //double max_grad = 0.0;
    //边界初始为NOTDEF
    for ( x = 0; x<imgx; x++)
    {
        (*angles)->data[x]=NOTDEF;
        (*angles)->data[(imgy-1)*imgx+x]=NOTDEF;
        (mod)->data[x]=NOTDEF;
        (mod)->data[(imgy-1)*imgx+x]=NOTDEF;
    }
    for ( y = 0; y<imgy; y++)
    {
        (*angles)->data[y*imgx] = NOTDEF;
        (*angles)->data[y*imgx+imgx-1] = NOTDEF;
        (mod)->data[y*imgx] = NOTDEF;
        (mod)->data[y*imgx+imgx-1] = NOTDEF;
    }
    /* compute gradient on the remaining pixels */
    for(x=1;x<imgx-1;x++)
        for(y=1;y<imgy-1;y++)
        {
            adr = y*imgx+x;
            /*
             Norm 2 computation using 2x2 pixel window:
             A B C
             D E F
             G H I
             and
             com1 = C-G,  com2 = I-A.
             Then
             gx = C+2F+I - (A+2D+G)=com1+com2+2(F-D)   horizontal difference
             gy = G+2H+I - (A+2B+C)=-com1+com2+2(H-B)   vertical difference
             com1 and com2 are just to avoid 2 additions.
             */
            com1 = img_in[adr-imgx+1] - img_in[adr+imgx-1];
            com2 = img_in[adr+imgx+1] - img_in[adr-imgx-1];
            
            gx = (com1+com2+2*(img_in[adr+1] - img_in[adr-1]))/(8.0*255); /* gradient x component */
            gy = (-com1+com2+2*(img_in[adr+imgx] - img_in[adr-imgx]))/(8.0*255); /* gradient y component */
            norm_square = gx*gx+gy*gy;
            sum+=norm_square;
            
            norm = sqrt( norm_square); /* gradient norm */
            
            (mod)->data[adr] = norm; /* store gradient norm */
            /* gradient angle computation */
            (*angles)->data[adr] = atan2(gy,gx);
        }
    threshold = 2*sqrt(sum/(imgx*imgy));//自动阈值
    //non maximum suppression
    for(x=1;x<imgx-1;x++)
        for(y=1;y<imgy-1;y++)
        {
            adr = y*imgx+x;
            value = (*angles)->data[adr];
            if((mod)->data[adr] < threshold )
            {
                (*angles)->data[adr] = NOTDEF;
                continue;
            }
            if( (value > -M_1_8_PI && value<=M_1_8_PI) || (value <= -M_7_8_PI ) || (value > M_7_8_PI))
            {
                if((mod)->data[adr] <= (mod)->data[adr+1] || (mod)->data[adr] <= (mod)->data[adr-1])
                    (*angles)->data[adr] = NOTDEF;
            }
            else if( (value> M_1_8_PI && value<= M_3_8_PI) || (value> -M_7_8_PI && value<= -M_5_8_PI) )
            {
                if((mod)->data[adr] <= (mod)->data[adr-imgx-1] || (mod)->data[adr] <= (mod)->data[adr+imgx+1])
                    (*angles)->data[adr] = NOTDEF;
            }
            else if((value> M_3_8_PI && value<= M_5_8_PI) || (value> -M_5_8_PI && value<= -M_3_8_PI))
            {
                if((mod)->data[adr] <= (mod)->data[adr-imgx] || (mod)->data[adr] <= (mod)->data[adr+imgx])
                    (*angles)->data[adr] = NOTDEF;
            }
            else
            {
                if((mod)->data[adr] <= (mod)->data[adr-imgx+1] || (mod)->data[adr] <= (mod)->data[adr+imgx-1])
                    (*angles)->data[adr] = NOTDEF;
            }
        }
    //也标记到mod图上面
    //for(x=1;x<imgx-1;x++)
    //  for(y=1;y<imgy-1;y++)
    //  {
    //    if((*angles)->data[y*imgx+x] == NOTDEF)
    //      (mod)->data[y*imgx+x] = NOTDEF;
    //  }
    free_image_double(mod);
}



//canny
void calculateGradient3( double * img_in, unsigned int imgx, unsigned int imgy, image_double * angles)
{
    Mat1b edge;
    Mat1s DX,DY;
    Mat gray = cv::Mat::zeros((int)(imgy),(int)(imgx),CV_8UC1);
    unsigned int x,y,addr;
    (*angles) = new_image_double(imgx,imgy);
    //copy to gray image
    for ( y = 0; y<imgy; y++)
        for ( x = 0; x<imgx; x++)
        {
            addr = y*imgx+x;
            gray.data[addr] = (uchar)(img_in[addr]);
        }
    //canny, return edge points as well as dx dy for computing gradient orientation
    Canny3(gray,edge,DX,DY,3,false);
    for ( y = 0; y<imgy; y++)
    {
        short* _dx = DX.ptr<short>(y);
        short* _dy = DY.ptr<short>(y);
        uchar* _e = edge.ptr<uchar>(y);
        for ( x = 0; x<imgx; x++)
        {
            if(_e[x] > 0)//0 or 255
            {
                (*angles)->data[y*imgx+x]  = atan2((double)_dy[x],(double)_dx[x]);//calculate gradient
            }
            else
                (*angles)->data[y*imgx+x] = NOTDEF;
        }
    }
    edge.release();
    DX.release();
    DY.release();
    gray.release();
}
