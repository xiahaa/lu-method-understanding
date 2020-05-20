//
//  canny.hpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#ifndef canny_hpp
#define canny_hpp

#include "opencv2/core/core.hpp"
//#include "opencv2/features2d/features2d.hpp"
//#include "opencv2/highgui/highgui.hpp"
#include "opencv2/opencv.hpp"


void Canny3(  cv::InputArray image, cv::OutputArray _edges,
            cv::OutputArray _sobel_x, cv::OutputArray _sobel_y,
            int apertureSize, bool L2gradient );

#endif /* canny_hpp */
