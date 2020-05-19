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
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/opencv.hpp"

using namespace cv;

void Canny3(  InputArray image, OutputArray _edges,
            OutputArray _sobel_x, OutputArray _sobel_y,
            int apertureSize, bool L2gradient );

#endif /* canny_hpp */
