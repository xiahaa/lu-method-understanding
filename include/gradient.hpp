//
//  gradient.hpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#ifndef gradient_hpp
#define gradient_hpp

#include "datatype.h"
#include <stdio.h>
#include "image.hpp"

// this one compute sobel gradient
void calculateGradient2( double * img_in, unsigned int imgx, unsigned int imgy, image_double * angles);

// this one compute canny
void calculateGradient3( double * img_in, unsigned int imgx, unsigned int imgy, image_double * angles);

#endif /* gradient_hpp */
