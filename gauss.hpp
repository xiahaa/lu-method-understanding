//
//  gauss.hpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#ifndef gauss_hpp
#define gauss_hpp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tuple.hpp"
#include "datatype.h"
#include "utils.hpp"
#include "image.hpp"

// this does gaussian smoothing and downsampling if scale if not 1
image_double gaussian_sampler( image_double in, double scale,
                              double sigma_scale );

#endif /* gauss_hpp */
