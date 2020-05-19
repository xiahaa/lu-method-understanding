//
//  elsd.hpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#ifndef elsd_hpp
#define elsd_hpp

#include <stdio.h>
#include <stdlib.h>
#include "image.hpp"
#include "datatype.h"

double * LineSegmentDetection( int * n_out,
                              double * img, int X, int Y,
                              double scale, double sigma_scale, double quant,
                              double ang_th, double log_eps, double density_th,
                              int n_bins,
                              int ** reg_img, int * reg_x, int * reg_y );

#endif /* elsd_hpp */
