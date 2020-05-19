//
//  utils.hpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#ifndef utils_hpp
#define utils_hpp

#include <stdio.h>
#include <stdlib.h>

inline double min(double v1,double v2);
inline double max(double v1,double v2);
int double_equal(double a, double b);
double angle_diff(double a, double b);
double angle_diff_signed(double a, double b);
void error(const char * msg);
double dist(double x1, double y1, double x2, double y2);
double dotProduct(point2d vec1, point2d vec2);

int isaligned( int x, int y, image_double angles, double theta,
              double prec );

#endif /* utils_hpp */
