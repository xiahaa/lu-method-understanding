//
//  ellipse_geometry.hpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#ifndef ellipse_geometry_hpp
#define ellipse_geometry_hpp

#include <stdio.h>
#include "datatype.h"

int fitEllipse(point2d* dataxy, int datanum, double* ellipara);
void addFitMatrix(double * S1, double * S2, double * S_out);
int fitEllipse2(double * S, double* ellicoeff);
int ellipse2Param(double *p,double param[]);
void calcuFitMatrix(point2d* dataxy, int datanum, double * S);

#endif /* ellipse_geometry_hpp */
