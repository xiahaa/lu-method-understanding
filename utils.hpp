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
#include <math.h>
#include <float.h>
#include "datatype.h"

inline double min(double v1,double v2)
{
    return (v1<v2?v1:v2);
}
inline double max(double v1,double v2)
{
    return (v1>v2?v1:v2);
}
int double_equal(double a, double b);
double angle_diff(double a, double b);
double angle_diff_signed(double a, double b);
void error(const char * msg);
double dist(double x1, double y1, double x2, double y2);
double dotProduct(point2d vec1, point2d vec2);

int isaligned( int x, int y, image_double angles, double theta,
              double prec );

inline PairGroupList * pairGroupListInit( int length)
{
    if(length <= 0)
        error("paired groups length less equal than 0");
    PairGroupList * pairGroupList = (PairGroupList*)malloc(sizeof(PairGroupList));
    pairGroupList->length = length;
    pairGroupList->pairGroup = (PairGroup*)malloc(sizeof(PairGroup)*length);
    if(pairGroupList->pairGroup == NULL)
        error("pairGroupListInit,not enough memory");
    return pairGroupList;
}

inline void freePairGroupList( PairGroupList * list)
{
    if(list == NULL || list->pairGroup == NULL)
        error("freePairGroupList,invalidate free");
    free(list->pairGroup);
    free(list);
    list->pairGroup = NULL;
    list = NULL;
}

#endif /* utils_hpp */
