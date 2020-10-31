//
//  clustering.hpp
//  lu
//
//  Created by xiao.hu on 2020/5/20.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#ifndef clustering_hpp
#define clustering_hpp

#include <stdio.h>
#include "datatype.h"

// not sure whether or not this is useful
void generateEllipseCandidates( PairGroupList * pairGroupList, double distance_tolerance, double * & ellipse_candidates, int * candidates_num);

int  cluster1DDatas(double * datas, int datas_num, double distance_tolerance, double * & centers, int * centers_num);

void clusterByDistance(double * points, int nPoints, int nDims, double distance_threshold, int number_control, double * & outPoints, int * nOutPoints);

void fastgenerateEllipseCandidates(PairGroupList * pairGroupList, double distance_tolerance, double * & ellipse_candidates, int * candidates_num);

#endif /* clustering_hpp */
