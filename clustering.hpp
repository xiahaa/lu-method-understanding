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

void generateEllipseCandidates( PairGroupList * pairGroupList, double distance_tolerance, double * & ellipse_candidates, int * candidates_num);

#endif /* clustering_hpp */
