//
//  gen_init_set.hpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#ifndef gen_init_set_hpp
#define gen_init_set_hpp

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "datatype.h"
// contribution from lu
PairGroupList * getValidInitialEllipseSet( double * lines, int line_num, std::vector<std::vector<int>> * groups, double * coverages, 
	image_double angles, double distance_tolerance, int specified_polarity);

#endif /* gen_init_set_hpp */
