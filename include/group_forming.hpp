//
//  group_forming.hpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#ifndef group_forming_hpp
#define group_forming_hpp

#include <stdio.h>
#include <vector>

// group ls
void groupLSs(double *lines, int line_num, int * region, int imgx, int imgy, std::vector<std::vector<int>> * groups);

// compute group coverage angle
void calcuGroupCoverage(double * lines, int line_num, std::vector<std::vector<int>> groups, double * &coverages);

#endif /* group_forming_hpp */
