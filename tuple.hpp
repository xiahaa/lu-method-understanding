//
//  tuple.hpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#ifndef tuple_hpp
#define tuple_hpp

#include <stdio.h>
#include "datatype.h"
#include "utils.hpp"

void free_ntuple_list(ntuple_list in);
ntuple_list new_ntuple_list(int dim);
void enlarge_ntuple_list(ntuple_list n_tuple);
void add_8tuple( ntuple_list out, double v1, double v2, double v3,
                double v4, double v5, double v6, double v7, int v8);

#endif /* tuple_hpp */
