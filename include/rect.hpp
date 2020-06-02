//
//  rect.hpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#ifndef rect_hpp
#define rect_hpp

#include <stdio.h>
#include <stdlib.h>
#include "datatype.h"
#include "utils.hpp"

// copied from elsd
void rect_copy(struct rect * in, struct rect * out);
void ri_del(rect_iter * iter);
int ri_end(rect_iter * i);
void ri_inc(rect_iter * i);
rect_iter * ri_ini(struct rect * r);


#endif /* rect_hpp */
