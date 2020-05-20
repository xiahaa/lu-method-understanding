//
//  image.hpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#ifndef image_hpp
#define image_hpp

#include <stdio.h>
#include "datatype.h"
#include <stdlib.h>
#include "utils.hpp"

void free_image_double(image_double i);
image_double new_image_double(int xsize, int ysize);
image_double new_image_double_ptr( int xsize,
                                  int ysize, double * data );

void free_image_char(image_char i);
image_char new_image_char(unsigned int xsize, unsigned int ysize);
image_char new_image_char_ini( unsigned int xsize, unsigned int ysize,
                              unsigned char fill_value );

image_int new_image_int(unsigned int xsize, unsigned int ysize);
image_int new_image_int_ini( unsigned int xsize, unsigned int ysize,
                                   int fill_value );

#endif /* image_hpp */
