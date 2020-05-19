//
//  image.cpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#include "image.hpp"


/*----------------------------------------------------------------------------*/
/** Free memory used in image_double 'i'.
 */
void free_image_double(image_double i)
{
    if( i == NULL || i->data == NULL )
        error("free_image_double: invalid input image.");
    free( (void *) i->data );
    free( (void *) i );
}

/*----------------------------------------------------------------------------*/
/** Create a new image_double of size 'xsize' times 'ysize'.
 */
image_double new_image_double(int xsize, int ysize)
{
    image_double image;
    
    /* check parameters */
    if( xsize == 0 || ysize == 0 ) error("new_image_double: invalid image size.");
    
    /* get memory */
    image = (image_double) malloc( sizeof(struct image_double_s) );
    if( image == NULL ) error("not enough memory.");
    image->data = (double *) calloc( (size_t) (xsize*ysize), sizeof(double) );
    if( image->data == NULL ) error("not enough memory.");
    
    /* set image size */
    image->xsize = xsize;
    image->ysize = ysize;
    
    return image;
}

/*----------------------------------------------------------------------------*/
/** Create a new image_double of size 'xsize' times 'ysize'
 with the data pointed by 'data'.
 */
image_double new_image_double_ptr( int xsize,
                                  int ysize, double * data )
{
    image_double image;
    
    /* check parameters */
    if( xsize == 0 || ysize == 0 )
        error("new_image_double_ptr: invalid image size.");
    if( data == NULL ) error("new_image_double_ptr: NULL data pointer.");
    
    /* get memory */
    image = (image_double) malloc( sizeof(struct image_double_s) );
    if( image == NULL ) error("not enough memory.");
    
    /* set image */
    image->xsize = xsize;
    image->ysize = ysize;
    image->data = data;
    
    return image;
}

/*----------------------------------------------------------------------------*/
/** Free memory used in image_char 'i'.
 */
void free_image_char(image_char i)
{
    if( i == NULL || i->data == NULL )
        error("free_image_char: invalid input image.");
    free( (void *) i->data );
    free( (void *) i );
}

/*----------------------------------------------------------------------------*/
/** Create a new image_char of size 'xsize' times 'ysize'.
 */
image_char new_image_char(unsigned int xsize, unsigned int ysize)
{
    image_char image;
    
    /* check parameters */
    if( xsize == 0 || ysize == 0 ) error("new_image_char: invalid image size.");
    
    /* get memory */
    image = (image_char) malloc( sizeof(struct image_char_s) );
    if( image == NULL ) error("not enough memory.");
    image->data = (unsigned char *) calloc( (size_t) (xsize*ysize),
                                           sizeof(unsigned char) );
    if( image->data == NULL ) error("not enough memory.");
    
    /* set image size */
    image->xsize = xsize;
    image->ysize = ysize;
    
    return image;
}

/*----------------------------------------------------------------------------*/
/** Create a new image_char of size 'xsize' times 'ysize',
 initialized to the value 'fill_value'.
 */
image_char new_image_char_ini( unsigned int xsize, unsigned int ysize,
                                     unsigned char fill_value )
{
    image_char image = new_image_char(xsize,ysize); /* create image */
    unsigned int N = xsize*ysize;
    unsigned int i;
    
    /* check parameters */
    if( image == NULL || image->data == NULL )
        error("new_image_char_ini: invalid image.");
    
    /* initialize */
    for(i=0; i<N; i++) image->data[i] = fill_value;
    
    return image;
}

/*----------------------------------------------------------------------------*/
/** Create a new image_int of size 'xsize' times 'ysize'.
 */
image_int new_image_int(unsigned int xsize, unsigned int ysize)
{
    image_int image;
    
    /* check parameters */
    if( xsize == 0 || ysize == 0 ) error("new_image_int: invalid image size.");
    
    /* get memory */
    image = (image_int) malloc( sizeof(struct image_int_s) );
    if( image == NULL ) error("not enough memory.");
    image->data = (int *) calloc( (size_t) (xsize*ysize), sizeof(int) );
    if( image->data == NULL ) error("not enough memory.");
    
    /* set image size */
    image->xsize = xsize;
    image->ysize = ysize;
    
    return image;
}

/*----------------------------------------------------------------------------*/
/** Create a new image_int of size 'xsize' times 'ysize',
 initialized to the value 'fill_value'.
 */
static image_int new_image_int_ini( unsigned int xsize, unsigned int ysize,
                                   int fill_value )
{
    image_int image = new_image_int(xsize,ysize); /* create image */
    unsigned int N = xsize*ysize;
    unsigned int i;
    
    /* initialize */
    for(i=0; i<N; i++) image->data[i] = fill_value;
    
    return image;
}
