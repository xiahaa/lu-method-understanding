//
//  rect.cpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
// done

#include "rect.hpp"

/*----------------------------------------------------------------------------*/
/** Copy one rectangle structure to another.
 */
void rect_copy(struct rect * in, struct rect * out)//in is the src, out is the dst
{
    /* check parameters */
    if( in == NULL || out == NULL ) error("rect_copy: invalid 'in' or 'out'.");
    
    /* copy values */
    out->x1 = in->x1;
    out->y1 = in->y1;
    out->x2 = in->x2;
    out->y2 = in->y2;
    out->width = in->width;
    out->x = in->x;
    out->y = in->y;
    out->theta = in->theta;
    out->dx = in->dx;
    out->dy = in->dy;
    out->polarity = in->polarity;
    out->prec = in->prec;
    out->p = in->p;
}

/*----------------------------------------------------------------------------*/
/** Interpolate y value corresponding to 'x' value given, in
 the line 'x1,y1' to 'x2,y2'; if 'x1=x2' return the smaller
 of 'y1' and 'y2'.
 
 The following restrictions are required:
 - x1 <= x2
 - x1 <= x
 - x  <= x2
 */
double inter_low(double x, double x1, double y1, double x2, double y2)
{
    /* check parameters */
    if( x1 > x2 || x < x1 || x > x2 )
        error("inter_low: unsuitable input, 'x1>x2' or 'x<x1' or 'x>x2'.");
    
    /* interpolation */
    if( double_equal(x1,x2) && y1<y2 ) return y1;
    if( double_equal(x1,x2) && y1>y2 ) return y2;
    return y1 + (x-x1) * (y2-y1) / (x2-x1);
}

/*----------------------------------------------------------------------------*/
/** Interpolate y value corresponding to 'x' value given, in
 the line 'x1,y1' to 'x2,y2'; if 'x1=x2' return the larger
 of 'y1' and 'y2'.
 
 The following restrictions are required:
 - x1 <= x2
 - x1 <= x
 - x  <= x2
 */
double inter_hi(double x, double x1, double y1, double x2, double y2)
{
    /* check parameters */
    if( x1 > x2 || x < x1 || x > x2 )
        error("inter_hi: unsuitable input, 'x1>x2' or 'x<x1' or 'x>x2'.");
    
    /* interpolation */
    if( double_equal(x1,x2) && y1<y2 ) return y2;
    if( double_equal(x1,x2) && y1>y2 ) return y1;
    return y1 + (x-x1) * (y2-y1) / (x2-x1);
}

/*----------------------------------------------------------------------------*/
/** Free memory used by a rectangle iterator.
 */
void ri_del(rect_iter * iter)
{
    if( iter == NULL ) error("ri_del: NULL iterator.");
    free( (void *) iter );
}

/*----------------------------------------------------------------------------*/
/** Check if the iterator finished the full iteration.
 
 See details in \ref rect_iter
 */
int ri_end(rect_iter * i)
{
    /* check input */
    if( i == NULL ) error("ri_end: NULL iterator.");
    
    /* if the current x value is larger than the largest
     x value in the rectangle (vx[2]), we know the full
     exploration of the rectangle is finished. */
    return (double)(i->x) > i->vx[2];
}

/*----------------------------------------------------------------------------*/
/** Increment a rectangle iterator.
 
 See details in \ref rect_iter
 */
void ri_inc(rect_iter * i)
{
    /* check input */
    if( i == NULL ) error("ri_inc: NULL iterator.");
    
    /* if not at end of exploration,
     increase y value for next pixel in the 'column' */
    if( !ri_end(i) ) i->y++;
    
    /* if the end of the current 'column' is reached,
     and it is not the end of exploration,
     advance to the next 'column' */
    while( (double) (i->y) > i->ye && !ri_end(i) )
    {
        /* increase x, next 'column' */
        i->x++;
        
        /* if end of exploration, return */
        if( ri_end(i) ) return;
        
        /* update lower y limit (start) for the new 'column'.
         
         We need to interpolate the y value that corresponds to the
         lower side of the rectangle. The first thing is to decide if
         the corresponding side is
         
         vx[0],vy[0] to vx[3],vy[3] or
         vx[3],vy[3] to vx[2],vy[2]
         
         Then, the side is interpolated for the x value of the
         'column'. But, if the side is vertical (as it could happen if
         the rectangle is vertical and we are dealing with the first
         or last 'columns') then we pick the lower value of the side
         by using 'inter_low'.
         */
        if( (double) i->x < i->vx[3] )
            i->ys = inter_low((double)i->x,i->vx[0],i->vy[0],i->vx[3],i->vy[3]);
        else
            i->ys = inter_low((double)i->x,i->vx[3],i->vy[3],i->vx[2],i->vy[2]);
        
        /* update upper y limit (end) for the new 'column'.
         
         We need to interpolate the y value that corresponds to the
         upper side of the rectangle. The first thing is to decide if
         the corresponding side is
         
         vx[0],vy[0] to vx[1],vy[1] or
         vx[1],vy[1] to vx[2],vy[2]
         
         Then, the side is interpolated for the x value of the
         'column'. But, if the side is vertical (as it could happen if
         the rectangle is vertical and we are dealing with the first
         or last 'columns') then we pick the lower value of the side
         by using 'inter_low'.
         */
        if( (double)i->x < i->vx[1] )
            i->ye = inter_hi((double)i->x,i->vx[0],i->vy[0],i->vx[1],i->vy[1]);
        else
            i->ye = inter_hi((double)i->x,i->vx[1],i->vy[1],i->vx[2],i->vy[2]);
        
        /* new y */
        i->y = (int) ceil(i->ys);
    }
}

/*----------------------------------------------------------------------------*/
/** Create and initialize a rectangle iterator.
 
 See details in \ref rect_iter
 */
rect_iter * ri_ini(struct rect * r)
{
    double vx[4],vy[4];
    int n,offset;
    rect_iter * i;
    
    /* check parameters */
    if( r == NULL ) error("ri_ini: invalid rectangle.");
    
    /* get memory */
    i = (rect_iter *) malloc(sizeof(rect_iter));
    if( i == NULL ) error("ri_ini: Not enough memory.");
    
    /* build list of rectangle corners ordered
     in a circular way around the rectangle */
    //从线段的起点(x1,y1)处的一端开始按照逆时针重构出矩形的四个定点
    vx[0] = r->x1 - r->dy * r->width / 2.0;
    vy[0] = r->y1 + r->dx * r->width / 2.0;
    vx[1] = r->x2 - r->dy * r->width / 2.0;
    vy[1] = r->y2 + r->dx * r->width / 2.0;
    vx[2] = r->x2 + r->dy * r->width / 2.0;
    vy[2] = r->y2 - r->dx * r->width / 2.0;
    vx[3] = r->x1 + r->dy * r->width / 2.0;
    vy[3] = r->y1 - r->dx * r->width / 2.0;
    
    /* compute rotation of index of corners needed so that the first
     point has the smaller x.
     
     if one side is vertical, thus two corners have the same smaller x
     value, the one with the largest y value is selected as the first.
     */
    if( r->x1 < r->x2 && r->y1 <= r->y2 ) offset = 0;
    else if( r->x1 >= r->x2 && r->y1 < r->y2 ) offset = 1;
    else if( r->x1 > r->x2 && r->y1 >= r->y2 ) offset = 2;
    else offset = 3;
    
    /* apply rotation of index. */
    for(n=0; n<4; n++)
    {
        i->vx[n] = vx[(offset+n)%4];
        i->vy[n] = vy[(offset+n)%4];
    }
    
    /* Set an initial condition.
     
     The values are set to values that will cause 'ri_inc' (that will
     be called immediately) to initialize correctly the first 'column'
     and compute the limits 'ys' and 'ye'.
     
     'y' is set to the integer value of vy[0], the starting corner.
     
     'ys' and 'ye' are set to very small values, so 'ri_inc' will
     notice that it needs to start a new 'column'.
     
     The smallest integer coordinate inside of the rectangle is
     'ceil(vx[0])'. The current 'x' value is set to that value minus
     one, so 'ri_inc' (that will increase x by one) will advance to
     the first 'column'.
     */
    i->x = (int) ceil(i->vx[0]) - 1;
    i->y = (int) ceil(i->vy[0]);
    i->ys = i->ye = -DBL_MAX;
    
    /* advance to the first pixel */
    ri_inc(i);
    
    return i;
}
