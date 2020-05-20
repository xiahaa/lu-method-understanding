//
//  utils.cpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
// done

#include "utils.hpp"
#include "datatype.h"

//==================================================================================================
//=============================miscellaneous functions==============================================

/** Compare doubles by relative error.
 
 The resulting rounding error after floating point computations
 depend on the specific operations done. The same number computed by
 different algorithms could present different rounding errors. For a
 useful comparison, an estimation of the relative rounding error
 should be considered and compared to a factor times EPS. The factor
 should be related to the cumulated rounding error in the chain of
 computation. Here, as a simplification, a fixed factor is used.
 */
int double_equal(double a, double b)
{
    double abs_diff,aa,bb,abs_max;
    
    /* trivial case */
    if( a == b ) return TRUE;
    
    abs_diff = fabs(a-b);
    aa = fabs(a);
    bb = fabs(b);
    abs_max = aa > bb ? aa : bb;
    
    /* DBL_MIN is the smallest normalized number, thus, the smallest
     number whose relative error is bounded by DBL_EPSILON. For
     smaller numbers, the same quantization steps as for DBL_MIN
     are used. Then, for smaller numbers, a meaningful "relative"
     error should be computed by dividing the difference by DBL_MIN. */
    if( abs_max < DBL_MIN ) abs_max = DBL_MIN;
    
    /* equal if relative error <= factor x eps */
    return (abs_diff / abs_max) <= (RELATIVE_ERROR_FACTOR * DBL_EPSILON); //RELATIVE_ERROR_FACTOR=100.0,
}

/*----------------------------------------------------------------------------*/
/** Absolute value angle difference.
 */
//得到2个弧度制角度的夹角的绝对值
double angle_diff(double a, double b)
{
    a -= b;
    while( a <= -M_PI ) a += M_2__PI;
    while( a >   M_PI ) a -= M_2__PI;
    if( a < 0.0 ) a = -a;
    return a;
}
/*----------------------------------------------------------------------------*/
/** Signed angle difference.
 */
double angle_diff_signed(double a, double b)
{
    a -= b;
    while( a <= -M_PI ) a += M_2__PI;
    while( a >   M_PI ) a -= M_2__PI;
    return a;
}

/*----------------------------------------------------------------------------*/
/** Fatal error, print a message to standard-error output and exit.
 */
void error(const char * msg)
{
    // fprintf(stderr,"circleDetection Error: %s\n",msg);
    //   exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*/
/** Computes Euclidean distance between point (x1,y1) and point (x2,y2).
 */
double dist(double x1, double y1, double x2, double y2)
{
    return sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
}

//向量内积
double dotProduct(point2d vec1, point2d vec2)
{
    return (vec1.x*vec2.x+vec1.y*vec2.y);
}


/*----------------------------------------------------------------------------*/
/** Is point2i (x,y) aligned to angle theta, up to precision 'prec'?
 */
int isaligned( int x, int y, image_double angles, double theta,
                     double prec )
{
    double a;
    
    /* check parameters */
    if( angles == NULL || angles->data == NULL )
        error("isaligned: invalid image 'angles'.");
    if( x < 0 || y < 0 || x >= (int) angles->xsize || y >= (int) angles->ysize )
        error("isaligned: (x,y) out of the image.");
    if( prec < 0.0 ) error("isaligned: 'prec' must be positive.");
    
    /* angle at pixel (x,y) */
    a = angles->data[ x + y * angles->xsize ];
    
    /* pixels whose level-line angle is not defined
     are considered as NON-aligned */
    if( a == NOTDEF ) return FALSE;  /* there is no need to call the function
                                      'double_equal' here because there is
                                      no risk of problems related to the
                                      comparison doubles, we are only
                                      interested in the exact NOTDEF value */
    
    /* it is assumed that 'theta' and 'a' are in the range [-pi,pi] */
    theta -= a;
    if( theta < 0.0 ) theta = -theta;
    if( theta > M_3_2_PI )
    {
        //--------------------------------------
        //origin code
        /* theta -= M_2__PI;
         if( theta < 0.0 ) theta = -theta;*/
        //--------------------------------------
        //-------------------------------------
        //mycode
        theta = M_2__PI-theta;
        if(theta < 0.0)
            theta = -theta;
        //--------------------------------------
    }
    
    return theta <= prec;
}

