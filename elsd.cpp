//
//  elsd.cpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#include "elsd.hpp"
#include "utils.hpp"
#include <math.h>
#include <limits.h>
#include "rect.hpp"
#include "nfa.hpp"
#include "gauss.hpp"

/*----------------------------------------------------------------------------*/
/*--------------------------------- Gradient ---------------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** Computes the direction of the level line of 'in' at each point2i.
 
 The result is:
 - an image_double with the angle at each pixel, or NOTDEF if not defined.
 - the image_double 'modgrad' (a point2ier is passed as argument)
 with the gradient magnitude at each point2i.
 - a list of pixels 'list_p' roughly ordered by decreasing
 gradient magnitude. (The order is made by classifying point2is
 into bins by gradient magnitude. The parameters 'n_bins' and
 'max_grad' specify the number of bins and the gradient modulus
 at the highest bin. The pixels in the list would be in
 decreasing gradient magnitude, up to a precision of the size of
 the bins.)
 - a point2ier 'mem_p' to the memory used by 'list_p' to be able to
 free the memory when it is not used anymore.
 */
//返回一张梯度角度顺时针旋转90°后的align角度图angles，如果梯度角度是(gx,gy)->(-gy,gx)，
//和梯度的模的图modgrad,然后按照n_bins进行伪排序返回链表的头指针list_p,里面存的是坐标
static image_double ll_angle( image_double in, double threshold,
                             struct coorlist ** list_p,
                             image_double * modgrad, unsigned int n_bins )
{
    image_double g;
    unsigned int n,p,x,y,adr,i;
    double com1,com2,gx,gy,norm,norm2;
    /* the rest of the variables are used for pseudo-ordering
     the gradient magnitude values */
    int list_count = 0;
    //struct coorlist * list;
    struct coorlist *temp;
    struct coorlist ** range_l_s; /* array of point2iers to start of bin list,表示1024个bin的头指针的指针数组 */
    struct coorlist ** range_l_e; /* array of point2iers to end of bin list，表示1024个bin的尾指针的指针数组*/
    struct coorlist * start;
    struct coorlist * end;
    double max_grad = 0.0;
    
    /* check parameters */
    if( in == NULL || in->data == NULL || in->xsize == 0 || in->ysize == 0 )
        error("ll_angle: invalid image.");
    if( threshold < 0.0 ) error("ll_angle: 'threshold' must be positive.");
    if( list_p == NULL ) error("ll_angle: NULL point2ier 'list_p'.");
    // if( mem_p == NULL ) error("ll_angle: NULL point2ier 'mem_p'.");
    if( modgrad == NULL ) error("ll_angle: NULL point2ier 'modgrad'.");
    if( n_bins == 0 ) error("ll_angle: 'n_bins' must be positive.");
    
    /* image size shortcuts */
    n = in->ysize;
    p = in->xsize;
    
    /* allocate output image */
    g = new_image_double(in->xsize,in->ysize);
    
    /* get memory for the image of gradient modulus */
    *modgrad = new_image_double(in->xsize,in->ysize);
    
    /* get memory for "ordered" list of pixels */
    //list = (struct coorlist *) calloc( (size_t) (n*p), sizeof(struct coorlist) );
    //*mem_p = (void *) list;
    range_l_s = (struct coorlist **) calloc( (size_t) n_bins,
                                            sizeof(struct coorlist *) );
    range_l_e = (struct coorlist **) calloc( (size_t) n_bins,
                                            sizeof(struct coorlist *) );
    // if( list == NULL || range_l_s == NULL || range_l_e == NULL )
    if( range_l_s == NULL || range_l_e == NULL )
        error("not enough memory.");
    for(i=0;i<n_bins;i++) range_l_s[i] = range_l_e[i] = NULL;
    
    /* 'undefined' on the down and right boundaries */
    for(x=0;x<p;x++) g->data[(n-1)*p+x] = NOTDEF;// p = in->xsize
    for(y=0;y<n;y++) g->data[p*y+p-1]   = NOTDEF;// n = in->ysize;
    
    /* compute gradient on the remaining pixels */
    for(x=0;x<p-1;x++)
        for(y=0;y<n-1;y++)
        {
            adr = y*p+x;
            /*
             Norm 2 computation using 2x2 pixel window:
             A B
             C D
             and
             com1 = D-A,  com2 = B-C.
             Then
             gx = B+D - (A+C)   horizontal difference
             gy = C+D - (A+B)   vertical difference
             com1 and com2 are just to avoid 2 additions.
             */
            com1 = in->data[adr+p+1] - in->data[adr];
            com2 = in->data[adr+1]   - in->data[adr+p];
            
            gx = com1+com2; /* gradient x component */
            gy = com1-com2; /* gradient y component */
            norm2 = gx*gx+gy*gy;
            norm = sqrt( norm2 / 4.0 ); /* gradient norm */
            
            (*modgrad)->data[adr] = norm; /* store gradient norm */
            
            if( norm <= threshold ) /* norm too small, gradient no defined */
                g->data[adr] = NOTDEF; /* gradient angle not defined */
            else
            {
                /* gradient angle computation */
                g->data[adr] = atan2(gx,-gy);
                /* look for the maximum of the gradient */
                if( norm > max_grad ) max_grad = norm;
            }
        }
    
    /* compute histogram of gradient values */
    for(x=0;x<p-1;x++)
        for(y=0;y<n-1;y++)
        {
            temp = new coorlist();
            if(temp == NULL)
            {
                printf("not enough memory");
                system("pause");
            }
            norm = (*modgrad)->data[y*p+x];
            /* store the point2i in the right bin according to its norm */
            i = (unsigned int) (norm * (double) n_bins / max_grad);
            if( i >= n_bins ) i = n_bins-1;
            if( range_l_e[i] == NULL )
                range_l_s[i] = range_l_e[i] = temp;//记录第i个区域的头指针到range_l_s[i]
            else
            {
                range_l_e[i]->next = temp;//第i个区域由尾指针range_l_e[i]完成勾链
                range_l_e[i] = temp;
            }
            range_l_e[i]->x = (int) x;//将坐标(x,y)记录到第i个分区
            range_l_e[i]->y = (int) y;
            range_l_e[i]->next = NULL;
        }
    
    /* Make the list of pixels (almost) ordered by norm value.
     It starts by the larger bin, so the list starts by the
     pixels with the highest gradient value. Pixels would be ordered
     by norm value, up to a precision given by max_grad/n_bins.
     */
    for(i=n_bins-1; i>0 && range_l_s[i]==NULL; i--);//找到第一个不为空的分区bin
    start = range_l_s[i];
    end = range_l_e[i];
    if( start != NULL )
        while(i>0)
        {
            --i;
            if( range_l_s[i] != NULL )
            {
                end->next = range_l_s[i];
                end = range_l_e[i];
            }
        }
    *list_p = start;
    // *mem_p  = start;
    /* free memory */
    free( (void *) range_l_s );
    free( (void *) range_l_e );
    
    return g;
}


/*---------------------------------- Regions ---------------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** Compute region's angle as the principal inertia axis of the region.
 
 The following is the region inertia matrix A:
 @f[
 
 A = \left(\begin{array}{cc}
 Ixx & Ixy \\
 Ixy & Iyy \\
 \end{array}\right)
 
 @f]
 where
 
 Ixx =   sum_i G(i).(y_i - cx)^2
 
 Iyy =   sum_i G(i).(x_i - cy)^2
 
 Ixy = - sum_i G(i).(x_i - cx).(y_i - cy)
 
 and
 - G(i) is the gradient norm at pixel i, used as pixel's weight.
 - x_i and y_i are the coordinates of pixel i.
 - cx and cy are the coordinates of the center of th region.
 
 lambda1 and lambda2 are the eigenvalues of matrix A,
 with lambda1 >= lambda2. They are found by solving the
 characteristic polynomial:
 
 det( lambda I - A) = 0
 
 that gives:
 
 lambda1 = ( Ixx + Iyy + sqrt( (Ixx-Iyy)^2 + 4.0*Ixy*Ixy) ) / 2
 
 lambda2 = ( Ixx + Iyy - sqrt( (Ixx-Iyy)^2 + 4.0*Ixy*Ixy) ) / 2
 
 To get the line segment direction we want to get the angle the
 eigenvector associated to the smallest eigenvalue. We have
 to solve for a,b in:
 
 a.Ixx + b.Ixy = a.lambda2
 
 a.Ixy + b.Iyy = b.lambda2
 
 We want the angle theta = atan(b/a). It can be computed with
 any of the two equations:
 
 theta = atan( (lambda2-Ixx) / Ixy )
 
 or
 
 theta = atan( Ixy / (lambda2-Iyy) )
 
 When |Ixx| > |Iyy| we use the first, otherwise the second (just to
 get better numeric precision).
 */
static double get_theta( point2i * reg, int reg_size, double x, double y,
                        image_double modgrad, double reg_angle, double prec )
{
    double lambda,theta,weight;
    double Ixx = 0.0;
    double Iyy = 0.0;
    double Ixy = 0.0;
    double temp1,temp2;
    int i;
    
    /* check parameters */
    if( reg == NULL ) error("get_theta: invalid region.");
    if( reg_size <= 1 ) error("get_theta: region size <= 1.");
    if( modgrad == NULL || modgrad->data == NULL )
        error("get_theta: invalid 'modgrad'.");
    if( prec < 0.0 ) error("get_theta: 'prec' must be positive.");
    
    /* compute inertia matrix */
    for(i=0; i<reg_size; i++)
    {
        weight = modgrad->data[ reg[i].x + reg[i].y * modgrad->xsize ];
        Ixx += ( (double) reg[i].y - y ) * ( (double) reg[i].y - y ) * weight;
        Iyy += ( (double) reg[i].x - x ) * ( (double) reg[i].x - x ) * weight;
        Ixy -= ( (double) reg[i].x - x ) * ( (double) reg[i].y - y ) * weight;
    }
    if( double_equal(Ixx,0.0) && double_equal(Iyy,0.0) && double_equal(Ixy,0.0) )//判断Ixx、Iyy、Ixy与0是否非常接近，由于它们为double类型，故需要专门的函数判断
        error("get_theta: null inertia matrix.");
    
    /* compute smallest eigenvalue */
    lambda = 0.5 * ( Ixx + Iyy - sqrt( (Ixx-Iyy)*(Ixx-Iyy) + 4.0*Ixy*Ixy ) );
    
    /* compute angle */
    theta = fabs(Ixx)>fabs(Iyy) ? atan2(lambda-Ixx,Ixy) : atan2(Ixy,lambda-Iyy);
    /* The previous procedure doesn't cares about orientation,
     so it could be wrong by 180 degrees. Here is corrected if necessary. */
    temp1 = angle_diff(theta,reg_angle);
    if( temp1 > prec )//这是由于用惯性矩阵算出的两个正交轴的较小特征值对应的角度和该区域的角度可能相差180°
    {
        //------------------------------------------
        //theta += M_PI;   //origin code
        //------------------------------------------
        //------------------------------------------
        //my code,增加该段代码，限制theta在 (-pi,pi)之间
        //int flag = 0;
        temp2 = angle_diff(theta+M_PI,reg_angle);// add xiaohu, been solved in elsd
        if(temp2 < prec)
        {
            theta += M_PI;
            if(theta > M_PI)
            {
                theta -= M_2__PI;
                //flag = 1;
                //if(angle_diff(theta,reg_angle) > prec)
                //{
                //   //flag = 2;
                //   theta = reg_angle;
                // }
            }
        }
        else
        {
            theta = (temp2 <= temp1) ? (theta+M_PI) : theta;
            while( theta <= -M_PI ) theta += M_2__PI;
            while( theta >   M_PI ) theta -= M_2__PI;
        }
        
        //--------------------------------------------
    }
    return theta;
}


/*----------------------------------------------------------------------------*/
/** Computes a rectangle that covers a region of point2is.
 */
static void region2rect( point2i * reg, int reg_size,
                        image_double modgrad, double reg_angle,
                        double prec, double p, struct rect * rec )
{
    double x,y,dx,dy,l,w,theta,weight,sum,l_min,l_max,w_min,w_max;
    int i;
    
    /* check parameters */
    if( reg == NULL ) error("region2rect: invalid region.");
    if( reg_size <= 1 ) error("region2rect: region size <= 1.");
    if( modgrad == NULL || modgrad->data == NULL )
        error("region2rect: invalid image 'modgrad'.");
    if( rec == NULL ) error("region2rect: invalid 'rec'.");
    
    /* center of the region:
     
     It is computed as the weighted sum of the coordinates
     of all the pixels in the region. The norm of the gradient
     is used as the weight of a pixel. The sum is as follows:
     cx = \sum_i G(i).x_i
     cy = \sum_i G(i).y_i
     where G(i) is the norm of the gradient of pixel i
     and x_i,y_i are its coordinates.
     */
    //获得质心 x,y
    x = y = sum = 0.0;
    for(i=0; i<reg_size; i++)
    {
        weight = modgrad->data[ reg[i].x + reg[i].y * modgrad->xsize ];
        x += (double) reg[i].x * weight;
        y += (double) reg[i].y * weight;
        sum += weight;
    }
    if( sum <= 0.0 ) error("region2rect: weights sum equal to zero.");
    x /= sum;
    y /= sum;
    
    /* theta */
    //运用惯性矩阵获得更为精确的角度估计
    theta = get_theta(reg,reg_size,x,y,modgrad,reg_angle,prec);
    dx = cos(theta);
    dy = sin(theta);
    
    /* length and width:
     
     'l' and 'w' are computed as the distance from the center of the
     region to pixel i, projected along the rectangle axis (dx,dy) and
     to the orthogonal axis (-dy,dx), respectively.
     
     The length of the rectangle goes from l_min to l_max, where l_min
     and l_max are the minimum and maximum values of l in the region.
     Analogously, the width is selected from w_min to w_max, where
     w_min and w_max are the minimum and maximum of w for the pixels
     in the region.
     */
    //因为区域的方向向量为 (dx,dy)
    /*
     ------------------->x
     |\
     | \
     |  \(dx,dy)
     |
     \|/
     y
     因此顺时针旋转90°是 (-dy,dx)
     */
    l_min = l_max = w_min = w_max = 0.0;
    for(i=0; i<reg_size; i++)//用向量内积求在线段方向和与线段方向垂直方向的投影求l,w
    {
        l =  ( (double) reg[i].x - x) * dx + ( (double) reg[i].y - y) * dy;
        w = -( (double) reg[i].x - x) * dy + ( (double) reg[i].y - y) * dx;
        
        if( l > l_max ) l_max = l;
        if( l < l_min ) l_min = l;
        if( w > w_max ) w_max = w;
        if( w < w_min ) w_min = w;
    }
    
    /* store values */
    rec->x1 = x + l_min * dx;
    rec->y1 = y + l_min * dy;
    rec->x2 = x + l_max * dx;
    rec->y2 = y + l_max * dy;
    rec->width = w_max - w_min;
    rec->x = x;
    rec->y = y;
    rec->theta = theta;
    rec->dx = dx;
    rec->dy = dy;
    rec->prec = prec;
    rec->p = p;
    
    /* we impose a minimal width of one pixel
     
     A sharp horizontal or vertical step would produce a perfectly
     horizontal or vertical region. The width computed would be
     zero. But that corresponds to a one pixels width transition in
     the image.
     */
    if( rec->width < 1.0 )
        rec->width = 1.0;
}

//区域质心和角度已经计算好了，因此只进行矩形近似。而region2rect此外还进行了质心和角度计算。
static void region2rect2(point2i * reg, int reg_size,double reg_center_x,double reg_center_y,
                         double reg_theta,double prec, double p, struct rect * rec )
{
    double dx,dy,l,w,l_min,l_max,w_min,w_max;
    int i;
    /* check parameters */
    if( reg == NULL ) error("region2rect: invalid region.");
    if( reg_size <= 1 ) error("region2rect: region size <= 1.");
    if( rec == NULL ) error("region2rect: invalid 'rec'.");
    
    //获得区域的方向向量(dx,dy)
    dx = cos(reg_theta);
    dy = sin(reg_theta);
    l_min = l_max = w_min = w_max = 0.0;
    for(i=0; i<reg_size; i++)//用向量内积求在线段方向和与线段方向垂直方向的投影求l,w
    {
        l =  ( (double) reg[i].x - reg_center_x) * dx + ( (double) reg[i].y - reg_center_y) * dy;
        w = -( (double) reg[i].x - reg_center_x) * dy + ( (double) reg[i].y - reg_center_y) * dx;
        
        if( l > l_max ) l_max = l;
        if( l < l_min ) l_min = l;
        if( w > w_max ) w_max = w;
        if( w < w_min ) w_min = w;
    }
    
    /* store values */
    rec->x1 = reg_center_x + l_min * dx;
    rec->y1 = reg_center_y + l_min * dy;
    rec->x2 = reg_center_x + l_max * dx;
    rec->y2 = reg_center_y + l_max * dy;
    rec->width = w_max - w_min;
    rec->x = reg_center_x;
    rec->y = reg_center_y;
    rec->theta = reg_theta;
    rec->dx = dx;
    rec->dy = dy;
    rec->prec = prec;
    rec->p = p;
    
    /* we impose a minimal width of one pixel
     
     A sharp horizontal or vertical step would produce a perfectly
     horizontal or vertical region. The width computed would be
     zero. But that corresponds to a one pixels width transition in
     the image.
     */
    if( rec->width < 1.0 )
        rec->width = 1.0;
}


/*----------------------------------------------------------------------------*/
/** Build a region of pixels that share the same angle, up to a
 tolerance 'prec', starting at point2i (x,y).
 */
static void region_grow( int x, int y, image_double angles, struct point2i * reg,
                        int * reg_size, double * reg_angle, image_char used,
                        double prec )
{
    double sumdx,sumdy;
    int xx,yy,i;
    
    /* check parameters */
    if( x < 0 || y < 0 || x >= (int) angles->xsize || y >= (int) angles->ysize )
        error("region_grow: (x,y) out of the image.");
    if( angles == NULL || angles->data == NULL )
        error("region_grow: invalid image 'angles'.");
    if( reg == NULL ) error("region_grow: invalid 'reg'.");
    if( reg_size == NULL ) error("region_grow: invalid point2ier 'reg_size'.");
    if( reg_angle == NULL ) error("region_grow: invalid point2ier 'reg_angle'.");
    if( used == NULL || used->data == NULL )
        error("region_grow: invalid image 'used'.");
    
    /* first point2i of the region */
    *reg_size = 1;
    reg[0].x = x;
    reg[0].y = y;
    *reg_angle = angles->data[x+y*angles->xsize];  /* region's angle */
    sumdx = cos(*reg_angle);
    sumdy = sin(*reg_angle);
    used->data[x+y*used->xsize] = USED;
    
    /* try neighbors as new region point2is */
    for(i=0; i<*reg_size; i++)
        for(xx=reg[i].x-1; xx<=reg[i].x+1; xx++)
            for(yy=reg[i].y-1; yy<=reg[i].y+1; yy++)
                if( xx>=0 && yy>=0 && xx<(int)used->xsize && yy<(int)used->ysize &&
                   used->data[xx+yy*used->xsize] != USED &&
                   isaligned(xx,yy,angles,*reg_angle,prec) )
                {
                    /* add point2i */
                    used->data[xx+yy*used->xsize] = USED;
                    reg[*reg_size].x = xx;
                    reg[*reg_size].y = yy;
                    ++(*reg_size);
                    
                    /* update region's angle */
                    sumdx += cos( angles->data[xx+yy*angles->xsize] );
                    sumdy += sin( angles->data[xx+yy*angles->xsize] );
                    *reg_angle = atan2(sumdy,sumdx);
                }
}



/*----------------------------------------------------------------------------*/
/** Try some rectangles variations to improve NFA value. Only if the
 rectangle is not meaningful (i.e., log_nfa <= log_eps).
 */
static double rect_improve( struct rect * rec, image_double angles,
                           double logNT, double log_eps )
{
    struct rect r;
    double log_nfa,log_nfa_new;
    double delta = 0.5;
    double delta_2 = delta / 2.0;
    int n;
    
    log_nfa = rect_nfa(rec,angles,logNT);
    
    if( log_nfa > log_eps ) return log_nfa;
    
    /* try finer precisions */
    rect_copy(rec,&r);
    for(n=0; n<5; n++)
    {
        r.p /= 2.0;
        r.prec = r.p * M_PI;
        log_nfa_new = rect_nfa(&r,angles,logNT);
        if( log_nfa_new > log_nfa )
        {
            log_nfa = log_nfa_new;
            rect_copy(&r,rec);
        }
    }
    
    if( log_nfa > log_eps ) return log_nfa;
    
    /* try to reduce width */
    rect_copy(rec,&r);
    for(n=0; n<5; n++)
    {
        if( (r.width - delta) >= 0.5 )
        {
            r.width -= delta;
            log_nfa_new = rect_nfa(&r,angles,logNT);
            if( log_nfa_new > log_nfa )
            {
                rect_copy(&r,rec);
                log_nfa = log_nfa_new;
            }
        }
    }
    
    if( log_nfa > log_eps ) return log_nfa;
    
    /* try to reduce one side of the rectangle */
    rect_copy(rec,&r);
    for(n=0; n<5; n++)
    {
        if( (r.width - delta) >= 0.5 )
        {
            r.x1 += -r.dy * delta_2;
            r.y1 +=  r.dx * delta_2;
            r.x2 += -r.dy * delta_2;
            r.y2 +=  r.dx * delta_2;
            r.width -= delta;
            log_nfa_new = rect_nfa(&r,angles,logNT);
            if( log_nfa_new > log_nfa )
            {
                rect_copy(&r,rec);
                log_nfa = log_nfa_new;
            }
        }
    }
    
    if( log_nfa > log_eps ) return log_nfa;
    
    /* try to reduce the other side of the rectangle */
    rect_copy(rec,&r);
    for(n=0; n<5; n++)
    {
        if( (r.width - delta) >= 0.5 )
        {
            r.x1 -= -r.dy * delta_2;
            r.y1 -=  r.dx * delta_2;
            r.x2 -= -r.dy * delta_2;
            r.y2 -=  r.dx * delta_2;
            r.width -= delta;
            log_nfa_new = rect_nfa(&r,angles,logNT);
            if( log_nfa_new > log_nfa )
            {
                rect_copy(&r,rec);
                log_nfa = log_nfa_new;
            }
        }
    }
    
    if( log_nfa > log_eps ) return log_nfa;
    
    /* try even finer precisions */
    rect_copy(rec,&r);
    for(n=0; n<5; n++)
    {
        r.p /= 2.0;
        r.prec = r.p * M_PI;
        log_nfa_new = rect_nfa(&r,angles,logNT);
        if( log_nfa_new > log_nfa )
        {
            log_nfa = log_nfa_new;
            rect_copy(&r,rec);
        }
    }
    
    return log_nfa;
}



/*----------------------------------------------------------------------------*/
/** Reduce the region size, by elimination the point2is far from the
 starting point2i, until that leads to rectangle with the right
 density of region point2is or to discard the region if too small.
 */
static int reduce_region_radius( struct point2i * reg, int * reg_size,
                                image_double modgrad, double reg_angle,
                                double prec, double p, struct rect * rec,
                                image_char used, image_double angles,
                                double density_th )
{
    double density,rad1,rad2,rad,xc,yc;
    int i;
    
    /* check parameters */
    if( reg == NULL ) error("reduce_region_radius: invalid point2ier 'reg'.");
    if( reg_size == NULL )
        error("reduce_region_radius: invalid point2ier 'reg_size'.");
    if( prec < 0.0 ) error("reduce_region_radius: 'prec' must be positive.");
    if( rec == NULL ) error("reduce_region_radius: invalid point2ier 'rec'.");
    if( used == NULL || used->data == NULL )
        error("reduce_region_radius: invalid image 'used'.");
    if( angles == NULL || angles->data == NULL )
        error("reduce_region_radius: invalid image 'angles'.");
    
    /* compute region point2is density */ //该密度判断已经在函数外判断过，应该可以不用在判断了吧
    density = (double) *reg_size / ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );
    
    // if the density criterion is satisfied there is nothing to do
    if( density >= density_th ) return TRUE;
    
    /* compute region's radius */
    xc = (double) reg[0].x;
    yc = (double) reg[0].y;
    rad1 = dist( xc, yc, rec->x1, rec->y1 );
    rad2 = dist( xc, yc, rec->x2, rec->y2 );
    rad = rad1 > rad2 ? rad1 : rad2;
    
    /* while the density criterion is not satisfied, remove farther pixels */
    while( density < density_th )
    {
        rad *= 0.75; /* reduce region's radius to 75% of its value */
        
        /* remove point2is from the region and update 'used' map */
        for(i=0; i<*reg_size; i++)
            if( dist( xc, yc, (double) reg[i].x, (double) reg[i].y ) > rad )
            {
                /* point2i not kept, mark it as NOTUSED */
                used->data[ reg[i].x + reg[i].y * used->xsize ] = NOTUSED;
                /* remove point2i from the region */
                reg[i].x = reg[*reg_size-1].x; /* if i==*reg_size-1 copy itself */
                reg[i].y = reg[*reg_size-1].y;
                --(*reg_size);
                --i; /* to avoid skipping one point2i */
            }
        
        /* reject if the region is too small.
         2 is the minimal region size for 'region2rect' to work. */
        if( *reg_size < 2 ) return FALSE;
        
        /* re-compute rectangle */
        region2rect(reg,*reg_size,modgrad,reg_angle,prec,p,rec);
        
        /* re-compute region point2is density */
        density = (double) *reg_size /
        ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );
    }
    
    /* if this point2i is reached, the density criterion is satisfied */
    return TRUE;
}

/*----------------------------------------------------------------------------*/
/** Refine a rectangle.
 
 For that, an estimation of the angle tolerance is performed by the
 standard deviation of the angle at point2is near the region's
 starting point2i. Then, a new region is grown starting from the same
 point2i, but using the estimated angle tolerance. If this fails to
 produce a rectangle with the right density of region point2is,
 'reduce_region_radius' is called to try to satisfy this condition.
 */
static int refine( struct point2i * reg, int * reg_size, image_double modgrad,
                  double reg_angle, double prec, double p, struct rect * rec,
                  image_char used, image_double angles, double density_th )
{
    double angle,ang_d,mean_angle,tau,density,xc,yc,ang_c,sum,s_sum;
    int i,n;
    
    /* check parameters */
    if( reg == NULL ) error("refine: invalid point2ier 'reg'.");
    if( reg_size == NULL ) error("refine: invalid point2ier 'reg_size'.");
    if( prec < 0.0 ) error("refine: 'prec' must be positive.");
    if( rec == NULL ) error("refine: invalid point2ier 'rec'.");
    if( used == NULL || used->data == NULL )
        error("refine: invalid image 'used'.");
    if( angles == NULL || angles->data == NULL )
        error("refine: invalid image 'angles'.");
    
    /* compute region point2is density */
    density = (double) *reg_size /
    ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );
    
    /* if the density criterion is satisfied there is nothing to do */
    if( density >= density_th ) return TRUE;
    
    /*------ First try: reduce angle tolerance ------*/
    
    /* compute the new mean angle and tolerance */
    xc = (double) reg[0].x;
    yc = (double) reg[0].y;
    ang_c = angles->data[ reg[0].x + reg[0].y * angles->xsize ];
    sum = s_sum = 0.0;
    n = 0;
    for(i=0; i<*reg_size; i++)
    {
        used->data[ reg[i].x + reg[i].y * used->xsize ] = NOTUSED;
        if( dist( xc, yc, (double) reg[i].x, (double) reg[i].y ) < rec->width )
        {
            angle = angles->data[ reg[i].x + reg[i].y * angles->xsize ];
            ang_d = angle_diff_signed(angle,ang_c);
            sum += ang_d;//加上角度差
            s_sum += ang_d * ang_d;//加上角度差的平方
            ++n;
        }
    }
    mean_angle = sum / (double) n;
    //以2倍标准差作为新的角度容忍度，最开始为22.5°*pi/180
    tau = 2.0 * sqrt( (s_sum - 2.0 * mean_angle * sum) / (double) n  +  mean_angle*mean_angle ); /* 2 * standard deviation */
    //以新的角度容忍度重新进行区域生长
    /* find a new region from the same starting point2i and new angle tolerance */
    region_grow(reg[0].x,reg[0].y,angles,reg,reg_size,&reg_angle,used,tau);
    
    /* if the region is too small, reject */
    if( *reg_size < 2 ) return FALSE;
    
    /* re-compute rectangle */
    region2rect(reg,*reg_size,modgrad,reg_angle,prec,p,rec);
    
    /* re-compute region point2is density */
    density = (double) *reg_size /
    ( dist(rec->x1,rec->y1,rec->x2,rec->y2) * rec->width );
    
    /*------ Second try: reduce region radius ------*/
    if( density < density_th )
        return reduce_region_radius( reg, reg_size, modgrad, reg_angle, prec, p,
                                    rec, used, angles, density_th );
    
    /* if this point2i is reached, the density criterion is satisfied */
    return TRUE;
}

//--------------------------------------------------------
//my code, read this, check the paper. add xiaohu done. too tricky
static bool isArcSegment(point2i * reg, int reg_size, struct rect * main_rect, image_double ll_angles,image_char used,image_char pol,
                  double prec, double p, rect * rect_up, rect * rect_down)
{
    point2i * reg_up = (point2i*)malloc(reg_size*sizeof(point2i));
    point2i * reg_down = (point2i*)malloc(reg_size*sizeof(point2i));
    int   reg_up_size,reg_down_size;
    double reg_up_theta,reg_down_theta, main_theta;
    double reg_up_sin_s,reg_up_cos_s,reg_down_sin_s,reg_down_cos_s;
    double reg_up_x,reg_up_y,reg_down_x,reg_down_y;
    //double weight,sum;
    double temp1,temp2;
    int same_pol_cnt,opp_pol_cnt;
    int i;
    
    same_pol_cnt = opp_pol_cnt = 0;
    reg_up_size = reg_down_size = 0;
    
    for ( i = 0; i < reg_size; i++)
    {
        switch(pol->data[reg[i].y*pol->xsize+reg[i].x])
        {
            case SAME_POL: same_pol_cnt++;break;//统计同极性的pixel数量
            case OPP_POL : opp_pol_cnt++; break;//统计反极性的pixel数量
            default:break;
        }
        //选与theta角度为法线方向，过质心的直线方程为 dx*(x-xi)+dy*(y-yi)=0,则与方向相同的点代入方程得到距离d,d>=0归入reg_up,d<0归入reg_down
        if( main_rect->dx*( reg[i].x - main_rect->x ) + main_rect->dy*( reg[i].y - main_rect->y ) >= 0)
            reg_up[reg_up_size++] = reg[i];
        else
            reg_down[reg_down_size++] = reg[i];
    }
    //对于已经被标记过极性的区域，我们没必要再进行极性分析
    if( (same_pol_cnt + opp_pol_cnt) > reg_size/2)
    {
        if(same_pol_cnt > opp_pol_cnt )
        {
            main_rect->polarity = 1;
            rect_up->polarity = 1;
            rect_down->polarity = 1;
        }
        else
        {
            main_rect->polarity = -1;
            rect_up->polarity = -1;
            rect_down->polarity = -1;
        }
        return TRUE;
    }
    //计算与主方向相同的上半部分区域质心, this is not center of mass, but a mean of sin cos values
    reg_up_x = reg_up_y = 0;
    //sum = 0;
    reg_up_sin_s = reg_up_cos_s = 0;
    for ( i = 0; i< reg_up_size; i++)
    {
        //weight = modgrad->data[ reg_up[i].x + reg_up[i].y * modgrad->xsize ];
        //reg_up_x += (double)weight*reg_up[i].x;
        //reg_up_y += (double)weight*reg_up[i].y;
        //sum += weight;
        reg_up_sin_s += sin(ll_angles->data[ reg_up[i].x + reg_up[i].y * ll_angles->xsize ]);
        reg_up_cos_s += cos(ll_angles->data[ reg_up[i].x + reg_up[i].y * ll_angles->xsize ]);
    }
    //reg_up_x /= sum;
    //reg_up_y /= sum;
    reg_up_theta = atan2(reg_up_sin_s,reg_up_cos_s);
    //计算主方向上的下半部分区域质心
    reg_down_x = reg_down_y = 0;
    //sum = 0;
    reg_down_sin_s = reg_down_cos_s = 0;
    for ( i = 0; i< reg_down_size; i++)
    {
        //weight = modgrad->data[ reg_down[i].x + reg_down[i].y * modgrad->xsize ];
        //reg_down_x += (double)weight*reg_down[i].x;
        //reg_down_y += (double)weight*reg_down[i].y;
        //sum += weight;
        reg_down_sin_s += sin(ll_angles->data[ reg_down[i].x + reg_down[i].y * ll_angles->xsize ]);
        reg_down_cos_s += cos(ll_angles->data[ reg_down[i].x + reg_down[i].y * ll_angles->xsize ]);
    }
    //reg_down_x /= sum;
    //reg_down_y /= sum;
    reg_down_theta = atan2(reg_down_sin_s,reg_down_cos_s);
    main_theta  = atan2(reg_up_sin_s+reg_down_sin_s,reg_up_cos_s+reg_down_cos_s);
    //估计两个区域方向
    //reg_up_theta = get_theta(reg_up,reg_up_size,reg_up_x,reg_up_y,modgrad,main_rect->theta,prec);
    //reg_down_theta = get_theta(reg_down,reg_down_size,reg_down_x,reg_down_y,modgrad,main_rect->theta,prec);
    //旋转到0°进行比较theta,reg_up_theta,reg_down_theta
    temp1 = angle_diff_signed(reg_up_theta,main_theta);
    temp2 = angle_diff_signed(reg_down_theta,main_theta);
    /*if(temp1>= M_PI/2 || temp1 <= -M_PI/2)
     temp1 += 0;
     if(temp2>= M_PI/2 || temp2 <= -M_PI/2)
     temp2 += 0;*/
    //if(temp1 >= prec/10 && temp2 <= -prec/10)//顺时针,边缘的梯度方向与弧的指向圆心方向相反，polarity = -1
    if(temp1 >= M_1_8_PI/10 && temp2 <= -M_1_8_PI/10)//实验证明取定值效果更好
    {// full of trick, no theoratical sound.
        main_rect->polarity = -1;
        rect_up->polarity = -1;
        rect_down->polarity = -1;
        //标记极性
        for ( i = 0; i < reg_size; i++)
        {
            pol->data[reg[i].y*pol->xsize+reg[i].x] = OPP_POL;//-1
        }
    }
    //else if(temp1 <= -prec/10 && temp2 >= prec/10)//逆时针，边缘的梯度方向与弧的指向圆心方向相同，polarity = 1
    else if(temp1 <= -M_1_8_PI/10 && temp2 >= M_1_8_PI/10)//实验证明取定值效果更好
    {
        main_rect->polarity = 1;
        rect_up->polarity = 1;
        rect_down->polarity = 1;
        //标记极性
        for ( i = 0; i < reg_size; i++)
        {
            pol->data[reg[i].y*pol->xsize+reg[i].x] = SAME_POL;//1
        }
    }
    else
    {
        //在region_grow中已经置为USED了
        //for ( i = 0; i< reg_size; i++)
        //  used->data[reg[i].y*used->xsize+reg[i].x] = USED;
        return FALSE;
    }
    
    //region2rect2(reg_up,reg_up_size,reg_up_x,reg_up_y,reg_up_theta,prec,p,rect_up);
    //region2rect2(reg_down,reg_down_size,reg_down_x,reg_down_y,reg_down_theta,prec,p,rect_down);
    
    free(reg_up);
    free(reg_down);
    return TRUE;
}

/*----------------------------------------------------------------------------*/
/*-------------------------- Line Segment Detector ---------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** LSD full interface.
 */
double * LineSegmentDetection( int * n_out,
                              double * img, int X, int Y,
                              double scale, double sigma_scale, double quant,
                              double ang_th, double log_eps, double density_th,
                              int n_bins,
                              int ** reg_img, int * reg_x, int * reg_y )
{
    image_double image;
    ntuple_list out = new_ntuple_list(8);
    double * return_value;
    image_double scaled_image,angles,modgrad;
    image_char used;
    image_char pol;  //对于构成圆弧的像素标记极性，如果梯度的方向和弧的方向指向一致，则为SAME_POLE,否则为OPP_POLE,该标记初始是为0
    image_int region = NULL;
    struct coorlist * list_p;
    struct coorlist * list_p_temp;
    //  struct coorlist * mem_p;
    struct rect main_rect;//main rect
    struct rect rect_up,rect_down;//divide the rect into 2 rects:rect_up and rect_down
    struct point2i * reg;
    int reg_size,min_reg_size,i;
    unsigned int xsize,ysize;
    double rho,reg_angle,prec,p;
    double log_nfa = -1,logNT;
    //  double log_nfa1,log_nfa2;
    int ls_count = 0;                   /* line segments are numbered 1,2,3,... */
    int seed_cnt = 0;
    int refine_cnt = 0;
    int reg_size_toosmall_cnt=0;
    
    /* check parameters */
    if( img == NULL || X <= 0 || Y <= 0 ) error("invalid image input.");
    if( scale <= 0.0 ) error("'scale' value must be positive.");
    if( sigma_scale <= 0.0 ) error("'sigma_scale' value must be positive.");
    if( quant < 0.0 ) error("'quant' value must be positive.");
    if( ang_th <= 0.0 || ang_th >= 180.0 )
        error("'ang_th' value must be in the range (0,180).");
    if( density_th < 0.0 || density_th > 1.0 )
        error("'density_th' value must be in the range [0,1].");
    if( n_bins <= 0 ) error("'n_bins' value must be positive.");
    
    /* angle tolerance */
    prec = M_PI * ang_th / 180.0;
    p = ang_th / 180.0;
    
    rho = quant / sin(prec); /* gradient magnitude threshold */
    
    /* load and scale image (if necessary) and compute angle at each pixel */
    image = new_image_double_ptr( (unsigned int) X, (unsigned int) Y, img );
    if( scale != 1.0 )
    {
        //按照scale进行高斯降采样的图像，注意宽高是上取整，设采样后高宽为imgN*imgM
        scaled_image = gaussian_sampler( image, scale, sigma_scale );
        //返回一张梯度角度顺时针旋转90°后的align角度图angles，如果梯度角度是(gx,gy)->(-gy,gx)，
        //和梯度的模的图modgrad,然后按照n_bins进行伪排序返回链表的头指针list_p,里面存的是坐标
        angles = ll_angle( scaled_image, rho, &list_p,&modgrad, (unsigned int) n_bins );
        free_image_double(scaled_image);
    }
    else
        angles = ll_angle( image, rho, &list_p,&modgrad,(unsigned int) n_bins );
    xsize = angles->xsize;//降采样后的图像的x size，宽度imgM
    ysize = angles->ysize;//降采样后的图像的y size，高度imgN
    
    /* Number of Tests - NT
     
     The theoretical number of tests is Np.(XY)^(5/2)
     where X and Y are number of columns and rows of the image.
     Np corresponds to the number of angle precisions considered.
     As the procedure 'rect_improve' tests 5 times to halve the
     angle precision, and 5 more times after improving other factors,
     11 different precision values are potentially tested. Thus,
     the number of tests is
     11 * (X*Y)^(5/2)
     whose logarithm value is
     log10(11) + 5/2 * (log10(X) + log10(Y)).
     */
    logNT = 5.0 * ( log10( (double) xsize ) + log10( (double) ysize ) ) / 2.0
    + log10(11.0);
    min_reg_size = (int) (-logNT/log10(p)); /* minimal number of point2is in region that can give a meaningful event，每个矩形区域内align point2i最小数量*/
    /* initialize some structures */
    if( reg_img != NULL && reg_x != NULL && reg_y != NULL ) /* save region data */
        region = new_image_int_ini(angles->xsize,angles->ysize,0);//申请与降采样后图像一样大小的int类型的内存，该内存的作用是将检测到的线段序号标到相应的图像格子里，该部分可有可无
    used = new_image_char_ini(xsize,ysize,NOTUSED);//申请与降采样后图像一样大小的char类型的内存
    pol  = new_image_char_ini(xsize,ysize,NOTDEF_POL);//像素点处的梯度和弧指向的方向的极性标记
    reg = (struct point2i *) calloc( (size_t) (xsize*ysize), sizeof(struct point2i) );
    if( reg == NULL ) error("not enough memory!");
    
    list_p_temp = list_p;//记录头链表的头指针，后面需要利用该头指针进行内存释放
    /* search for line segments */
    for(; list_p_temp != NULL; list_p_temp = list_p_temp->next )
        if( used->data[ list_p_temp->x + list_p_temp->y * used->xsize ] == NOTUSED &&
           angles->data[ list_p_temp->x + list_p_temp->y * angles->xsize ] != NOTDEF )
        /* there is no risk of double comparison problems here
         because we are only interested in the exact NOTDEF value */
        {
            /* find the region of connected point2i and ~equal angle */
            //reg是长度为imgN*imgM的一维point2i型数组，有足够大的空间存储生长的区域，reg_size是里面存储了数据的数量，记录的是区域的point2i
            //reg_angle是该区域的主方向的double型变量，存的角度是弧度制
            seed_cnt ++;
            region_grow( list_p_temp->x, list_p_temp->y, angles, reg, &reg_size,&reg_angle, used, prec );
            
            /* reject small regions */
            if( reg_size < min_reg_size )
            {
                reg_size_toosmall_cnt++;
                continue;
            }
            
            /* construct rectangular approximation for the region */
            //根据生长的区域得到近似外接矩阵的参数，矩形参数包括:起点，终点，方向theta，宽度等
            region2rect(reg,reg_size,modgrad,reg_angle,prec,p,&main_rect);
            // very tricky since what you get is a rect (ls) you want to justify it
            // is a arc.... added xiaohu
            if( FALSE == isArcSegment(reg,reg_size,&main_rect,angles,used,pol,prec,p,&rect_up,&rect_down))
                continue;
            /* Check if the rectangle exceeds the minimal density of
             region point2is. If not, try to improve the region.
             The rectangle will be rejected if the final one does
             not fulfill the minimal density condition.
             This is an addition to the original LSD algorithm published in
             "LSD: A Fast Line Segment Detector with a False Detection Control"
             by R. Grompone von Gioi, J. Jakubowicz, J.M. Morel, and G. Randall.
             The original algorithm is obtained with density_th = 0.0.
             */
            
            //提纯，通过重新生长区域来达到期望的密度阈值
            if( !refine( reg, &reg_size, modgrad, reg_angle,
                        prec, p, &main_rect, used, angles, density_th ) ) continue;
            
            refine_cnt++;
            // compute NFA value
            log_nfa = rect_improve(&main_rect,angles,logNT,log_eps);//通过改善矩形区域以尝试得到期望的nfa值
            if( log_nfa <= log_eps ) //错误控制
                continue;
            // A New Line Segment was found!
            ++ls_count;  // increase line segment counter
            
            //
            //  The gradient was computed with a 2x2 mask, its value corresponds to
            //  point2is with an offset of (0.5,0.5), that should be added to output.
            //  The coordinates origin is at the center of pixel (0,0).
            //
            main_rect.x1 += 0.5; main_rect.y1 += 0.5;
            main_rect.x2 += 0.5; main_rect.y2 += 0.5;
            
            // scale the result values if a subsampling was performed */
            if( scale != 1.0 )
            {
                main_rect.x1 /= scale; main_rect.y1 /= scale;
                main_rect.x2 /= scale; main_rect.y2 /= scale;
                //  main_rect.width /= scale;
            }
            
            /* add line segment found to output */
            add_8tuple( out, main_rect.x1, main_rect.y1, main_rect.x2, main_rect.y2,main_rect.dx,main_rect.dy,
                       dist(main_rect.x1, main_rect.y1, main_rect.x2, main_rect.y2), main_rect.polarity);
            
            //-------------------------------------------------------------------------------------------------
            /*
             cout<<ls_count<<'\t'<<main_rect.theta<<'\t'<<main_rect.theta*180/M_PI<<"\t polarity:"<<main_rect.polarity<<endl;//打印theta
             
             fstream file1,file2;
             if(ls_count == 1)//清空内容
             {
             file1.open("D:\\Graduate Design\\picture\\sp\\coor.txt",ios::out | ios::trunc);
             file1.close();
             file2.open("D:\\Graduate Design\\picture\\sp\\reg.txt",ios::out | ios::trunc);
             file2.close();
             }
             
             file1.open("D:\\Graduate Design\\picture\\sp\\coor.txt",ios::app);
             file1<<main_rect.x1<<'\t'<<main_rect.y1<<'\t'<<main_rect.x2<<'\t'<<main_rect.y2<<'\t'<<(main_rect.theta*180/M_PI)<<endl;
             file1.close();
             
             if(ls_count == 1)//保持第1根线段的区域
             {
             file2.open("D:\\Graduate Design\\picture\\sp\\reg.txt",ios::app);
             for(i=0; i<reg_size; i++)
             file2<<angles->data[ reg[i].x + reg[i].y * angles->xsize ]*180/M_PI<<endl;
             file2.close();
             }
             */
            //-------------------------------------------------------------------------------------------------------
            /* add region number to 'region' image if needed */ //将检测到的线段序号标到相应的图像格子里，该部分可有可无
            if( region != NULL )
                for(i=0; i<reg_size; i++)
                    region->data[ reg[i].x + reg[i].y * region->xsize ] = ls_count;
        }
    
    
    /* free memory */
    free( (void *) image );   /* only the double_image structure should be freed,
                               the data point2ier was provided to this functions
                               and should not be destroyed.                 */
    free_image_double(angles);
    free_image_double(modgrad);
    free_image_char(used);
    free_image_char(pol);
    free( (void *) reg );
    //  free( (void *) mem_p );
    //释放分成1024区的存储梯度从大到小的链表,mycode
    //---------------------------------------
    list_p_temp = list_p->next;
    while(list_p_temp != NULL)
    {
        free(list_p);
        list_p = list_p_temp;
        list_p_temp = list_p->next;
    }
    free(list_p);
    
    //cout<<"seed cnt:"<<seed_cnt<<endl;
    //cout<<"refine cnt:"<<refine_cnt<<endl;
    //cout<<"reg_size_toosmall cnt:"<<reg_size_toosmall_cnt<<endl;
    //----------------------------------------
    /* return the result */
    if( reg_img != NULL && reg_x != NULL && reg_y != NULL )
    {
        if( region == NULL ) error("'region' should be a valid image.");
        *reg_img = region->data;
        if( region->xsize > (unsigned int) INT_MAX ||
           region->xsize > (unsigned int) INT_MAX )
            error("region image to big to fit in INT sizes.");
        *reg_x = (int) (region->xsize);
        *reg_y = (int) (region->ysize);
        
        /* free the 'region' structure.
         we cannot use the function 'free_image_int' because we need to keep
         the memory with the image data to be returned by this function. */
        free( (void *) region );
    }
    if( out->size > (unsigned int) INT_MAX )
        error("too many detections to fit in an INT.");
    *n_out = (int) (out->size);
    
    return_value = out->values;
    free( (void *) out );  /* only the 'ntuple_list' structure must be freed,
                            but the 'values' point2ier must be keep to return
                            as a result. */
    return return_value;
}
