//
//  tuple.cpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#include "tuple.hpp"

/*----------------------------------------------------------------------------*/
/** Free memory used in n-tuple 'in'.
 */
void free_ntuple_list(ntuple_list in)
{
    if( in == NULL || in->values == NULL )
        error("free_ntuple_list: invalid n-tuple input.");
    free( (void *) in->values );
    free( (void *) in );
}

/*----------------------------------------------------------------------------*/
/** Create an n-tuple list and allocate memory for one element.
 @param dim the dimension (n) of the n-tuple.
 */
ntuple_list new_ntuple_list(int dim)
{
    ntuple_list n_tuple;
    
    /* check parameters */
    if( dim == 0 ) error("new_ntuple_list: 'dim' must be positive.");
    
    /* get memory for list structure */
    n_tuple = (ntuple_list) malloc( sizeof(struct ntuple_list_s) );
    if( n_tuple == NULL ) error("not enough memory.");
    
    /* initialize list */
    n_tuple->size = 0;
    n_tuple->max_size = 1;
    n_tuple->dim = dim;
    
    /* get memory for tuples */
    n_tuple->values = (double *) malloc( dim*n_tuple->max_size * sizeof(double) );
    if( n_tuple->values == NULL ) error("not enough memory.");
    
    return n_tuple;
}

/*----------------------------------------------------------------------------*/
/** Enlarge the allocated memory of an n-tuple list.
 */
void enlarge_ntuple_list(ntuple_list n_tuple)
{
    /* check parameters */
    if( n_tuple == NULL || n_tuple->values == NULL || n_tuple->max_size == 0 )
        error("enlarge_ntuple_list: invalid n-tuple.");
    
    /* duplicate number of tuples */
    n_tuple->max_size *= 2;
    
    /* realloc memory */
    n_tuple->values = (double *) realloc( (void *) n_tuple->values,
                                         n_tuple->dim * n_tuple->max_size * sizeof(double) );
    if( n_tuple->values == NULL ) error("not enough memory.");
}

/*----------------------------------------------------------------------------*/
/** Add a 8-tuple to an n-tuple list.
 */
void add_8tuple( ntuple_list out, double v1, double v2, double v3,
                       double v4, double v5, double v6, double v7, int v8)
{
    /* check parameters */
    if( out == NULL ) error("add_8tuple: invalid n-tuple input.");
    if( out->dim != 8 ) error("add_8tuple: the n-tuple must be a 8-tuple.");
    
    /* if needed, alloc more tuples to 'out' */
    if( out->size == out->max_size ) enlarge_ntuple_list(out);
    if( out->values == NULL ) error("add_8tuple: invalid n-tuple input.");
    
    /* add new 8-tuple */
    out->values[ out->size * out->dim + 0 ] = v1;
    out->values[ out->size * out->dim + 1 ] = v2;
    out->values[ out->size * out->dim + 2 ] = v3;
    out->values[ out->size * out->dim + 3 ] = v4;
    out->values[ out->size * out->dim + 4 ] = v5;
    out->values[ out->size * out->dim + 5 ] = v6;
    out->values[ out->size * out->dim + 6 ] = v7;
    out->values[ out->size * out->dim + 7 ] = v8;
    
    /* update number of tuples counter */
    out->size++;
}
