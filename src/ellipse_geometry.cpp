//
//  ellipse_geometry.cpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
// done, ellipse fitting

#include "ellipse_geometry.hpp"
#include "utils.hpp"
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include "datatype.h"
#include <float.h>
#include <iostream>
#include <stdio.h>

#ifdef MEX_COMPILE
    #include "lapack.h"  //matlab
#else
    #ifdef __APPLE__
//#include <Accelerate/Accelerate.h>
        #include <clapack.h>
    #else
        #include "lapack.h"  // make sure you have lapack
    #endif
#endif
//#include "mex.h"
//#endif

/*----------------------------------------------------------------------------*/
/** Approximate the distance between a point and an ellipse using Rosin distance.
 */
inline double d_rosin (double *param, double x, double y)
{
    double ae2 = param[2]*param[2];
    double be2 = param[3]*param[3];
    x = x - param[0];
    y = y - param[1];
    double xp = x*cos(-param[4])-y*sin(-param[4]);
    double yp = x*sin(-param[4])+y*cos(-param[4]);
    double fe2;
    fe2 = ae2-be2;
    double X = xp*xp;
    double Y = yp*yp;
    double delta = (X+Y+fe2)*(X+Y+fe2)-4*X*fe2;
    double A = (X + Y + fe2 - sqrt(delta))/2.0;
    double ah = sqrt(A);
    double bh2 = fe2-A;
    double term = (A*be2+ae2*bh2);
    double xi = ah*sqrt(ae2*(be2+bh2)/term);
    double yi = param[3]*sqrt(bh2*(ae2-A)/term);
    double d[4],dmin;
    
    
    d[0] = dist(xp,yp,xi,yi);
    d[1] = dist(xp,yp,xi,-yi);
    d[2] = dist(xp,yp,-xi,yi);
    d[3] = dist(xp,yp,-xi,-yi);
    dmin = DBL_MAX;
    for ( int i = 0; i<4; i++)
    {
        if( d[i] <= dmin)
            dmin = d[i];
    }
    //  if (X+Y>xi*xi+yi*yi)
    //    return dmin;
    //  else return -dmin;
    return dmin;
}
/*----------------------------------------------------------------------------*/


//=============================================================================
/** Convert ellipse from matrix form to common form:
 ellipse = (centrex,centrey,ax,ay,orientation).
 */
int ellipse2Param(double *p,double param[])
{
    // ax^2 + bxy + cy^2 + dx + ey + f = 0
    double a,b,c,d,e,f;
    double thetarad,cost,sint,cos_squared,sin_squared,cos_sin,Ao,Au,Av,Auu,Avv,tuCentre,tvCentre,wCentre,uCentre,vCentre,Ru,Rv;
    a = p[0];
    b = p[1];
    c = p[2];
    d = p[3];
    e = p[4];
    f = p[5];
    
    thetarad=0.5*atan2(b,a-c);
    cost=cos(thetarad);
    sint=sin(thetarad);
    sin_squared=sint*sint;
    cos_squared=cost*cost;
    cos_sin=sint*cost;
    Ao=f;
    Au=d*cost+e* sint;
    Av=-d*sint+e* cost;
    Auu=a*cos_squared+c*sin_squared+b*cos_sin;
    Avv=a*sin_squared+c*cos_squared-b*cos_sin;
    
    if(Auu==0 || Avv==0){ param[0]=0;param[1]=0;param[2]=0;param[3]=0;param[4]=0;return 0;}
    else
    {
        tuCentre=-Au/(2.*Auu);
        tvCentre=-Av/(2.*Avv);
        wCentre=Ao-Auu*tuCentre*tuCentre-Avv*tvCentre*tvCentre;
        uCentre=tuCentre*cost-tvCentre*sint;
        vCentre=tuCentre*sint+tvCentre*cost;
        Ru=-wCentre/Auu;
        Rv=-wCentre/Avv;
        //     if (Ru>0) Ru=pow(Ru,0.5);
        //     else Ru=-pow(-Ru,0.5);
        //     if (Rv>0) Rv=pow(Rv,0.5);
        //     else Rv=-pow(-Rv,0.5);
        if (Ru <= 0 || Rv <= 0)//长短轴小于0的情况？？？
            return 0;
        Ru = sqrt(Ru);
        Rv = sqrt(Rv);
        param[0]=uCentre;param[1]=vCentre;
        param[2]=Ru;param[3]=Rv;param[4]=thetarad;
        //会出现Ru < Rv情况，对调一下
        if(Ru < Rv )
        {
            param[2] = Rv;
            param[3] = Ru;
            if(thetarad < 0)//调换长短轴，使得第三个参数为长轴，第四个为短轴
                param[4] += M_1_2_PI;
            else
                param[4] -= M_1_2_PI;
            if(thetarad < - M_1_2_PI)//长轴倾角限定在-pi/2 ~ pi/2，具备唯一性
                param[4] += M_PI;
            if(thetarad > M_1_2_PI)
                param[4] -= M_PI;
        }
    }
    return 1;
}
//input : (xi,yi)
//output: x0,y0,a,b,phi,ellipara需要事先申请内存
//successfull, return 1; else return 0
int fitEllipse(point2d* dataxy, int datanum, double* ellipara)
{
    double* D = (double*)malloc(datanum*6*sizeof(double));
    double S[36];
    double C[36];
    memset(D,0,sizeof(double)*datanum);
    memset(S,0,sizeof(double)*36);
    memset(C,0,sizeof(double)*36);
    for ( int i = 0; i<datanum; i++)
    {
        D[i*6]  = dataxy[i].x*dataxy[i].x;
        D[i*6+1]= dataxy[i].x*dataxy[i].y;
        D[i*6+2]= dataxy[i].y*dataxy[i].y;
        D[i*6+3]= dataxy[i].x;
        D[i*6+4]= dataxy[i].y;
        D[i*6+5]= 1;
    }
    for ( int i = 0; i<6; i++)
        for ( int j = i; j<6; j++)
        {
            //S[i*6+j]
            for ( int k = 0; k<datanum; k++)
                S[i*6+j] += D[k*6+i]*D[k*6+j];
        }
    free(D);//释放内存
    //对称矩阵赋值
    for ( int i = 0; i<6; i++)
        for ( int j = 0; j<i; j++)
            S[i*6+j]=S[j*6+i];
    C[0*6+2] = 2;
    C[1*6+1] = -1;
    C[2*6+0] = 2;
    // eig(S,C) eig(inv(S)*C)
    double alphar[6],alphai[6],beta[6];
    double vl[36] = {0};//此处不用
    double vr[36] = {0};
    char JOBVL = 'N';
    char JOBVR = 'V';
    double fitWork[64];
    
#ifdef MEX_COMPILE
    ptrdiff_t fitN = 6;
    ptrdiff_t workLen = 64;
    ptrdiff_t info;
    //info = LAPACKE_dggev(LAPACK_ROW_MAJOR,'N','V',6,S,6,C,6,alphar,alphai,beta,vl,6,vr,6);
    //注意S为对称矩阵，故转置后等于本身，变成列优先，S可以不变
    dggev(&JOBVL,&JOBVR,&fitN,S,&fitN,C,&fitN,alphar,\
          alphai,beta,vl,&fitN,vr,&fitN,fitWork,&workLen,&info);
#else
    #ifdef __APPLE__
        __CLPK_integer fitN = 6;
        __CLPK_integer workLen = 64;
        __CLPK_integer info;
        dggev_(&JOBVL, &JOBVR, &fitN, S, &fitN, C, &fitN, alphar, alphai, beta, vl, &fitN, vr, &fitN, fitWork, &workLen, &info);
    #else
        ptrdiff_t fitN = 6;
        ptrdiff_t workLen = 64;
        ptrdiff_t info;
        //info = LAPACKE_dggev(LAPACK_ROW_MAJOR,'N','V',6,S,6,C,6,alphar,alphai,beta,vl,6,vr,6);
        //注意S为对称矩阵，故转置后等于本身，变成列优先，S可以不变
        dggev(&JOBVL,&JOBVR,&fitN,S,&fitN,C,&fitN,alphar,\
              alphai,beta,vl,&fitN,vr,&fitN,fitWork,&workLen,&info);
    #endif
#endif
    
    if(info == 0)
    {
        int index = -1;
        for ( int i = 0; i<6; i++)
            if( (alphar[i]>=-(2.2204460492503131e-014)) && (alphai[i] == 0) && (beta[i] != 0)) // 100*DBL_EPSILON, eigenvalue = (alphar + i*alphai)/beta
                index = i;//vr[:,i],vr第i列对应的特征向量则为拟合参数
        if(index == -1)//再试一次，放宽对实部>0的约束，放宽到>-0.005
        {
            double temp = -0.005;//这个参数很关键
            for ( int i = 0; i<6; i++)
                if( (alphar[i]>=temp) && (alphai[i] == 0) && (beta[i] != 0)) // 100*DBL_EPSILON, eigenvalue = (alphar + i*alphai)/beta
                {
                    temp = alphar[i];
                    index = i;//vr[:,i],vr第i列对应的特征向量则为拟合参数
                }
        }
        if(index != -1)
        {
            //此处借用beta来传递参数
            //beta[0] = vr[6*0+index];
            //beta[1] = vr[6*1+index];
            //beta[2] = vr[6*2+index];
            //beta[3] = vr[6*3+index];
            //beta[4] = vr[6*4+index];
            //beta[5] = vr[6*5+index];
            beta[0] = vr[6*index+0];
            beta[1] = vr[6*index+1];
            beta[2] = vr[6*index+2];
            beta[3] = vr[6*index+3];
            beta[4] = vr[6*index+4];
            beta[5] = vr[6*index+5];
            ellipse2Param(beta,ellipara);//ax^2 + bxy + cy^2 + dx + ey + f = 0, transform to (x0,y0,a,b,phi)
            return 1;
        }
    }
    return 0;
}

//input: dataxy为数据点(xi,yi)d,总共有datanum个
//output: 拟合矩阵S. 注意：S需要事先申请内存，double S[36].
void calcuFitMatrix(point2d* dataxy, int datanum, double * S)
{
    double* D = (double*)malloc(datanum*6*sizeof(double));
    memset(D,0,sizeof(double)*datanum);
    for ( int i = 0; i<datanum; i++)
    {
        D[i*6]  = dataxy[i].x*dataxy[i].x;
        D[i*6+1]= dataxy[i].x*dataxy[i].y;
        D[i*6+2]= dataxy[i].y*dataxy[i].y;
        D[i*6+3]= dataxy[i].x;
        D[i*6+4]= dataxy[i].y;
        D[i*6+5]= 1;
    }
    for ( int i = 0; i<6; i++)
    {
        for ( int j = i; j<6; j++)
        {
            //S[i*6+j]
            for ( int k = 0; k<datanum; k++)
                S[i*6+j] += D[k*6+i]*D[k*6+j];
        }
    }
    free(D);//释放内存
    //对称矩阵赋值
    for ( int i = 0; i<6; i++)
        for ( int j = 0; j<i; j++)
            S[i*6+j]=S[j*6+i];
}
//input: fit matrixes S1,S2. length is 36.
//output: fit matrix S_out. S_out = S1 + S2.
//S_out事先需要申请内存
void addFitMatrix(double * S1, double * S2, double * S_out)
{
    int ind;
    for ( int i = 0; i<6; i++ )
        for ( int j = i; j<6; j++)
        {
            ind = i*6+j;
            S_out[ind] = S1[ind]+S2[ind];
        }
    //对称矩阵赋值
    for ( int i = 0; i<6; i++)
        for ( int j = 0; j<i; j++)
            S_out[i*6+j]=S_out[j*6+i];
}


//input : S矩阵，6 x 6 = 36
//output: (A,B,C,D,E,F)且A>0, ellicoeff需要事先申请内存. 当要转换成(x0,y0,a,b,phi)时，则要用
//ellipse2Param(ellicoeff,ellipara); ax^2 + bxy + cy^2 + dx + ey + f = 0, transform to (x0,y0,a,b,phi)
//successfull, return 1; else return 0
int fitEllipse2(double * S, double* ellicoeff)
{
    double C[36];
    memset(C,0,sizeof(double)*36);
    
    C[0*6+2] = 2;
    C[1*6+1] = -1;
    C[2*6+0] = 2;
    // eig(S,C) eig(inv(S)*C)
    double alphar[6],alphai[6],beta[6];
    double vl[36] = {0};//此处不用
    double vr[36] = {0};
    char JOBVL = 'N';
    char JOBVR = 'V';
    double fitWork[64];
    
#ifdef MEX_COMPILE
    ptrdiff_t fitN = 6;
    ptrdiff_t workLen = 64;
    ptrdiff_t info;
    //info = LAPACKE_dggev(LAPACK_ROW_MAJOR,'N','V',6,S,6,C,6,alphar,alphai,beta,vl,6,vr,6);
    dggev(&JOBVL,&JOBVR,&fitN,S,&fitN,C,&fitN,alphar,alphai,\
          beta,vl,&fitN,vr,&fitN,fitWork,&workLen,&info);
#else
    #ifdef __APPLE__
        __CLPK_integer fitN = 6;
        __CLPK_integer workLen = 64;
        __CLPK_integer info;
        dggev_(&JOBVL, &JOBVR, &fitN, S, &fitN, C, &fitN, alphar, alphai, beta, vl, &fitN, vr, &fitN, fitWork, &workLen, &info);
    #else
        ptrdiff_t fitN = 6;
        ptrdiff_t workLen = 64;
        ptrdiff_t info;
        //info = LAPACKE_dggev(LAPACK_ROW_MAJOR,'N','V',6,S,6,C,6,alphar,alphai,beta,vl,6,vr,6);
        dggev(&JOBVL,&JOBVR,&fitN,S,&fitN,C,&fitN,alphar,alphai,\
               beta,vl,&fitN,vr,&fitN,fitWork,&workLen,&info);
    #endif
#endif
    
    if(info == 0)
    {
        int index = -1;
        for ( int i = 0; i<6; i++)
            if( (alphar[i]>=-(2.2204460492503131e-014)) && (alphai[i] == 0) && (beta[i] != 0)) // 100*DBL_EPSILON, eigenvalue = (alphar + i*alphai)/beta
                index = i;//vr[:,i],vr第i列对应的特征向量则为拟合参数
        if(index == -1)//再试一次，放宽对实部>0的约束，放宽到>-0.005
        {
            double temp = -0.005;//这个参数很关键
            for ( int i = 0; i<6; i++)
                if( (alphar[i]>=temp) && (alphai[i] == 0) && (beta[i] != 0)) // 100*DBL_EPSILON, eigenvalue = (alphar + i*alphai)/beta
                {
                    temp = alphar[i];
                    index = i;//vr[:,i],vr第i列对应的特征向量则为拟合参数
                }
        }
        if(index != -1)
        {
            //此处借用beta来传递参数
            if(vr[6*index+0] < 0)//注意列优先
            {
                ellicoeff[0] = -vr[6*index+0]; //-vr[6*0+index];
                ellicoeff[1] = -vr[6*index+1]; //-vr[6*1+index];
                ellicoeff[2] = -vr[6*index+2]; //-vr[6*2+index];
                ellicoeff[3] = -vr[6*index+3]; //-vr[6*3+index];
                ellicoeff[4] = -vr[6*index+4]; //-vr[6*4+index];
                ellicoeff[5] = -vr[6*index+5]; //-vr[6*5+index];
            }
            else
            {
                ellicoeff[0] = vr[6*index+0];//vr[6*0+index];
                ellicoeff[1] = vr[6*index+1];//vr[6*1+index];
                ellicoeff[2] = vr[6*index+2];//vr[6*2+index];
                ellicoeff[3] = vr[6*index+3];//vr[6*3+index];
                ellicoeff[4] = vr[6*index+4];//vr[6*4+index];
                ellicoeff[5] = vr[6*index+5];//vr[6*5+index];
            }
            return 1;
        }
    }
    return 0;
}
