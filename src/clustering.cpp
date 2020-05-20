//
//  clustering.cpp
//  lu
//
//  Created by xiao.hu on 2020/5/20.
//  Copyright © 2020 xiao.hu. All rights reserved.
// done

#include "clustering.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <vector>
#include "utils.hpp"

//==============================================================================
//====================================================================================================
//================================clustering==========================================================
//聚类
//求points中第i行与initializations中第j行里每个元素的平方差总和,每行维度都为nDims
// done, add xiaohu
inline double squaredDifference(int & nDims, double *& points, int & i, double *& initializations, int & j)
{
    double result = 0;
    for (int k = 0; k < nDims; ++k)
        result += pow(points[i*nDims+k] - initializations[j*nDims+k], 2);
    return result;
}
/**
 *输入
 *points: 待均值漂移的点集，总共有nPoints个点，每个点有nDims维度，是一维数组
 *initPoints: 均值漂移初始化位置，在nxd空间中找均值漂移初始时开始搜索的位置，总共有initLength个点，每个点有nDims维度
 *sigma = 1
 *window_size: window parameter = distance_tolerance或者window parameter = distance_tolerance/2
 *accuracy_tolerance: 收敛容忍误差1e-6
 *iter_times: 迭代次数50
 *输出
 *收敛的位置，位置个数与初始化搜索位置个数一样,我们将结果更新到initPoints,也就是它既是输入参数，也是输出参数，节省内存
 done, add xiaohu
 */
void meanShift( double * points, int nPoints, int nDims, double * initPoints, int initLength, double sigma, double window_size, double accuracy_tolerance, int iter_times )
{
    //  for (int i = 0; i<initLength; i++)
    //    cout<<initPoints[2*i]<<'\t'<<initPoints[2*i+1]<<endl;
    int nQuerries = initLength;
    double * initializations = (double*)malloc(nQuerries * nDims * sizeof(double));
    memcpy(initializations, initPoints , nQuerries * nDims * sizeof(double));//copy
    
    double sigma2 = sigma*sigma;//sigma平方
    double radius2 = window_size *window_size;//平方
    double tolerance = accuracy_tolerance;
    int iters, maxiters = iter_times;//最大迭代次数
    //返回与初始搜索点集一样大小的最终定位点集
    double * finals = (double*)malloc(nQuerries * nDims * sizeof(double));;//最终定位点集的指针
    memcpy(finals, initializations, nQuerries * nDims * sizeof(double));
    double * distances = (double*)malloc(nPoints*sizeof(double));
    //printf("meanShift: nPoints:%d \tnDims: %d \tnQuerries:%d \n",nPoints,nDims,nQuerries);//打印
    for (int loop = 0; loop < nQuerries; ++loop)
    {
        iters = 0;
        while (iters < maxiters)
        {
            bool flag = false;
            double denominator = 0;//分母
            for (int i = 0; i < nPoints; ++i)//对所有的点集进行遍历，找到落在搜索圆域内的点
            {
                distances[i] = squaredDifference(nDims, points, i, initializations, loop);//求距离的平方
                if (distances[i] <= radius2)//在第loop个搜索中心的以sqrt(radius2)为半径的圆域内
                {
                    flag = true;
                    denominator += exp(-distances[i] / sigma2);
                }
            }
            if (!flag)
                break;
            for (int j = 0; j < nDims; ++j)
                finals[loop*nDims+j] = 0;//对最终定位点集中的第loop个点的向量赋值为0
            for (int i = 0; i < nPoints; ++i)
                if (distances[i] <= radius2)
                {
                    for (int j = 0; j < nDims; ++j)//每个内点向量的以一定权值累加
                        finals[loop*nDims+j] += exp(-distances[i] / sigma2) * points[i*nDims+j];
                }
            for (int j = 0; j < nDims; ++j)//权值归一化
                finals[loop*nDims+j] /= denominator;
            if (sqrt(squaredDifference(nDims, finals, loop, initializations, loop)) < tolerance)//相继两次的迭代中心在误差内了，则认为已经收敛，没必要再继续迭代
                break;
            iters = iters + 1;
            for (int j = 0; j < nDims; ++j)//更新迭代的搜索中心
                initializations[loop*nDims+j] = finals[loop*nDims+j];
        }
    }
    memcpy(initPoints, finals, nQuerries * nDims * sizeof(double));
    free(distances);
    free(initializations);
    free(finals);
}

/***
 *输入
 *points,待聚类的点集,为一维数组,nPoints个点，每个点维度是nDims
 *distance_threshold 决定聚类的距离阈值
 *输出 outPoints
 *聚类后的点集 nOutPoints x nDims
 *该函数要千万注意，当被调用后，函数内部会多申请nOutPoints个double型的数组内存，在外边使用完毕后，切记free(outPoints).
 done, add xiaohu
 */
void clusterByDistance(double * points, int nPoints, int nDims, double distance_threshold,int number_control, double * & outPoints, int * nOutPoints)
{
    double threshold2 = distance_threshold*distance_threshold;
    std::vector<double*> centers;
    std::vector<int> counts;
    centers.clear();
    counts.clear();
    char * labeled = (char*)malloc(sizeof(char)*nPoints);
    memset(labeled, 0, nPoints * sizeof(char));//初始化bool型标签为0
    if(nPoints == 1)
    {
        centers.push_back((double*)malloc(sizeof(double)*nDims));
        for (int k = 0; k < nDims; ++k)
            centers[centers.size() - 1][k] = points[k];
        counts.push_back(1);
    }
    else
    {
        for (int i = 0; i < nPoints; ++i)
        {
            if (!labeled[i])
            {
                labeled[i] = 1;
                centers.push_back((double*)malloc(sizeof(double)*nDims));
                counts.push_back(1);
                for (int k = 0; k < nDims; ++k)
                {
                    centers[centers.size() - 1][k] = points[i*nDims+k];
                }
                for (int j = i+1; j < nPoints; ++j)
                {
                    if (!labeled[j])
                    {
                        double d = 0;
                        for (int k = 0; k < nDims; ++k)
                            d += pow(centers[centers.size() - 1][k] / counts[centers.size() - 1] - points[j*nDims+k], 2);
                        if (d <= threshold2)
                        {
                            ++counts[centers.size() - 1];
                            for (int k = 0; k < nDims; ++k)
                                centers[centers.size() - 1][k] += points[j*nDims+k];
                            labeled[j] = 1;
                            if(counts[centers.size() - 1] >= number_control)//聚类数量控制，防止均值中心漂的太远  圆心聚类时20  半径聚类时10
                                break;
                        }
                    }
                }
            }
        }
    }
    free(labeled);
    centers.shrink_to_fit();
    counts.shrink_to_fit();
    int m = (int) centers.size();
    outPoints = (double*)malloc(sizeof(double)*m*nDims);
    (*nOutPoints) = m;
    for (unsigned int i = 0; i < centers.size(); ++i)
    {
        for (int j = 0; j < nDims; ++j)
        {
            outPoints[i*nDims+j] = centers[i][j] / counts[i];
            //      cout<<out[i*nDims+j]<<'\t';
        }
        //    cout<<endl;
        free(centers[i]);
    }
    centers.resize(0);
    counts.resize(0);
    //  vector<double*>().swap(centers);//释放回收vector内存
    //  vector<int>().swap(counts);
}

//聚类算法，均值漂移
//三个步骤，一是选取初始迭代点，二是均值漂移，三是去除重复点，从而得到聚类中心
//获得候选圆心的聚类中心(xi,yi)
//输入：
//points，一维点数据,长度为points_num x 2
//distance_tolerance,数据点聚类的半径
//输出：
//二维数据点的聚类中心 centers是一维double数组， 大小为 centers_num x 2
//正确返回值为1，出现错误为0. 例如points为空
//切记切记！！！ centers为函数内部申请的内存，用来返回centers_num个点的聚类中心，使用完后一定要释放，记住free(centers)！！！
// done
int  cluster2DPoints( double * points, int points_num, double distance_tolerance, double * & centers, int * centers_num)
{
    double xmax,xmin,ymax,ymin,xdelta,ydelta;
    int nbins_x,nbins_y;
    int x,y;
    int i;
    unsigned int addr,addr2;
    xmax = ymax = 0;
    xmin = ymin = DBL_MAX;
    for( i = 0; i< points_num; i++ )
    {
        addr = 2*i;
        if( points[addr] > xmax)
            xmax = points[addr];
        if( points[addr] < xmin)
            xmin = points[addr];
        if( points[addr+1] > ymax)
            ymax = points[addr+1];
        if( points[addr+1] < ymin)
            ymin = points[addr+1];
    }
    xmax += xmax*0.02;//避免xdelta、ydelta为0
    xmin -= xmin*0.02;
    ymax += ymax*0.02;
    ymin -= ymin*0.02;
    xdelta = (xmax-xmin);
    ydelta = (ymax-ymin);//有问题，假设所有数据一样大，此处为0
    nbins_x = (int)ceil(xdelta/distance_tolerance);
    nbins_y = (int)ceil(ydelta/distance_tolerance);
    if(nbins_x <= 0 )
    {
        nbins_x = 1;//至少保留1个bin
        //error("generateCircleCandidates,nbins_x,nbins_y error");
    }
    if(nbins_y <= 0)
    {
        nbins_y = 1;
    }
    point2d1i * center_bins;
    center_bins = (point2d1i *)calloc(nbins_y*nbins_x, sizeof(point2d1i));//(x,y,z),x用来记sum(xi),y用来记sum(yi),z用来记落在格子里的数量
    memset(center_bins,0,sizeof(point2d1i)*nbins_y*nbins_x);
    if(center_bins == NULL)
        error("cluster2DPoints, not enough memory");
    //  cout<<"2D原始数据:"<<points_num<<endl;
    for ( i = 0; i< points_num; i++ )//将圆心记录到格子里面，同时落在相应格子里面的数量++
    {
        addr = 2*i;
        
        //    cout<<points[addr]<<'\t'<<points[addr+1]<<endl;
        
        x = (int)((points[addr]   - xmin)/xdelta*nbins_x+0.5);//四舍五入
        y = (int)((points[addr+1] - ymin)/ydelta*nbins_y+0.5);
        if( x >= nbins_x)
            x = nbins_x-1;
        if( y >= nbins_y)
            y = nbins_y-1;
        addr2 = y*nbins_x+x;
        center_bins[addr2].x += points[addr];
        center_bins[addr2].y += points[addr+1];
        center_bins[addr2].z ++;
    }
    int initCentersLength = 0;
    
    //将vote后非0的格子里面的point取均值，并按照顺序重写到center_bins里面，无内存消耗
    for ( y = 0; y<nbins_y;y++)
        for ( x = 0; x<nbins_x; x++)
        {
            addr = y*nbins_x+x;
            if(center_bins[addr].z > 0)
            {
                center_bins[initCentersLength].x = center_bins[addr].x/center_bins[addr].z;
                center_bins[initCentersLength].y = center_bins[addr].y/center_bins[addr].z;
                initCentersLength++;
            }
        }
    if(initCentersLength == 0)
    {
        (*centers_num) = 0;
        centers = NULL;
        //cout<<"cluster2DPoints,points number:"<<points_num<<endl;
        //cout<<"cluster2DPoints,initCentersLength equals 0"<<endl;
        return 0;
        //error("generateCircleCandidates,initCentersLength equals 0");
    }
    double * initCenters; //initCentersLength x 2
    initCenters = (double*)malloc(sizeof(double)*initCentersLength*2);
    //将记录在链表里面的分区后的圆心均值记录到数组里，便于作为初始点进行均值漂移
    for ( i = 0; i<initCentersLength; i++ )// initCenters 大小是 initCentersLength*2
    {
        int addr = 2*i;
        initCenters[addr]   = center_bins[i].x;
        initCenters[addr+1] = center_bins[i].y;
    }
    free((void*)center_bins);//赶紧释放该内存
    
    //  cout<<"2D均值漂移前初始迭代点："<<endl;
    //  for (int  i = 0; i<initCentersLength; i++)
    //    cout<<initCenters[2*i]<<'\t'<<initCenters[2*i+1]<<endl;
    
    //均值漂移的结果会更新到initCenters里面
    meanShift(points,points_num,\
              2,initCenters,initCentersLength,1,distance_tolerance,1e-6,50);//迭代20次
    
    //  cout<<"2D均值漂移后的聚类中心:"<<endl;
    //  for (int  i = 0; i<initCentersLength; i++)
    //    cout<<initCenters[2*i]<<'\t'<<initCenters[2*i+1]<<endl;
    
    //聚类
    //千万要注意centers_num是int型指针，++--时要(*centers_num).
    clusterByDistance(initCenters,initCentersLength,2,distance_tolerance/2,40,centers, centers_num);//此处控制参数要改，要调节
    
    //  cout<<"2D距离聚类，去除重复点后的点集:"<<endl;
    //  for (int  i = 0; i<(*centers_num); i++)
    //    cout<<centers[2*i]<<'\t'<<centers[2*i+1]<<endl;
    
    if((*centers_num) <= 0)//可无
    {
        return 0;  //不懂为什么，聚类中心的周围确没有最靠近它的点
        //system("pause");
        //error("cluster2DPoints,(*centers_num)<=0");
    }
    free(initCenters);
    //cout<<"2D聚类后数量:"<<(*centers_num)<<endl;
    return 1;
}

//聚类算法，均值漂移
//三个步骤，一是选取初始迭代点，二是均值漂移，三是去除重复点，从而得到聚类中心
//获得候选圆心的聚类中心(xi,yi)
//输入：
//datas，一维点数据,长度为datas_num x 1
//distance_tolerance,数据点聚类的半径
//输出：
//一维数据点的聚类中心 centers是一维double数组， 大小为 centers_num x 1
//正确返回值为1，出现错误为0. 例如points为空
//切记切记！！！ centers为函数内部申请的内存，用来返回centers_num个点的聚类中心，使用完后一定要释放，记住free(centers)！！！
// done
int  cluster1DDatas( double * datas, int datas_num, double distance_tolerance, double * & centers, int * centers_num)
{
    double rmin,rmax,rdelta;
    int r;
    int i;
    rmin = DBL_MAX;
    rmax = 0;
    for( i  = 0; i < datas_num; i++)//将链表里的r集合复制到数组
    {
        if(datas[i] < rmin)//在这一次遍历中，记录最大最小值
            rmin = datas[i];
        if(datas[i] > rmax)
            rmax = datas[i];
    }
    int nbins_r = 0;
    point1d1i * center_bins;
    rmax += rmin*0.02;//避免rmax-rmin = 0
    rmin -= rmin*0.02;
    rdelta = rmax - rmin;
    nbins_r = (int)ceil((rdelta)/distance_tolerance);
    if(nbins_r <= 0)//至少有一个bin
        nbins_r = 1;
    center_bins = (point1d1i *)malloc(sizeof(point1d1i)*nbins_r);
    memset(center_bins,0,sizeof(point1d1i)*nbins_r);//初始化为0
    //  cout<<"1D原始数据:"<<datas_num<<endl;
    for( i = 0; i<datas_num; i++)//对分区vote
    {
        //    cout<<datas[i]<<endl;
        r = int((datas[i]-rmin)/rdelta*nbins_r+0.5);
        if(r>=nbins_r)
            r = nbins_r-1;
        center_bins[r].data += datas[i];
        center_bins[r].cnt  ++;
    }
    int init_r_length = 0;
    for( i = 0; i<nbins_r; i++)
    {
        if(center_bins[i].cnt > 0)//统计非0分区,并且对每一个bin取均值，按照顺序重写到center_bins里面，无内存消耗
        {
            center_bins[init_r_length].data = center_bins[i].data/center_bins[i].cnt;
            init_r_length++;
        }
    }
    if(init_r_length == 0)
    {
        (*centers_num) = 0;
        centers = NULL;
        //cout<<"cluster1DDatas,points number:"<<datas_num<<endl;
        //cout<<"cluster2DDatas,init_r_length equals 0"<<endl;
        return 0;
        //error("generateCircleCandidates,initCentersLength equals 0");
    }
    double * initCenters; //init_r_length x 1
    initCenters = (double*)malloc(sizeof(double)*init_r_length);
    //将记录在链表里面的分区后的圆心均值记录到数组里，便于作为初始点进行均值漂移
    for ( i = 0; i<init_r_length; i++ )// initCenters 大小是 init_r_length x 1
    {
        initCenters[i] = center_bins[i].data;
    }
    free(center_bins);//赶紧释放该内存
    
    //  cout<<"1D均值漂移前初始迭代点："<<endl;
    //  for (int  i = 0; i<init_r_length; i++)
    //    cout<<initCenters[i]<<'\t';
    //  cout<<endl;
    
    //至此，得到了均值漂移初始的initCenters，为一维double数组，长度是init_r_length
    meanShift(datas, datas_num, 1, initCenters, init_r_length, 1, distance_tolerance, 1e-6, 20);//迭代20次
    
    //  cout<<"1D均值漂移后的聚类中心:"<<endl;
    //  for (int  i = 0; i<init_r_length; i++)
    //    cout<<initCenters[i]<<'\t';
    //  cout<<endl;
    
    //聚类
    //千万要注意centers_num是int型指针，++--时要(*centers_num).
    clusterByDistance(initCenters, init_r_length, 1, distance_tolerance/2, 40, centers, centers_num);//控制参数40，最多40个点合成1个点
    
    //  cout<<"1D距离聚类，去除重复点后的点集:"<<endl;
    //  for (int  i = 0; i<(*centers_num); i++)
    //    cout<<centers[i]<<'\t';
    //  cout<<endl;
    
    if((*centers_num) <= 0)//可无
    {
        return 0;  //不懂为什么，聚类中心的周围确没有最靠近它的点
        //system("pause");
        //error("cluster1DDatas,(*centers_num)<=0");
    }
    free(initCenters);
    //  cout<<"1D聚类后数量::"<<(*centers_num)<<endl;
    return 1;
}





/*
 void drawEllipse(Mat img, double * ellipara)
 {
 Point peliicenter(ellipara[0],ellipara[1]);
 Size  saxis(ellipara[2],ellipara[3]);
 //Mat ellimat = Mat::zeros(img.rows,img.cols,CV_8UC3);
 //ellimat.setTo(255);
 static int ccc = 0;
 static unsigned int cnt = 0;
 if(cnt % 2 == 0 )
 ccc = 0;
 else
 {
 ccc = 255;
 cout<<cnt/2<<'\t'<<ellipara[0]<<'\t'<<ellipara[1]<<"\t"<<ellipara[2]<<'\t'<<ellipara[3]<<'\t'<<ellipara[4]<<endl;
 }
 cnt++;
 
 Mat imgtemp = img.clone();
 ellipse(imgtemp,peliicenter,saxis,ellipara[4]*180/M_PI,0,360,(Scalar(0,255,ccc)),2);
 namedWindow("w1");
 imshow("w1",imgtemp);
 //waitKey(0);
 }
 void drawEdge(Mat img, point2d * dataxy, int num)
 {
 static int ccc = 0;
 static int cnt = 0;
 cnt++;
 if(cnt % 2 == 0 )
 ccc = 0;
 else
 ccc = 255;
 Mat imgtemp = img.clone();
 for (int i = 0; i<num; i++)
 {
 imgtemp.at<Vec3b>(dataxy[i].y,dataxy[i].x) = (Vec3b(ccc,255,0));
 }
 namedWindow("w2");
 imshow("w2",imgtemp);
 }
 */








void generateEllipseCandidates( PairGroupList * pairGroupList, double distance_tolerance, double * & ellipse_candidates, int * candidates_num)
{
    if( pairGroupList->length <= 0 )//检测，至少要有1个样本用来产生候选
    {
        ellipse_candidates = NULL;
        (*candidates_num) = 0;
        return;
    }
    double * centers;
    int center_num; //椭圆中心(xi,yi)的聚类数量
    double * phis;
    int phi_num;    //针对每一个椭圆中心(xi,yi)，倾斜角度phi的聚类数量
    double * axises;
    int axis_num;   //针对每一个椭圆中心和倾角(xi,yi,phi),长短半轴(a,b)的聚类数量
    double * bufferXY = (double*)calloc(pairGroupList->length*2,sizeof(double));
    double * bufferPhi = (double*)calloc(pairGroupList->length,sizeof(double));
    double * bufferAB = (double*)calloc(pairGroupList->length*2,sizeof(double));
    point2i * bufferIndexes = (point2i *)calloc(pairGroupList->length,sizeof(point2i));//point[i].x记录第i个分类在bufferXX中的起始索引位置，point[i].y记录第i个分类在bufferXX中的长度
    double  * buffer2AB = (double*)calloc(pairGroupList->length*2,sizeof(double));
    point2i * buffer2Indexes = (point2i *)calloc(pairGroupList->length,sizeof(point2i));//point[i].x记录第i个分类在bufferXX中的起始索引位置，point[i].y记录第i个分类在bufferXX中的长度
    int     * buffer_temp = (int*)calloc(pairGroupList->length,sizeof(int));
    int addr,addr2,info,ind;
    double dis_min,dis_temp;
    if ( bufferXY == NULL || bufferPhi == NULL || bufferAB == NULL || bufferIndexes == NULL ||
        buffer2AB == NULL || buffer2Indexes == NULL || buffer_temp == NULL
        )
    {
        ellipse_candidates = NULL;
        (*candidates_num) = 0;
        error("generateEllipseCandidates, not enough memory");
    }
    (*candidates_num) = 0; //候选椭圆数量，初始化为0,非常重要
    //copy
    for ( int i = 0; i<pairGroupList->length; i++)
    {
        addr = 2*i;
        bufferXY[addr] = pairGroupList->pairGroup[i].center.x;
        bufferXY[addr+1] = pairGroupList->pairGroup[i].center.y;
    }
    //cluster the ellipses' centers
    info = cluster2DPoints(bufferXY,pairGroupList->length,distance_tolerance,centers,&center_num);
    if( info == 0)
    {
        ellipse_candidates = NULL;
        (*candidates_num) = 0;
        error("generateEllipseCandidates, cluster2DPoints, error in clustering elliptic centers");
    }
    //classification,寻找每个点归属的聚类中心
    for ( int i = 0; i<pairGroupList->length; i++)
    {
        dis_min = DBL_MAX;
        ind = -1;
        for ( int j = 0; j<center_num; j++)
        {
            addr = 2*j;
            dis_temp = (pairGroupList->pairGroup[i].center.x - centers[addr])*(pairGroupList->pairGroup[i].center.x - centers[addr]) + (pairGroupList->pairGroup[i].center.y - centers[addr+1])*(pairGroupList->pairGroup[i].center.y - centers[addr+1]);
            if(dis_temp < dis_min)
            {
                dis_min = dis_temp;
                ind = j; //record the nearest center's index
            }
        }
        buffer_temp[i] = ind; //此处借用buffer2来记下第i个初始椭圆对应第ind个椭圆聚类中心
    }
    //将分类结果按顺序存到bufferXY,bufferPhi,bufferAB中，且bufferIndexes[i]存着第i个聚类中心的起始索引位置和长度
    memset(bufferIndexes,0,sizeof(point2i)*pairGroupList->length);
    ind = 0;//清零，样本点起始位置，索引位置是ind*2,分区的基址
    for ( int i = 0; i<center_num; i++)
    {
        bufferIndexes[i].x = ind;
        for ( int j = 0; j<pairGroupList->length; j++)
        {
            if ( buffer_temp[j] == i)
            {
                addr = ind*2;//切记长短半轴是一组一组寸储的，需要 x 2
                addr2 = bufferIndexes[i].y*2;
                bufferPhi[ind+bufferIndexes[i].y] = pairGroupList->pairGroup[j].phi;
                bufferAB[addr+addr2] = pairGroupList->pairGroup[j].axis.x;
                bufferAB[addr+addr2+1] = pairGroupList->pairGroup[j].axis.y;
                bufferIndexes[i].y++;//第i个聚类中心周围的点数量加1
            }
        }
        if(bufferIndexes[i].y == 0)//聚类中心周围没有靠近的点
        {
            error("generateEllipseCandidates, no XY points near to the clustering center");
        }
        ind += bufferIndexes[i].y;
    }
    //cout<<"2D cluster centers over"<<endl;
    //对每一个椭圆中心的周围的点进行倾角聚类
    //第i个椭圆聚类中心，其邻近点的索引范围是：bufferIndexs[i].x ~ (bufferIndex[i].x + bufferIndex[i].y-1)
    for ( int i = 0; i<center_num; i++)
    {
        double * phi_pointer_temp = bufferPhi+bufferIndexes[i].x;//倾角指针
        double * ab_pointer_temp = bufferAB+bufferIndexes[i].x*2;//长短半轴的指针,记住 x 2
        info = cluster1DDatas(phi_pointer_temp, bufferIndexes[i].y, 0.0873, phis, &phi_num);//对phi聚类, pi/180*5 = 0.0873, 5°误差
        if (info == 0) //不懂为什么，聚类中心centers[i]的周围可能没有最靠近它的点,数量bufferIndexes[i].y = 0
        {
            //cout<<"generateEllipseCandidates, cluster2DPoints, error in clustering elliptic phis"<<endl;
            continue;
            //error("generateEllipseCandidates, cluster2DPoints, error in clustering elliptic phis");
        }
        //classification,寻找每个点归属的聚类中心
        for ( int j = 0; j<bufferIndexes[i].y; j++ )
        {
            dis_min = DBL_MAX;
            ind = -1;
            for ( int k = 0; k<phi_num; k++)
            {
                dis_temp = (*(phi_pointer_temp+j)-phis[k]) * (*(phi_pointer_temp+j)-phis[k]);
                if(dis_temp < dis_min)
                {
                    dis_min = dis_temp;
                    ind = k;//record the nearest phi's index
                }
            }
            buffer_temp[j] = ind;
        }
        //将分类结果按顺序存储到buffer2AB中，且buffer2Indexes[j].x对应第i个phi的聚类中心起始点，buffer2Indexes[j].y对应数量(长度)
        memset(buffer2Indexes,0,sizeof(point2i)*bufferIndexes[i].y);
        ind = 0;
        for ( int j = 0; j<phi_num; j++)
        {
            buffer2Indexes[j].x = ind;//起始点
            for ( int k = 0; k<bufferIndexes[i].y; k++)
            {
                if ( buffer_temp[k] == j)
                {
                    addr = ind*2;
                    addr2 = buffer2Indexes[j].y*2;
                    buffer2AB[addr+addr2] = *(ab_pointer_temp+k*2);
                    buffer2AB[addr+addr2+1] = *(ab_pointer_temp+k*2+1);
                    buffer2Indexes[j].y++;//长度加1
                }
            }
            ind += buffer2Indexes[j].y;
        }
        for ( int j = 0; j<phi_num; j++ )
        {
            double * ab_pointer_temp2 = buffer2AB+buffer2Indexes[j].x*2; //长短半轴的指针,记住 x 2
            info = cluster2DPoints(ab_pointer_temp2, buffer2Indexes[j].y, distance_tolerance, axises, &axis_num);
            if (info == 0) //不懂为什么，聚类中心phi_j的周围可能没有最靠近它的点,数量buffer2Indexes[j].y = 0
            {
                //cout<<"generateEllipseCandidates, cluster2DPoints, error in clustering elliptic axises"<<endl;
                continue;
                //error("generateEllipseCandidates, cluster2DPoints, error in clustering elliptic axises");
            }
            //将候选椭圆重写到bufferXY,bufferPhi,bufferAB里面, 候选椭圆数量(*candidates_num)++
            for ( int k = 0; k<axis_num; k++)
            {
                addr = (*candidates_num)*2;
                bufferXY[addr] = centers[i*2];
                bufferXY[addr+1] = centers[i*2+1];
                bufferPhi[(*candidates_num)] = phis[j];
                bufferAB[addr] = axises[k*2];
                bufferAB[addr+1] = axises[k*2+1];
                (*candidates_num)++;
            }
            free(axises);//cluster2DPoints严格要求，用完axises后，需要释放函数内部申请的内存
        }
        free(phis);//cluster1DDatas严格要求，用完phis后，需要释放函数内部申请的内存
    }
    free(centers);//cluster2DPoints严格要求，用完centers后，需要释放函数内部申请的内存
    //释放在函数开头申请的部分内存
    free(buffer_temp); //此处释放出问题
    free(buffer2Indexes);
    free(buffer2AB);
    free(bufferIndexes);
    ellipse_candidates = (double*)malloc(sizeof(double)*(*candidates_num)*5);
    for ( int i = 0; i < (*candidates_num); i++ )
    {
        addr = 2*i;
        ellipse_candidates[i*5]  = bufferXY[addr];
        ellipse_candidates[i*5+1]= bufferXY[addr+1];
        ellipse_candidates[i*5+2]= bufferAB[addr];
        ellipse_candidates[i*5+3]= bufferAB[addr+1];
        ellipse_candidates[i*5+4]= bufferPhi[i];
    }
    //释放在函数开头申请的内存
    free(bufferAB);
    free(bufferPhi);
    free(bufferXY);
    if((*candidates_num)<= 0)
    {
        *candidates_num = 0;
        ellipse_candidates = NULL;
        //cout<<"no any candidates generated!"<<endl;
    }
}

