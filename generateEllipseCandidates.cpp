#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <limits.h>
#include <float.h>
#include <iostream>
#include "lapack.h"  //matlab 
//#include "include/lapacke_config.h"  //lapack手动，未成功
//#include "include/lapacke.h"
#include "opencv2/core/core.hpp" 
#include "opencv2/features2d/features2d.hpp"
//#include "opencv2/nonfree/features2d.hpp"
#include "opencv2/xfeatures2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/opencv.hpp"

#include "canny.hpp"
#include "datatype.h"
#include "utils.hpp"
#include "image.hpp"
#include "rect.hpp"
#include "tuple.hpp"
#include "elsd.hpp"

using namespace cv;

/*------------------------------------------------------------------------------------------------*/
/**
my code,Alan Lu
输入
img  : 输入图像的一维double型数组,大小为Y*X，按照行优先存储，传入前需要拥有内存
X    : 输入图像的columns
Y    ：输入图像的rows
输出
n_out: lsd算法检测得到的线段的数量n，return的返回值是n条线段，为一维double型数组，长度为8*n，每8个为一组，存着x1,y1,x2,y2,dx,dy,width,polarity
reg_img: 输出标记区域，是一维的int型数组，大小reg_y*reg_x,在相应的像素位置标记着它属于的线段(1,2,3,...n),如果值为0表示不属于任何线段.
         假如外部是int * region_img,则只需要 &region_img,就可以得到标记区域的返回，不需要时直接NULL传入
reg_x  : 输出标记区域的columns,不需要时直接NULL传入
reg_y  : 输出标记区域的rows,不需要时直接NULL传入
*/
double * mylsd(int * n_out, double * img, int X, int Y, int ** reg_img, int * reg_x, int * reg_y)
{
   /* LSD parameters */
  double scale = 0.8;       /* Scale the image by Gaussian filter to 'scale'. */
  double sigma_scale = 0.6; /* Sigma for Gaussian filter is computed as
                                sigma = sigma_scale/scale.                    */
  double quant = 2.0;       /* Bound to the quantization error on the
                                gradient norm.                                */
  double ang_th = 22.5;     /* Gradient angle tolerance in degrees.           */
  double log_eps = 0.0;     /* Detection threshold: -log10(NFA) > log_eps     */
  double density_th = 0.7;  /* Minimal density of region point2is in rectangle. */
  int n_bins = 1024;        /* Number of bins in pseudo-ordering of gradient
                               modulus.                                       */ 

  return LineSegmentDetection( n_out, img, X, Y, scale, sigma_scale, quant,
                               ang_th, log_eps, density_th, n_bins,
                               reg_img, reg_x, reg_y );
}

/*----------------------------------------------------------------------------*/
//输入：
//start_angle,end_angle, 角度方位是(-pi,pi).  
//  pi    ------->x  0
//        |
//        |
//       y\/ pi/2
//polarity: 当polarity为1时，表示的是从start_angle按照逆时针方向旋转到end_angle的角度;当polarity为-1时，表示的是从start_angle按照顺时针方向旋转到end_angle的角度;
//返回值： 旋转角度coverage
inline double rotateAngle(double start_angle, double end_angle, int polarity)
{
  double coverage;
  //首先需要将angle1和angle2转换到 0 ~ 2pi
  if(start_angle < 0) start_angle += M_2__PI;//限制角度在0~2pi之间
  if(end_angle < 0 ) end_angle += M_2__PI;
  if(polarity == 1)//极性为1
  {
    coverage = start_angle - end_angle;
  }
  else //极性为-1
  { 
    coverage = end_angle - start_angle;
  }
  if(coverage < 0) coverage += M_2__PI;
  return coverage;
}
//对线段按照凸性和距离进行分组
//lines: 输入的lines_num条线段，每条线段8个值，存着x1,y1,x2,y2,dx,dy,length,polarity
//lines_num:
//输出分组groups. 每个组是一个vector<int>
//注意：切记用完region,需要在函数外面手动释放region
// todo, read this and check paper
void groupLSs(double *lines, int line_num, int * region, int imgx, int imgy, std::vector<std::vector<int>> * groups)
{
  if(line_num == 0)
  {
    groups = NULL;
    return;
  }
  unsigned char isEnd = 0;//是否还可以继续搜寻
  int currentLine; //当前线段
  char * label = (char*)calloc(line_num, sizeof(char));
  memset(label,0,sizeof(char)*line_num); //init the label all to be zero
  int * group_up = (int*)malloc(sizeof(int)*line_num);//申请足够内存，存储延线段方向得到的分组的线段
  int * group_down = (int*)malloc(sizeof(int)*line_num);//存储线段反方向分组的线段
  int group_up_cnt,group_down_cnt;
  //coorlist * head,*tail;
  std::vector<int> group_temp;
  point2d dir_vec1,dir_vec2;
  point2i *votebin = (point2i*)calloc(line_num,sizeof(point2i));//申请足够内存，用来投票. x记录线段索引，y记录票数
  int bincnt = 0;
  int xx,yy,temp;
  double start_angle,end_angle,angle_delta;
  for ( int i = 0; i<line_num; i++)
  {
    if( label[i] == 0)//未被分组过
    {
      group_up_cnt = group_down_cnt = 0;//每开始寻找一组，需要置零
      //先从第i条线段的头部开始搜索，进行分组,结果存在group_up里面
      group_up[group_up_cnt++] = i;//记录线段i,注意线段是0~line_num-1
      isEnd = 0;//置零，表示还可以从当前线段开始搜索，还未结束
        currentLine = i;
      while(isEnd == 0)
      {
        label[currentLine] = 1; //标记该线段已经被分组
        //head = tail = NULL;
            bincnt = 0;
        dir_vec1.x = lines[currentLine*8+4];
        dir_vec1.y = lines[currentLine*8+5];
        if ( lines[currentLine*8+7] == 1)//极性为正
        {
            //将dir_vec1逆时针旋转45°
            dir_vec2.x = (dir_vec1.x + dir_vec1.y)*0.707106781186548; // sqrt(2)/2 = 0.707106781186548
            dir_vec2.y = (-dir_vec1.x + dir_vec1.y)*0.707106781186548;
        }
        else
        {
            //将dir_vec1顺时针旋转45°
            dir_vec2.x = (dir_vec1.x - dir_vec1.y)*0.707106781186548; // sqrt(2)/2 = 0.707106781186548
            dir_vec2.y = (dir_vec1.x + dir_vec1.y)*0.707106781186548;
        }
        for ( int j = 1; j<=4; j++)
          for ( int k = 1; k<=4; k++)//在4x4邻域内搜索
          {
            xx = (int)(lines[currentLine*8+2]*0.8+j*dir_vec1.x+k*dir_vec2.x);
            yy = (int)(lines[currentLine*8+3]*0.8+j*dir_vec1.y+k*dir_vec2.y);
            if(xx < 0 || xx >= imgx || yy < 0 || yy >= imgy)//越界
              continue;
            temp = region[yy*imgx+xx];
            if(temp>0)//表示有线段的支持区域，在1~line_num
            {
              region[yy*imgx+xx] = -temp;//取负数标记
              for (xx = 0; xx<bincnt; xx++)
              {
                if(votebin[xx].x == temp - 1)//如果以前投票过，直接在相应的bin的票数上加1
                {
                  votebin[xx].y++;
                  break;
                }
              }
              if(xx == bincnt)//如果以前没有投票过，增加该线段，并记录票数为1
              {
                if(bincnt == line_num)
                  error("group ls error1");
                votebin[bincnt].x = temp - 1;
                votebin[bincnt].y = 1;
                bincnt++; //bin的总数加1
              }
            }
          }
          //寻找投票最多的线段，并且需要满足数量大于一定值
          temp = 0;
          for ( int j = 0; j<bincnt; j++)
          {
            if(votebin[j].y>temp)
            {
              temp = votebin[j].y;
              xx = votebin[j].x;//借用xx变量
            }
          }
          if ( temp >= 5 && label[xx] == 0 && lines[8*xx+7] == lines[8*i+7] )//待实验调整参数值
         {
            if(group_up_cnt == line_num)
             error("group ls error2");
            yy = group_up_cnt-1;//借用yy变量
            start_angle = atan2(lines[8*group_up[yy]+5],lines[8*group_up[yy]+4]);
            end_angle = atan2(lines[8*xx+5],lines[8*xx+4]);
            angle_delta = rotateAngle(start_angle,end_angle,(int)lines[8*i+7]);
            if(angle_delta <= M_3_8_PI)//相邻两线段的旋转夹角也需要满足在pi/4内
            {
              group_up[group_up_cnt++] = xx;//压入线段
              currentLine = xx; //更新当前搜索线段
            }
            else
              isEnd = 1;
          }
          else
            isEnd = 1;//结束，已经找不到可以分组的线段了
        }
      //先从第i条线段的尾部开始搜索，进行分组,结果存在group_down里面。记住，第i条线段在group_up和group_down中的0索引处都储存了
      group_down[group_down_cnt++] = i; 
      isEnd = 0;//置零，表示还可以从当前线段开始搜索，还未结束
        currentLine = i;
      while(isEnd == 0)
      {
        label[currentLine] = 1; //标记该线段已经被分组
        //head = tail = NULL;
            bincnt = 0;
        dir_vec1.x = -lines[currentLine*8+4];
        dir_vec1.y = -lines[currentLine*8+5];
        if ( lines[currentLine*8+7] == 1)//极性相同
        {
          //将dir_vec1顺时针旋转45°
          dir_vec2.x = (dir_vec1.x - dir_vec1.y)*0.707106781186548; // sqrt(2)/2 = 0.707106781186548
            dir_vec2.y = (dir_vec1.x + dir_vec1.y)*0.707106781186548;
        }
        else
        {
          //将dir_vec1顺时针旋转45°
          dir_vec2.x = (dir_vec1.x + dir_vec1.y)*0.707106781186548; // sqrt(2)/2 = 0.707106781186548
            dir_vec2.y = (-dir_vec1.x + dir_vec1.y)*0.707106781186548;
        }
        for ( int j = 1; j<=4; j++)
          for ( int k = 1; k<=4; k++)//在4x4邻域内搜索
          {
            xx = (int)(lines[currentLine*8+0]*0.8+j*dir_vec1.x+k*dir_vec2.x);
            yy = (int)(lines[currentLine*8+1]*0.8+j*dir_vec1.y+k*dir_vec2.y);
            if(xx < 0 || xx >= imgx || yy < 0 || yy >= imgy)//越界
              continue;
            temp = region[yy*imgx+xx];
            if(temp>0)//表示有线段的支持区域，在1~line_num
            {
              region[yy*imgx+xx] = -temp;//取负数标记
              for (xx = 0; xx<bincnt; xx++)
              {
                if(votebin[xx].x == temp - 1)//如果以前投票过，直接在相应的bin的票数上加1
                {
                  votebin[xx].y++;
                  break;
                }
              }
              if(xx == bincnt)//如果以前没有投票过，增加该线段，并记录票数为1
              {
                if(bincnt == line_num)
                  error("group ls error3");
                votebin[bincnt].x = temp - 1;
                votebin[bincnt].y = 1;
                bincnt++; //bin的总数加1
              }
            }
          }
          //寻找投票最多的线段，并且需要满足数量大于一定值
          temp = 0;
        for ( int j = 0; j<bincnt; j++)
        {
          if(votebin[j].y>temp)
          {
            temp = votebin[j].y;
            xx = votebin[j].x;//借用xx变量
          }
        }
        if ( temp >= 5 && label[xx] == 0 && lines[8*xx+7] == lines[8*i+7])//待实验调整参数值
        {
          if(group_down_cnt == line_num)
             error("group ls error2");
          yy = group_down_cnt-1;//借用yy变量
          start_angle = atan2(lines[8*group_down[yy]+5],lines[8*group_down[yy]+4]);
          end_angle = atan2(lines[8*xx+5],lines[8*xx+4]);
          angle_delta = rotateAngle(end_angle,start_angle,(int)lines[8*i+7]);//注意此时需要调换一下，因为是从尾部开始搜索
          if(angle_delta < M_3_8_PI)//相邻两线段的旋转夹角也需要满足在pi/4内,pi*3/8 = 66.5°
          {
            group_down[group_down_cnt++] = xx; //压入线段
            currentLine = xx; //更新当前搜索线段
          }
          else
            isEnd = 1;
        }
        else
          isEnd = 1;//结束，已经找不到可以分组的线段了
      }
      (*groups).push_back(group_temp); //添加线段分组
      temp = (*groups).size()-1;
      for (int j = group_down_cnt-1; j>= 0; j--)
      {
        (*groups)[temp].push_back(group_down[j]);
      }
      for (int j = 1; j<group_up_cnt; j++)//由于第i条线段在group_up和group_down都储存了，所以就从索引1开始
      {
        (*groups)[temp].push_back(group_up[j]);
      }
    }
  }
  free(label);
  free(group_up);
  free(group_down);
  free(votebin);
}
//计算groups中每个组的跨度
//输入：
//lines: 输入的lines_num条线段，每条线段8个值，存着x1,y1,x2,y2,dx,dy,length,polarity
//lines_num:
//groups: 分组，每个分组都存着线段的索引
//输出:
//coverages: 每个组的跨度，当组内线段只有1条时，跨度为0. coverages的长度等于组的数量 = groups.size()
//注意，coverages用前不需要申请内存，coverages用完后，需要在函数外手动释放内存，长度等于分组数量
void calcuGroupCoverage(double * lines, int line_num, std::vector<std::vector<int>> groups, double * &coverages)
{
  int groups_num = groups.size();
  int temp;
  double start_angle,end_angle;
  coverages = (double*)malloc(sizeof(double)*groups_num);
  for ( int i = 0; i<groups_num; i++)
  {
    temp = groups[i].size()-1;
    if(groups[i].size() == 0)//第i个分组只有1条线段，则跨度为0
    {
      coverages[i] = 0;
    }
    else
    {
      start_angle = atan2(lines[8*groups[i][0]+5],lines[8*groups[i][0]+4]);
      end_angle = atan2(lines[8*groups[i][temp]+5],lines[8*groups[i][temp]+4]);
      coverages[i] = rotateAngle(start_angle,end_angle,(int)lines[8*groups[i][0]+7]);
    }
  }
}

//==============================================================================
//====================================================================================================
//================================clustering==========================================================
//聚类
//求points中第i行与initializations中第j行里每个元素的平方差总和,每行维度都为nDims
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
 */
void meanShift( double * points, int nPoints, int nDims, double * & initPoints, int initLength, double sigma, double window_size, double accuracy_tolerance, int iter_times )
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
  for ( y = 0; y<nbins_y; y++)//将vote后非0的格子里面的point取均值，并按照顺序重写到center_bins里面，无内存消耗
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
  meanShift(points,points_num,2,initCenters,initCentersLength,1,distance_tolerance,1e-6,50);//迭代20次

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

//================================Generate Ellipse Candidates=========================================
//匹配组对，组对的索引参数，椭圆参数
typedef struct PairGroup_s
{
  point2i pairGroupInd;
  point2d center;  //(x0,y0)
  point2d axis;    //(a,b)
  double  phi;     //angle of orientation  
}PairGroup;

//匹配组对节点
typedef struct PairGroupNode_s
{
  point2i pairGroupInd;
  point2d center;  //(x0,y0)
  point2d axis;    //(a,b)
  double  phi;     //angle of orientation  
  PairGroupNode_s* next;
}PairGroupNode;

typedef struct  PairGroupList_s
{
  int length;
  PairGroup * pairGroup;
}PairGroupList;

typedef struct Point2dNode_s
{
  point2d point;
  Point2dNode_s * next;
}Point2dNode;

typedef struct Point3dNode_s
{
  point3d point;
  Point3dNode_s * next;
}Point3dNode;

typedef struct Point5dNode_s
{
  point2d center;
  point2d axis;
  double  phi;
  Point5dNode_s * next;
}Point5dNode;

typedef struct Point1dNode_s
{
  double data;
  Point1dNode_s * next;
}Point1dNode;

PairGroupList * pairGroupListInit( int length)
{
  if(length <= 0)
    error("paired groups length less equal than 0");
  PairGroupList * pairGroupList = (PairGroupList*)malloc(sizeof(PairGroupList));
  pairGroupList->length = length;
  pairGroupList->pairGroup = (PairGroup*)malloc(sizeof(PairGroup)*length);
  if(pairGroupList->pairGroup == NULL)
    error("pairGroupListInit,not enough memory");
  return pairGroupList;
}

void freePairGroupList( PairGroupList * list)
{
  if(list == NULL || list->pairGroup == NULL)
    error("freePairGroupList,invalidate free");
  free(list->pairGroup);
  free(list);
  list->pairGroup = NULL;
  list = NULL;
}

void calculateGradient2( double * img_in, unsigned int imgx, unsigned int imgy, image_double * angles)
{
  if(img_in == NULL || imgx == 0 || imgy == 0)
    error("calculateGradient error!");
  image_double mod = new_image_double(imgx,imgy);
  (*angles) = new_image_double(imgx,imgy);
  unsigned int x,y,adr;
  double com1,com2;
  double gx,gy;
  double norm,norm_square;
  double threshold;
  double sum = 0;
  double value;  
  //double max_grad = 0.0;
  //边界初始为NOTDEF
  for ( x = 0; x<imgx; x++) 
  {
    (*angles)->data[x]=NOTDEF;
    (*angles)->data[(imgy-1)*imgx+x]=NOTDEF;
    (mod)->data[x]=NOTDEF;
    (mod)->data[(imgy-1)*imgx+x]=NOTDEF;
  }
  for ( y = 0; y<imgy; y++) 
  {
    (*angles)->data[y*imgx] = NOTDEF;
    (*angles)->data[y*imgx+imgx-1] = NOTDEF;
    (mod)->data[y*imgx] = NOTDEF;
    (mod)->data[y*imgx+imgx-1] = NOTDEF;
  }
   /* compute gradient on the remaining pixels */
  for(x=1;x<imgx-1;x++)
    for(y=1;y<imgy-1;y++)
    {
      adr = y*imgx+x;
      /*
         Norm 2 computation using 2x2 pixel window:
           A B C
           D E F
         G H I
         and
           com1 = C-G,  com2 = I-A.
         Then
           gx = C+2F+I - (A+2D+G)=com1+com2+2(F-D)   horizontal difference
           gy = G+2H+I - (A+2B+C)=-com1+com2+2(H-B)   vertical difference
         com1 and com2 are just to avoid 2 additions.
       */
      com1 = img_in[adr-imgx+1] - img_in[adr+imgx-1];
      com2 = img_in[adr+imgx+1] - img_in[adr-imgx-1];

      gx = (com1+com2+2*(img_in[adr+1] - img_in[adr-1]))/(8.0*255); /* gradient x component */
      gy = (-com1+com2+2*(img_in[adr+imgx] - img_in[adr-imgx]))/(8.0*255); /* gradient y component */
      norm_square = gx*gx+gy*gy;
      sum+=norm_square;

      norm = sqrt( norm_square); /* gradient norm */

      (mod)->data[adr] = norm; /* store gradient norm */
       /* gradient angle computation */
       (*angles)->data[adr] = atan2(gy,gx);
    }
  threshold = 2*sqrt(sum/(imgx*imgy));//自动阈值
  //non maximum suppression
  for(x=1;x<imgx-1;x++)
    for(y=1;y<imgy-1;y++)
    {
      adr = y*imgx+x;
      value = (*angles)->data[adr];
      if((mod)->data[adr] < threshold )
      {
        (*angles)->data[adr] = NOTDEF;
        continue;
      }
      if( (value > -M_1_8_PI && value<=M_1_8_PI) || (value <= -M_7_8_PI ) || (value > M_7_8_PI))
      {
        if((mod)->data[adr] <= (mod)->data[adr+1] || (mod)->data[adr] <= (mod)->data[adr-1])
          (*angles)->data[adr] = NOTDEF;
      }
      else if( (value> M_1_8_PI && value<= M_3_8_PI) || (value> -M_7_8_PI && value<= -M_5_8_PI) )
      {
        if((mod)->data[adr] <= (mod)->data[adr-imgx-1] || (mod)->data[adr] <= (mod)->data[adr+imgx+1])
          (*angles)->data[adr] = NOTDEF;
      }
      else if((value> M_3_8_PI && value<= M_5_8_PI) || (value> -M_5_8_PI && value<= -M_3_8_PI))
      {
        if((mod)->data[adr] <= (mod)->data[adr-imgx] || (mod)->data[adr] <= (mod)->data[adr+imgx])
          (*angles)->data[adr] = NOTDEF;
      }
      else 
      {
        if((mod)->data[adr] <= (mod)->data[adr-imgx+1] || (mod)->data[adr] <= (mod)->data[adr+imgx-1])
          (*angles)->data[adr] = NOTDEF;
      }
    }
    //也标记到mod图上面
  //for(x=1;x<imgx-1;x++)
  //  for(y=1;y<imgy-1;y++)
  //  {
  //    if((*angles)->data[y*imgx+x] == NOTDEF)
  //      (mod)->data[y*imgx+x] = NOTDEF;
  //  }
    free_image_double(mod);
}



//canny
void calculateGradient3( double * img_in, unsigned int imgx, unsigned int imgy, image_double * angles)
{
  Mat1b edge;
  Mat1s DX,DY;
  Mat1b gray = Mat::zeros(imgy,imgx,CV_8UC1);
  unsigned int x,y,addr;
  (*angles) = new_image_double(imgx,imgy);
  //copy to gray image
  for ( y = 0; y<imgy; y++)
    for ( x = 0; x<imgx; x++)
    {
      addr = y*imgx+x;
      gray.data[addr] = (uchar)(img_in[addr]);
    }
  //canny
   Canny3(gray,edge,DX,DY,3,false);
   for ( y = 0; y<imgy; y++)
   {
      short* _dx = DX.ptr<short>(y);
    short* _dy = DY.ptr<short>(y);
    uchar* _e = edge.ptr<uchar>(y);
    for ( x = 0; x<imgx; x++)
    {
      if(_e[x] > 0)//0 or 255
      {
        (*angles)->data[y*imgx+x]  = atan2((double)_dy[x],(double)_dx[x]);//calculate gradient 
      }
      else
        (*angles)->data[y*imgx+x] = NOTDEF;
    }
   }
   edge.release();
   DX.release();
   DY.release();
   gray.release();
}


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
  ptrdiff_t fitN = 6;
  double fitWork[64];
  ptrdiff_t workLen = 64;
  ptrdiff_t info;
  //info = LAPACKE_dggev(LAPACK_ROW_MAJOR,'N','V',6,S,6,C,6,alphar,alphai,beta,vl,6,vr,6);
  //注意S为对称矩阵，故转置后等于本身，变成列优先，S可以不变
  dggev(&JOBVL,&JOBVR,&fitN,S,&fitN,C,&fitN,alphar,alphai,beta,vl,&fitN,vr,&fitN,fitWork,&workLen,&info);
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

//input: dataxy为数据点(xi,yi),总共有datanum个
//output: 拟合矩阵S. 注意：S需要事先申请内存，double S[36].
inline void calcuFitMatrix(point2d* dataxy, int datanum, double * S)
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
inline void addFitMatrix(double * S1, double * S2, double * S_out)
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
  ptrdiff_t fitN = 6;
  double fitWork[64];
  ptrdiff_t workLen = 64;
  ptrdiff_t info;
  //info = LAPACKE_dggev(LAPACK_ROW_MAJOR,'N','V',6,S,6,C,6,alphar,alphai,beta,vl,6,vr,6);
  dggev(&JOBVL,&JOBVR,&fitN,S,&fitN,C,&fitN,alphar,alphai,beta,vl,&fitN,vr,&fitN,fitWork,&workLen,&info);
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

//入参：e1 = (x1,y1,a1,b1,phi1), e2 = (x2,y2,a2,b2,phi2)
//输出：相等为1，否则为0
inline bool isEllipseEqual(double * ellipse1, double * ellipse2, double centers_distance_threshold, double semimajor_errorratio, double semiminor_errorratio, double angle_errorratio, double iscircle_ratio)
{
  bool con1 = ( abs(ellipse1[0] - ellipse2[0]) < centers_distance_threshold && abs(ellipse1[1] - ellipse2[1]) < centers_distance_threshold &&
    abs(ellipse1[2] - ellipse2[2])/MAX(ellipse1[2],ellipse2[2]) < semimajor_errorratio && abs(ellipse1[3] - ellipse2[3])/MIN(ellipse1[3],ellipse2[3]) < semiminor_errorratio );
  bool con2 = ( ellipse1[3]/ellipse1[2] >= iscircle_ratio );//0.9 0.85
  bool con3 = ( ellipse2[3]/ellipse2[2] >= iscircle_ratio );
  bool con4 = ( (con2 && con3) || (con2 == false && con3 == false && abs(ellipse1[4]-ellipse2[4])<= angle_errorratio*M_PI) );
  return (con1 && con4);
}

inline bool regionLimitation( point2d point_g1s, point2d g1s_ls_dir, point2d point_g1e, point2d g1e_ls_dir, point2d point_g2s, point2d g2s_ls_dir, point2d point_g2e, point2d g2e_ls_dir, double polarity, double region_limitation_dis_tolerance)
{
  point2d g1m_ls_dir, g2m_ls_dir;
  point2d g1s_arc_dir,g1e_arc_dir,g1m_arc_dir,g2s_arc_dir,g2e_arc_dir,g2m_arc_dir;
  point2d test_vec1,test_vec2,test_vec3; //弧指向圆心的向量和测试向量
  //组的pend<-pstart构成的向量为gim_arc_dir
  double xdelta, ydelta, theta;
  xdelta = point_g1e.x - point_g1s.x;
  ydelta = point_g1e.y - point_g1s.y;
  theta = atan2(ydelta,xdelta);
  g1m_ls_dir.x = cos(theta);
  g1m_ls_dir.y = sin(theta);
  xdelta = point_g2e.x - point_g2s.x;
  ydelta = point_g2e.y - point_g2s.y;
  theta = atan2(ydelta,xdelta);
  g2m_ls_dir.x = cos(theta);
  g2m_ls_dir.y = sin(theta);

  if( polarity == 1)// polarity is equal 1, arc vector = (dy,-dx)
  {
    g1s_arc_dir.x = g1s_ls_dir.y;
    g1s_arc_dir.y = -g1s_ls_dir.x;
    g1e_arc_dir.x = g1e_ls_dir.y;
    g1e_arc_dir.y = -g1e_ls_dir.x;
    g1m_arc_dir.x = g1m_ls_dir.y;
    g1m_arc_dir.y = -g1m_ls_dir.x;
    g2s_arc_dir.x = g2s_ls_dir.y;
    g2s_arc_dir.y = -g2s_ls_dir.x;
    g2e_arc_dir.x = g2e_ls_dir.y;
    g2e_arc_dir.y = -g2e_ls_dir.x;
    g2m_arc_dir.x = g2m_ls_dir.y;
    g2m_arc_dir.y = -g2m_ls_dir.x;
  }
  else// polarity is equal -1, arc vector = (-dy,dx)
  {
    g1s_arc_dir.x = -g1s_ls_dir.y;
    g1s_arc_dir.y = g1s_ls_dir.x;
    g1e_arc_dir.x = -g1e_ls_dir.y;
    g1e_arc_dir.y = g1e_ls_dir.x;
    g1m_arc_dir.x = -g1m_ls_dir.y;
    g1m_arc_dir.y = g1m_ls_dir.x;
    g2s_arc_dir.x = -g2s_ls_dir.y;
    g2s_arc_dir.y = g2s_ls_dir.x;
    g2e_arc_dir.x = -g2e_ls_dir.y;
    g2e_arc_dir.y = g2e_ls_dir.x;
    g2m_arc_dir.x = -g2m_ls_dir.y;
    g2m_arc_dir.y = g2m_ls_dir.x;
  }
  test_vec1.x = (point_g2e.x - point_g1s.x);
  test_vec1.y = (point_g2e.y - point_g1s.y);
  test_vec2.x = (point_g2s.x - point_g1e.x);
  test_vec2.y = (point_g2s.y - point_g1e.y);
  test_vec3.x = (test_vec1.x + test_vec2.x)/2;
  test_vec3.y = (test_vec1.y + test_vec2.y)/2;
  double t1,t2,t3,t4,t5,t6;
  t1 = dotProduct(g1s_arc_dir,test_vec1);
  t2 = dotProduct(g1e_arc_dir,test_vec2);
  t3 = dotProduct(g1m_arc_dir,test_vec3);
  t4 = -dotProduct(g2e_arc_dir,test_vec1);
  t5 = -dotProduct(g2s_arc_dir,test_vec2);
  t6 = -dotProduct(g2m_arc_dir,test_vec3);

  if(  dotProduct(g1s_arc_dir,test_vec1)  >= region_limitation_dis_tolerance && \
     dotProduct(g1e_arc_dir,test_vec2)  >= region_limitation_dis_tolerance && \
     dotProduct(g1m_arc_dir,test_vec3)  >= region_limitation_dis_tolerance && \
    -dotProduct(g2e_arc_dir,test_vec1) >= region_limitation_dis_tolerance && \
      -dotProduct(g2s_arc_dir,test_vec2) >= region_limitation_dis_tolerance && \
    -dotProduct(g2m_arc_dir,test_vec3) >= region_limitation_dis_tolerance
    )
    return TRUE;
  return FALSE;
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

//输入
//lsd算法检测得到的线段集合lines的数量line_num，return的返回值是line_nums条线段，为一维double型数组lines，长度为8*n，每8个为一组
//存着x1,y1,x2,y2,dx,dy,length,polarity
//groups: 线段分组，每个组存按照几何分布顺序顺时针或者逆时针存储着线段索引，线段索引范围是0~line_num-1. 这里由于是指针，使用时要注意(*group)
//first_group_ind、second_group_ind是匹配组队的索引，当提取salient hypothesis时，second_group_ind = -1, fit_matrix2 = NULL.
//fit_matrix1, fit_matrix2, 分别是组队的对应的拟合矩阵
//angles, 是边缘点图+梯度方向。 无边缘点时是NODEF
//distance_tolerance:
//group_inliers_num:记录着各个组的支持内点数量的数组，实时更新，初始时为0
//输出
//ellipara
bool calcEllipseParametersAndValidate( double * lines, int line_num, std::vector<std::vector<int>> * groups, int first_group_ind,int second_group_ind, double * fit_matrix1, double * fit_matrix2, image_double angles, double distance_tolerance, unsigned int * group_inliers_num, point5d *ellipara)
{
  double S[36]; //拟合矩阵S
  double Coefficients[6] = {0,0,0,0,0,0};// ax^2 + bxy + cy^2 + dx + ey + f = 0 
  double param[5], param2[5];
  int info,addr;
  rect rec;
  rect_iter* iter;
  int rec_support_cnt,rec_inliers_cnt;
  bool flag1 = TRUE, flag2 = TRUE;
  double point_normalx, point_normaly, point_normal, temp;
  std::vector<point2i> first_group_inliers, second_group_inliers;
  point2i pixel_temp;
  double semimajor_errorratio,semiminor_errorratio,iscircle_ratio;
  if( fit_matrix2 == NULL || second_group_ind == -1)//只对一个覆盖度较大的组进行拟合
  {
    for ( int i  = 0; i < 36; i++)
      S[i] = fit_matrix1[i];
  }
  else
  {
    addFitMatrix(fit_matrix1,fit_matrix2,S);//对组对进行拟合， S = fit_matrix1 + fit_matrix2
  }
  info = fitEllipse2(S, Coefficients);// ax^2 + bxy + cy^2 + dx + ey + f = 0, a > 0
  if ( info == 0 )//拟合失败
  {
    ellipara = NULL;
    return FALSE;
  }
  ellipse2Param(Coefficients, param);// (x0,y0,a,b,phi)
  if ( min(param[2],param[3]) < 3*distance_tolerance || max(param[2],param[3]) > min(angles->xsize,angles->ysize) ||  param[0] < 0 || param[0] > angles->xsize || param[1] < 0 || param[1] > angles->ysize )
  {
    ellipara = NULL;
    return FALSE;
  }
  //if ( first_group_ind == 2 && second_group_ind == 8)
  //drawEllipse(img,param);
  //组队中的 first group先进行内点准则验证，并且更新组的支持内点数量
  for ( unsigned int i = 0; i<(*groups)[first_group_ind].size(); i++)
  {
    addr = (*groups)[first_group_ind][i] * 8; //第first_group_ind分组的第i条线段索引*8
    rec.x1 = lines[addr];
    rec.y1 = lines[addr+1];
    rec.x2 = lines[addr+2];
    rec.y2 = lines[addr+3];
    rec.x  = (rec.x1 + rec.x2)/2;
    rec.y  = (rec.y1 + rec.y2)/2;
    rec.dx = lines[addr+4];
    rec.dy = lines[addr+5];
    rec.width = 3*distance_tolerance;
    //line_length[i] = (int)lines[addr+6];//记录线段长度到数组line_length[i]
    rec_support_cnt = rec_inliers_cnt = 0;//清零很重要
    if ( lines[addr+7] == 1) //极性一致
    {
      for(iter = ri_ini(&rec);!ri_end(iter);ri_inc(iter))//线段1
      {
        //外接矩形可能会越界
        if(iter->x >= 0 && iter->y >= 0 && iter->x < angles->xsize && iter->y < angles->ysize)
        {
          temp  = angles->data[iter->y*angles->xsize+iter->x] ;//内点的梯度方向
          if(temp!= NOTDEF )
          {
            //test point's normal is (ax0+by0/2+d/2, cy0+bx0/2+e/2)
            point_normalx = Coefficients[0]*iter->x + (Coefficients[1]*iter->y + Coefficients[3])/2;
            point_normaly = Coefficients[2]*iter->y + (Coefficients[1]*iter->x + Coefficients[4])/2;
            point_normal = atan2(-point_normaly,-point_normalx); //边缘点的法线方向,指向椭圆内侧
            rec_inliers_cnt++;
            if(angle_diff(point_normal,temp) <= M_1_8_PI ) //+- 22.5°内 且 || d - r || < 3 dis_t
            {
              rec_support_cnt++;
              pixel_temp.x = iter->x; pixel_temp.y = iter->y;
              first_group_inliers.push_back(pixel_temp);//添加该线段对应的内点
            }
          } 
        }
      }
    }
    else//极性相反
    {
      for(iter = ri_ini(&rec);!ri_end(iter);ri_inc(iter))//线段1
      {
        //外接矩形可能会越界
        if(iter->x >= 0 && iter->y >= 0 && iter->x < angles->xsize && iter->y < angles->ysize)
        {
          temp  = angles->data[iter->y*angles->xsize+iter->x] ;//内点的梯度方向
          if(temp!= NOTDEF )
          {
            //test point's normal is (ax0+by0/2+d/2, cy0+bx0/2+e/2)
            point_normalx = Coefficients[0]*iter->x + (Coefficients[1]*iter->y + Coefficients[3])/2;
            point_normaly = Coefficients[2]*iter->y + (Coefficients[1]*iter->x + Coefficients[4])/2;
            point_normal = atan2(point_normaly,point_normalx); //边缘点的法线方向,指向椭圆外侧
            rec_inliers_cnt++;
            if(angle_diff(point_normal,temp) <= M_1_8_PI ) //+- 22.5°内 且 || d - r || < 3 dis_t
            {
              rec_support_cnt++;
              pixel_temp.x = iter->x; pixel_temp.y = iter->y;
              first_group_inliers.push_back(pixel_temp);//添加该线段对应的内点
            }
          } 
        }
      }
    }
    if( !( rec_support_cnt > 0 && ( rec_support_cnt >= 0.8*lines[addr+6] || rec_support_cnt*1.0/rec_inliers_cnt >= 0.6) ) )
    {
      flag1 = FALSE; //flag1 初始化时为TRUE, 一旦组内有一条线段不满足要求，直接false, 内点准则验证不通过
      break;
    }
  }
  if ( flag1 == TRUE && first_group_inliers.size() >= 0.8*group_inliers_num[first_group_ind] )//靠近最大统计过的内点,通过验证
  {
    if( first_group_inliers.size() >= group_inliers_num[first_group_ind])//更新组出现过的最大内点数
      group_inliers_num[first_group_ind] =  first_group_inliers.size();
  }
  else 
    flag1 = FALSE;
  //第一个组完成验证
  if ( second_group_ind == -1 || fit_matrix2 == NULL)//只对一个覆盖度较大的组进行拟合
  {
    ellipara->x = param[0];//因为无论如何，都需要返回显著性强的椭圆
      ellipara->y = param[1];
      ellipara->a = param[2];
      ellipara->b = param[3];
      ellipara->phi = param[4];
    if ( flag1 == TRUE)//通过内点再次拟合，提高质量
    {
      point2d * dataxy = (point2d*)malloc(sizeof(point2d)*first_group_inliers.size());
      for ( unsigned int i = 0; i<first_group_inliers.size(); i++)
      {
        dataxy[i].x = first_group_inliers[i].x;
        dataxy[i].y = first_group_inliers[i].y;
      }
      info = fitEllipse(dataxy,first_group_inliers.size(), param2);
      free(dataxy); //释放内存
      if ( info == 1  && isEllipseEqual(param2,param,3*distance_tolerance,0.1,0.1,0.1,0.9) )
      {
        ellipara->x = param2[0];//更新椭圆，提高品质
          ellipara->y = param2[1];
          ellipara->a = param2[2];
          ellipara->b = param2[3];
          ellipara->phi = param2[4];
          //drawEllipse(img,param2);
      }
    }
    return TRUE;//对于只有一个组的提取椭圆，此时直接返回
  }
  //接下来，对组队中的 second group进行内点准则验证，并且更新组的支持内点数量
  if (flag1 == FALSE)//在组队运算中，如果第一个组都无法满足内点要求，直接返回false
    return FALSE;
  for ( unsigned int i = 0; i<(*groups)[second_group_ind].size(); i++)
  {
    addr = (*groups)[second_group_ind][i] * 8; //第first_group_ind分组的第i条线段索引*8
    rec.x1 = lines[addr];
    rec.y1 = lines[addr+1];
    rec.x2 = lines[addr+2];
    rec.y2 = lines[addr+3];
    rec.x  = (rec.x1 + rec.x2)/2;
    rec.y  = (rec.y1 + rec.y2)/2;
    rec.dx = lines[addr+4];
    rec.dy = lines[addr+5];
    rec.width = 3*distance_tolerance;
    //line_length[i] = (int)lines[addr+6];//记录线段长度到数组line_length[i]
    rec_support_cnt = rec_inliers_cnt = 0;//清零很重要
    if ( lines[addr+7] == 1) //极性一致
    {
      for(iter = ri_ini(&rec);!ri_end(iter);ri_inc(iter))//线段1
      {
        //外接矩形可能会越界
        if(iter->x >= 0 && iter->y >= 0 && iter->x < angles->xsize && iter->y < angles->ysize)
        {
          temp  = angles->data[iter->y*angles->xsize+iter->x] ;//内点的梯度方向
          if(temp!= NOTDEF )
          {
            //test point's normal is (ax0+by0/2+d/2, cy0+bx0/2+e/2)
            point_normalx = Coefficients[0]*iter->x + (Coefficients[1]*iter->y + Coefficients[3])/2;
            point_normaly = Coefficients[2]*iter->y + (Coefficients[1]*iter->x + Coefficients[4])/2;
            point_normal = atan2(-point_normaly,-point_normalx); //边缘点的法线方向,指向椭圆内侧
            rec_inliers_cnt++;
            if(angle_diff(point_normal,temp) <= M_1_8_PI ) //+- 22.5°内 且 || d - r || < 3 dis_t
            {
              rec_support_cnt++;
              pixel_temp.x = iter->x; pixel_temp.y = iter->y;
              second_group_inliers.push_back(pixel_temp);//添加该线段对应的内点
            }
          } 
        }
      }
    }
    else//极性相反
    {
      for(iter = ri_ini(&rec);!ri_end(iter);ri_inc(iter))//线段1
      {
        //外接矩形可能会越界
        if(iter->x >= 0 && iter->y >= 0 && iter->x < angles->xsize && iter->y < angles->ysize)
        {
          temp  = angles->data[iter->y*angles->xsize+iter->x] ;//内点的梯度方向
          if(temp!= NOTDEF )
          {
            //test point's normal is (ax0+by0/2+d/2, cy0+bx0/2+e/2)
            point_normalx = Coefficients[0]*iter->x + (Coefficients[1]*iter->y + Coefficients[3])/2;
            point_normaly = Coefficients[2]*iter->y + (Coefficients[1]*iter->x + Coefficients[4])/2;
            point_normal = atan2(point_normaly,point_normalx); //边缘点的法线方向,指向椭圆外侧
            rec_inliers_cnt++;
            if(angle_diff(point_normal,temp) <= M_1_8_PI ) //+- 22.5°内 且 || d - r || < 3 dis_t
            {
              rec_support_cnt++;
              pixel_temp.x = iter->x; pixel_temp.y = iter->y;
              second_group_inliers.push_back(pixel_temp);//添加该线段对应的内点
            }
          } 
        }
      }
    }
    if( !(rec_support_cnt > 0 && ( rec_support_cnt >= 0.8*lines[addr+6] || rec_support_cnt*1.0/rec_inliers_cnt >= 0.6) ) )
    {
      flag2 = FALSE; //flag1 初始化时为TRUE, 一旦组内有一条线段不满足要求，直接false, 内点准则验证不通过
      break;
    }
  }
  if ( flag2 == TRUE && second_group_inliers.size() >= 0.8*group_inliers_num[second_group_ind] )//靠近最大统计过的内点,通过验证
  {
    if(second_group_inliers.size() >= group_inliers_num[second_group_ind])//更新组出现过的最大内点数
      group_inliers_num[second_group_ind] = second_group_inliers.size();
  }
  else 
    flag2 = FALSE;
  //第二个组完成验证
  if ( flag1 == TRUE && flag2 == TRUE)
  {
    point2d * dataxy = (point2d*)malloc(sizeof(point2d)*(first_group_inliers.size() + second_group_inliers.size()));
    for ( unsigned int i = 0; i<first_group_inliers.size(); i++)
    {
      dataxy[i].x = first_group_inliers[i].x;
      dataxy[i].y = first_group_inliers[i].y;
    }
    addr = first_group_inliers.size();
    for ( unsigned int i = 0; i<second_group_inliers.size(); i++)//连接两个数组时一定要注意索引范围
    {
      dataxy[addr+i].x = second_group_inliers[i].x;
      dataxy[addr+i].y = second_group_inliers[i].y;
    }
//    drawEdge(img,dataxy,(first_group_inliers.size() + second_group_inliers.size()));
    info = fitEllipse(dataxy,(first_group_inliers.size() + second_group_inliers.size()), param2);
    free(dataxy); //释放内存
    //小长短轴的椭圆需要放宽参数
    if ( param[2] <= 50 )
      semimajor_errorratio = 0.25;
    else if (param[2] <= 100 )
      semimajor_errorratio = 0.15;
    else
      semimajor_errorratio = 0.1;
    if ( param[3] <= 50 )
      semiminor_errorratio = 0.25;
    else if ( param[3] <= 100)
      semiminor_errorratio = 0.15;
    else
      semiminor_errorratio = 0.1;
    if (param[2] <= 50 && param[3] <= 50 )
      iscircle_ratio = 0.75;
    else if (param[2] >= 50 && param[2] <= 100 &&  param[3] >= 50 && param[3] <= 100 )
      iscircle_ratio = 0.85;
    else
      iscircle_ratio = 0.9;
    if ( info == 1  && isEllipseEqual(param2,param,3*distance_tolerance,semimajor_errorratio,semiminor_errorratio,0.1, iscircle_ratio) )
    {
      ellipara->x = param2[0];//更新椭圆，提高品质
        ellipara->y = param2[1];
        ellipara->a = param2[2];
        ellipara->b = param2[3];
        ellipara->phi = param2[4];
        //drawEllipse(img,param2);
      return TRUE;
    }
  }
  return FALSE;
}


//输入
//lsd算法检测得到的线段集合lines的数量line_num，return的返回值是line_nums条线段，为一维double型数组lines，长度为8*n，每8个为一组
//存着x1,y1,x2,y2,dx,dy,length,polarity
//groups: 线段分组，每个组存按照几何分布顺序顺时针或者逆时针存储着线段索引，线段索引范围是0~line_num-1
//coverages: 每个分组的角度覆盖范围0~2pi，如果组里只有1条线段，覆盖角度为0。数组长度等于分组的数量。
//angles 存边缘点的梯度方向gradient direction, 无边缘点位NOTDEF
//返回值 PairedGroupList* list 返回的是初始椭圆集合的数组，长度list->length. 
//切记，该内存在函数内申请，用完该函数记得释放内存，调用函数freePairedSegmentList()进行释放

PairGroupList * getValidInitialEllipseSet( double * lines, int line_num, std::vector<std::vector<int>> * groups, double * coverages, image_double angles, double distance_tolerance, int specified_polarity)
{
  //加速计算
  //int* lineInliersIndex = (int*)malloc(sizeof(int)*line_num);//如果第i条线段找到了内点，则记录其索引为j = length(supportInliers),即supportInliers.at(j)存着该线段的支持内点,没找到内点的线段对应索引为初始值-1.
    //vector<vector<point2d>> supportInliers;//保存相应线段的支持内点
  //memset(lineInliersIndex,-1,sizeof(int)*line_num);//此处要实践确实可行，对于整数可以初始化为0，-1.对于浮点数则只可以为0.

  PairGroupList * pairGroupList = NULL;
  PairGroupNode *head, *tail;
  int pairlength = 0;
  point2d pointG1s,pointG1e,pointG2s,pointG2e,g1s_ls_dir,g1e_ls_dir,g2s_ls_dir,g2e_ls_dir;
  double polarity;
  point5d ellipara;
    int groupsNum = (*groups).size();//组的数量
  double * fitMatrixes = (double*)malloc(sizeof(double)*groupsNum*36);//定义拟合矩阵S_{6 x 6}. 每个组都有一个拟合矩阵
  unsigned int * supportInliersNum = (unsigned int*)malloc(sizeof(int)*groupsNum);//用于存储每个组曾经最大出现的支持内点数量
  memset(fitMatrixes,0,sizeof(double)*groupsNum*36);
  memset(supportInliersNum, 0, sizeof(unsigned int)*groupsNum);//初始化为0.
  //double distance_tolerance = max( 2.0, 0.005*min(angles->xsize,angles->ysize) ); // 0.005%*min(xsize,ysize)
    int i,j;
  int cnt_temp,ind_start,ind_end;
  bool info;
    
  //实例化拟合矩阵Si
  point2d * dataxy = (point2d*)malloc(sizeof(point2d)*line_num*2);//申请足够大内存, line_num条线段，共有2line_num个端点
  for ( i = 0; i<groupsNum; i++)
  {
    cnt_temp = 0;//千万注意要清0
    for ( j = 0; j<(*groups)[i].size(); j++)
    {
      //每一条线段有2个端点
      dataxy[cnt_temp].x = lines[(*groups)[i][j]*8];
      dataxy[cnt_temp++].y = lines[(*groups)[i][j]*8+1];
      dataxy[cnt_temp].x = lines[(*groups)[i][j]*8+2];
      dataxy[cnt_temp++].y = lines[(*groups)[i][j]*8+3];
    }
    calcuFitMatrix(dataxy,cnt_temp, fitMatrixes+i*36);
  }
  free(dataxy);//释放内存

  head = tail = NULL;//将初始椭圆集合存储到链表中
  //selection of salient elliptic hypothesis
  for ( i = 0; i<groupsNum; i++)
  {
    if(coverages[i] >= M_4_9_PI )//当组的覆盖角度>= 4pi/9 = 80°, 我们认为具有很大的显著性，可直接拟合提取
    {
      //加入极性判断,只提取指定极性的椭圆
      if (specified_polarity == 0 || (lines[(*groups)[i][0]*8+7] == specified_polarity))
      {
        //显著性大的初始椭圆提取，一定会返回TRUE，因此没必要再判断
        info = calcEllipseParametersAndValidate(lines,line_num,groups,i,-1,(fitMatrixes+i*36),NULL,angles,distance_tolerance,supportInliersNum,&ellipara);
        if (info == FALSE) 
        {
          continue;
          error("getValidInitialEllipseSet, selection of salient ellipses failed!");//这种情况会出现？？,跑54.jpg出现该问题
        }
        PairGroupNode * node = (PairGroupNode*)malloc(sizeof(PairGroupNode));
        node->center.x = ellipara.x;
        node->center.y = ellipara.y;
        node->axis.x   = ellipara.a;
        node->axis.y   = ellipara.b;
        node->phi      = ellipara.phi;
        node->pairGroupInd.x = i;
        node->pairGroupInd.y = -1;//无
        if(head != NULL)
        {
          tail->next = node;
          tail = node;
        }
        else
        {
          head = tail = node;
        }
        pairlength++;
      }
    }
  }
    //selection of pair group hypothesis
  for ( i = 0; i<groupsNum-1; i++)
    for ( j = i+1; j<groupsNum; j++)
      {
        //加入极性判断,只提取指定极性的椭圆
         if (specified_polarity == 0 || (lines[(*groups)[i][0]*8+7] == specified_polarity))
          {
          //group i 's polarity is the same as group j; and the number of two paired groups should be >= 3.
          if( lines[(*groups)[i][0]*8+7] == lines[(*groups)[j][0]*8+7] && ((*groups)[i].size() + (*groups)[j].size()) >= 3)
          {
            ind_start = (*groups)[i][0];//第i组的最开始一条线段索引
            ind_end   = (*groups)[i][(*groups)[i].size()-1];//第i组的最后一条线段索引
            pointG1s.x = lines[ind_start*8];
            pointG1s.y = lines[ind_start*8+1];
            g1s_ls_dir.x = lines[ind_start*8+4];
            g1s_ls_dir.y = lines[ind_start*8+5];
            pointG1e.x = lines[ind_end*8+2];
            pointG1e.y = lines[ind_end*8+3];
            g1e_ls_dir.x = lines[ind_end*8+4];
            g1e_ls_dir.y = lines[ind_end*8+5];
            ind_start = (*groups)[j][0];//第j组的最开始一条线段索引
            ind_end   = (*groups)[j][(*groups)[j].size()-1];//第j组的最后一条线段索引
            pointG2s.x = lines[ind_start*8];
            pointG2s.y = lines[ind_start*8+1];
            g2s_ls_dir.x = lines[ind_start*8+4];
            g2s_ls_dir.y = lines[ind_start*8+5];
            pointG2e.x = lines[ind_end*8+2];
            pointG2e.y = lines[ind_end*8+3];
            g2e_ls_dir.x = lines[ind_end*8+4];
            g2e_ls_dir.y = lines[ind_end*8+5];
            polarity = lines[ind_start*8+7]; //i,j两组的极性
            if(regionLimitation(pointG1s,g1s_ls_dir,pointG1e,g1e_ls_dir,pointG2s,g2s_ls_dir,pointG2e,g2e_ls_dir,polarity,-3*distance_tolerance))//都在彼此的线性区域内
            {
              //if ( i == 2)
              //  drawPairGroup(img,lines,(*groups),i,j);

              if(calcEllipseParametersAndValidate(lines,line_num,groups,i,j,(fitMatrixes+i*36),(fitMatrixes+j*36),angles,distance_tolerance,supportInliersNum,&ellipara))//二次一般方程线性求解，线段的内点支持比例
              {
                PairGroupNode * node = (PairGroupNode*)malloc(sizeof(PairGroupNode));
                node->center.x = ellipara.x;
                node->center.y = ellipara.y;
                node->axis.x   = ellipara.a;
                node->axis.y   = ellipara.b;
                node->phi      = ellipara.phi;
                node->pairGroupInd.x = i;
                node->pairGroupInd.y = -1;//无
                if(head != NULL)
                {
                  tail->next = node;
                  tail = node;
                }
                else
                {
                  head = tail = node;
                }
                pairlength++;
              }
            }
            
          }
         }
      }
  if(pairlength > 0)
  {
    PairGroupNode *p;
    p = head;
    pairGroupList = pairGroupListInit(pairlength);
    for( int i = 0; i<pairGroupList->length; i++)
    {
      pairGroupList->pairGroup[i].center.x = p->center.x;
      pairGroupList->pairGroup[i].center.y = p->center.y;
      pairGroupList->pairGroup[i].axis.x = p->axis.x;
      pairGroupList->pairGroup[i].axis.y = p->axis.y;
      pairGroupList->pairGroup[i].phi = p->phi;
      pairGroupList->pairGroup[i].pairGroupInd.x = p->pairGroupInd.x;//记录组对(i,j),由groups中的第i个组和第j个组构成的匹配组产生该有效椭圆参数
      pairGroupList->pairGroup[i].pairGroupInd.y = p->pairGroupInd.y;
      p = p->next;
    }
    tail->next = NULL;
    while (head != NULL)
    {
      p = head;
      head = head->next;
      free(p);
    }
  }
  //supportInliers.resize(0);
  //free(lineInliersIndex);//释放线段内点的索引
  free(supportInliersNum);//释放存储各个组的支持内点数量的数组
  free(fitMatrixes);//释放存储各个组的拟合矩阵
  return pairGroupList;
}


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







//==========================================END=======================================================================
/**
输入：
prhs[0]: 输入的灰度图像，单通道，大小是imgy x imgx
prhs[1]: 边缘提取选择，1 canny; 2 sobel
prhs[2]: 检测指定的椭圆极性
输出：
plhs[0]: 候选椭圆组合(xi,yi,ai,bi,phi_i)', 5 x m
plhs[1]: 边缘图，大小是imgy x imgx，设边缘点总数为 edgepix_n. 二值化，0 或者 255
plhs[2]: 边缘点的梯度向量矩阵，大小是 2 x edgepix_n, (cos(theta_rad),sin(theta_rad))'...
plhs[3]: 线段图，大小是imgy x imgx 
*/
/*
compile：
mex generateEllipseCandidates.cpp -IF:\OpenCV\opencv2.4.9\build\include -IF:\OpenCV\opencv2.4.9\build\include\opencv -IF:\OpenCV\opencv2.4.9\build\include\opencv2 -LF:\OpenCV\opencv2.4.9\build\x64\vc11\lib -IF:\Matlab\settlein\extern\include -LF:\Matlab\settlein\extern\lib\win64\microsoft -lopencv_core249 -lopencv_highgui249 -lopencv_imgproc249 -llibmwlapack.lib
*/
//======================================MEX function==================================================================

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs!=3)
        mexErrMsgIdAndTxt( "MATLAB:revord:invalidNumInputs","One input required.");
    else if(nlhs > 4)
        mexErrMsgIdAndTxt( "MATLAB:revord:maxlhs","Too many output arguments.");
    uchar * inputimg = (uchar*)mxGetData(prhs[0]);
    int imgy,imgx;
    int edge_process_select = (int)mxGetScalar(prhs[1]);//边缘提取选择，1 canny; 2 sobel
    int specified_polarity  = (int)mxGetScalar(prhs[2]);//1,指定检测的椭圆极性要为正; -1指定极性为负; 0表示两种极性椭圆都检测
    imgy = (int)mxGetM(prhs[0]);
    imgx = (int)mxGetN(prhs[0]);
    double *data=(double*)malloc(imgy*imgx*sizeof(double));//将输入矩阵中的图像数据转存到一维数组中
    for(int c=0;c<imgx;c++)
    {
        for(int r=0;r<imgy;r++)
        {
            data[c+r*imgx]=inputimg[r+c*imgy];
        }
    }
    int n;//线段数量
    //int new_n;
    std::vector<std::vector<int>> groups;
    double * coverages;
    int * reg;
    int reg_x;
    int reg_y;
    double* out=mylsd(&n, data,imgx,imgy,&reg,&reg_x,&reg_y);
    groupLSs(out,n,reg,reg_x,reg_y,&groups);//分组
    free(reg); //释放内存
    calcuGroupCoverage(out,n,groups,coverages);//计算每个组的覆盖角度

    printf("The number of output arc-support line segments: %i\n",n);
    printf("The number of arc-support groups:%i\n",groups.size());
  /*int groups_t = 0;
  for (int i = 0; i<groups.size(); i++)
  { 
    groups_t+= groups[i].size();
  }
  printf("Groups' total ls num:%i\n",groups_t);*/

    image_double angles;
    if(edge_process_select == 1)
        calculateGradient2(data,imgx,imgy,&angles); //version2, sobel; version 3 canny
    else
        calculateGradient3(data,imgx,imgy,&angles); //version2, sobel; version 3 canny
    PairGroupList * pairGroupList;
    double distance_tolerance = 2;//max( 2.0, 0.005*min(angles->xsize,angles->ysize) ); // 0.005%*min(xsize,ysize)
    double * candidates; //候选椭圆
    double * candidates_out;//输出候选椭圆指针
    int  candidates_num = 0;//候选椭圆数量
    //rejectShortLines(out,n,&new_n);
    pairGroupList = getValidInitialEllipseSet(out,n,&groups,coverages,angles,distance_tolerance,specified_polarity);
    if(pairGroupList != NULL)
    {
        printf("The number of initial ellipses：%i \n",pairGroupList->length);
        generateEllipseCandidates(pairGroupList, distance_tolerance, candidates, &candidates_num);
        printf("The number of ellipse candidates: %i \n",candidates_num);
    
        plhs[0] = mxCreateDoubleMatrix(5,candidates_num,mxREAL);
        candidates_out = (double*)mxGetPr(plhs[0]);
        //候选圆组合(xi,yi,ai,bi,phi_i)', 5 x candidates_num, 复制到矩阵candidates_out中
        memcpy(candidates_out,candidates,sizeof(double)*5*candidates_num);

        freePairGroupList(pairGroupList);
        free(candidates);
   }
   else
   {
        printf("The number of initial ellipses：%i \n",0);
        double *candidates_out;
        plhs[0] = mxCreateDoubleMatrix(5,1,mxREAL);
        candidates_out = (double*)mxGetPr(plhs[0]);
        candidates_out[0] = candidates_out[1] = candidates_out[2] = candidates_out[3] = candidates_out[4] = 0;
   }
   uchar *edgeimg_out;
   unsigned long edge_pixels_total_num = 0;//边缘总像素
   double *gradient_vec_out;
   plhs[1] = mxCreateNumericMatrix(imgy,imgx,mxUINT8_CLASS,mxREAL);
   edgeimg_out = (uchar*)mxGetData(plhs[1]);
   //将边缘图复制到矩阵edgeimg_out中
   //将梯度向量存到矩阵gradient_vec_out中
   unsigned long addr,g_cnt = 0;
   for ( int c = 0; c < imgx; c++ )
     for ( int r = 0; r < imgy; r++)
     {
       addr = r*imgx+c;
       if(angles->data[addr] == NOTDEF)
         edgeimg_out[c*imgy+r] = 0;
       else
       {
         edgeimg_out[c*imgy+r] = 255;//为边缘点，赋值为白色
         //------------------------------------------------
         edge_pixels_total_num++;
       }
     }
   printf("edge pixel number: %i\n",edge_pixels_total_num);
  //申请edge_pixels_total_num x 2 来保存每一个边缘点的梯度向量，以列为优先，符合matlab的习惯
   plhs[2] = mxCreateDoubleMatrix(2,edge_pixels_total_num,mxREAL);
   gradient_vec_out = (double*)mxGetPr(plhs[2]);
    for ( int c = 0; c < imgx; c++ )
     for ( int r = 0; r < imgy; r++)
     {
       addr = r*imgx+c;
       if(angles->data[addr] != NOTDEF)
       {
         gradient_vec_out[g_cnt++] = cos(angles->data[addr]);
         gradient_vec_out[g_cnt++] = sin(angles->data[addr]);
       }
     }
   //---------------------------------------------------------------------
  //输出线段检测的图像
  if(nlhs == 4)
  {
    Mat ls_mat = Mat::zeros(imgy,imgx,CV_8UC1);
    for ( int i = 0; i<n ; i++)//draw lines
    {
      Point2d p1(out[8*i],out[8*i+1]),p2(out[8*i+2],out[8*i+3]);
      line(ls_mat,p1,p2,Scalar(255,0,0));
    }
    if(candidates_num > 0)//draw ellipses
    {
      for ( int i = 0; i<candidates_num; i++)
        ellipse(ls_mat,cv::Point((int)candidates_out[i*5],(int)candidates_out[i*5+1]),cv::Size(candidates_out[i*5+2],candidates_out[i*5+3]),candidates_out[i*5+4]*180/M_PI,0,360,(Scalar(255,0,0)),1);
    }
    plhs[3] = mxCreateDoubleMatrix(imgy,imgx,mxREAL);
    double * ls_img_out = (double*)mxGetPr(plhs[3]);
    //memcpy(ls_out_mat,ls_mat.data ,sizeof(unsigned char)*M*N);
    for (int i = 0; i<imgx; i++)
      for (int j = 0; j<imgy;j++)
        ls_img_out[i*imgy+j]=ls_mat.data[j*imgx+i];
  }
  //---------------------------------------------------------------------
  //这里的free是释放程序中用于产生候选圆所用到的一系列内存
  free(data);
  free(coverages);
  free(out);
  free_image_double(angles);

}
