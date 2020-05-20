//
//  group_forming.cpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
//

#include "group_forming.hpp"
#include "datatype.h"
#include "utils.hpp"
#include <math.h>

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
// region: a img with each pixel indicate whether the pixel is located on a line (i) or not
// section III-A, done
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
                dir_vec1.x = lines[currentLine*8+4];// tuple-8
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
                    if(angle_delta < M_3_8_PI)//相邻两线段的旋转夹角也需要满足在pi/4内,pi*3/8 =  66.5°
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
