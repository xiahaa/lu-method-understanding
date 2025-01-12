//
//  gen_init_set.cpp
//  lu
//
//  Created by xiao.hu on 2020/5/19.
//  Copyright © 2020 xiao.hu. All rights reserved.
// another contribution from lu

#include "gen_init_set.hpp"
#include <math.h>
#include "utils.hpp"
#include "opencv2/core/core.hpp"
#include "rect.hpp"
#include "ellipse_geometry.hpp"
#include <iostream>

using namespace cv;

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


//输入
//lsd算法检测得到的线段集合lines的数量line_num，return的返回值是line_nums条线段，为一维double型数组lines，长度为TUPLELENGTH*n，每8个为一组
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
        addr = (*groups)[first_group_ind][i] * TUPLELENGTH; //第first_group_ind分组的第i条线段索引*TUPLELENGTH
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
        addr = (*groups)[second_group_ind][i] * TUPLELENGTH; //第first_group_ind分组的第i条线段索引*TUPLELENGTH
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


bool find_P(const point2d &Q1, const point2d &Q2, const point2d &Q3, const point2d &Q4, point2d& Pout)
{
	cv::Vec3d Q1h = cv::Vec3d(Q1.x, Q1.y, 1);
	cv::Vec3d Q2h = cv::Vec3d(Q2.x, Q2.y, 1);
	cv::Vec3d Q3h = cv::Vec3d(Q3.x, Q3.y, 1);
	cv::Vec3d Q4h = cv::Vec3d(Q4.x, Q4.y, 1);

	cv::Vec3d l1 = Q1h.cross(Q2h);
	cv::Vec3d l2 = Q3h.cross(Q4h);

	cv::Vec3d P = l1.cross(l2);

	if (fabs(P[2]) < 1e-3)
	{
		Pout.x = 0;
		Pout.y = 0;
		return false;
	}
	else
	{
		Pout.x = P[0] / P[2];
		Pout.y = P[1] / P[2];
		return true;
	}
}

void find_ab(const point2d &P1, const point2d &P2, const point2d &Q1, double &a, double &b)
{
	cv::Mat M = (cv::Mat_<double>(2, 2) << P1.x, P2.x, P1.y, P2.y);
	cv::Mat m = (cv::Mat_<double>(2, 1) << Q1.x, Q1.y);
	cv::Mat x;
	cv::solve(M, m, x);
	a = x.at<double>(0, 0);
	b = x.at<double>(1, 0);
}

// st, mid, end
char checkCN(const point2d &Q11, const point2d &Q12, const point2d &Q21,
	const point2d &Q22, const point2d &Q31, const point2d &Q32)
{
	point2d P1, P2, P3;
	if (!find_P(Q11, Q12, Q31, Q32, P1)) return 0;
	if (!find_P(Q11, Q12, Q21, Q22, P2)) return 0;
	if (!find_P(Q21, Q22, Q31, Q32, P3)) return 0;

	double a11, b11, a12, b12, a21, b21, a22, b22, a31, b31, a32, b32;
	find_ab(P1, P2, Q11, a11, b11);
	find_ab(P1, P2, Q12, a12, b12);
	find_ab(P2, P3, Q21, a21, b21);
	find_ab(P2, P3, Q22, a22, b22);
	find_ab(P3, P1, Q31, a31, b31);
	find_ab(P3, P1, Q32, a32, b32);

	double cn_score = b11 / a11 * b12 / a12 * b21 / a21 * b22 / a22 * b31 / a31 * b32 / a32;

	if (fabs(cn_score - 1) < 1.5)
		return 1;
	else
		return -1;

}

//输入
//lsd算法检测得到的线段集合lines的数量line_num，return的返回值是line_nums条线段，为一维double型数组lines，长度为TUPLELENGTH*n，每TUPLELENGTH个为一组
//存着x1,y1,x2,y2,dx,dy,length,polarity
//groups: 线段分组，每个组存按照几何分布顺序顺时针或者逆时针存储着线段索引，线段索引范围是0~line_num-1
//coverages: 每个分组的角度覆盖范围0~2pi，如果组里只有1条线段，覆盖角度为0。数组长度等于分组的数量。
//angles 存边缘点的梯度方向gradient direction, 无边缘点位NOTDEF
//返回值 PairedGroupList* list 返回的是初始椭圆集合的数组，长度list->length.
//切记，该内存在函数内申请，用完该函数记得释放内存，调用函数freePairedSegmentList()进行释放

PairGroupList * getValidInitialEllipseSet( double * lines, int line_num, std::vector<std::vector<int>> * groups, 
	double * coverages, image_double angles, double distance_tolerance, int specified_polarity)
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

	// label, if an arc has been connected
	/* the idea here is that: whenever an ellipse is fitted, find affine A, warp other arcs back to circle, vote to check
	if this arc belongs to the ellipse, if yes, then this arc does not need to be combined for ellipse fitting. 
	I hope this can reduce the number of initial candidate ellipse since the for for loop may produce lots of candidate ellipse.
	*/
#if USE_AFFINE_WARP
	unsigned char *isLabelled = (unsigned char*)calloc(groupsNum, sizeof(unsigned char));
#endif

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
            dataxy[cnt_temp].x = lines[(*groups)[i][j]* TUPLELENGTH];
            dataxy[cnt_temp++].y = lines[(*groups)[i][j]* TUPLELENGTH +1];
            dataxy[cnt_temp].x = lines[(*groups)[i][j]* TUPLELENGTH +2];
            dataxy[cnt_temp++].y = lines[(*groups)[i][j]* TUPLELENGTH +3];
        }
        calcuFitMatrix(dataxy,cnt_temp, fitMatrixes+i*36);
    }
    free(dataxy);//释放内存
    
    head = tail = NULL;//将初始椭圆集合存储到链表中
    //selection of salient elliptic hypothesis
    for ( i = 0; i<groupsNum; i++)
    {
#if USE_AFFINE_WARP
		if (isLabelled[i] != 0)
			continue;
#endif

        if(coverages[i] >= M_4_9_PI )//当组的覆盖角度>= 4pi/9 = 80°, 我们认为具有很大的显著性，可直接拟合提取
        {
            //加入极性判断,只提取指定极性的椭圆
            if (specified_polarity == 0 || (lines[(*groups)[i][0]* TUPLELENGTH +7] == specified_polarity))
            {
                //显著性大的初始椭圆提取，一定会返回TRUE，因此没必要再判断
                info = calcEllipseParametersAndValidate(lines,line_num,groups,i,-1,(fitMatrixes+i*36),NULL,angles,distance_tolerance,supportInliersNum,&ellipara);
                if (info == FALSE)
                {
                    continue;
                    error("getValidInitialEllipseSet, selection of salient ellipses failed!");//这种情况会出现？？,跑54.jpg出现该问题
                }

#if USE_AFFINE_WARP
				// check arcs
				float Phi = ellipara.phi;
				float cphi = cos(Phi);
				float sphi = sin(Phi);
				float a[3][3] = { { ellipara.a * cphi, -ellipara.b * sphi, ellipara.x}, \
									{ ellipara.a * sphi, ellipara.b * cphi, ellipara.y}, \
									{ 0, 0, 1} };
				// this is the affine matrix that warps points back to unit circle
				Mat Ainv = Mat(3, 3, CV_32F, a).inv();

#if USE_SPAN_CHECK
				// check coverage span
				// start point of first LSD
				float x1 = Ainv.at<float>(0, 0) * lines[(*groups)[i][0] * TUPLELENGTH] + Ainv.at<float>(0, 1) * lines[(*groups)[i][0] * TUPLELENGTH + 1] + Ainv.at<float>(0, 2);
				float y1 = Ainv.at<float>(1, 0) * lines[(*groups)[i][0] * TUPLELENGTH] + Ainv.at<float>(1, 1) * lines[(*groups)[i][0] * TUPLELENGTH + 1] + Ainv.at<float>(1, 2);
				int kk = (*groups)[i].size()-1;
				// end point of last LSD
				float x2 = Ainv.at<float>(0, 0) * lines[(*groups)[i][kk] * TUPLELENGTH + 2] + Ainv.at<float>(0, 1) * lines[(*groups)[i][kk] * TUPLELENGTH + 3] + Ainv.at<float>(0, 2);
				float y2 = Ainv.at<float>(1, 0) * lines[(*groups)[i][kk] * TUPLELENGTH + 2] + Ainv.at<float>(1, 1) * lines[(*groups)[i][kk] * TUPLELENGTH + 3] + Ainv.at<float>(1, 2);
				float ang1 = atan2(y1, x1);
				ang1 = ang1 < 0 ? ang1 + M_PI * 2 : ang1;
				float ang2 = atan2(y2, x2);
				ang2 = ang2 < 0 ? ang2 + M_PI * 2 : ang2;
				float span;
				if (lines[(*groups)[i][0] * TUPLELENGTH + 7] == -1)// counter-clockwise
				{
					ang2 = ang2 < ang1 ? ang2 + M_PI * 2 : ang2;
					span = ang2 - ang1;
				}
				else if (lines[(*groups)[i][0] * TUPLELENGTH + 7] == 1)// clockwise
				{
					ang1 = ang1 < ang2 ? ang1 + M_PI * 2 : ang1;
					span = ang1 - ang2;
				}

				if (coverages[i] / (M_PI * 2) < SPAN_CHECK_THRESH) continue;
#endif

				//std::cout << Ainv << std::endl;
				//Mat M = (Mat_<double>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);
#if USE_LABEL_MATCH
				std::vector<int> candidates;
				for (int j = 0; j < groupsNum; j++)
				{
					if (j == i) continue;
					if (lines[(*groups)[i][0] * TUPLELENGTH + 7] != lines[(*groups)[j][0] * TUPLELENGTH + 7]) continue;

					int validsum = 0;
					for (int k = 0; k<(*groups)[j].size(); k++)
					{
						//std::cout << lines[(*groups)[j][k] * TUPLELENGTH] << ";" << lines[(*groups)[j][k] * TUPLELENGTH + 1] << std::endl;
						//std::cout << lines[(*groups)[j][k] * TUPLELENGTH+2] << ";" << lines[(*groups)[j][k] * TUPLELENGTH + 3] << std::endl;

						//每一条线段有2个端点
						float x1 = Ainv.at<float>(0, 0) * lines[(*groups)[j][k] * TUPLELENGTH] + Ainv.at<float>(0, 1) * lines[(*groups)[j][k] * TUPLELENGTH + 1] + Ainv.at<float>(0, 2);
						float y1 = Ainv.at<float>(1, 0) * lines[(*groups)[j][k] * TUPLELENGTH] + Ainv.at<float>(1, 1) * lines[(*groups)[j][k] * TUPLELENGTH + 1] + Ainv.at<float>(1, 2);

						float x2 = Ainv.at<float>(0, 0) * lines[(*groups)[j][k] * TUPLELENGTH +2] + Ainv.at<float>(0, 1) * lines[(*groups)[j][k] * TUPLELENGTH + 3] + Ainv.at<float>(0, 2);
						float y2 = Ainv.at<float>(1, 0) * lines[(*groups)[j][k] * TUPLELENGTH +2] + Ainv.at<float>(1, 1) * lines[(*groups)[j][k] * TUPLELENGTH + 3] + Ainv.at<float>(1, 2);

						if (fabs(x1*x1 + y1*y1 - 1) < WARP_RADIUS_THRESH && fabs(x2*x2 + y2*y2 - 1) < WARP_RADIUS_THRESH)
							validsum++;
					}
					if (float(validsum) / (*groups)[j].size() > 0.5)
						candidates.push_back(j);
				}
				
				if (candidates.size() > 0)
				{
					// re fit ellipse using all combined arcs
					cv::Mat refitMat(6, 6, CV_64F, fitMatrixes + i * 36);
					double allcoverage = coverages[i];
					for (int j = 0; j < candidates.size(); j++)
					{
						cv::Mat refitMatj(6, 6, CV_64F, fitMatrixes + candidates[j] * 36);
						refitMat = refitMat + refitMatj;
						allcoverage += coverages[candidates[j]];
					}

					if (allcoverage >= M_PI)
					{
						for (int j = 0; j < candidates.size(); j++)
							isLabelled[j] = 1;
						isLabelled[i] = 1;
					}
					
					double S[36]; //拟合矩阵S
					memcpy(S, refitMat.data, sizeof(double) * 36);
					double Coefficients[6] = { 0,0,0,0,0,0 };// ax^2 + bxy + cy^2 + dx + ey + f = 0
					info = fitEllipse2(S, Coefficients);// ax^2 + bxy + cy^2 + dx + ey + f = 0, a > 0
					double param[5];
					if (info != 0)//拟合失败
					{
						ellipse2Param(Coefficients, param);// (x0,y0,a,b,phi)
						ellipara.x = param[0];//因为无论如何，都需要返回显著性强的椭圆
						ellipara.y = param[1];
						ellipara.a = param[2];
						ellipara.b = param[3];
						ellipara.phi = param[4];
					}
				}
#endif // USE_LABEL_MATCH
#endif

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
	for (i = 0; i < groupsNum - 1; i++)
	{
#if USE_AFFINE_WARP
		if (isLabelled[i] != 0)
			continue;
#endif
		for (j = i + 1; j < groupsNum; j++)
		{
#if USE_AFFINE_WARP
			if (isLabelled[j] != 0)
				continue;
#endif
			//加入极性判断,只提取指定极性的椭圆
			if (specified_polarity == 0 || (lines[(*groups)[i][0] * TUPLELENGTH + 7] == specified_polarity))
			{
				//group i 's polarity is the same as group j; and the number of two paired groups should be >= 3.
				if (lines[(*groups)[i][0] * TUPLELENGTH + 7] == lines[(*groups)[j][0] * TUPLELENGTH + 7] && ((*groups)[i].size() + (*groups)[j].size()) >= 3)
				{
					ind_start = (*groups)[i][0];//第i组的最开始一条线段索引
					ind_end = (*groups)[i][(*groups)[i].size() - 1];//第i组的最后一条线段索引
					pointG1s.x = lines[ind_start * TUPLELENGTH];
					pointG1s.y = lines[ind_start * TUPLELENGTH + 1];
					g1s_ls_dir.x = lines[ind_start * TUPLELENGTH + 4];
					g1s_ls_dir.y = lines[ind_start * TUPLELENGTH + 5];
					pointG1e.x = lines[ind_end * TUPLELENGTH + 2];
					pointG1e.y = lines[ind_end * TUPLELENGTH + 3];
					g1e_ls_dir.x = lines[ind_end * TUPLELENGTH + 4];
					g1e_ls_dir.y = lines[ind_end * TUPLELENGTH + 5];

					ind_start = (*groups)[j][0];//第j组的最开始一条线段索引
					ind_end = (*groups)[j][(*groups)[j].size() - 1];//第j组的最后一条线段索引
					pointG2s.x = lines[ind_start * TUPLELENGTH];
					pointG2s.y = lines[ind_start * TUPLELENGTH + 1];
					g2s_ls_dir.x = lines[ind_start * TUPLELENGTH + 4];
					g2s_ls_dir.y = lines[ind_start * TUPLELENGTH + 5];
					pointG2e.x = lines[ind_end * TUPLELENGTH + 2];
					pointG2e.y = lines[ind_end * TUPLELENGTH + 3];
					g2e_ls_dir.x = lines[ind_end * TUPLELENGTH + 4];
					g2e_ls_dir.y = lines[ind_end * TUPLELENGTH + 5];

					polarity = lines[ind_start * TUPLELENGTH + 7]; //i,j两组的极性
					if (regionLimitation(pointG1s, g1s_ls_dir, pointG1e, g1e_ls_dir, pointG2s, g2s_ls_dir, pointG2e, g2e_ls_dir, polarity, -3 * distance_tolerance))//都在彼此的线性区域内
					{
						//if ( i == 2)
						//  drawPairGroup(img,lines,(*groups),i,j);
#if USE_CNSCORE
						int ind_mid = (*groups)[i][(*groups)[i].size() / 2];
						point2d pointG1m = point2d(lines[ind_mid * TUPLELENGTH], lines[ind_mid * TUPLELENGTH]);// chose the mid line's start point
						ind_mid = (*groups)[j][(*groups)[j].size() / 2];
						point2d pointG2m = point2d(lines[ind_mid * TUPLELENGTH], lines[ind_mid * TUPLELENGTH]);
						// add cn score check
						if (checkCN(pointG1s, pointG1m, pointG1e, pointG2s, pointG2m, pointG2e) == -1)continue;
#endif

						if (calcEllipseParametersAndValidate(lines, line_num, groups, i, j, (fitMatrixes + i * 36), (fitMatrixes + j * 36), angles, distance_tolerance, supportInliersNum, &ellipara))//二次一般方程线性求解，线段的内点支持比例
						{

#if USE_AFFINE_WARP
							// check arcs
							float Phi = ellipara.phi;
							float cphi = cos(Phi);
							float sphi = sin(Phi);
							float a[3][3] = { { ellipara.a * cphi, -ellipara.b * sphi, ellipara.x }, \
							{ ellipara.a * sphi, ellipara.b * cphi, ellipara.y}, \
							{ 0, 0, 1} };
							// this is the affine matrix that warps points back to unit circle
							Mat Ainv = Mat(3, 3, CV_32F, a).inv();

#if USE_SPAN_CHECK
							// check coverage span
							// start point of first LSD
							float x1 = Ainv.at<float>(0, 0) * lines[(*groups)[i][0] * TUPLELENGTH] + Ainv.at<float>(0, 1) * lines[(*groups)[i][0] * TUPLELENGTH + 1] + Ainv.at<float>(0, 2);
							float y1 = Ainv.at<float>(1, 0) * lines[(*groups)[i][0] * TUPLELENGTH] + Ainv.at<float>(1, 1) * lines[(*groups)[i][0] * TUPLELENGTH + 1] + Ainv.at<float>(1, 2);
							int kk = (*groups)[i].size() - 1;
							// end point of last LSD
							float x2 = Ainv.at<float>(0, 0) * lines[(*groups)[i][kk] * TUPLELENGTH + 2] + Ainv.at<float>(0, 1) * lines[(*groups)[i][kk] * TUPLELENGTH + 3] + Ainv.at<float>(0, 2);
							float y2 = Ainv.at<float>(1, 0) * lines[(*groups)[i][kk] * TUPLELENGTH + 2] + Ainv.at<float>(1, 1) * lines[(*groups)[i][kk] * TUPLELENGTH + 3] + Ainv.at<float>(1, 2);
							float ang1 = atan2(y1, x1);
							ang1 = ang1 < 0 ? ang1 + M_PI * 2 : ang1;
							float ang2 = atan2(y2, x2);
							ang2 = ang2 < 0 ? ang2 + M_PI * 2 : ang2;
							float span1;
							if (lines[(*groups)[i][0] * TUPLELENGTH + 7] == -1)// counter-clockwise
							{
								ang2 = ang2 < ang1 ? ang2 + M_PI * 2 : ang2;
								span1 = ang2 - ang1;
							}
							else if (lines[(*groups)[i][0] * TUPLELENGTH + 7] == 1)// clockwise
							{
								ang1 = ang1 < ang2 ? ang1 + M_PI * 2 : ang1;
								span1 = ang1 - ang2;
							}

							x1 = Ainv.at<float>(0, 0) * lines[(*groups)[j][0] * TUPLELENGTH] + Ainv.at<float>(0, 1) * lines[(*groups)[j][0] * TUPLELENGTH + 1] + Ainv.at<float>(0, 2);
							y1 = Ainv.at<float>(1, 0) * lines[(*groups)[j][0] * TUPLELENGTH] + Ainv.at<float>(1, 1) * lines[(*groups)[j][0] * TUPLELENGTH + 1] + Ainv.at<float>(1, 2);
							kk = (*groups)[j].size() - 1;
							// end point of last LSD
							x2 = Ainv.at<float>(0, 0) * lines[(*groups)[j][kk] * TUPLELENGTH + 2] + Ainv.at<float>(0, 1) * lines[(*groups)[j][kk] * TUPLELENGTH + 3] + Ainv.at<float>(0, 2);
							y2 = Ainv.at<float>(1, 0) * lines[(*groups)[j][kk] * TUPLELENGTH + 2] + Ainv.at<float>(1, 1) * lines[(*groups)[j][kk] * TUPLELENGTH + 3] + Ainv.at<float>(1, 2);
							ang1 = atan2(y1, x1);
							ang1 = ang1 < 0 ? ang1 + M_PI * 2 : ang1;
							ang2 = atan2(y2, x2);
							ang2 = ang2 < 0 ? ang2 + M_PI * 2 : ang2;
							float span2;
							if (lines[(*groups)[j][0] * TUPLELENGTH + 7] == -1)// counter-clockwise
							{
								ang2 = ang2 < ang1 ? ang2 + M_PI * 2 : ang2;
								span2 = ang2 - ang1;
							}
							else if (lines[(*groups)[j][0] * TUPLELENGTH + 7] == 1)// clockwise
							{
								ang1 = ang1 < ang2 ? ang1 + M_PI * 2 : ang1;
								span2 = ang1 - ang2;
							}
							double allcoverage = coverages[i] + coverages[j];
							if ((allcoverage) / (M_PI * 2) < SPAN_CHECK_THRESH) continue;
							
#else
							double allcoverage = coverages[i] + coverages[j];
#endif

							

#if USE_LABEL_MATCH
							if (allcoverage >= M_4_9_PI)
							{
								//std::cout << Ainv << std::endl;
								//Mat M = (Mat_<double>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);
								std::vector<int> candidates;
								for (int k = 0; k < groupsNum; k++)
								{
									if (k == i || k == j) continue;
									if (lines[(*groups)[i][0] * TUPLELENGTH + 7] != lines[(*groups)[k][0] * TUPLELENGTH + 7]) continue;
									int validsum = 0;
									for (int kk = 0; kk<(*groups)[k].size(); kk++)
									{
										//std::cout << lines[(*groups)[j][k] * TUPLELENGTH] << ";" << lines[(*groups)[j][k] * TUPLELENGTH + 1] << std::endl;
										//std::cout << lines[(*groups)[j][k] * TUPLELENGTH+2] << ";" << lines[(*groups)[j][k] * TUPLELENGTH + 3] << std::endl;

										//每一条线段有2个端点
										float x1 = Ainv.at<float>(0, 0) * lines[(*groups)[k][kk] * TUPLELENGTH] + Ainv.at<float>(0, 1) * lines[(*groups)[k][kk] * TUPLELENGTH + 1] + Ainv.at<float>(0, 2);
										float y1 = Ainv.at<float>(1, 0) * lines[(*groups)[k][kk] * TUPLELENGTH] + Ainv.at<float>(1, 1) * lines[(*groups)[k][kk] * TUPLELENGTH + 1] + Ainv.at<float>(1, 2);

										float x2 = Ainv.at<float>(0, 0) * lines[(*groups)[k][kk] * TUPLELENGTH + 2] + Ainv.at<float>(0, 1) * lines[(*groups)[k][kk] * TUPLELENGTH + 3] + Ainv.at<float>(0, 2);
										float y2 = Ainv.at<float>(1, 0) * lines[(*groups)[k][kk] * TUPLELENGTH + 2] + Ainv.at<float>(1, 1) * lines[(*groups)[k][kk] * TUPLELENGTH + 3] + Ainv.at<float>(1, 2);

										if (fabs(x1*x1 + y1*y1 - 1) < WARP_RADIUS_THRESH && fabs(x2*x2 + y2*y2 - 1) < WARP_RADIUS_THRESH)
											validsum++;
									}
									if (float(validsum) / (*groups)[k].size() > 0.5)
										candidates.push_back(k);
								}
								if (candidates.size() > 0)
								{
									// re fit ellipse using all combined arcs
									cv::Mat refitMat1(6, 6, CV_64F, fitMatrixes + i * 36);
									cv::Mat refitMat2(6, 6, CV_64F, fitMatrixes + j * 36);
									cv::Mat refitMat = refitMat1 + refitMat2;

									for (int k = 0; k < candidates.size(); k++)
									{
										cv::Mat refitMatk(6, 6, CV_64F, fitMatrixes + candidates[k] * 36);
										refitMat = refitMat + refitMatk;
										allcoverage += coverages[candidates[k]];
									}

									if (allcoverage >= M_PI)
									{
										for (int k = 0; k < candidates.size(); k++)
											isLabelled[candidates[k]] = 1;
										isLabelled[i] = 1;
										isLabelled[j] = 1;
									}
									double S[36]; //拟合矩阵S
									memcpy(S, refitMat.data, sizeof(double) * 36);
									double Coefficients[6] = { 0,0,0,0,0,0 };// ax^2 + bxy + cy^2 + dx + ey + f = 0
									info = fitEllipse2(S, Coefficients);// ax^2 + bxy + cy^2 + dx + ey + f = 0, a > 0
									double param[5];
									if (info != 0)//拟合失败
									{
										ellipse2Param(Coefficients, param);// (x0,y0,a,b,phi)
										ellipara.x = param[0];//因为无论如何，都需要返回显著性强的椭圆
										ellipara.y = param[1];
										ellipara.a = param[2];
										ellipara.b = param[3];
										ellipara.phi = param[4];
									}
								}
							}
#endif
#endif
							PairGroupNode * node = (PairGroupNode*)malloc(sizeof(PairGroupNode));
							node->center.x = ellipara.x;
							node->center.y = ellipara.y;
							node->axis.x = ellipara.a;
							node->axis.y = ellipara.b;
							node->phi = ellipara.phi;
							node->pairGroupInd.x = i;
							node->pairGroupInd.y = -1;//无 // bug????
							if (head != NULL)
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
            pairGroupList->pairGroup[i].pairGroupInd.y = p->pairGroupInd.y;// never group pair ij, seems a bug???
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
