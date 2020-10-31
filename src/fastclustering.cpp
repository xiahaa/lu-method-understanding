//
//  fastclustering.cpp
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
#include "nanoflann.hpp"

template <typename T>
struct Point4
{
	T  x, y, a, b;

	Point4(void) :x(0), y(0), a(0), b(0) {}
	Point4(T a1, T a2, T a3, T a4) :x(a1), y(a2), a(a3), b(a4) {}

	T operator()(const int dim)
	{
		T ret_val = 0;
		switch (dim)
		{
		case 0: ret_val = x; break;
		case 1: ret_val = y; break;
		case 2: ret_val = a; break;
		case 3: ret_val = b; break;
		}
		return ret_val;
	}

	Point4& operator+=(const Point4& rhs)
	{
		this->a += rhs.a;
		this->b += rhs.b;
		this->x += rhs.x;
		this->y += rhs.y;

		return *this;
	}
	Point4& operator/=(size_t num)
	{
		this->a /= num;
		this->b /= num;
		this->x /= num;
		this->y /= num;

		return *this;
	}
	Point4& operator=(const Point4& rhs)
	{
		this->a = rhs.a;
		this->b = rhs.b;
		this->x = rhs.x;
		this->y = rhs.y;

		return *this;
	}
	Point4& operator*(const double s)
	{
		Point4 ret;
		ret.a = this->a*s;
		ret.b = this->b*s;
		ret.x = this->x*s;
		ret.y = this->y*s;

		return ret;
	}
};

template <typename T>
struct PointCloud_XYAB
{
	std::vector<Point4<T> >  pts;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline T kdtree_get_pt(const size_t idx, const size_t dim) const
	{
		if (dim == 0) return pts[idx].x;
		else if (dim == 1) return pts[idx].y;
		else if (dim == 2) return pts[idx].a;
		else return pts[idx].b;
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};

// construct a kd-tree index:
typedef nanoflann::KDTreeSingleIndexAdaptor<
	nanoflann::L2_Simple_Adaptor<double, PointCloud_XYAB<double> >,
	PointCloud_XYAB<double>,
	4 /* dim */
> my_4d_tree_t;

inline double squaredDifference(int & nDims, double *& points, int & i, double *& initializations, int & j)
{
	double result = 0;
	for (int k = 0; k < nDims; ++k)
		result += pow(points[i*nDims + k] - initializations[j*nDims + k], 2);
	return result;
}

int kdSearchRadius4(void *index,
	const double *query_pt, const double search_radius,
	const nanoflann::SearchParams &params,
	std::vector<std::pair<size_t, double> >   &ret_matches)
{
	my_4d_tree_t *index4 = (my_4d_tree_t *)index;
	int nMatches = index4->radiusSearch(query_pt, search_radius, ret_matches, params);
	return nMatches;
}

typedef int(*fkdpointer)(void *, const double *, const double, const nanoflann::SearchParams &, std::vector<std::pair<size_t, double> > &);
void kdmeanShift(double *points, void *index, fkdpointer fhandler, int nDims, double * initPoints, int initLength,
	double sigma, double window_size, double accuracy_tolerance, int iter_times)
{
	int nQuerries = initLength;
	double * initializations = (double*)malloc(nQuerries * nDims * sizeof(double));
	memcpy(initializations, initPoints, nQuerries * nDims * sizeof(double));//copy

	double sigma2 = sigma*sigma;//sigma平方
	double radius2 = window_size *window_size;//平方
	double tolerance = accuracy_tolerance;
	int iters, maxiters = iter_times;//最大迭代次数
									 //返回与初始搜索点集一样大小的最终定位点集
	double * finals = (double*)malloc(nQuerries * nDims * sizeof(double));;//最终定位点集的指针
	memcpy(finals, initializations, nQuerries * nDims * sizeof(double));

	const double search_radius = static_cast<double>(window_size);
	nanoflann::SearchParams params;

	double *query_pt = (double *)malloc(sizeof(double) * nDims);

	for (int loop = 0; loop < nQuerries; ++loop)
	{
		iters = 0;
		while (iters < maxiters)
		{
			bool flag = false;
			double denominator = 0;//分母

			std::vector<std::pair<size_t, double> >   ret_matches;

			for (int i = 0; i < nDims; i++)
				query_pt[i] = initializations[loop*nDims + i];
			//const double query_pt[4] = { initializations[loop*nDims + 0], initializations[loop*nDims + 1], 
			// initializations[loop*nDims + 2], initializations[loop*nDims + 3] };
			
			const int nMatches = (*fhandler)(index, query_pt, search_radius, params, ret_matches);		
			
			//const size_t nMatches = index.radiusSearch(&query_pt[0], search_radius, ret_matches, params);
			if (nMatches > 0)
			{
				flag = true;
				for (size_t j = 0; j < nMatches;j++)
					denominator += exp(-(ret_matches[j].second) / sigma2);// todo, check if the radius is squared or square root.
			}
				
			if (!flag)
				break;

			for (int j = 0; j < nDims; ++j)
				finals[loop*nDims + j] = 0;//对最终定位点集中的第loop个点的向量赋值为0
			for (size_t j = 0; j < nMatches; j++)
			{
				int i = ret_matches[j].first;
				double distance = ret_matches[j].second;
				for (int j = 0; j < nDims; ++j)//每个内点向量的以一定权值累加
					finals[loop*nDims + j] += exp(-distance / sigma2) * points[i*nDims + j];
			}
			for (int j = 0; j < nDims; ++j)//权值归一化
				finals[loop*nDims + j] /= denominator;

			if (sqrt(squaredDifference(nDims, finals, loop, initializations, loop)) < tolerance)//相继两次的迭代中心在误差内了，则认为已经收敛，没必要再继续迭代
				break;
			iters = iters + 1;
			for (int j = 0; j < nDims; ++j)//更新迭代的搜索中心
				initializations[loop*nDims + j] = finals[loop*nDims + j];
		}
	}
	memcpy(initPoints, finals, nQuerries * nDims * sizeof(double));
	free(initializations);
	free(finals);
}

int  cluster4DPoints(double *points,
	my_4d_tree_t   &index, int points_num, double distance_tolerance, double * & centers, int * centers_num)
{
	int i;
	const double search_radius = static_cast<double>(distance_tolerance);
	nanoflann::SearchParams params;

	std::vector<unsigned char> label(points_num,0);
	
	std::vector<Point4<double> > xyabcenters(points_num);
	int valid_centers = 0;
	for (i = 0; i < points_num; i++)
	{
		if (label[i] == 1)continue;
		label[i] = 1;
		std::vector<std::pair<size_t, double> >   ret_matches;
		
		const double query_pt[4] = { points[i*4], points[i * 4 + 1], points[i * 4 + 2], points[i * 4 + 3] };
		const size_t nMatches = index.radiusSearch(&query_pt[0], search_radius, ret_matches, params);

		if (nMatches == 0) continue;
		Point4<double> xyabcenter;
		for (size_t j = 0; j < nMatches; j++)
		{
			label[ret_matches[j].first] = 1;
			int ii = ret_matches[j].first;
			xyabcenter += Point4<double>(points[ii * 4], points[ii * 4 + 1], points[ii * 4 + 2], points[ii * 4 + 3]);
		}
		xyabcenter /= nMatches;
		xyabcenters[valid_centers] = xyabcenter;
		valid_centers++;
	}
	int initCentersLength = valid_centers;
	if (initCentersLength == 0)
	{
		(*centers_num) = 0;
		centers = NULL;
		return 0;
	}

	double * initCenters; //initCentersLength x 4
	initCenters = (double*)malloc(sizeof(double)*initCentersLength * 4);
	for (i = 0; i<initCentersLength; i++)// initCenters 大小是 initCentersLength*2
	{
		int addr = 4 * i;
		initCenters[addr] = xyabcenters[i].x;
		initCenters[addr + 1] = xyabcenters[i].y;
		initCenters[addr + 2] = xyabcenters[i].a;
		initCenters[addr + 3] = xyabcenters[i].b;
	}

	fkdpointer f4dsearch = &kdSearchRadius4;
	kdmeanShift(points, (void *)(&index), f4dsearch, 4, initCenters, initCentersLength, 1, distance_tolerance, 1e-6, 50);//迭代20次
						
	//聚类
	//千万要注意centers_num是int型指针，++--时要(*centers_num).
	clusterByDistance(initCenters, initCentersLength, 4, distance_tolerance / 2, 40, centers, centers_num);//此处控制参数要改，要调节

	if ((*centers_num) <= 0)//可无
	{
		return 0;  //不懂为什么，聚类中心的周围确没有最靠近它的点
				   //system("pause");
				   //error("cluster2DPoints,(*centers_num)<=0");
	}
	free(initCenters);
	return 1;
}

void fastgenerateEllipseCandidates(PairGroupList * pairGroupList, double distance_tolerance, double * & ellipse_candidates, int * candidates_num)
{
	if (pairGroupList->length <= 0)//检测，至少要有1个样本用来产生候选
	{
		ellipse_candidates = NULL;
		(*candidates_num) = 0;
		return;
	}

	double *xyabes; // dims that use euclidean distance
	int xyab_num;
	double *phis; // dim that cannot use euclidean distance
	int phi_num;

	double * buffer_xyab = (double*)calloc(pairGroupList->length * 4, sizeof(double));
	double * bufferPhi = (double*)calloc(pairGroupList->length, sizeof(double));

	point2i * bufferIndexes = (point2i *)calloc(pairGroupList->length, sizeof(point2i));//point[i].x记录第i个分类在bufferXX中的起始索引位置，point[i].y记录第i个分类在bufferXX中的长度
	int     * buffer_temp = (int*)calloc(pairGroupList->length, sizeof(int));
	
	int addr, info, ind;

	if (buffer_xyab == NULL || bufferPhi == NULL || bufferIndexes == NULL || buffer_temp == NULL)
	{
		ellipse_candidates = NULL;
		(*candidates_num) = 0;
		error("generateEllipseCandidates, not enough memory");
	}
	(*candidates_num) = 0; //候选椭圆数量，初始化为0,非常重要

	PointCloud_XYAB<double> cloud;
	cloud.pts.resize(pairGroupList->length);
	for (int i = 0; i<pairGroupList->length; i++)
	{
		cloud.pts[i].x = pairGroupList->pairGroup[i].center.x;
		cloud.pts[i].y = pairGroupList->pairGroup[i].center.y;
		cloud.pts[i].a = pairGroupList->pairGroup[i].axis.x;
		cloud.pts[i].b = pairGroupList->pairGroup[i].axis.y;
		addr = 4 * i;
		buffer_xyab[addr] = pairGroupList->pairGroup[i].center.x;
		buffer_xyab[addr + 1] = pairGroupList->pairGroup[i].center.y;
		buffer_xyab[addr + 2] = pairGroupList->pairGroup[i].axis.x;
		buffer_xyab[addr + 3] = pairGroupList->pairGroup[i].axis.y;
	}
	
	my_4d_tree_t   index(4 /*dim*/, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	index.buildIndex();

	//cluster the ellipses' centers
	info = cluster4DPoints(buffer_xyab, index, pairGroupList->length, distance_tolerance, xyabes, &xyab_num);

	if (info == 0)
	{
		ellipse_candidates = NULL;
		(*candidates_num) = 0;
		error("generateEllipseCandidates, cluster2DPoints, error in clustering elliptic centers");
	}
	
	// cloud with xyabs
	PointCloud_XYAB<double> newcloud;
	newcloud.pts.resize(xyab_num);
	for (int i = 0; i<xyab_num; i++)
	{
		int addr = 4 * i;
		newcloud.pts[i].x = xyabes[addr];
		newcloud.pts[i].y = xyabes[addr+1];
		newcloud.pts[i].a = xyabes[addr+2];
		newcloud.pts[i].b = xyabes[addr+3];
	}

	my_4d_tree_t   newindex(4 /*dim*/, newcloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	newindex.buildIndex();
	
	//classification,寻找每个点归属的聚类中心
	for (int i = 0; i < pairGroupList->length; i++)
	{
		size_t num_results = 1;
		std::vector<size_t>   ret_index(num_results);
		std::vector<double> out_dist_sqr(num_results);
		const double *query_pt = &buffer_xyab[4*i];
		num_results = newindex.knnSearch(query_pt, num_results, &ret_index[0], &out_dist_sqr[0]);
		buffer_temp[i] = ret_index[0]; 
	}
	
	//将分类结果按顺序存到bufferXY,bufferPhi,bufferAB中，且bufferIndexes[i]存着第i个聚类中心的起始索引位置和长度
	memset(bufferIndexes, 0, sizeof(point2i)*pairGroupList->length);
	ind = 0;//清零，样本点起始位置，索引位置是ind*2,分区的基址
	for (int i = 0; i<xyab_num; i++)
	{
		bufferIndexes[i].x = ind;
		for (int j = 0; j<pairGroupList->length; j++)
		{
			if (buffer_temp[j] == i)// if this ellipse belongs to this center
			{
				bufferPhi[ind + bufferIndexes[i].y] = pairGroupList->pairGroup[j].phi;
				bufferIndexes[i].y++;//第i个聚类中心周围的点数量加1
			}
		}
		if (bufferIndexes[i].y == 0)//聚类中心周围没有靠近的点
		{
			error("generateEllipseCandidates, no XY points near to the clustering center");
		}
		ind += bufferIndexes[i].y;
	}
	
	//对每一个椭圆中心的周围的点进行倾角聚类
	//第i个椭圆聚类中心，其邻近点的索引范围是：bufferIndexs[i].x ~ (bufferIndex[i].x + bufferIndex[i].y-1)
	for (int i = 0; i<xyab_num; i++)
	{
		double * phi_pointer_temp = bufferPhi + bufferIndexes[i].x;//倾角指针
		info = cluster1DDatas(phi_pointer_temp, bufferIndexes[i].y, 0.0873, phis, &phi_num);//对phi聚类, pi/180*5 = 0.0873, 5°误差
		if (info == 0) //不懂为什么，聚类中心centers[i]的周围可能没有最靠近它的点,数量bufferIndexes[i].y = 0
		{
			continue;
		}

		for (int j = 0; j<phi_num; j++)
		{
			addr = (*candidates_num) * 4;
			buffer_xyab[addr] = xyabes[i * 4];
			buffer_xyab[addr+1] = xyabes[i * 4+1];
			buffer_xyab[addr+2] = xyabes[i * 4+2];
			buffer_xyab[addr+3] = xyabes[i * 4+3];
			bufferPhi[(*candidates_num)] = phis[j];
			(*candidates_num)++;
		}
		free(phis);//cluster1DDatas严格要求，用完phis后，需要释放函数内部申请的内存
	}
	free(xyabes);//cluster2DPoints严格要求，用完centers后，需要释放函数内部申请的内存
	//释放在函数开头申请的部分内存
	free(buffer_temp); //此处释放出问题
	free(bufferIndexes);
	ellipse_candidates = (double*)malloc(sizeof(double)*(*candidates_num) * 5);
	for (int i = 0; i < (*candidates_num); i++)
	{
		addr = 4 * i;
		ellipse_candidates[i * 5] =		buffer_xyab[addr];
		ellipse_candidates[i * 5 + 1] = buffer_xyab[addr + 1];
		ellipse_candidates[i * 5 + 2] = buffer_xyab[addr + 2];
		ellipse_candidates[i * 5 + 3] = buffer_xyab[addr + 3];
		ellipse_candidates[i * 5 + 4] = bufferPhi[i];
	}
	//释放在函数开头申请的内存
	free(bufferPhi);
	free(buffer_xyab);
	if ((*candidates_num) <= 0)
	{
		*candidates_num = 0;
		ellipse_candidates = NULL;
		//cout<<"no any candidates generated!"<<endl;
	}
}

