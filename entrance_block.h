#pragma once
#include "geometry.h"
#include "identification.h"
#include "outer_poly.h"
#include <vector>


namespace BI
{

class outer_bound;

struct angle_edge
{
	struct sub_angle
	{
		//每个子角表示成一系列半空间的交集。子角都是凸集
		std::vector <vector_3> vnormal;		//其中的每个向量表示一个半空间,每个向量vnormal[i]指向子角的外部
		//返回false这个子角存在（内点的集合不为空），返回true表示这个子角不存在（内点集为空）
		bool is_degenerate () const 
		{
			std::vector <plane_3> vpl;
			vpl.reserve (vnormal.size ());
			for (int i(0);i!=(int)vnormal.size ();++i)
				vpl.push_back (plane_3(vnormal[i].x(), vnormal[i].y(), vnormal[i].z(),0));
			switch (BI::intersection(vpl))
			{
			case Intersection_Status::Inter_Boundary:
				return true;
			case Intersection_Status::Inter_Unbounded:
				return false;
			default:
				throw std::logic_error ("bool angle_edge::sub_angle::is_degenerate () const, LP results error");
			}
		}
		//计算这个子角和v所表示的半空间是否相交，其中v指向半空间外部！外部！！
		//判断这个子角和给定的半空间是否内部相交(顶点相同时)
		bool is_internal_intersected (const vector_3 &v) const;
		//判断这个子角和给定的子角sa是否内部相交(顶点相同时)
		bool is_internal_intersected (const sub_angle &sa) const;
		
		//判断这个角和给定的角sa是否相交
		//其中p_this表示this这个角的顶点，p_sa表示给定角sa的顶点
		Intersection_Status intersection (const sub_angle&sa, const point_3 &p_this, const point_3 &p_sa) const;
		
	};

	struct angle
	{
		int e;						//这个角的顶点，所有点用索引表示
		std::vector <sub_angle> vsa;	//每个角被表示为若干子角的交集

		//判断这个角和给定的角a是否内部相交(顶点相同时)
		bool is_internal_intersected (const angle &a) const;
		//判断这个角和给定的半空间是否内部相交,v指向半空间外部
		bool is_internal_intersected (const vector_3 &v) const;

		//判断这个角和给定的角a是否相交
		//其中p_this表示this这个角的顶点，p_a表示给定角a的顶点
		Intersection_Status intersection (const angle &a, const point_3&p_this, const point_3&p_a) const;
	};

	struct sub_edge
	{
		vector_3 n1, n2;		//表示一个凸子边，由n1和n2表示的半空间的交集所确定。其中n1和n2指向子边的外部，注意，是外部！！
								//如果n1==n2，则表示这个子边是一个半空间
		//判断两个子边是否内部相交，让这两个子边中的半空间都过原点，然后看这些半空间是否有内部交集。
		//判断这个子边和给定的子边se是否内部相交(顶点相同时)
		bool is_internal_intersected (const sub_edge& se) const;

		//判断这个边和给定的边se是否相交
		//其中p_this表示this这个边棱上的某个点，p_se表示给定边se棱上的某个点
		Intersection_Status intersection (const sub_edge&se, const point_3 &p_this, const point_3 &p_se) const;
	};

	struct edge
	{
		int e1,e2;		//表示这条边的棱，点都由索引表示。并且要保证e1<e2
		std::vector<sub_edge> vse;		//表示这条边的子边，注意vse中的法向向量一定垂直于向量e1e2
		edge (const int&e1_, const int&e2_)
		{
			if (e1_<e2_)
			{e1 = e1_; e2 = e2_;}
			else
			{e1 = e2_; e2 = e1_;}
		}

		//判断这个边和给定的边e是否内部相交(顶点相同时)
		bool is_internal_intersected (const edge &e) const;

		//判断两个边是否内部相交。
		//判断这个边和给定的边e是否相交
		//其中p_this表示this这个边棱上的某个点，p_e表示给定边e棱上的某个点
		Intersection_Status intersection (const edge& e, const point_3& p_this, const point_3 &p_e) const;
	};

	std::vector<point_3> vp;
	std::vector<angle> va;
	std::vector<edge> ve;
};

//struct cover
//{
//	enum TYPE {Edge2Edge, Vertex2Poly};
//	int ia,ib;
//
//};

inline bool operator== (const angle_edge::edge &edge1, const angle_edge::edge &edge2)
{
	if (edge1.e1 == edge2.e1 && edge1.e2 == edge2.e2)
		return true;
	else
		return false;
}

class entrance_block
{
public:
	std::vector<poly_frac> vcover;		//每个裂隙是一个覆盖而已 vector of cover
	outer_bound outer_b;						//生成外表面多边形作为进入块体
	block_system bs;
	cblock inner;		//进入块体内部的点被表示为一个单元块体的组合
private:
	int gen_cover (const block_system &bs, const int &cbi, const int &cbj);		//读取复合块体的形状，初始化覆盖vcov
	int gen_cover (const outer_bound &oba, const outer_bound &obb, const point_3 &a0);				//用outer_bound表示A和B,要求A和B都是已经regulated之后的
	int gen_cover (const outer_bound &ob_a, const angle_edge &ae_a, const outer_bound &ob_b, const angle_edge &ae_b, const point_3 &a0 = point_3 (0,0,0));
public:
	int gen_entrance_block (const block_system &bs_, const cblock& cba, const cblock& cbb, const point_3 &a0, bool debug_info = false);		//接口设计的不是很好，有待改进
	int gen_entrance_block (const block_system &bs_, const cblock& cba, const cblock& cbb, const point_3 &a0,
							const outer_bound &ob_a, const angle_edge &ae_a, const outer_bound &ob_b, const angle_edge &ae_b, bool debug_info = false);
	
	friend std::ostream& operator<< (std::ostream& os, const entrance_block& eb);
	friend std::istream& operator>> (std::istream& is, entrance_block& eb);
	friend bool operator== (const entrance_block& eb1, const entrance_block& eb2);
};

int outer_bound2angle_edge (const block_system & bs, const cblock & cb,const outer_bound &ob, angle_edge &ae);

inline bool operator== (const entrance_block& eb1, const entrance_block& eb2)
{
	return (eb1.vcover == eb2.vcover &&
		eb1.outer_b == eb2.outer_b &&
		eb1.bs == eb2.bs &&
		eb1.inner == eb2.inner );
}

//判断一个点和一个单元块体的位置关系,vfi是一个返回值，表示如果这个点在该单元块体上，则记录这个点位于哪几个平面上
Oriented_side point_eblock_location (const block_system &bs, const eblock &eb, const point_3&tar, std::vector<int>& vfi);


}