#pragma once
#include "geometry.h"
#include <CGAL\Triangle_3.h>
#include <algorithm>
#include <CGAL\Constrained_Delaunay_triangulation_2.h>
#include <CGAL\Triangulation_face_base_with_info_2.h>

namespace BI
{

struct tri_facet
{
	int pn[3];
	FT tag;
	tri_facet (const int &n1, const int &n2, const int &n3,const FT& tag_)
		:tag(tag_)
	{pn[0] = n1; pn[1] = n2; pn[2] = n3; std::sort(pn,pn+3);}
};

inline bool operator== (const tri_facet &tf1, const tri_facet &tf2)
{
	if (tf1.pn[0] == tf2.pn[0] && tf1.pn[1] == tf2.pn[1] && tf1.pn[2] == tf2.pn[2])
		return true;
	if (tf1.pn[0] == tf2.pn[1] && tf1.pn[1] == tf2.pn[2] && tf1.pn[2] == tf2.pn[0])
		return true;
	if (tf1.pn[0] == tf2.pn[2] && tf1.pn[1] == tf2.pn[0] && tf1.pn[2] == tf2.pn[1])
		return true;
	if (tf1.pn[0] == tf2.pn[0] && tf1.pn[1] == tf2.pn[2] && tf1.pn[2] == tf2.pn[1])
		return true;
	if (tf1.pn[0] == tf2.pn[2] && tf1.pn[1] == tf2.pn[1] && tf1.pn[2] == tf2.pn[0])
		return true;
	if (tf1.pn[0] == tf2.pn[1] && tf1.pn[1] == tf2.pn[0] && tf1.pn[2] == tf2.pn[2])
		return true;
	return false;
}

class curve_frac
{
	struct FaceInfo2
	{
		int nesting_level;
		FaceInfo2 () {nesting_level=0;}
		bool in_domain ()	//nesting_level等于0的时候表示初始状态，等于1的时候表示被mark_domain这个函数标记为了外围三角形
		{
			if (nesting_level !=1) return true;
			else return false;
		}
	};

	typedef CGAL::Exact_intersections_tag					EXACT_TAG;
	typedef CGAL::Triangulation_vertex_base_2<K>			VB;
	typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,K> Fbb;
	typedef CGAL::Constrained_triangulation_face_base_2<K,Fbb>	FB;
	typedef CGAL::Triangulation_data_structure_2<VB,FB>		TDS;
	typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, EXACT_TAG>	CDT;

	void mark_domains (CDT &cdt);

	typedef K::Triangle_3 triangle_3;
	std::vector <point_3> vp;
	std::vector <tri_facet> vtri;
public:
	void swap (curve_frac & c) 
	{vp.swap (c.vp); vtri.swap (c.vtri);}
	void clear ()
	{vp.clear(); vtri.clear();}
	int read_STL (std::istream &infile,const FT & tag);

	//在new_vtri中添加三角形，不进行查重。新的顶点直接添加到this->vp中。
	void insert_tri (std::vector <tri_facet> &new_vtri, const triangle_3 &t, const FT& tag);		
	void insert_tri (std::vector <tri_facet> &new_vtri, const point_3 &p1, const point_3 &p2, const point_3 &p3, const FT& tag);
	void insert_tri (std::vector <tri_facet> &new_vtri, const CDT &cdt, const FT& tag,const plane_3 &pl);

	//在new_vtri中添加三角形，进行查重。新的顶点直接添加到this->vp中。
	void join_tri (std::vector <tri_facet> &new_vtri, const triangle_3 &t, const FT &tag);			//添加三角形的时候会检查是否和已存在的三角形重复
	void join_tri (std::vector <tri_facet> &new_vtri, const point_3 &p1, const point_3 &p2, const point_3 &p3, const FT & tag);
	void join_tri (std::vector <tri_facet> &new_vtri, const CDT &cdt, const FT& tag,const plane_3 &pl);

	//将this和给定的c做相交运算。相交的三角形进行remesh
	//返回0才算正确返回
	int join (const curve_frac & c);

	int output_STL (std::ostream &outfile);

private:
	typedef CGAL::cpp11::result_of <K::Intersect_3(triangle_3,triangle_3)>::type intersect_type;

	bool is_overlap (const intersect_type &it) const;		//其实不用实现成成员函数的

	//给cdt添加inter表示的约束，inter表示的元素可能是point_3, segment_3, triangle_3, vector<point_3>。这些元素在pl上投影后添加到cdt中
	//如果函数inter表示的元素是triangle_3或者是 vector<point_3>，那么就返回true，表示有共面三角形。否则返回false，表示没有共面三角形
	bool add_constraint (CDT &cdt, const intersect_type &inter, const plane_3 &pl);
	//这个函数的意义同上，只不过是由一个约束变成了一组约束.
	//只要这一组相交元素中有一个相交元素是triangle_3或者vector<point_3>，那么就返回true，表示有共面三角形。否则返回false，表示没有共面三角形
	bool add_constraint (CDT &cdt, const std::vector <intersect_type> &vinter, const plane_3 &pl);
	bool add_constraint (CDT &cdt, const std::vector <intersect_type> &vinter, const std::vector <int> &vindex, const plane_3 &pl);


	//添加t这个三角形到cdt中，作为其中的约束，其中pl是t所在的平面，要利用pl.to_2d进行转换
	void add_constraint (CDT &cdt, const point_3 &p1, const point_3 &p2, const point_3 &p3, const plane_3 &pl);
	void add_constraint (CDT &cdt, const triangle_3 &t, const plane_3 &pl);
	void add_constraint (CDT &cdt, const std::vector<point_3> &vp, const tri_facet &f, const plane_3 &pl);
	
};

}