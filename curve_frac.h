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
		bool in_domain ()	//nesting_level����0��ʱ���ʾ��ʼ״̬������1��ʱ���ʾ��mark_domain����������Ϊ����Χ������
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

	//��new_vtri����������Σ������в��ء��µĶ���ֱ����ӵ�this->vp�С�
	void insert_tri (std::vector <tri_facet> &new_vtri, const triangle_3 &t, const FT& tag);		
	void insert_tri (std::vector <tri_facet> &new_vtri, const point_3 &p1, const point_3 &p2, const point_3 &p3, const FT& tag);
	void insert_tri (std::vector <tri_facet> &new_vtri, const CDT &cdt, const FT& tag,const plane_3 &pl);

	//��new_vtri����������Σ����в��ء��µĶ���ֱ����ӵ�this->vp�С�
	void join_tri (std::vector <tri_facet> &new_vtri, const triangle_3 &t, const FT &tag);			//��������ε�ʱ������Ƿ���Ѵ��ڵ��������ظ�
	void join_tri (std::vector <tri_facet> &new_vtri, const point_3 &p1, const point_3 &p2, const point_3 &p3, const FT & tag);
	void join_tri (std::vector <tri_facet> &new_vtri, const CDT &cdt, const FT& tag,const plane_3 &pl);

	//��this�͸�����c���ཻ���㡣�ཻ�������ν���remesh
	//����0������ȷ����
	int join (const curve_frac & c);

	int output_STL (std::ostream &outfile);

private:
	typedef CGAL::cpp11::result_of <K::Intersect_3(triangle_3,triangle_3)>::type intersect_type;

	bool is_overlap (const intersect_type &it) const;		//��ʵ����ʵ�ֳɳ�Ա������

	//��cdt���inter��ʾ��Լ����inter��ʾ��Ԫ�ؿ�����point_3, segment_3, triangle_3, vector<point_3>����ЩԪ����pl��ͶӰ����ӵ�cdt��
	//�������inter��ʾ��Ԫ����triangle_3������ vector<point_3>����ô�ͷ���true����ʾ�й��������Ρ����򷵻�false����ʾû�й���������
	bool add_constraint (CDT &cdt, const intersect_type &inter, const plane_3 &pl);
	//�������������ͬ�ϣ�ֻ��������һ��Լ�������һ��Լ��.
	//ֻҪ��һ���ཻԪ������һ���ཻԪ����triangle_3����vector<point_3>����ô�ͷ���true����ʾ�й��������Ρ����򷵻�false����ʾû�й���������
	bool add_constraint (CDT &cdt, const std::vector <intersect_type> &vinter, const plane_3 &pl);
	bool add_constraint (CDT &cdt, const std::vector <intersect_type> &vinter, const std::vector <int> &vindex, const plane_3 &pl);


	//���t��������ε�cdt�У���Ϊ���е�Լ��������pl��t���ڵ�ƽ�棬Ҫ����pl.to_2d����ת��
	void add_constraint (CDT &cdt, const point_3 &p1, const point_3 &p2, const point_3 &p3, const plane_3 &pl);
	void add_constraint (CDT &cdt, const triangle_3 &t, const plane_3 &pl);
	void add_constraint (CDT &cdt, const std::vector<point_3> &vp, const tri_facet &f, const plane_3 &pl);
	
};

}