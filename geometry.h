#pragma once
#include <cmath>
#include <vector>
#include <list>
#include <CGAL\Exact_predicates_exact_constructions_kernel.h>
#include <CGAL\Boolean_set_operations_2.h>
#include <CGAL\Polygon_2.h>
#include <CGAL\Polygon_set_2.h>
#include <CGAL\Polygon_with_holes_2.h>
#include <CGAL\Bbox_2.h>
#include <CGAL\Gmpq.h>
#include <CGAL\Triangle_3.h>
#include <CGAL\Tetrahedron_3.h>
#include <iostream>
#define M_PI_2 6.28318530717958647693
#define M_PI 3.14159265358979323846

namespace BI	//block identification的简写
{
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2				point_2;
typedef K::Point_3				point_3;
typedef K::Plane_3				plane_3;
typedef K::Circle_3				circle_3;
typedef K::Vector_3				vector_3;
typedef K::Segment_3			segment_3;
typedef K::Triangle_3			triangle_3;
typedef K::Tetrahedron_3		tetrahedron_3;
typedef CGAL::Polygon_2<K>		polygon_2;
typedef CGAL::Polygon_set_2<K>	polygon_set_2;
typedef CGAL::Polygon_with_holes_2<K>	polygon_with_holes_2;
typedef CGAL::Bbox_2	mbox_2;
typedef CGAL::Bbox_3	mbox_3;
typedef K::FT		FT;
typedef int			INT;
typedef K::Aff_transformation_3		aff_tm_3;

vector_3 dip2normal (const FT& dip_dir_, const FT& dip_,const FT& dx_);
template <typename K>
inline int add_set (std::vector<K> &p, const K &po);
enum Oriented_side {IN_, OUT_, ON_};		//标记
enum Intersection_Status {Separate, Inter_Boundary, Inter_Internal, Inter_Unbounded};
class block_system;

class frac_base
{
protected:
	frac_base ()
		: frac_id(-1),
		  plane_id(-1),
		  cohesion(0),
		  f_angle(0),
		  aperture(0){}

public:
	INT		plane_id;			//裂隙所在平面的ID
	INT		frac_id;			//裂隙的ID,正数
	FT		cohesion;
	FT		f_angle;
	FT		aperture;
};

class disc_frac: public frac_base
{
public:
	point_3 center;
	FT		r;

	//序列化输入输出函数
	friend std::ostream& operator<< (std::ostream &os, const disc_frac &df);
	friend std::istream& operator>> (std::istream &is, disc_frac& df);
	friend bool operator== (const disc_frac &df1, const disc_frac&df2);
};

class poly_frac: public frac_base
{
public:
	typedef std::vector <point_3> v_loop;
	v_loop vp_ob;		//表示多边形裂隙外表面。构造的时候要为每一个点取一个projection.识别块体的时候旋转方向并不重要
	std::vector <v_loop> vp_ib;		//表示多边形裂隙内表面
	
	void proj_plane (const plane_3& pl)		//将这个裂隙上的点都投影到pl上
	{
		for (int i(0);i!=(int)vp_ob.size ();++i)
			vp_ob[i] = pl.projection(vp_ob[i]);
		for (int i(0);i!=(int)vp_ib.size ();++i)
			for (int j(0);j!=(int)vp_ib[i].size ();++j)
				vp_ib[i][j] = pl.projection (vp_ib[i][j]);
	}

	//序列化输入输出函数
	friend std::ostream& operator<< (std::ostream &os, const poly_frac &pf);
	friend std::istream& operator>> (std::istream &is, poly_frac& pf);
	friend bool operator== (const poly_frac &pf1, const poly_frac &pf2);
};

class eblock
{
public:
	struct face
	{
		std::vector<INT> np;	//顶点序列应该逆时针指向外
		INT frac_id;			//形成这个面的裂隙的ID,如果frac_id是负数，则表示该面是输入的初始子块体的某个面，并且如果frac_id<=-4，则该面不用和其他面合并。如果frac_id>=0时，则该面是由某个裂隙生成的，合并块体的时候应该合并。如果-3<=frac_id<0时，该面应该合并。
		INT plane_id;			//为每一个面记录一个平面方程
		bool outer_normal;		//true：该面的外法线方向为plane_id所确定面的正方向（positive size）
		void swap (face &f)
		{
			INT frac_id_t = frac_id;
			INT plane_id_t = plane_id;
			bool outer_normal_t = outer_normal;
			frac_id = f.frac_id; plane_id = f.plane_id; outer_normal = f.outer_normal;
			np.swap(f.np);
			f.frac_id = frac_id_t; f.plane_id = plane_id_t; f.outer_normal = outer_normal_t;
		}
		//返回这个平面上的一个点, 也就是各个顶点的平均值
		point_3 get_rand_point (const std::vector<point_3> &vp) const
		{
			FT x(0), y(0), z(0);
			for (int i(0);i!=(int) np.size ();++i)
			{
				x = x + vp[np[i]].x(); y = y + vp[np[i]].y(); z = z + vp[np[i]].z();
			}
			return point_3(x/FT(np.size ()), y/FT(np.size ()), z/FT(np.size ()));
		}

		//序列化输入输出函数
		friend std::ostream& operator<< (std::ostream &os, const face &face);
		friend std::istream& operator>> (std::istream &is, face& face);
		friend bool operator== (const face &f1, const face &f2);
	};

	std::vector<face> vf;
	FT	attribute;		//可能有的性质，比如密度之类的
	int cblock_index;	//属于哪个复合块体

	void swap (eblock &tar)
	{
		FT temp_a(attribute);
		attribute = tar.attribute;
		tar.attribute = temp_a;
		int temp_i(cblock_index);
		cblock_index = tar.cblock_index;
		tar.cblock_index = temp_i;
		vf.swap(tar.vf);
	}
	
	//判断pi这个索引是否是这个单元的块体的一个顶点的索引。
	bool p_is_b_index (const int &pi) const
	{
		for (int i(0);i!=vf.size();++i)
			if (find (vf[i].np.begin(),vf[i].np.end (), pi) != vf[i].np.end ())
				return true;
		return false;
	}

	//返回一个点的索引，要求该索引表示的点在该块体的边界上，但是不在vf[fi]这个面上.返回-1表示出现错误
	//要求这个eblock已经添加完整了
	int p_inb_outf (const int &fi) const
	{
		int res(0);
		if (fi >= vf.size () || fi<0)
			throw std::logic_error("p_inb_outf, range error");
		for (int i(0);i!=vf.size ();++i)
		{
			if (i==fi)
				continue;
			for (int j(0);j!=vf[i].np.size ();++j)
				if (find (vf[fi].np.begin(),vf[fi].np.end (),vf[i].np[j]) == vf[fi].np.end())
					return vf[i].np[j];
		}
		return -1;
	}

	//计算这个eb的顶点的数量，总之是一个很昂贵的操作
	int get_point_count () const
	{
		std::vector<int> vpo;
		vpo.reserve (vf.size ()*6);
		for (int i(0);i!=vf.size ();++i)
			for (int j(0);j!=vf[i].np.size ();++j)
				add_set(vpo,vf[i].np[j]);
		return (int)vpo.size ();
	}

	//返回在这个单元块体内的一点，也就是各个顶点的平均值
	point_3 get_rand_point (const std::vector<point_3> &vp) const
	{
		FT x(0), y(0), z(0);
		for (int i(0);i!=(int) vf.size ();++i)
		{
			point_3 temp_p(vf[i].get_rand_point(vp));
			x = x + temp_p.x (); y = y + temp_p.y(); z = z + temp_p.z();
		}
		return point_3 (x/FT(vf.size ()), y/FT(vf.size ()), z/FT(vf.size ()));
	}

	//序列化输入输出函数
	friend std::ostream& operator<< (std::ostream &os, const eblock &eb);
	friend std::istream& operator>> (std::istream &is, eblock& eb);
	friend bool operator== (const eblock&eb1, const eblock&eb2);
};

class cblock
{
public:
	std::vector<INT> vb;
	//返回一个包含在这个复合块体中的点
	point_3 get_rand_p (const std::vector <eblock> &veb, const std::vector<point_3> &vp) const
	{
		return veb[vb[0]].get_rand_point(vp);
	}

	//序列化输入输出函数
	friend std::ostream& operator<< (std::ostream &os, const cblock &cb);
	friend std::istream& operator>> (std::istream &is, cblock& cb);
	friend bool operator== (const cblock&cb1, const cblock&cb2);
};

//将一个倾向倾角表示的方向计算出来，当然这种计算不是严格的，因为中间用到了sin和cos函数
inline vector_3 dip2normal (const FT& dip_dir_, const FT& dip_,const FT& dx_)
{
	double dip_dir = CGAL::to_double(dip_dir_*M_PI/180.0);
	double dip = CGAL::to_double(dip_*M_PI/180.0);
	double dx = CGAL::to_double(dx_*M_PI/180.0);
	return vector_3 (std::sin(dip)*std::cos(dip_dir-dx),
					-std::sin(dip)*std::sin(dip_dir-dx),
					std::cos(dip));
}

point_3 intersection_3planes (const plane_3&p1, const plane_3&p2, const plane_3&p3);

//这两个相交函数用线性规划的方法计算所给的两个凸集是否相交。
//将问题转化为求解最优化问题：min s  s.t. AX + D <= s
//计算两个单元块体是否相交，注意：也有可能返回Intersection_Status::Inter_Unbounded.
Intersection_Status intersection (const std::vector<plane_3> &vpl, const eblock& eb1, const eblock& eb2, const vector_3& offset1 = vector_3(0,0,0), const vector_3& offset2 = vector_3(0,0,0));		//判断两个单元块体是否相交,其中v1和v2分别是两个块体的偏移量
//计算两个符合块体是否相交，注意，如果两个符合块体中的单元块体无界相交，则报错。但是如果该函数返回Intersection_Status::Inter_Internal，不能说明不可能存在无界相交的单元块体
Intersection_Status intersection (const std::vector<plane_3> &vpl, const std::vector <eblock> &veb, const cblock& cb1, const cblock& cb2, const vector_3& offset1 = vector_3(0,0,0), const vector_3& offset2 = vector_3(0,0,0));		//判断两个复合块体是否相交,其中v1和v2分别是两个块体的偏移量
Intersection_Status intersection (const std::vector<plane_3> &vpl);		//判断vpl所表示的半空间(ax+by+cz+d<0)是否有交集,vpl中平面的正方向指向各个半空间外部

//返回值0表示成功，交点被放在最后一个变量中
//返回值其他值表示不成功
int intersection_segp (const point_3 &p1, const point_3 &p2, const plane_3 &pl, point_3 &res);

//注意，返回的不是精确结果，其中tol表示最终结果比原结果扩大了多少
//计算平面tar与一个圆盘（由pl，c，r表示）的交线段，由于无法用严格表示圆盘，所得到的结果并不是严格的，但是保证所得到的结果包含精确结果
//最终结果保存在res中，函数返回值是0时表示有交线段，此时res有意义，返回其他值时表示要么tar和pl平行要么没有交线段，此时res无意义
int intersect (const plane_3 &tar, const plane_3 &pl, const point_3 &c, const FT& r, segment_3 &res_s, const double &tol = 0.000001);

//计算np中索引所指向的点组成的法线方向。假设pn中所有点都在pl上（exactly）
vector_3 oriented_normal (const std::vector<point_3> &vp, const std::vector<int> &np, const plane_3 &pl);
//同上，计算np中索引所指向的点组成的法线方向。假设pn中所有点都（exactly）在同一个平面上
vector_3 oriented_normal (const std::vector<point_3> &vp, const std::vector<int> &np);
//同上，计算vp中所有点的组成边的法线方向。假设vp中所有点都（exactly）在同一个平面上
vector_3 oriented_normal (const std::vector<point_3> &vp);

//将空间圆盘裂隙转化为平面多边形，最后一个参数指定多边形边数，假设pl和df所在平面重合，res是返回值，坐标变换不保距离
//返回值始终是0
//近似表示，将圆盘表示成一个多边形本身就是一个近似
int proj_polygon (const plane_3 & pl, const disc_frac& df, polygon_2 &res, const int &n=30);

//将多边形裂隙转化为平面多边形，假设pl和pf所在平面重合，res是返回值,要求多边形裂隙是合法的，坐标变化不保距离
//返回0表示构建成功，返回其他值表示失败
int proj_polygon (const plane_3 & pl, const poly_frac& pf, polygon_with_holes_2 &res);

//将np投影到pl上，res是结果，不保证结果的形状，返回值始终是0，坐标变化不保距离
int proj_polygon (const plane_3 & pl, const std::vector<point_3>& np, polygon_2 &res);

//将索引vpi上的点投影到pl上，res是结果，不保证结果形状，返回值始终是0，坐标变化不保距离
int proj_polygon (const plane_3 & pl, const std::vector<point_3>& vpo, const std::vector<int>& vpi, polygon_2 &res);

//为多边形构建mbox，其中offset是外扩距离
mbox_2 build_mbox_2 (const polygon_2 &p,const double &offset = 0.001);

//构建三维的mbox。
mbox_3 build_mbox_3 (const point_3* vpo, const int & count, const double &offset = 0.001);
mbox_3 build_mbox_3 (const std::vector<point_3> &vpo, const double &offset = 0.001);
mbox_3 build_mbox_3 (const std::vector<point_3> &vpo, const std::vector<int> &vpi, const double &offset = 0.001);
mbox_3 build_mbox_3 (const poly_frac &df, const double &offset = 0.001);
mbox_3 build_mbox_3 (const point_3 &p, const plane_3&pl, const FT &r, const double &offset = 0.001);	//构建一个空间圆盘的三维最小边界盒
mbox_3 build_mbox_3 (const std::vector<point_3> &vpo, const eblock &eb, const double &offset = 0.001);
mbox_3 build_mbox_3 (const block_system &bs, const cblock &cb, const double offset = 0.001);
bool if_contained_mbox_3 (const mbox_3 &tar, const mbox_3 &b);	//判断b是否包含在tar中

//计算带洞的多边形的面积
inline FT cal_area_pwh (const polygon_with_holes_2 &pwh)
{
	FT area = 0;
	area += pwh.outer_boundary().area();
	for (polygon_with_holes_2::Hole_const_iterator i = pwh.holes_begin();i!=pwh.holes_end();++i)
		area += (i->area());
	return area;
}

//计算polyset的面积
inline FT cal_area_pset (const polygon_set_2 &pset)
{
	FT area(0.0);
	std::vector<polygon_with_holes_2> temp_v_polygons;
	pset.polygons_with_holes (std::back_inserter (temp_v_polygons));
	for (int i(0);i!=temp_v_polygons.size();++i)
		area += cal_area_pwh(temp_v_polygons[i]);
	return area;
}

template <typename K>
inline int add_set (std::vector<K> &p, const K &po)
{
	int i(0);
	for (i;i!=p.size ();++i)
		if (p[i] == po)
			break;
	if (i==p.size ())
		p.push_back (po);
	return i;
}

template <typename K>
inline typename std::list<K>::iterator add_set (std::list<K> &p, const K &po)
{
	std::list<K>::iterator i,last;
	for (i = last = p.begin ();i!=p.end ();last = i++)
		if ((*i) == po)
			break;
	if (i==p.end())
	{
		p.push_back (po);		//还没有添加
		return ++last;
	}
	else
		return i;				//已经添加
}

inline int add_set_pl (std::vector<plane_3> &vp, const plane_3 &p)
{
	int i(0);
	for (i;i!=vp.size ();++i)
		if (CGAL::parallel(vp[i],p) && vp[i].has_on(p.point()))
			break;

	if (i==vp.size ())
		vp.push_back (p);
	return i;
}

template <typename K, typename C>
inline int add_set (std::vector<K> &p, const K &po, const C&equal)
{
	int i(0);
	for (i;i!=p.size ();++i)
		if (equal(p[i],po))
			break;
	if (i==p.size ())
		p.push_back (po);
	return i;
}

template <typename K>
inline int find_n (const std::vector<K> &vp, const K &p)
{
	int i(0);
	for (;i!=vp.size ();++i)
		if (vp[i] == p)
			break;
	return i;
}

inline bool  radius_bigger (const disc_frac &f1, const disc_frac &f2)
{
	return (f1.r > f2.r);
}

//输出CAD脚本
template<typename P>
void out_line_3(std::ostream &os,const P& p1,const P &p2)
{
	os<<"line "<<CGAL::to_double(p1.x())<<','<<CGAL::to_double(p1.y())<<','<<CGAL::to_double(p1.z())<<' '
		<<CGAL::to_double(p2.x())<<','<<CGAL::to_double(p2.y())<<','<<CGAL::to_double(p2.z())<<"  "<<std::endl;
}

template<typename P>
void out_line_2(std::ostream &os,const P& p1,const P &p2)
{
	os<<"line "<<p1.x()<<','<<p1.y()<<' '
		<<p2.x()<<','<<p2.y()<<"  "<<std::endl;
}

inline void out_line_3 (std::ostream &os, const poly_frac::v_loop &vl)
{
	for (int i(0);i!=vl.size ();++i)
	{
		int i1(i+1);
		if (i1>=vl.size ())
			i1 = 0;
		out_line_3 (os, vl[i], vl[i1]);
	}
}

//inline void out_line_3_rhino (std::ostream &os, const poly_frac::v_loop &v1)
//{
//	os<<"polyline ";
//	for (int i(0);i!= (int)v1.size ();++i)
//		os<<v1[i].x()<<','<<v1[i].y()<<','<<v1[i].z()<<endl;
//	os<<"c"<<endl;
//	os<<"sellast"<<endl;
//	os<<"planarsrf"<<endl;
//	os<<"delete"<<endl<<endl;
//}

inline void out_line_3(std::ostream &os, const poly_frac &pf)
{
	out_line_3 (os,pf.vp_ob);
	for (int i(0);i!=pf.vp_ib.size ();++i)
		out_line_3 (os, pf.vp_ib[i]);
}

inline void out_line_3_rhino (std::ostream &os, const poly_frac &pf)
{
	os<<"polyline ";
	for (int i(0);i!=(int) pf.vp_ob.size ();++i)
		os<<pf.vp_ob[i].x()<<','<<pf.vp_ob[i].y()<<','<<pf.vp_ob[i].z()<<std::endl;
	os<<"c"<<std::endl;
	for (int i(0);i!= (int)pf.vp_ib.size ();++i)
	{
		os<<"polyline ";
		for (int j(0);j!= (int)pf.vp_ib[i].size ();++j)
			os<<pf.vp_ib[i][j].x()<<','<<pf.vp_ib[i][j].y()<<','<<pf.vp_ib[i][j].z()<<std::endl;
		os<<"c"<<std::endl;
	}
	os<<"selcrv"<<std::endl;
	os<<"planarsrf"<<std::endl;
	os<<"delete"<<std::endl;
}

inline void out_line_2(std::ostream &os, const polygon_2&poly)
{
	for (polygon_2::Edge_const_iterator i = poly.edges_begin(); i!=poly.edges_end();++i)
		out_line_2(os,i->start(), i->end());
}

inline void out_line_2(std::ostream &os, const polygon_with_holes_2&pwh)
{
	out_line_2(os,pwh.outer_boundary());
	for (polygon_with_holes_2::Hole_const_iterator i = pwh.holes_begin();
		i!=pwh.holes_end();++i)
		out_line_2(os, *i);
}

inline void out_line_2(std::ostream &os, const polygon_set_2 &ps)
{
	std::vector<polygon_with_holes_2> vpwh;
	vpwh.reserve (ps.number_of_polygons_with_holes());
	ps.polygons_with_holes(std::back_inserter(vpwh));

	for (int i(0);i!=(int)vpwh.size ();++i)
		out_line_2(os,vpwh[i]);

}

polygon_2 change_proj_plane (const plane_3 &oldp, const plane_3 &newp, const polygon_2 &tar);
polygon_with_holes_2 change_proj_plane (const plane_3 &oldp, const plane_3 &newp, const polygon_with_holes_2 &tar);
polygon_set_2 change_proj_plane (const plane_3 &oldp, const plane_3 &newp, const polygon_set_2 &tar);

//测试用函数，可以把一个圆盘裂隙转化为多边形裂隙(近似地)，测试多边形裂隙的代码,n表示转化为几边形
poly_frac disc2poly (const disc_frac & df, const plane_3 &pl,const int &n = 20, const FT &dx=0);

//将二维多边形投影到E3中，并且沿off_set偏移。pf是返回值.其中inverse表示是否需要在偏移前将生成的三维点取镜像（i.e. 如果inverse=true  则在偏移前先把(x,y,z)转化为(-x,-y,-z)）。
int construct_poly_frac (poly_frac & pf, const polygon_with_holes_2 &pwh, const plane_3& pl, const vector_3& off_set = vector_3(0,0,0), const bool &inverse = false);

inline bool operator== (const disc_frac&df1, const disc_frac&df2)
{
	return (df1.aperture == df2.aperture &&
		df1.center == df2.center &&
		df1.cohesion == df2.cohesion &&
		df1.frac_id == df2.frac_id &&
		df1.f_angle == df2.f_angle &&
		df1.plane_id == df2.plane_id &&
		df1.r == df2.r);
}

inline bool operator== (const poly_frac& pf1, const poly_frac& pf2)
{
	return (pf1.aperture == pf2.aperture &&
		pf1.cohesion == pf2.cohesion &&
		pf1.frac_id == pf2.frac_id &&
		pf1.f_angle == pf2.f_angle &&
		pf1.plane_id == pf2.plane_id &&
		pf1.vp_ib == pf2.vp_ib &&
		pf1.vp_ob == pf2.vp_ob);
}

inline bool operator== (const eblock::face& f1, const eblock::face& f2)
{
	return (f1.frac_id == f2.frac_id &&
		f1.np == f2.np &&
		f1.outer_normal == f2.outer_normal &&
		f1.plane_id == f2.plane_id);
}

inline bool operator== (const eblock& eb1, const eblock& eb2)
{
	return (eb1.attribute == eb2.attribute &&
		eb1.cblock_index == eb2.cblock_index &&
		eb1.vf == eb2.vf);
}

inline bool operator== (const cblock& cb1, const cblock& cb2)
{
	return (cb1.vb == cb2.vb);
}


//精简多边形的顶点，poly是输入多边形，res是输出。
//其中res和poly的形状应该相同，但是res和poly比起来少了一些多余的顶点
//vp_retain表示在res中需要得到保留的点，这些点不会被删掉
int simplify_polygon_2 (const polygon_2& poly, polygon_2 & res, const std::vector<point_2> &vp_retain = std::vector<point_2>());
//利用simplify_polygon_2这个函数精简一个带洞的多边形
int simplify_polygon_2 (const polygon_with_holes_2 &pwh, polygon_with_holes_2 & res);

//将relative simple 的多边形分割为若干simple的多边形，前提是poly一定不是simple的
int split_polygon_2 (const polygon_2 &poly, std::vector<polygon_2> &res);
//前提是pwh的外表面不是simple的
int split_polygon_2 (const polygon_with_holes_2 &pwh, std::vector<polygon_with_holes_2> &res);

//用来标记三角形面，在输出outer_boundary和inner_boundary的STL时有用到。只能用来标记嵌套过一次的面
template <typename CDT>
inline void mark_domains (CDT &cdt, typename CDT::Face_handle start, int index, std::list<typename CDT::Edge>& border)
{
	if (start->info().nesting_level !=-1)
		return ;
	std::list<CDT::Face_handle> queue;
	queue.push_back (start);
	while (!queue.empty())
	{
		CDT::Face_handle fh = queue.front();
		queue.pop_front();
		if (fh->info().nesting_level == -1)
		{
			fh->info().nesting_level = index;
			for (int i(0);i<3;++i)
			{
				CDT::Edge e(fh,i);
				CDT::Face_handle n = fh->neighbor(i);
				if (n->info().nesting_level == -1)
					if (cdt.is_constrained(e)) border.push_back (e);
					else queue.push_back(n);
			}
		}
	}
}

//和上面那个函数结合使用，（其实是上面那个函数和这个函数结合使用）
template <typename CDT>
inline void mark_domains (CDT &cdt)
{
	for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it!=cdt.all_faces_end();++it)
		it->info().nesting_level = -1;
	std::list<CDT::Edge> border;
	mark_domains(cdt,cdt.infinite_face(),0,border);
	while (!border.empty())
	{
		CDT::Edge e = border.front();
		border.pop_front();
		CDT::Face_handle n = e.first->neighbor(e.second);
		if (n->info().nesting_level == -1)
			mark_domains(cdt,n,e.first->info().nesting_level+1,border);
	}
}

//这个函数用来给一个Constrained Delaunay Triangulation (CDT)中添加多边形限制的
template <typename CDT, typename Polygon_2>
inline void insert_polygon(CDT& cdt, const Polygon_2& polygon)
{
	if (polygon.is_empty())
		return;
	CDT::Vertex_handle v_prev = cdt.insert(* CGAL::cpp11::prev(polygon.vertices_end()));
	for (Polygon_2::Vertex_iterator vit = polygon.vertices_begin();
		vit!=polygon.vertices_end();++vit)
	{
		CDT::Vertex_handle vh = cdt.insert(*vit);
		cdt.insert_constraint (vh,v_prev);
		v_prev = vh;
	}
}

//poly表示一个在平面pl上的多边形，normal_dir表示这个面的朝向（true这个面的外法线方向与pl这个平面相同）tranx,trany, tranz是偏移量，flag表示选用哪种输出格式（true表示将这个多边形输出成单独的solid，false表示将这个多边形输出成某个solid中的面片）
int output_STL_poly (std::ostream &outfile, const polygon_with_holes_2 &poly, const plane_3& pl, const bool &normal_dir, const FT&tranx, const FT&trany, const FT&tranz, const bool &flag = true);

}