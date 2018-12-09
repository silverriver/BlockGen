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

namespace BI	//block identification�ļ�д
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
enum Oriented_side {IN_, OUT_, ON_};		//���
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
	INT		plane_id;			//��϶����ƽ���ID
	INT		frac_id;			//��϶��ID,����
	FT		cohesion;
	FT		f_angle;
	FT		aperture;
};

class disc_frac: public frac_base
{
public:
	point_3 center;
	FT		r;

	//���л������������
	friend std::ostream& operator<< (std::ostream &os, const disc_frac &df);
	friend std::istream& operator>> (std::istream &is, disc_frac& df);
	friend bool operator== (const disc_frac &df1, const disc_frac&df2);
};

class poly_frac: public frac_base
{
public:
	typedef std::vector <point_3> v_loop;
	v_loop vp_ob;		//��ʾ�������϶����档�����ʱ��ҪΪÿһ����ȡһ��projection.ʶ������ʱ����ת���򲢲���Ҫ
	std::vector <v_loop> vp_ib;		//��ʾ�������϶�ڱ���
	
	void proj_plane (const plane_3& pl)		//�������϶�ϵĵ㶼ͶӰ��pl��
	{
		for (int i(0);i!=(int)vp_ob.size ();++i)
			vp_ob[i] = pl.projection(vp_ob[i]);
		for (int i(0);i!=(int)vp_ib.size ();++i)
			for (int j(0);j!=(int)vp_ib[i].size ();++j)
				vp_ib[i][j] = pl.projection (vp_ib[i][j]);
	}

	//���л������������
	friend std::ostream& operator<< (std::ostream &os, const poly_frac &pf);
	friend std::istream& operator>> (std::istream &is, poly_frac& pf);
	friend bool operator== (const poly_frac &pf1, const poly_frac &pf2);
};

class eblock
{
public:
	struct face
	{
		std::vector<INT> np;	//��������Ӧ����ʱ��ָ����
		INT frac_id;			//�γ���������϶��ID,���frac_id�Ǹ��������ʾ����������ĳ�ʼ�ӿ����ĳ���棬�������frac_id<=-4������治�ú�������ϲ������frac_id>=0ʱ�����������ĳ����϶���ɵģ��ϲ������ʱ��Ӧ�úϲ������-3<=frac_id<0ʱ������Ӧ�úϲ���
		INT plane_id;			//Ϊÿһ�����¼һ��ƽ�淽��
		bool outer_normal;		//true��������ⷨ�߷���Ϊplane_id��ȷ�����������positive size��
		void swap (face &f)
		{
			INT frac_id_t = frac_id;
			INT plane_id_t = plane_id;
			bool outer_normal_t = outer_normal;
			frac_id = f.frac_id; plane_id = f.plane_id; outer_normal = f.outer_normal;
			np.swap(f.np);
			f.frac_id = frac_id_t; f.plane_id = plane_id_t; f.outer_normal = outer_normal_t;
		}
		//�������ƽ���ϵ�һ����, Ҳ���Ǹ��������ƽ��ֵ
		point_3 get_rand_point (const std::vector<point_3> &vp) const
		{
			FT x(0), y(0), z(0);
			for (int i(0);i!=(int) np.size ();++i)
			{
				x = x + vp[np[i]].x(); y = y + vp[np[i]].y(); z = z + vp[np[i]].z();
			}
			return point_3(x/FT(np.size ()), y/FT(np.size ()), z/FT(np.size ()));
		}

		//���л������������
		friend std::ostream& operator<< (std::ostream &os, const face &face);
		friend std::istream& operator>> (std::istream &is, face& face);
		friend bool operator== (const face &f1, const face &f2);
	};

	std::vector<face> vf;
	FT	attribute;		//�����е����ʣ������ܶ�֮���
	int cblock_index;	//�����ĸ����Ͽ���

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
	
	//�ж�pi��������Ƿ��������Ԫ�Ŀ����һ�������������
	bool p_is_b_index (const int &pi) const
	{
		for (int i(0);i!=vf.size();++i)
			if (find (vf[i].np.begin(),vf[i].np.end (), pi) != vf[i].np.end ())
				return true;
		return false;
	}

	//����һ�����������Ҫ���������ʾ�ĵ��ڸÿ���ı߽��ϣ����ǲ���vf[fi]�������.����-1��ʾ���ִ���
	//Ҫ�����eblock�Ѿ����������
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

	//�������eb�Ķ������������֮��һ���ܰ���Ĳ���
	int get_point_count () const
	{
		std::vector<int> vpo;
		vpo.reserve (vf.size ()*6);
		for (int i(0);i!=vf.size ();++i)
			for (int j(0);j!=vf[i].np.size ();++j)
				add_set(vpo,vf[i].np[j]);
		return (int)vpo.size ();
	}

	//�����������Ԫ�����ڵ�һ�㣬Ҳ���Ǹ��������ƽ��ֵ
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

	//���л������������
	friend std::ostream& operator<< (std::ostream &os, const eblock &eb);
	friend std::istream& operator>> (std::istream &is, eblock& eb);
	friend bool operator== (const eblock&eb1, const eblock&eb2);
};

class cblock
{
public:
	std::vector<INT> vb;
	//����һ��������������Ͽ����еĵ�
	point_3 get_rand_p (const std::vector <eblock> &veb, const std::vector<point_3> &vp) const
	{
		return veb[vb[0]].get_rand_point(vp);
	}

	//���л������������
	friend std::ostream& operator<< (std::ostream &os, const cblock &cb);
	friend std::istream& operator>> (std::istream &is, cblock& cb);
	friend bool operator== (const cblock&cb1, const cblock&cb2);
};

//��һ��������Ǳ�ʾ�ķ�������������Ȼ���ּ��㲻���ϸ�ģ���Ϊ�м��õ���sin��cos����
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

//�������ཻ���������Թ滮�ķ�����������������͹���Ƿ��ཻ��
//������ת��Ϊ������Ż����⣺min s  s.t. AX + D <= s
//����������Ԫ�����Ƿ��ཻ��ע�⣺Ҳ�п��ܷ���Intersection_Status::Inter_Unbounded.
Intersection_Status intersection (const std::vector<plane_3> &vpl, const eblock& eb1, const eblock& eb2, const vector_3& offset1 = vector_3(0,0,0), const vector_3& offset2 = vector_3(0,0,0));		//�ж�������Ԫ�����Ƿ��ཻ,����v1��v2�ֱ������������ƫ����
//�����������Ͽ����Ƿ��ཻ��ע�⣬����������Ͽ����еĵ�Ԫ�����޽��ཻ���򱨴���������ú�������Intersection_Status::Inter_Internal������˵�������ܴ����޽��ཻ�ĵ�Ԫ����
Intersection_Status intersection (const std::vector<plane_3> &vpl, const std::vector <eblock> &veb, const cblock& cb1, const cblock& cb2, const vector_3& offset1 = vector_3(0,0,0), const vector_3& offset2 = vector_3(0,0,0));		//�ж��������Ͽ����Ƿ��ཻ,����v1��v2�ֱ������������ƫ����
Intersection_Status intersection (const std::vector<plane_3> &vpl);		//�ж�vpl����ʾ�İ�ռ�(ax+by+cz+d<0)�Ƿ��н���,vpl��ƽ���������ָ�������ռ��ⲿ

//����ֵ0��ʾ�ɹ������㱻�������һ��������
//����ֵ����ֵ��ʾ���ɹ�
int intersection_segp (const point_3 &p1, const point_3 &p2, const plane_3 &pl, point_3 &res);

//ע�⣬���صĲ��Ǿ�ȷ���������tol��ʾ���ս����ԭ��������˶���
//����ƽ��tar��һ��Բ�̣���pl��c��r��ʾ���Ľ��߶Σ������޷����ϸ��ʾԲ�̣����õ��Ľ���������ϸ�ģ����Ǳ�֤���õ��Ľ��������ȷ���
//���ս��������res�У���������ֵ��0ʱ��ʾ�н��߶Σ���ʱres�����壬��������ֵʱ��ʾҪôtar��plƽ��Ҫôû�н��߶Σ���ʱres������
int intersect (const plane_3 &tar, const plane_3 &pl, const point_3 &c, const FT& r, segment_3 &res_s, const double &tol = 0.000001);

//����np��������ָ��ĵ���ɵķ��߷��򡣼���pn�����е㶼��pl�ϣ�exactly��
vector_3 oriented_normal (const std::vector<point_3> &vp, const std::vector<int> &np, const plane_3 &pl);
//ͬ�ϣ�����np��������ָ��ĵ���ɵķ��߷��򡣼���pn�����е㶼��exactly����ͬһ��ƽ����
vector_3 oriented_normal (const std::vector<point_3> &vp, const std::vector<int> &np);
//ͬ�ϣ�����vp�����е����ɱߵķ��߷��򡣼���vp�����е㶼��exactly����ͬһ��ƽ����
vector_3 oriented_normal (const std::vector<point_3> &vp);

//���ռ�Բ����϶ת��Ϊƽ�����Σ����һ������ָ������α���������pl��df����ƽ���غϣ�res�Ƿ���ֵ������任��������
//����ֵʼ����0
//���Ʊ�ʾ����Բ�̱�ʾ��һ������α������һ������
int proj_polygon (const plane_3 & pl, const disc_frac& df, polygon_2 &res, const int &n=30);

//���������϶ת��Ϊƽ�����Σ�����pl��pf����ƽ���غϣ�res�Ƿ���ֵ,Ҫ��������϶�ǺϷ��ģ�����仯��������
//����0��ʾ�����ɹ�����������ֵ��ʾʧ��
int proj_polygon (const plane_3 & pl, const poly_frac& pf, polygon_with_holes_2 &res);

//��npͶӰ��pl�ϣ�res�ǽ��������֤�������״������ֵʼ����0������仯��������
int proj_polygon (const plane_3 & pl, const std::vector<point_3>& np, polygon_2 &res);

//������vpi�ϵĵ�ͶӰ��pl�ϣ�res�ǽ��������֤�����״������ֵʼ����0������仯��������
int proj_polygon (const plane_3 & pl, const std::vector<point_3>& vpo, const std::vector<int>& vpi, polygon_2 &res);

//Ϊ����ι���mbox������offset����������
mbox_2 build_mbox_2 (const polygon_2 &p,const double &offset = 0.001);

//������ά��mbox��
mbox_3 build_mbox_3 (const point_3* vpo, const int & count, const double &offset = 0.001);
mbox_3 build_mbox_3 (const std::vector<point_3> &vpo, const double &offset = 0.001);
mbox_3 build_mbox_3 (const std::vector<point_3> &vpo, const std::vector<int> &vpi, const double &offset = 0.001);
mbox_3 build_mbox_3 (const poly_frac &df, const double &offset = 0.001);
mbox_3 build_mbox_3 (const point_3 &p, const plane_3&pl, const FT &r, const double &offset = 0.001);	//����һ���ռ�Բ�̵���ά��С�߽��
mbox_3 build_mbox_3 (const std::vector<point_3> &vpo, const eblock &eb, const double &offset = 0.001);
mbox_3 build_mbox_3 (const block_system &bs, const cblock &cb, const double offset = 0.001);
bool if_contained_mbox_3 (const mbox_3 &tar, const mbox_3 &b);	//�ж�b�Ƿ������tar��

//��������Ķ���ε����
inline FT cal_area_pwh (const polygon_with_holes_2 &pwh)
{
	FT area = 0;
	area += pwh.outer_boundary().area();
	for (polygon_with_holes_2::Hole_const_iterator i = pwh.holes_begin();i!=pwh.holes_end();++i)
		area += (i->area());
	return area;
}

//����polyset�����
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
		p.push_back (po);		//��û�����
		return ++last;
	}
	else
		return i;				//�Ѿ����
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

//���CAD�ű�
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

//�����ú��������԰�һ��Բ����϶ת��Ϊ�������϶(���Ƶ�)�����Զ������϶�Ĵ���,n��ʾת��Ϊ������
poly_frac disc2poly (const disc_frac & df, const plane_3 &pl,const int &n = 20, const FT &dx=0);

//����ά�����ͶӰ��E3�У�������off_setƫ�ơ�pf�Ƿ���ֵ.����inverse��ʾ�Ƿ���Ҫ��ƫ��ǰ�����ɵ���ά��ȡ����i.e. ���inverse=true  ����ƫ��ǰ�Ȱ�(x,y,z)ת��Ϊ(-x,-y,-z)����
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


//�������εĶ��㣬poly���������Σ�res�������
//����res��poly����״Ӧ����ͬ������res��poly����������һЩ����Ķ���
//vp_retain��ʾ��res����Ҫ�õ������ĵ㣬��Щ�㲻�ᱻɾ��
int simplify_polygon_2 (const polygon_2& poly, polygon_2 & res, const std::vector<point_2> &vp_retain = std::vector<point_2>());
//����simplify_polygon_2�����������һ�������Ķ����
int simplify_polygon_2 (const polygon_with_holes_2 &pwh, polygon_with_holes_2 & res);

//��relative simple �Ķ���ηָ�Ϊ����simple�Ķ���Σ�ǰ����polyһ������simple��
int split_polygon_2 (const polygon_2 &poly, std::vector<polygon_2> &res);
//ǰ����pwh������治��simple��
int split_polygon_2 (const polygon_with_holes_2 &pwh, std::vector<polygon_with_holes_2> &res);

//��������������棬�����outer_boundary��inner_boundary��STLʱ���õ���ֻ���������Ƕ�׹�һ�ε���
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

//�������Ǹ��������ʹ�ã�����ʵ�������Ǹ�����������������ʹ�ã�
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

//�������������һ��Constrained Delaunay Triangulation (CDT)����Ӷ�������Ƶ�
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

//poly��ʾһ����ƽ��pl�ϵĶ���Σ�normal_dir��ʾ�����ĳ���true�������ⷨ�߷�����pl���ƽ����ͬ��tranx,trany, tranz��ƫ������flag��ʾѡ�����������ʽ��true��ʾ��������������ɵ�����solid��false��ʾ���������������ĳ��solid�е���Ƭ��
int output_STL_poly (std::ostream &outfile, const polygon_with_holes_2 &poly, const plane_3& pl, const bool &normal_dir, const FT&tranx, const FT&trany, const FT&tranz, const bool &flag = true);

}