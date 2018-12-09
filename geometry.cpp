#include "geometry.h"
#include <exception>
#include <CGAL\intersections.h>
#include <CGAL\Constrained_Delaunay_triangulation_2.h>
#include <CGAL\QP_functions.h>
#include <CGAL\Triangulation_face_base_with_info_2.h>
#include "identification.h"
using namespace std;

namespace BI
{

int simplify_polygon_2 (const polygon_2& poly, polygon_2 & res, const std::vector<point_2> &vp_retain)
{
	res.clear();
	if (poly.size () < 3)
		return 1;	//输入的多边形顶点数小于3

	for (int i(0);i!= (int)poly.size ();++i)
	{
		int i1 (i+1), i_1(i-1);
		if (i1 >= (int)poly.size ())
			i1 = 0;
		if (i_1 < 0)
			i_1 = (int)poly.size ()-1;
		if (find_n(vp_retain,poly[i]) != (int)vp_retain.size () || (!CGAL::collinear(poly[i_1],poly[i],poly[i1])))
			res.push_back (poly[i]);	//只有保留点和角点才加到新的polygon中，其他点不保留
	}
	return 0;
}

int simplify_polygon_2 (const polygon_with_holes_2 &pwh, polygon_with_holes_2 & res)
{
	polygon_2 outer;
	vector<polygon_2> vhole;
	vhole.reserve (pwh.number_of_holes());
	vector<point_2> vp_retained;		//记录哪些点应该被保留
	vp_retained.reserve (pwh.number_of_holes()*5);

	simplify_polygon_2(pwh.outer_boundary(),outer);
	for (int i(0);i!=(int)outer.size ();++i)
		vp_retained.push_back (outer[i]);		//外表面上的点也有可能在外表面上，所以也需要记录下

	for (polygon_with_holes_2::Hole_const_iterator hi = pwh.holes_begin(); hi!=pwh.holes_end();++hi)
	{
		polygon_2 temp;
		simplify_polygon_2(*hi,temp);
		for (int i(0);i!=(int)temp.size ();++i)
			vp_retained.push_back (temp[i]);	//这些点是hole的角点
	}

	simplify_polygon_2 (pwh.outer_boundary(),outer,vp_retained);
	res = polygon_with_holes_2(outer);

	for (polygon_with_holes_2::Hole_const_iterator hi = pwh.holes_begin(); hi!=pwh.holes_end();++hi)
	{
		polygon_2 temp;
		simplify_polygon_2(*hi,temp,vp_retained);
		res.add_hole (temp);
	}
	return 0;
}

//int simplify_polygon_2 (const polygon_set_2 &pset, polygon_set_2 &res)
//{
//	res.clear ();
//	vector<polygon_with_holes_2> temp;
//	temp.reserve(pset.number_of_polygons_with_holes());
//	pset.polygons_with_holes(back_inserter(temp));
//	for (int i(0);i!=temp.size ();++i)
//	{
//		polygon_with_holes_2 pwh_temp;
//		if (simplify_polygon_2(temp[i],pwh_temp)!=0)
//			return 1;		//简化多边形不成功
//		res.insert(pwh_temp);
//	}
//	return 0;
//}

int split_polygon_2 (const polygon_2 &poly, std::vector<polygon_2> &res)
{
	if (poly.is_simple())
	{
		res.push_back (poly);
		return 0;
	}
	for (int i(0);i!=poly.size ();++i)
		for (int j(i+1); j<poly.size ();++j)
		{
			if (poly[i] == poly[j])		//这里可以不用递归的，可以用迭代代替，并且也可以在函数开头的地方不用检查是否是simple的。这些等以后再优化吧
			{
				polygon_2 p1,p2;
				int temp_i(i);
				while (1)
				{
					if (i==j)
						break;
					p1.push_back (poly[i]);
					i++;
				}
				while (1)
				{
					if (j == temp_i)
						break;
					p2.push_back (poly[j]);
					j++;
					if (j >= (int)poly.size ())
						j = 0;
				}
				return split_polygon_2(p1,res) + split_polygon_2(p2,res) +1;
			}
		}
	return 1;		//理论上不应该在这里跳出
}

int split_polygon_2 (const polygon_with_holes_2 &pwh, std::vector<polygon_with_holes_2> &res)
{
	if (pwh.outer_boundary().is_simple())
		return 0;

	vector<polygon_2> vouter;
	vouter.reserve (2);
	split_polygon_2(pwh.outer_boundary(),vouter);
	res.clear();
	res.reserve (vouter.size ());
	for (int i(0);i!=vouter.size();++i)
		res.push_back (polygon_with_holes_2(vouter[i]));
	for (polygon_with_holes_2::Hole_const_iterator hi = pwh.holes_begin();
		hi!=pwh.holes_end();++hi)
	{
		polygon_2 p_temp(*hi);
		if (p_temp.is_clockwise_oriented())
			p_temp.reverse_orientation();
		int i(0);
		for (;i!=res.size ();++i)
			if (CGAL::do_intersect(res[i].outer_boundary(),p_temp))
				break;
		if (i == res.size ())
			throw logic_error ("split_polygon_2, can not find a poly for hole");
		res[i].add_hole(*hi);
	}
	return 1;
}

point_3 intersection_3planes (const plane_3&p1, const plane_3&p2, const plane_3&p3)
{
	CGAL::cpp11::result_of<K::Intersect_3(K::Plane_3,K::Plane_3)>::type
		res_line = CGAL::intersection(p1,p2);
	if (res_line)
	{
		if (const K::Line_3* l = boost::get<K::Line_3> (&*res_line))
		{
			CGAL::cpp11::result_of<K::Intersect_3(K::Line_3,K::Plane_3)>::type
				res_p = CGAL::intersection(*l,p3);
			if (res_p)
			{
				if (const K::Point_3* p = boost::get<K::Point_3> (&*res_p))
					return *p;
				else
					throw logic_error ("intersection_3planes, line and p3 parallel");
			}
			else
				throw logic_error ("intersection_3planes, line and p3 do not intersect");
		}
		else
			throw logic_error ("intersection_3planes, p1 and p2 coplanar");
	}
	else
		throw logic_error ("intersection_3planes, p1 and p2 do not intersect");
}

int intersection_segp (const point_3 &p1, const point_3 &p2, const plane_3 &pl, point_3 &res)
{
	CGAL::cpp11::result_of<K::Intersect_3(K::Segment_3, K::Plane_3)>::type
		result = CGAL::intersection(pl,segment_3(p1,p2));
	if (result)
	{
		if (const K::Point_3 *s = boost::get<K::Point_3> (&*result))
		{
			res = *s;
			return 0;
		}
	}
	return 1;
}

int intersect (const plane_3 &tar, const plane_3 &pl, const point_3 &c_, const FT& r, segment_3 &res, const double &tol)
{
	if (CGAL::parallel(tar,pl))
		return 1;

	point_3 c (pl.projection(c_));

	CGAL::cpp11::result_of <K::Intersect_3(K::Plane_3, K::Plane_3)>::type
		result = CGAL::intersection (tar,pl);
	if (result)
	{
		if (const K::Line_3 *l = boost::get<K::Line_3> (&*result))
		{
			vector_3 r(*l);
			double x (CGAL::to_double(r.x())), y(CGAL::to_double(r.y())), z(CGAL::to_double(r.z()));
			r = (1/sqrt(x*x+y*y+z*z))*r;

			point_3 pp (l->projection (c));
			FT squared_dis (CGAL::squared_distance(pp,c));

			if (squared_dis > r*r)
				return 2;
			
			double temp = sqrt(CGAL::to_double(r*r-squared_dis))+tol;
			res = segment_3(pp+r*temp, pp-r*temp);
		}
		else
			throw logic_error ("intersect (const plane_3 &tar, const plane_3 &pl, const point_3 &c_, const FT& r, segment_3 &res), tar parallel to pl");
	}
	else 
		throw logic_error ("intersect (const plane_3 &tar, const plane_3 &pl, const point_3 &c_, const FT& r, segment_3 &res), intersection error");
	return 0;
}

Intersection_Status intersection (const std::vector<plane_3> &vpl, const eblock& eb1, const eblock& eb2, const vector_3& v1, const vector_3& v2)
{
	typedef CGAL::Quadratic_program<K::FT> Program;
	typedef CGAL::Quadratic_program_solution<K::FT> Solution;
	Program lp (CGAL::SMALLER,false,0,false,0);

	for (int i(0);i!=(int)eb1.vf.size ();++i)
	{
		K::FT offset1 = vpl[eb1.vf[i].plane_id].a ()*v1.x() + vpl[eb1.vf[i].plane_id].b ()*v1.y() + vpl[eb1.vf[i].plane_id].c ()*v1.z();
		if (eb1.vf[i].outer_normal )
		{
			lp.set_a (0,i,vpl[eb1.vf[i].plane_id].a ());
			lp.set_a (1,i,vpl[eb1.vf[i].plane_id].b ());
			lp.set_a (2,i,vpl[eb1.vf[i].plane_id].c ());
			lp.set_b (i,-vpl[eb1.vf[i].plane_id].d()+offset1);
		}
		else
		{
			lp.set_a (0,i,-vpl[eb1.vf[i].plane_id].a ());
			lp.set_a (1,i,-vpl[eb1.vf[i].plane_id].b ());
			lp.set_a (2,i,-vpl[eb1.vf[i].plane_id].c ());
			lp.set_b (i,vpl[eb1.vf[i].plane_id].d()-offset1);
		}
		lp.set_a (3,i,-1);
	}
	for (int i(0);i!=(int)eb2.vf.size ();++i)
	{
		K::FT offset2 = vpl[eb2.vf[i].plane_id].a ()*v2.x() + vpl[eb2.vf[i].plane_id].b ()*v2.y() + vpl[eb2.vf[i].plane_id].c ()*v2.z();
		if (eb2.vf[i].outer_normal )
		{
			lp.set_a (0,i+ (int) eb1.vf.size (),vpl[eb2.vf[i].plane_id].a ());
			lp.set_a (1,i+ (int) eb1.vf.size (),vpl[eb2.vf[i].plane_id].b ());
			lp.set_a (2,i+ (int) eb1.vf.size (),vpl[eb2.vf[i].plane_id].c ());
			lp.set_b (i+ (int) eb1.vf.size (),-vpl[eb2.vf[i].plane_id].d()+offset2);
		}
		else
		{
			lp.set_a (0,i+ (int) eb1.vf.size (),-vpl[eb2.vf[i].plane_id].a ());
			lp.set_a (1,i+ (int) eb1.vf.size (),-vpl[eb2.vf[i].plane_id].b ());
			lp.set_a (2,i+ (int) eb1.vf.size (),-vpl[eb2.vf[i].plane_id].c ());
			lp.set_b (i+ (int) eb1.vf.size (),vpl[eb2.vf[i].plane_id].d()-offset2);
		}
		lp.set_a (3,i,-1);
	}
	lp.set_c(3,1);
	Solution s = CGAL::solve_linear_program(lp,K::FT());
	if (s.is_infeasible())
		throw logic_error ("intersection (eblock, eblock), LP infeasible");
	if (s.is_unbounded())
		return Intersection_Status::Inter_Unbounded;
	if (s.objective_value() <0)					//内部相交
		return Intersection_Status::Inter_Internal;
	else if (s.objective_value() ==0)			//边界处相交
		return Intersection_Status::Inter_Boundary;
	else
		return Intersection_Status::Separate;	//不相交
}

Intersection_Status intersection (const std::vector<plane_3> &vpl, const std::vector <eblock> &veb, const cblock& cb1, const cblock& cb2, const vector_3& offset1, const vector_3& offset2)
{
	int n_bound(0), n_inter(0);
	for (int i(0);i!=(int)cb1.vb.size ();++i)
	{
		for (int j(0);j!=(int)cb2.vb.size ();++j)
		{
			switch (intersection (vpl,veb[cb1.vb[i]],veb[cb2.vb[j]],offset1,offset2))
			{
			case Intersection_Status::Inter_Boundary:
				++n_bound; break;
			case Intersection_Status::Inter_Internal:
				return Intersection_Status::Inter_Internal;
			case Intersection_Status::Inter_Unbounded:
				throw logic_error("intersection (cblock, cblock), intersection unbounded");
			}
		}
	}
	if (n_bound > 0)
		return Intersection_Status::Inter_Boundary;
	else
		return Intersection_Status::Separate;
}

Intersection_Status intersection (const std::vector<plane_3> &vpl)
{
	typedef CGAL::Quadratic_program<K::FT> Program;
	typedef CGAL::Quadratic_program_solution<K::FT> Solution;
	Program lp (CGAL::SMALLER,false,0,false,0);

	for (int i(0);i!=(int)vpl.size ();++i)
	{
		lp.set_a (0,i,vpl[i].a());
		lp.set_a (1,i,vpl[i].b());
		lp.set_a (2,i,vpl[i].c());
		lp.set_a (3,i,-1);
		lp.set_b (i,-vpl[i].d());
	}
	lp.set_c (3,1);
	Solution s = CGAL::solve_linear_program(lp,K::FT());
	if (s.is_infeasible())
		throw logic_error ("intersection (vector<plane_3>, LP infeasible");		//不可能没有可行域，可行域一定存在
	if (s.is_unbounded())		//最优解无限小的情况返回3
		return Intersection_Status::Inter_Unbounded;
	if (s.objective_value() < 0)		//内部相交
		return Intersection_Status::Inter_Internal;
	else if (s.objective_value() == 0)		//边界处相交
		return Intersection_Status::Inter_Boundary;
	else 
		return Intersection_Status::Separate;		//不相交，应为最优解大于零
}


vector_3 oriented_normal (const std::vector<point_3> &vp, const std::vector<int> &np, const plane_3 &pl)
{
	if (np.size ()<3)
		throw logic_error("oriented_normal, point count less than 3");
	point_3 ori(pl.point());
	vector_3 res(0,0,0);
	for (int i(0);i!=np.size ();++i)
	{
		 int i1(i+1);
		 if (i1>=(int)np.size ())
			 i1 = 0;
		 res = res + CGAL::cross_product(vector_3(ori,vp[np[i]]), vector_3(ori,vp[np[i1]]));
	}
	return res;
}

vector_3 oriented_normal (const std::vector<point_3> &vp, const std::vector<int> &np)
{
	if (np.size ()<3)
		throw logic_error("oriented_normal, point count less than 3");
	const point_3 &ori =vp[np[0]];
	vector_3 res(0,0,0);
	for (int i(1);i!=np.size ()-1;++i)
	{
		 int i1(i+1);
		 res = res + CGAL::cross_product(vector_3(ori,vp[np[i]]), vector_3(ori,vp[np[i1]]));
	}
	return res;
}

vector_3 oriented_normal (const std::vector<point_3> &vp)
{
	if(vp.size ()<3)
		throw logic_error ("oriented_normal, point count less than 3");
	const point_3 &ori = vp[0];
	vector_3 res(0,0,0);
	for (int i(1);i!=vp.size ()-1;++i)
	{
		 int i1(i+1);
		 res = res + CGAL::cross_product(vector_3(ori,vp[i]), vector_3(ori,vp[i1]));
	}
	return res;
}

int proj_polygon (const plane_3 & pl, const disc_frac& df, polygon_2 &res, const int &n)
{
	res.clear();
	vector_3 b1 (pl.base1());
	vector_3 b2 (pl.base2());
	b1 = (1/sqrt(CGAL::to_double(b1.squared_length()))) * b1;	//既然将圆盘近似为多边形本身就是一个近似，那么就用近似函数吧
	b2 = (1/sqrt(CGAL::to_double(b2.squared_length()))) * b2;
	
	double r = CGAL::to_double(df.r);
	double da = M_PI_2/n;
	double a = 0;
	point_3 center (pl.projection(df.center));
	for (int i(0);i!=n;++i)
	{
		res.push_back (pl.to_2d(center + df.r*b1*cos(a) + df.r*b2*sin(a)));
		a += da;
	}
	if (res.is_clockwise_oriented())
		res.reverse_orientation();
	return 0;
}

int proj_polygon (const plane_3 & pl, const std::vector<point_3>& np, polygon_2 &res)
{
	res.clear();
	for (int i(0);i!=np.size ();++i)
		res.push_back (pl.to_2d(np[i]));
	return 0;
}

int proj_polygon (const plane_3 & pl, const std::vector<point_3>& vpo, const std::vector<int>& vpi, polygon_2 &res)
{
	res.clear();
	for (int i(0);i!=vpi.size ();++i)
		res.push_back (pl.to_2d(vpo[vpi[i]]));
	return 0;
}

int proj_polygon (const plane_3 & pl, const poly_frac& pf, polygon_with_holes_2 &res)
{
	polygon_2 outer_b;
	vector<polygon_2> inner_b (pf.vp_ib.size ());

	proj_polygon(pl,pf.vp_ob,outer_b);
	for (int i(0);i != (int)pf.vp_ib.size ();++i)
		proj_polygon(pl,pf.vp_ib[i],inner_b[i]);

	if (outer_b.is_collinear_oriented())
		return 2;
	if (outer_b.is_clockwise_oriented())
		outer_b.reverse_orientation();

	for (int i(0);i!= (int) inner_b.size ();++i)
	{
		if (!inner_b[i].is_simple())
			return 3;
		if (inner_b[i].is_collinear_oriented())
			return 4;
		if (inner_b[i].is_counterclockwise_oriented())
			inner_b[i].reverse_orientation();
	}

	res = polygon_with_holes_2(outer_b,inner_b.begin (),inner_b.end ());
	
	if (res.is_unbounded())
		return 5;
	if (!polygon_set_2(res).is_valid())
		return 6;
	return 0;
}

mbox_2 build_mbox_2 (const polygon_2 &p,const double &offset)
{
	if (p.is_empty())
		throw logic_error ("build_bbox_2, poly is empty");
	double minx(CGAL::to_double(p[0].x())),maxx(CGAL::to_double(p[0].x()));
	double miny(CGAL::to_double(p[0].y())),maxy(CGAL::to_double(p[0].y()));
	for (polygon_2::Vertex_const_iterator i = p.vertices_begin();
		i!=p.vertices_end();++i)
	{
		if (i->x() > maxx) maxx = CGAL::to_double(i->x());
		if (i->x() < minx) minx = CGAL::to_double(i->x());
		if (i->y() > maxy) maxy = CGAL::to_double(i->y());
		if (i->y() < miny) miny = CGAL::to_double(i->y());
	}
	return mbox_2(minx-offset, miny-offset, maxx+offset, maxy+offset);
}

mbox_3 build_mbox_3 (const point_3* vpo, const int & count, const double &offset)
{
	if (count<1)
		throw logic_error ("build_mbox_3, vpo is empty");
	double minx(CGAL::to_double(vpo[0].x())), maxx(CGAL::to_double(vpo[0].x()));
	double miny(CGAL::to_double(vpo[0].y())), maxy(CGAL::to_double(vpo[0].y()));
	double minz(CGAL::to_double(vpo[0].z())), maxz(CGAL::to_double(vpo[0].z()));
	for (int i(0);i!=count;++i)
	{
		if (vpo[i].x() > maxx) maxx = CGAL::to_double(vpo[i].x());
		if (vpo[i].x() < minx) minx = CGAL::to_double(vpo[i].x());
		if (vpo[i].y() > maxy) maxy = CGAL::to_double(vpo[i].y());
		if (vpo[i].y() < miny) miny = CGAL::to_double(vpo[i].y());
		if (vpo[i].z() > maxz) maxz = CGAL::to_double(vpo[i].z());
		if (vpo[i].z() < minz) minz = CGAL::to_double(vpo[i].z());
	}
	return mbox_3(minx-offset, miny-offset, minz-offset, maxx+offset, maxy+offset, maxz+offset);
}

mbox_3 build_mbox_3 (const std::vector<point_3> &vpo, const double &offset)
{
	if (vpo.size ()<1)
		throw logic_error ("build_mbox_3, vpo is empty");
	double minx(CGAL::to_double(vpo[0].x())), maxx(CGAL::to_double(vpo[0].x()));
	double miny(CGAL::to_double(vpo[0].y())), maxy(CGAL::to_double(vpo[0].y()));
	double minz(CGAL::to_double(vpo[0].z())), maxz(CGAL::to_double(vpo[0].z()));
	for (int i(0);i!=vpo.size ();++i)
	{
		if (vpo[i].x() > maxx) maxx = CGAL::to_double(vpo[i].x());
		if (vpo[i].x() < minx) minx = CGAL::to_double(vpo[i].x());
		if (vpo[i].y() > maxy) maxy = CGAL::to_double(vpo[i].y());
		if (vpo[i].y() < miny) miny = CGAL::to_double(vpo[i].y());
		if (vpo[i].z() > maxz) maxz = CGAL::to_double(vpo[i].z());
		if (vpo[i].z() < minz) minz = CGAL::to_double(vpo[i].z());
	}
	return mbox_3(minx-offset, miny-offset, minz-offset, maxx+offset, maxy+offset, maxz+offset);
}

mbox_3 build_mbox_3 (const std::vector<point_3> &vpo, const std::vector<int> &vpi, const double &offset)
{
	if (vpi.size ()<1)
		throw logic_error ("build_mbox_3, vpo is empty");
	double minx(CGAL::to_double(vpo[vpi[0]].x())), maxx(CGAL::to_double(vpo[vpi[0]].x()));
	double miny(CGAL::to_double(vpo[vpi[0]].y())), maxy(CGAL::to_double(vpo[vpi[0]].y()));
	double minz(CGAL::to_double(vpo[vpi[0]].z())), maxz(CGAL::to_double(vpo[vpi[0]].z()));
	for (int i(0);i!=vpi.size ();++i)
	{
		if (vpo[vpi[i]].x() > maxx) maxx = CGAL::to_double(vpo[vpi[i]].x());
		if (vpo[vpi[i]].x() < minx) minx = CGAL::to_double(vpo[vpi[i]].x());
		if (vpo[vpi[i]].y() > maxy) maxy = CGAL::to_double(vpo[vpi[i]].y());
		if (vpo[vpi[i]].y() < miny) miny = CGAL::to_double(vpo[vpi[i]].y());
		if (vpo[vpi[i]].z() > maxz) maxz = CGAL::to_double(vpo[vpi[i]].z());
		if (vpo[vpi[i]].z() < minz) minz = CGAL::to_double(vpo[vpi[i]].z());
	}
	return mbox_3(minx-offset, miny-offset, minz-offset, maxx+offset, maxy+offset, maxz+offset);
}

mbox_3 build_mbox_3 (const poly_frac &df, const double &offset)
{
	return build_mbox_3(df.vp_ob,offset);
}

mbox_3 build_mbox_3 (const point_3 &p, const plane_3&pl, const FT &r_, const double &offset)
{
	using CGAL::to_double;
	double a(to_double(pl.orthogonal_vector().x()));
	double b(to_double(pl.orthogonal_vector().y()));
	double c(to_double(pl.orthogonal_vector().z()));
	double res = sqrt(a*a+b*b+c*c);
	a/=res; b/=res; c/=res;
	double r(to_double(r_)), x(to_double(p.x())), y(to_double(p.y())), z(to_double(p.z()));
	return mbox_3(x-r*sqrt(1-a*a)-offset, y-r*sqrt(1-b*b)-offset, z-r*sqrt(1-c*c)-offset,
		x+r*sqrt(1-a*a)+offset, y+r*sqrt(1-b*b)+offset, z+r*sqrt(1-c*c)+offset);
}

mbox_3 build_mbox_3 (const std::vector<point_3> &vpo, const eblock &eb, const double &offset)
{
	double minx(CGAL::to_double(vpo[eb.vf[0].np[0]].x())), maxx(CGAL::to_double(vpo[eb.vf[0].np[0]].x()));
	double miny(CGAL::to_double(vpo[eb.vf[0].np[0]].y())), maxy(CGAL::to_double(vpo[eb.vf[0].np[0]].y()));
	double minz(CGAL::to_double(vpo[eb.vf[0].np[0]].z())), maxz(CGAL::to_double(vpo[eb.vf[0].np[0]].z()));
	for (int i(0);i!=eb.vf.size ();++i)
		for (int j(0);j!=eb.vf[i].np.size ();++j)
		{
			if (vpo[eb.vf[i].np[j]].x() > maxx) maxx = CGAL::to_double(vpo[eb.vf[i].np[j]].x());
			if (vpo[eb.vf[i].np[j]].x() < minx) minx = CGAL::to_double(vpo[eb.vf[i].np[j]].x());
			if (vpo[eb.vf[i].np[j]].y() > maxy) maxy = CGAL::to_double(vpo[eb.vf[i].np[j]].y());
			if (vpo[eb.vf[i].np[j]].y() < miny) miny = CGAL::to_double(vpo[eb.vf[i].np[j]].y());
			if (vpo[eb.vf[i].np[j]].z() > maxz) maxz = CGAL::to_double(vpo[eb.vf[i].np[j]].z());
			if (vpo[eb.vf[i].np[j]].z() < minz) minz = CGAL::to_double(vpo[eb.vf[i].np[j]].z());
		}
	return mbox_3(minx-offset, miny-offset, minz-offset, maxx+offset, maxy+offset, maxz+offset);
}

mbox_3 build_mbox_3 (const block_system &bs, const cblock &cb, const double offset)
{
	//慎用，如果cb包含的单元块体过多的话会很慢
	mbox_3 res;
	if (cb.vb.empty ())
		return res;
	res = build_mbox_3(bs.vpo,cb.vb[0]);
	for (int i(0);i!=cb.vb .size ();++i)
	{
		res += build_mbox_3(bs.vpo, cb.vb[i]);
	}
	return res;
}

bool if_contained_mbox_3 (const mbox_3 &tar, const mbox_3 &b)
{
	for (int i(0);i!=3;++i)
	{
		if (tar.max (i) < b.max (i) || tar.min (i) > b.min(i))
			return false;
	}
	return true;
}

polygon_2 change_proj_plane (const plane_3 &oldp, const plane_3 &newp, const polygon_2 &tar)
{
	polygon_2 temp_p;
	for (int i(0);i!=tar.size ();++i)
		temp_p.push_back (newp.to_2d(oldp.to_3d(tar[i])));
	return temp_p;
}

polygon_with_holes_2 change_proj_plane (const plane_3 &oldp, const plane_3 &newp, const polygon_with_holes_2 &tar)
{
	polygon_with_holes_2 temp_p(change_proj_plane(oldp,newp,tar.outer_boundary ()));
	for (polygon_with_holes_2::Hole_const_iterator hi = tar.holes_begin(); hi!=tar.holes_end();++hi)
		temp_p.add_hole(change_proj_plane(oldp, newp,*hi));
	return temp_p;
}

polygon_set_2 change_proj_plane (const plane_3 &oldp, const plane_3 &newp, const polygon_set_2 &tar)
{
	vector<polygon_with_holes_2> temp_vp;
	temp_vp.reserve (tar.number_of_polygons_with_holes());
	tar.polygons_with_holes(back_inserter(temp_vp));
	for (int i(0);i!=temp_vp.size ();++i)
		temp_vp[i] = change_proj_plane(oldp,newp,temp_vp[i]);
	polygon_set_2 temp_set;
	for (int i(0);i!=temp_vp.size ();++i)
	{
		temp_set.insert (temp_vp[i]);
	}
	return temp_set;
}

poly_frac disc2poly (const disc_frac & df,  const plane_3 &pl, const int &n, const FT &dx)
{
	vector_3 b1 (pl.base1 ());
	vector_3 b2 (pl.base2());
	b1 = (1/sqrt(CGAL::to_double(b1.squared_length()))) * b1;	//既然将圆盘近似为多边形本身就是一个近似，那么就用近似函数吧
	b2 = (1/sqrt(CGAL::to_double(b2.squared_length()))) * b2;

	double r = CGAL::to_double(df.r);
	double da = M_PI_2/n;
	double a = 0;

	point_3 center (pl.projection(df.center));
	poly_frac res;
	for (int i(0);i!=n;++i)
	{
		res.vp_ob.push_back(center + df.r*b1*cos(a) + df.r*b2*sin(a));
		a += da;
	}
	res.proj_plane(pl);
	return res;
}

int construct_poly_frac (poly_frac & pf, const polygon_with_holes_2 &pwh, const plane_3& pl, const vector_3& off_set, const bool &inverse)
{
	pf.vp_ob.clear();
	pf.vp_ib.clear();
	pf.vp_ob.reserve (pwh.outer_boundary ().size ());
	pf.vp_ib.reserve (pwh.number_of_holes());
	const polygon_2 &outer = pwh.outer_boundary ();
	for (int i(0);i!=(int) outer.size ();++i)
	{
		point_3 p_temp(pl.to_3d(outer[i]));
		if (inverse)
			pf.vp_ob.push_back (point_3(-p_temp.x(),-p_temp.y(),-p_temp.z()) + off_set);
		else
			pf.vp_ob.push_back (p_temp + off_set);
	}

	for (polygon_with_holes_2::Hole_const_iterator hi = pwh.holes_begin ();
		hi!=pwh.holes_end();++hi)
	{
		pf.vp_ib.push_back (poly_frac::v_loop());
		pf.vp_ib.back ().reserve (hi->size ());
		for (int i(0);i!=hi->size ();++i)
		{
			point_3 p_temp(pl.to_3d((*hi)[i]));
			if(inverse)
				pf.vp_ib.back ().push_back (point_3(-p_temp.x(),-p_temp.y(),-p_temp.z()) + off_set);
			else
				pf.vp_ib.back ().push_back (p_temp + off_set);
		}
	}
	return 0;
}

int output_STL_poly (std::ostream &outfile, const polygon_with_holes_2 &poly, const plane_3& pl, const bool &normal_dir, const FT&tranx, const FT&trany, const FT&tranz, const bool &flag)
{
	if (flag)
		outfile<<"solid Polygon"<<endl;
	typedef struct
	{
		int nesting_level;
	} face_info;

	using CGAL::to_double;
	typedef CGAL::Triangulation_vertex_base_2<K>					Vb;
	typedef CGAL::Triangulation_face_base_with_info_2<face_info,K>	Fbb;
	typedef CGAL::Constrained_triangulation_face_base_2 <K, Fbb>	Fb;
	typedef CGAL::Triangulation_data_structure_2<Vb,Fb>				TDS;
	typedef CGAL::Exact_intersections_tag							TAG;
	typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS,TAG>	CDT;

	outfile.precision (20);
	CDT cdt;
	vector_3 normal (pl.orthogonal_vector());
	if (!normal_dir)
		normal = -normal;
	double a(to_double (normal.x ()));
	double b(to_double (normal.y ()));
	double c(to_double (normal.z ()));
	double d(sqrt(a*a+b*b+c*c));
	a/=d; b/=d; c/=d;

	if (poly.is_unbounded())
		throw logic_error ("outer_bound::output_STL, poly is unbounded");
	insert_polygon(cdt,poly.outer_boundary());
	for (polygon_with_holes_2::Hole_const_iterator j =poly.holes_begin();j!=poly.holes_end();++j)
		insert_polygon(cdt,*j);
	mark_domains(cdt);
	for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit!=cdt.finite_faces_end();++fit)
	{
		if (fit->info().nesting_level%2 !=1)
			continue;
		vector<point_3> vp;
		vp.push_back (pl.to_3d(fit->vertex(0)->point()));
		vp.push_back (pl.to_3d(fit->vertex(1)->point()));
		vp.push_back (pl.to_3d(fit->vertex(2)->point()));
		if (normal*oriented_normal(vp) < 0)
			reverse(vp.begin (),vp.end());
		outfile<<"  facet normal "<<a<<' '<<b<<' '<<c<<endl;
		outfile<<"    outer loop"<<endl;
		outfile<<"      vertex "<<vp[0].x() + tranx<<' '<<vp[0].y() + trany<<' '<<vp[0].z() + tranz<<endl;
		outfile<<"      vertex "<<vp[1].x() + tranx<<' '<<vp[1].y() + trany<<' '<<vp[1].z() + tranz<<endl;
		outfile<<"      vertex "<<vp[2].x() + tranx<<' '<<vp[2].y() + trany<<' '<<vp[2].z() + tranz<<endl;
		outfile<<"    endloop"<<endl;
		outfile<<"  endfacet"<<endl;
	}
	if (flag)
		outfile<<"endsolid Polygon"<<endl;
	return 0;
}

}