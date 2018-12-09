#include "inner_poly.h"
#include <iostream>
#include <CGAL\Constrained_Delaunay_triangulation_2.h>
#include <CGAL\Triangulation_face_base_with_info_2.h>
#include <CGAL\Segment_3.h>
using namespace std;

namespace BI
{

int inner_bound::construct_inner_bound (const block_system &bs, const std::vector<int> &vebi, const bool& debug_info)
{
	if (debug_info)
		cout<<"constructing common planes"<<endl;
	vpoly.clear();
	vplane.clear();
	vector<int> plane_index;
	vector<polygon_set_2> union_eb;		//记录于plane_index[i]相邻的eb的面的交集
	vector<polygon_set_2> outer_b;		//记录plane_index[i]上的外表面
	vector<polygon_set_2> frac_poly;	//用于记录在plane_index[i]平面上裂隙的形状

	for (int i(0);i!=(int)vebi.size ();++i)
	{
		const eblock& eb = bs.veb[vebi[i]];
		for (int j(0);j!=eb.vf.size ();++j)
		{
			int index = add_set (plane_index, eb.vf[j].plane_id);
			if (plane_index.size () > union_eb.size ())
			{
				//表明在plane_index中添加了新的平面
				union_eb.push_back (polygon_set_2());
				outer_b.push_back (polygon_set_2());
				frac_poly.push_back (polygon_set_2());
			}
			polygon_2 poly_temp;
			if (proj_polygon(bs.vpl[eb.vf[j].plane_id], bs.vpo, eb.vf[j].np,poly_temp) !=0)
				throw logic_error ("construct_inner_boundary, projection polygon false");
			if (poly_temp.is_clockwise_oriented())
				poly_temp.reverse_orientation ();
			union_eb[index].join (poly_temp);
			outer_b[index].symmetric_difference(poly_temp);
		}
	}
	if (debug_info)
		cout<<"Generating fractures"<<endl;
	for (int i(0);i!=(int)bs.vdisc.size ();++i)
	{
		int index = find_n(plane_index,bs.vdisc[i].plane_id);
		if (index == plane_index.size ())
			continue;		//这个裂隙不在plane_index中，意思是这个裂隙不可能参与这个块体的形成
		polygon_2 poly_temp;
		proj_polygon(bs.vpl[bs.vdisc[i].plane_id],bs.vdisc[i],poly_temp);
		frac_poly[index].join (poly_temp);
	}
	for (int i(0);i!=(int)bs.vpoly.size ();++i)
	{
		int index = find_n(plane_index,bs.vpoly[i].plane_id);
		if (index == plane_index.size ())
			continue;
		polygon_with_holes_2 poly_temp;
		proj_polygon(bs.vpl[bs.vpoly[i].plane_id],bs.vpoly[i],poly_temp);
		frac_poly[index].join (poly_temp);
	}

	if (debug_info)
		cout<<"Generating inner boundaries"<<endl;
	for (int i(0);i!=plane_index.size ();++i)
	{
		frac_poly[i].intersection(union_eb[i]);
		frac_poly[i].difference (outer_b[i]);
		if (!frac_poly[i].is_valid())
			throw logic_error ("construct_inner_boundary, result poly set is invaild");
		if (frac_poly[i].is_empty ())
			continue;		//如果是空的直接跳出就好
		
		vector<polygon_with_holes_2> temp_pwh;
		int vp_size = (int) vpoly.size ();
		frac_poly[i].polygons_with_holes(std::back_inserter(temp_pwh));
		for (int j(0);j!=(int) temp_pwh.size ();++j)
		{
			vector<polygon_with_holes_2> split_pwh;
			if (split_polygon_2(temp_pwh[j],split_pwh) == 0)
			{
				vpoly.push_back (polygon_with_holes_2());
				simplify_polygon_2(temp_pwh[j],vpoly.back ());
			}
			else
			{
				for (int k(0);k!=(int)  split_pwh.size ();++k)
				{
					vpoly.push_back (polygon_with_holes_2());
					simplify_polygon_2(split_pwh[k],vpoly.back ());
				}
			}
		}
		for (; vp_size<vpoly.size ();++vp_size)
			vplane.push_back (bs.vpl[plane_index[i]]);
	}

	if (debug_info)
		cout<<"Generating inner boundaries complete"<<endl;
	return 0;
}

int inner_bound::construct_inner_bound (const block_system &bs, const int &cbi, const bool& debug_info)
{
	return construct_inner_bound (bs,bs.vcb[cbi].vb,debug_info);
}

int inner_bound::output_STL_poly (std::ostream &outfile, const int &pi, const bool& flag, const FT& tranx, const FT& trany, const FT& tranz) const
{
	if (flag)
		outfile<<"solid Polygon"<<pi<<endl;
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

	outfile.precision(25);
	CDT cdt;
	vector_3 normal (vplane[pi].orthogonal_vector());

	double a(to_double (normal.x ()));
	double b(to_double (normal.y ()));
	double c(to_double (normal.z ()));
	double d(sqrt(a*a+b*b+c*c));
	a/=d; b/=d; c/=d;
	const polygon_with_holes_2 &poly = vpoly[pi];
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
		vp.push_back (vplane[pi].to_3d(fit->vertex(0)->point()));
		vp.push_back (vplane[pi].to_3d(fit->vertex(1)->point()));
		vp.push_back (vplane[pi].to_3d(fit->vertex(2)->point()));
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

int inner_bound::output_STL ( std::ostream &outfile, const FT& tranx, const FT& trany, const FT& tranz,const std::string& ob_name) const
{
	if (!outfile)
		return 1;
	outfile<<"solid "<<ob_name<<endl;
	for (int i(0);i!=vpoly.size ();++i)
		output_STL_poly(outfile, i,false,tranx,trany,tranz);
	outfile<<"endsolid "<<ob_name<<endl;
	return 0;
}

int inner_bound::output_exp_STL (std::ostream &outfile, const double &ratio , const FT& tranx, const FT& trany, const FT& tranz,const std::string& ob_name) const
{
	if (!outfile)
		return 1;

	vector<point_3> vp;
	vp.reserve (vpoly.size ());
	for (int i(0);i!=vpoly.size ();++i)
		vp.push_back (get_poly_center(i));
	double x(0), y(0), z(0);
	for (int i(0);i!=vp.size ();++i)
	{
		x += CGAL::to_double(vp[i].x());
		y += CGAL::to_double(vp[i].y());
		z += CGAL::to_double(vp[i].z());
	}
	x /= double(vp.size ());
	y /= double(vp.size ());
	z /= double(vp.size ());
	for (int i(0);i!=vpoly.size ();++i)
		output_STL_poly(outfile,i,true,(vp[i].x()-x)*ratio+tranx, (vp[i].y()-y)*ratio+trany, (vp[i].z()-z)*ratio+tranz);
	
	return 0;
}

point_3 inner_bound::get_poly_center (const int &pi) const
{
	FT x(0), y(0), z(0);
	for (polygon_2::Vertex_const_iterator vi = vpoly[pi].outer_boundary().vertices_begin();
		vi != vpoly[pi].outer_boundary().vertices_end();++vi)
	{
		point_3 p_temp (vplane[pi].to_3d(*vi));
		x = x + p_temp.x();
		y = y + p_temp.y();
		z = z + p_temp.z();
	}
	return point_3 (x/FT(vpoly[pi].outer_boundary ().size ()), y/FT(vpoly[pi].outer_boundary ().size ()), z/FT(vpoly[pi].outer_boundary ().size ()));
}
}