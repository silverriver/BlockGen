#include "outer_poly.h"
#include <iostream>
#include <list>
#include <string>
#include <CGAL\Constrained_Delaunay_triangulation_2.h>
#include <CGAL\Triangulation_face_base_with_info_2.h>
#include <CGAL\Segment_3.h>
using namespace std;

namespace BI
{
int outer_bound::construct_outer_bound (const block_system &bs, const std::vector <int> &vebi, const bool& debug_info)
{
	if (debug_info)
		cout<<"Begin to construct outer boundaries"<<endl<<"Constructing common planes"<<endl;
	vpoly.clear();	vplane.clear();	outer_normal.clear();
	vector<int> plane_index;			//计算每个平面在bs.vpl中的索引
	vector<polygon_set_2> p_poly;		//表示索引为plane_index[i]的平面正侧的多边形
	vector<polygon_set_2> n_poly;		//表示索引为plane_index[i]的平面负侧的多边形
	
	for (int i(0);i!=vebi.size ();++i)
	{
		const eblock& eb = bs.veb[vebi[i]];
		for (int j(0);j!=eb.vf.size ();++j)
		{
			int index = add_set(plane_index,eb.vf[j].plane_id);
			if (plane_index.size () > p_poly.size())
			{
				//表明在plane_index中添加了新的平面
				p_poly.push_back (polygon_set_2()); 
				n_poly.push_back (polygon_set_2());
			}
			polygon_2 poly_temp;
			if (proj_polygon(bs.vpl[eb.vf[j].plane_id],bs.vpo,eb.vf[j].np,poly_temp) != 0)
				throw logic_error ("construct_outer_boundary, projection polygon false");
			if (poly_temp.is_clockwise_oriented())
				poly_temp.reverse_orientation();
			if (eb.vf[j].outer_normal)
				n_poly[index].join (poly_temp);
			else
				p_poly[index].join (poly_temp);
		}
	}
	if (debug_info)
		cout<<"Generating polygons"<<endl;
	for (int i(0);i!=plane_index.size ();++i)
	{
		polygon_set_2 temp_set (n_poly[i]);		//这里可以直接用symmetric_difference这个运算，但是这样的话输出外法线方向的时候还需要单独判断各个多边形的外法线方向
		n_poly[i].difference(p_poly[i]);
		p_poly[i].difference(temp_set);

		if ((!n_poly[i].is_valid()) || (!p_poly[i].is_valid()))
			throw logic_error ("construct_outer_boundary, result poly set is invaild");
		if (n_poly[i].is_empty() && p_poly[i].is_empty())
			continue;	//如果是空的就直接跳出
		else
		{
			//polygon_set_2 n_poly_simplified, p_poly_simplified;
			vector<polygon_with_holes_2> v_n_poly, v_p_poly;
			int vp_size = (int)vpoly.size ();
			n_poly[i].polygons_with_holes(std::back_inserter(v_n_poly));
			for (int j(0);j!=(int)v_n_poly.size ();++j)
			{
				vector<polygon_with_holes_2> vpwh;
				if (split_polygon_2(v_n_poly[j],vpwh) == 0)		//先试着将每个多边形剖分一下，剖分成简单多边形
				{
					vpoly.push_back (polygon_with_holes_2());
					simplify_polygon_2(v_n_poly[j],vpoly.back ());		//添加前一定要记得要简化一下边界
				}
				else
				{
					for (int k(0);k!=(int)vpwh.size ();++k)
					{
						vpoly.push_back (polygon_with_holes_2());
						simplify_polygon_2(vpwh[k],vpoly.back ());
					}
				}
			}
			for (;vp_size<vpoly.size ();++vp_size)
			{
				vplane.push_back (bs.vpl[plane_index[i]]);
				outer_normal.push_back (true);
			}

			p_poly[i].polygons_with_holes(back_inserter(v_p_poly));
			for (int j(0);j!=(int) v_p_poly.size ();++j)
			{
				vector<polygon_with_holes_2> vpwh;
				if (split_polygon_2(v_p_poly[j],vpwh) == 0)
				{
					vpoly.push_back (polygon_with_holes_2());
					simplify_polygon_2(v_p_poly[j],vpoly.back ());
				}
				else
				{
					for (int k(0);k!=(int) vpwh.size ();++k)
					{
						vpoly.push_back (polygon_with_holes_2());
						simplify_polygon_2(vpwh[k],vpoly.back ());
					}
				}
			}
			//下面是以前的版本，在之前的版本里面没有split各个多边形
			//simplify_polygon_2(p_poly[i],p_poly_simplified);		//简化一下生成的多边形
			//p_poly_simplified.polygons_with_holes(std::back_inserter(vpoly));
			for (;vp_size<vpoly.size ();++vp_size)
			{
				vplane.push_back (bs.vpl[plane_index[i]]);
				outer_normal.push_back (false);
			}
		}
	}
	
	if (debug_info)
		cout<<"Constructing coner points"<<endl;
	
	integrate_boundary();
	if (debug_info)
		cout<<"Generating polygons complete"<<endl;

	return 0;
}

int outer_bound::construct_outer_bound(const block_system &bs, const int &cbi, const bool& debug_info)
{
	return construct_outer_bound (bs,bs.vcb[cbi].vb,debug_info);
}

int outer_bound::construct_outer_bound_using_vcb (const block_system &bs, const std::vector <int> &vcbi, const bool& debug_info)
{
	vector<int> veb;
	for (int i(0); i!=vcbi.size ();++i)
	{
		veb.insert (veb.end (), bs.vcb[vcbi[i]].vb.begin (), bs.vcb[vcbi[i]].vb.end ());
	}
	return construct_outer_bound (bs, veb,debug_info);
}

int outer_bound::output_STL_poly (std::ostream &outfile, const int &pi, const bool& flag,const FT& tranx, const FT& trany, const FT& tranz) const
{
	return BI::output_STL_poly(outfile,vpoly[pi],vplane[pi],outer_normal[pi],tranx,trany,tranz,flag);
}

int outer_bound::output_STL (std::ostream &outfile,const FT& tranx, const FT& trany, const FT& tranz, const std::string& ob_name) const
{	
	if (!outfile)
		return 1;
	if (ob_name.empty())
		outfile<<"solid OuterBoundary"<<endl;
	else
		outfile<<"solid "<<ob_name<<endl;
	for (int i(0);i!=vpoly.size ();++i)
		output_STL_poly(outfile, i,false,tranx,trany,tranz);

	if (ob_name.empty())
		outfile<<"endsolid OuterBoundary"<<endl;
	else
		outfile<<"endsolid "<<ob_name<<endl;
	return 0;
}

int outer_bound::output_exp_STL (std::ostream &outfile, const double &ratio,const FT& tranx, const FT& trany, const FT& tranz, const std::string& ob_name) const
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

int outer_bound::output_CAD_script (std::ostream &outfile,const bool& split_ploy,const bool &output_normal) const
{
	using CGAL::to_double;
	if(!outfile)
		return 1;
	outfile.precision(20);
	for (int i(0);i!=vpoly.size ();++i)
	{
		vector_3 normal (vplane[i].orthogonal_vector());
		if (!outer_normal[i])
			normal = -normal;
		double a(to_double (normal.x ()));
		double b(to_double (normal.y ()));
		double c(to_double (normal.z ()));
		double d(sqrt(a*a+b*b+c*c));
		a/=d; b/=d; c/=d;
		const polygon_2 &outer_p (vpoly[i].outer_boundary());
		for (polygon_2::Edge_const_iterator j=outer_p.edges_begin();j!=outer_p.edges_end();++j)
			out_line_3(outfile,vplane[i].to_3d(j->source ()), vplane[i].to_3d(j->end ()));

		for (polygon_with_holes_2::Hole_const_iterator j = vpoly[i].holes_begin(); j!=vpoly[i].holes_end();++j)
		{
			for (polygon_2::Edge_const_iterator k = j->edges_begin();k!=j->edges_end();++k)
				out_line_3 (outfile,vplane[i].to_3d(k->source()),vplane[i].to_3d(k->end ()));
		}
		if (split_ploy)
			outfile<<endl;
		if (output_normal)
		{
			out_line_3(outfile,vplane[i].to_3d(outer_p[0]),vplane[i].to_3d(outer_p[0])+vector_3(a,b,c));
			outfile<<endl;
		}
	}
	return 0;

}

int outer_bound::output_BlockAni_format (std::ostream &outfile, const int &id) const
{
	if (!outfile)
		return 1;
	outfile<<"BLOCK "<<id<<endl;
	outfile.precision (15);
	for (int i(0);i!=vplane.size ();++i)
	{
		vector_3 normal (vplane[i].orthogonal_vector());
		if (!outer_normal[i])
			normal = -normal;
		outfile<<"R "<<normal.x()<<" "<<normal.y()<<" "<<normal.z()<<endl;
		const polygon_2 &outer(vpoly[i].outer_boundary());
		outfile<<"L ";
		for (polygon_2::Vertex_const_iterator vi = outer.vertices_begin();
			vi!=outer.vertices_end();++vi)
		{
			point_3 p (vplane[i].to_3d(*vi));
			outfile<<p.x()<<" "<<p.y()<<" "<<p.z()<<endl;
		}
		for (polygon_with_holes_2::Hole_const_iterator hi = vpoly[i].holes_begin();
			hi!=vpoly[i].holes_end();++hi)
		{
			outfile<<"L ";
			for (polygon_2::Vertex_const_iterator vi = hi->vertices_begin();
				vi!=hi->vertices_end();++vi)
			{
				point_3 p(vplane[i].to_3d(*vi));
				outfile<<p.x()<<" "<<p.y()<<" "<<p.z()<<endl;
			}
		}
	}
	return 0;
}

point_3 outer_bound::get_center () const
{
	FT minx,miny,minz,maxx,maxy,maxz;
	for (int i(0);i!=vpoly.size ();++i)
	{
		const polygon_2 &poly(vpoly[i].outer_boundary());
		if (i==0)
		{
			point_3 temp_p(vplane[i].to_3d(poly[0]));
			minx = maxx = temp_p.x();
			miny = maxy = temp_p.y();
			minz = maxz = temp_p.z();
		}
		for (polygon_2::Vertex_iterator j(poly.vertices_begin());
			j!=poly.vertices_end();++j)
		{
			point_3 temp_p(vplane[i].to_3d(*j));
			if (temp_p.x() < minx)
				minx = temp_p.x();
			if (temp_p.x() > maxx)
				maxx = temp_p.x();
			if (temp_p.y() < miny)
				miny = temp_p.y();
			if (temp_p.y() > maxy)
				maxy = temp_p.y();
			if (temp_p.z() < minz)
				minz = temp_p.z();
			if (temp_p.z() > maxz)
				maxz = temp_p.z();
		}
	}
	return point_3((minx+maxx)/2, (miny+maxy)/2, (minz+maxz)/2);
}

int outer_bound::add_coner_point (const int &poly_n, point_3 **vpoint,  int  *p_count,  mbox_3 *vbox)
{
	//寻找位于平面vplane[poly_n]上的点，为了方便直接记录这些点的坐标，不记录索引
	vector<point_3> v_coner_p;
	v_coner_p.reserve(vpoly.size ()/3);
	for (int i(0);i!=(int)vpoly.size ();++i)
	{
		if (!CGAL::do_overlap(vbox[i],vbox[poly_n]))		//如果两个多边形都不相交，则跳过。注意如果此时这两个多边形是同一个，也需要添加,因为可能处理的是relative simple的多边形
			continue;
		for (int j(0);j!=p_count[i];++j)
			if (vplane[poly_n].has_on(vpoint[i][j]))
				add_set(v_coner_p,vpoint[i][j]);		//记录这个点的坐标
	}

	polygon_2 poly_temp;
	int point_added (0);		//记录在多边形边界上添加了几个点
	point_added = add_point_on_poly(vpoly[poly_n].outer_boundary (), vplane[poly_n], v_coner_p, poly_temp);
	if (point_added != 0)
		vpoly[poly_n].outer_boundary() = poly_temp;
	
	for (polygon_with_holes_2::Hole_iterator holei = vpoly[poly_n].holes_begin(); holei != vpoly[poly_n].holes_end(); ++holei)
	{
		int pa_temp(0);
		pa_temp += add_point_on_poly(*holei,vplane[poly_n],v_coner_p,poly_temp);
		if (pa_temp != 0)
			*holei = poly_temp;
		point_added += pa_temp;
	}
	return point_added;
}

point_3 outer_bound::get_poly_center (const int &pi) const
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

int outer_bound::add_point_on_poly (const polygon_2 &original, const plane_3& plane, const vector<point_3> & vconer, polygon_2 &new_poly)
{
	new_poly.clear();
	int point_added(0);
	vector<point_3> poly_3d;		//转化到三维空间
	poly_3d.reserve (original.size ());
	for (int i(0);i!=original.size ();++i)
		poly_3d.push_back (plane.to_3d(original[i]));

	for (int i(0); i!=original.size ();++i)
	{
		new_poly.push_back (original[i]);
		int i1 (i+1);
		if (i1>=(int)original.size ())
			i1 = 0;
		vector<int> vp;		//记录在线段i,i1上的需要添加的点在vconer上的索引
		vp.reserve (vconer.size ());
		segment_3 seg(poly_3d[i], poly_3d[i1]);
		for (int j(0);j!=(int)vconer.size ();++j)
			if (seg.has_on(vconer[j]) && vconer[j]!=seg.start() && vconer[j]!=seg.end ())	//给定点在线段上，并且不在端点
				vp.push_back (j);
		
		if (vp.empty ()) continue;
		point_added += (int)vp.size ();
		vector<FT> dis;		//计算vp中各个点离poly_3d[i]的距离
		dis.reserve (vp.size ());
		for (int j(0);j!=(int)vp.size ();++j)
			dis.push_back (CGAL::squared_distance(poly_3d[i],vconer[vp[j]]));
		vector<int> sorted_index;
		sort_index(dis,sorted_index);	//给需要添加的点排一下序，按照离poly_3d[i]从近到远排序
		for (int j(0);j!=sorted_index.size ();++j)		//依次添加各个点
			new_poly.push_back (plane.to_2d(vconer[vp[sorted_index[j]]]));
	}
	return point_added;		//构建完成
}

int outer_bound::join_veb (const block_system &bs, const std::vector <int> &vebi, const bool& debug_info)
{
	//啊，fuck，本来想要实现一个更加简便的算法，但是后来实现的时候发现想要实现这个算法会很麻烦，
	//并且比较容易出bug，现在为了省事，还是直接用相对笨一点的方法，比较安全省事一点
	return join(outer_bound(bs,vebi),debug_info);
}

int outer_bound::join_cb (const block_system &bs, const int cbi, const bool& debug_info)
{
	return join_veb(bs,bs.vcb[cbi].vb,debug_info);
}

int outer_bound::join (const outer_bound &ob, const bool &debug_info)
{
	//输入函数中的ob和*this可能来自于不同额块体系统，因此其所包含的vplane中的平面可能会存在反向的情况，需要注意
	vector<vector<int>> common_pl;	//找出this和ob的vplane中一共有多少不重合的平面。common_pl[i]中的平面都互相重合，这些平面用索引表示，其中ob中的平面索引在原来值的基础上加上this->vplane.size().如果ob中的某个平面和this中的某个平面相反，则将其索引标记为负数
	common_pl.reserve (vplane.size () + ob.vplane.size ());		//common_pl[i][0]是新的平面vector
	for (int i(0);i!=(int)vplane.size ();++i)
	{
		int count(0);
		for (;count!=(int)common_pl.size ();++count)
		{
			if (vplane[i] == vplane[common_pl[count][0]])
			{
				common_pl[count].push_back (i);
				break;
			}
		}
		if (count == common_pl.size ())
		{
			common_pl.push_back(vector<int>());
			common_pl.back ().push_back (i);
		}
	}
	for (int i(0);i!=(int)ob.vplane.size ();++i)
	{
		int count (0);
		for (;count!=(int)common_pl.size ();++count)
		{
			const plane_3 * p_temp;
			if (common_pl[count][0]<vplane.size ())		//注意common_pl[i][0]不可能小于0，在this->vpoly为空的时候才可能等于0
				p_temp = &vplane[common_pl[count][0]];	//如果common_pl[count]中有this中的平面
			else 
				p_temp = &ob.vplane [common_pl[count][0] - vplane.size ()];		//如果common_pl中没有this中的平面
				
			if(ob.vplane [i] == *p_temp)
			{
				common_pl[count].push_back (i+(int)vplane.size ()); break;
			}
			if (ob.vplane [i].opposite() == *p_temp)
			{
				common_pl[count].push_back (-i-(int)vplane.size ()); break;
			}
		}
		if (count == common_pl.size ())
		{
			common_pl.push_back (vector<int> ());
			common_pl.back ().push_back (i+(int)vplane.size ());		//如果this->vplane中没有这个平面的话就添加一下。
		}
	}

	vector<polygon_set_2> n_poly(common_pl.size ());		//表示在common_pl[i]所表示平面负侧的多边形
	vector<polygon_set_2> p_poly(common_pl.size ());		//表示在common_pl[i]所表示平面正侧的多边形 （外法线方向与该平面相反）
	for (int i(0);i!=(int)common_pl.size ();++i)
	{
		for (int j(0);j!=(int)common_pl[i].size ();++j)
		{
			if (abs(common_pl[i][j]) < vplane.size ())			//该平面在this中
				if (outer_normal[common_pl[i][j]])
					n_poly[i].join(vpoly[common_pl[i][j]]);		//将外法线方向与common_pl[i]所表示平面同向的多边形添加到n_poly中
				else 
					p_poly[i].join(vpoly[common_pl[i][j]]);
			else												//该平面在ob中。
				if (common_pl[i][j] >= 0)						//这种情况下表示这个平面和common_pl[i][0]相同。只有当this->vpoly.size() == 0 时才有可能等于0
					if (ob.outer_normal[common_pl[i][j]-vplane.size ()])
						n_poly[i].join (ob.vpoly[common_pl[i][j]-vplane.size ()]);
					else
						p_poly[i].join (ob.vpoly[common_pl[i][j]-vplane.size ()]);
				else											//这种情况下表示这个平面和common_pl[i][0]相反。这种情况下common_pl[i][0]只能是this中的平面，不可能是ob中的平面
					if (ob.outer_normal[-common_pl[i][j]-vplane.size ()])
						p_poly[i].join (change_proj_plane(ob.vplane[-common_pl[i][j]-vplane.size ()],vplane[common_pl[i][0]], ob.vpoly[-common_pl[i][j]-vplane.size ()]));		//总之这种情况指的是ob中相应的平面和this中方向相反，这个时候反过来添加.注意，这个时候应该重新生成一下目标多边形，待完成2015-10-28 11:59
					else
						n_poly[i].join (change_proj_plane(ob.vplane[-common_pl[i][j]-vplane.size ()],vplane[common_pl[i][0]], ob.vpoly[-common_pl[i][j]-vplane.size ()]));
		}
	}

	if (debug_info)
		cout<<"Generating polygon"<<endl;
	vector<plane_3> new_vplane;		//新的外表面边界对应的平面
	new_vplane.reserve (vplane.size ()+ob.vplane.size ());
	vector<polygon_with_holes_2> new_vpoly;		//新的外表面多边形数组
	new_vpoly.reserve (vplane.size ()+ob.vplane.size());
	vector<bool> new_outernormal;
	new_outernormal.reserve (vplane.size() + ob.vplane.size ());
	for (int i(0);i!=common_pl.size ();++i)
	{
		polygon_set_2 temp_set(n_poly[i]);
		n_poly[i].difference (p_poly[i]);
		p_poly[i].difference (temp_set);
		
		if ((!n_poly[i].is_valid()) || (!p_poly[i].is_valid ()))
			throw logic_error ("join, result poly set is invaild");
		if (n_poly[i].is_empty () && p_poly[i].is_empty ())
			continue;
		else
		{
			vector<polygon_with_holes_2> v_n_poly, v_p_poly;
			int vp_size = (int)new_vpoly.size ();
			n_poly[i].polygons_with_holes(std::back_inserter(v_n_poly));
			for (int j(0);j!=(int)v_n_poly.size ();++j)
			{
				vector<polygon_with_holes_2> vpwh;
				if (split_polygon_2(v_n_poly[j],vpwh) == 0)		//先试着将每个多边形剖分一下，剖分成简单多边形
				{
					new_vpoly.push_back (polygon_with_holes_2());
					simplify_polygon_2(v_n_poly[j],vpoly.back ());		//添加前一定要记得要简化一下边界
				}
				else
				{
					for (int k(0);k!=(int)vpwh.size ();++k)
					{
						new_vpoly.push_back (polygon_with_holes_2());
						simplify_polygon_2(vpwh[k],vpoly.back ());
					}
				}
			}
			for (;vp_size<new_vpoly.size ();++vp_size)
			{
				if (common_pl[i][0]<vpoly.size ())
					new_vplane.push_back (vplane[common_pl[i][0]]);
				else
					new_vplane.push_back (ob.vplane[common_pl[i][0]-vplane.size ()]);
				new_outernormal.push_back (true);
			}

			p_poly[i].polygons_with_holes(back_inserter(v_p_poly));
			for (int j(0);j!=(int) v_p_poly.size ();++j)
			{
				vector<polygon_with_holes_2> vpwh;
				if (split_polygon_2(v_p_poly[j],vpwh) == 0)
				{
					new_vpoly.push_back (polygon_with_holes_2());
					simplify_polygon_2(v_p_poly[j],vpoly.back ());
				}
				else
				{
					for (int k(0);k!=(int) vpwh.size ();++k)
					{
						new_vpoly.push_back (polygon_with_holes_2());
						simplify_polygon_2(vpwh[k],vpoly.back ());
					}
				}
			}
			//simplify_polygon_2(p_poly[i],p_poly_simplified);
			//p_poly_simplified.polygons_with_holes(std::back_inserter(new_vpoly));
			for (;vp_size<new_vpoly.size ();++vp_size)
			{
				if(common_pl[i][0] < vpoly.size ())
					new_vplane.push_back (vplane[common_pl[i][0]]);
				else
					new_vplane.push_back (ob.vplane[common_pl[i][0]-vplane.size ()]);
				new_outernormal.push_back (false);
			}
		}
	}
	vpoly.swap (new_vpoly);
	vplane.swap (new_vplane);
	outer_normal.swap (new_outernormal);
	if (debug_info)
		cout<<"Constructing coner points"<<endl;
	
	integrate_boundary();
	if (debug_info)
		cout<<"Generating polygons complete"<<endl;

	return 0;
}

int outer_bound::integrate_boundary ()
{
	//给vpoly中的多边形添加一下其他点，以保证不会存在多边形a的顶点落在多边形b边上的情况
	point_3 **vpoint = new point_3*[vpoly.size()];		//记录每个多边形顶点在三维平面上的坐标
	int *p_count = new int[vpoly.size ()];				//记录每个ｐｗｈ的拥有的顶点的数量
	//记录外表面边界中每个多边形的每个顶点
	for (int i(0);i!=(int)vpoly.size ();++i)
	{
		p_count[i] = 0;
		const polygon_2 &outer (vpoly[i].outer_boundary());
		p_count[i] += (int)outer.size ();
		for (polygon_with_holes_2::Hole_const_iterator hi = vpoly[i].holes_begin();	hi!=vpoly[i].holes_end();++hi)
			p_count[i] += (int)hi->size ();

		vpoint[i] = new point_3[p_count[i]];
		int p_num = 0;
		for (polygon_2::Vertex_const_iterator pi = outer.vertices_begin(); pi!=outer.vertices_end();++pi)
		{
			vpoint[i][p_num] = vplane[i].to_3d(*pi);
			++p_num;
		}
		for (polygon_with_holes_2::Hole_const_iterator hi = vpoly[i].holes_begin(); hi!=vpoly[i].holes_end();++hi)
		{
			for (polygon_2::Vertex_const_iterator pi = hi->vertices_begin(); pi!=hi->vertices_end();++pi)
			{
				vpoint[i][p_num] = vplane[i].to_3d(*pi);
				++p_num;
			}
		}
	}
	mbox_3 *vbox = new mbox_3[vpoly.size ()];		//为每个多边形都初始化一个最小包围盒
	for (int i(0);i!=(int)vpoly.size ();++i)
		vbox[i] = build_mbox_3(vpoint[i],p_count[i]);

	//下面这部分做的就是给vpoly中的每个多边形都添加一下角点
	for (int i(0);i!=(int)vpoly.size ();++i)
		add_coner_point(i,vpoint,p_count,vbox);		//给vpoly中的每个多边形都添加角点

	//释放内存
	for (int i(0);i!=(int)vpoly.size ();++i)
		delete [] vpoint[i];
	delete [] vpoint;
	delete [] p_count;
	delete [] vbox;
	return 0;
}

void sort_index (const std::vector<FT> &tar, std::vector<int> &index)
{
	index.reserve (tar.size ());
	vector<bool> used_tag(tar.size (),false);
	int count(0);
	while (count<tar.size ())
	{
		FT maxval;
		int maxi;
		for (int i(0);i!=tar.size ();++i)
		{
			if (used_tag[i]) continue;
			maxval = tar[i];		//初始化一下
			maxi = i;
			break;
		}
		for (int i(0);i!=tar.size ();++i)
		{
			if (used_tag[i]) continue;
			if (tar[i]>maxval)
			{
				maxi = i;
				maxval = tar[i];
			}
		}
		index.push_back (maxi);
		used_tag[maxi] = true;
		++count;
	}
	std::reverse(index.begin (),index.end ());
}


}