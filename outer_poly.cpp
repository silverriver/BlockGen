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
	vector<int> plane_index;			//����ÿ��ƽ����bs.vpl�е�����
	vector<polygon_set_2> p_poly;		//��ʾ����Ϊplane_index[i]��ƽ������Ķ����
	vector<polygon_set_2> n_poly;		//��ʾ����Ϊplane_index[i]��ƽ�渺��Ķ����
	
	for (int i(0);i!=vebi.size ();++i)
	{
		const eblock& eb = bs.veb[vebi[i]];
		for (int j(0);j!=eb.vf.size ();++j)
		{
			int index = add_set(plane_index,eb.vf[j].plane_id);
			if (plane_index.size () > p_poly.size())
			{
				//������plane_index��������µ�ƽ��
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
		polygon_set_2 temp_set (n_poly[i]);		//�������ֱ����symmetric_difference������㣬���������Ļ�����ⷨ�߷����ʱ����Ҫ�����жϸ�������ε��ⷨ�߷���
		n_poly[i].difference(p_poly[i]);
		p_poly[i].difference(temp_set);

		if ((!n_poly[i].is_valid()) || (!p_poly[i].is_valid()))
			throw logic_error ("construct_outer_boundary, result poly set is invaild");
		if (n_poly[i].is_empty() && p_poly[i].is_empty())
			continue;	//����ǿյľ�ֱ������
		else
		{
			//polygon_set_2 n_poly_simplified, p_poly_simplified;
			vector<polygon_with_holes_2> v_n_poly, v_p_poly;
			int vp_size = (int)vpoly.size ();
			n_poly[i].polygons_with_holes(std::back_inserter(v_n_poly));
			for (int j(0);j!=(int)v_n_poly.size ();++j)
			{
				vector<polygon_with_holes_2> vpwh;
				if (split_polygon_2(v_n_poly[j],vpwh) == 0)		//�����Ž�ÿ��������ʷ�һ�£��ʷֳɼ򵥶����
				{
					vpoly.push_back (polygon_with_holes_2());
					simplify_polygon_2(v_n_poly[j],vpoly.back ());		//���ǰһ��Ҫ�ǵ�Ҫ��һ�±߽�
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
			//��������ǰ�İ汾����֮ǰ�İ汾����û��split���������
			//simplify_polygon_2(p_poly[i],p_poly_simplified);		//��һ�����ɵĶ����
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
	//Ѱ��λ��ƽ��vplane[poly_n]�ϵĵ㣬Ϊ�˷���ֱ�Ӽ�¼��Щ������꣬����¼����
	vector<point_3> v_coner_p;
	v_coner_p.reserve(vpoly.size ()/3);
	for (int i(0);i!=(int)vpoly.size ();++i)
	{
		if (!CGAL::do_overlap(vbox[i],vbox[poly_n]))		//�����������ζ����ཻ����������ע�������ʱ�������������ͬһ����Ҳ��Ҫ���,��Ϊ���ܴ������relative simple�Ķ����
			continue;
		for (int j(0);j!=p_count[i];++j)
			if (vplane[poly_n].has_on(vpoint[i][j]))
				add_set(v_coner_p,vpoint[i][j]);		//��¼����������
	}

	polygon_2 poly_temp;
	int point_added (0);		//��¼�ڶ���α߽�������˼�����
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
	vector<point_3> poly_3d;		//ת������ά�ռ�
	poly_3d.reserve (original.size ());
	for (int i(0);i!=original.size ();++i)
		poly_3d.push_back (plane.to_3d(original[i]));

	for (int i(0); i!=original.size ();++i)
	{
		new_poly.push_back (original[i]);
		int i1 (i+1);
		if (i1>=(int)original.size ())
			i1 = 0;
		vector<int> vp;		//��¼���߶�i,i1�ϵ���Ҫ��ӵĵ���vconer�ϵ�����
		vp.reserve (vconer.size ());
		segment_3 seg(poly_3d[i], poly_3d[i1]);
		for (int j(0);j!=(int)vconer.size ();++j)
			if (seg.has_on(vconer[j]) && vconer[j]!=seg.start() && vconer[j]!=seg.end ())	//���������߶��ϣ����Ҳ��ڶ˵�
				vp.push_back (j);
		
		if (vp.empty ()) continue;
		point_added += (int)vp.size ();
		vector<FT> dis;		//����vp�и�������poly_3d[i]�ľ���
		dis.reserve (vp.size ());
		for (int j(0);j!=(int)vp.size ();++j)
			dis.push_back (CGAL::squared_distance(poly_3d[i],vconer[vp[j]]));
		vector<int> sorted_index;
		sort_index(dis,sorted_index);	//����Ҫ��ӵĵ���һ���򣬰�����poly_3d[i]�ӽ���Զ����
		for (int j(0);j!=sorted_index.size ();++j)		//������Ӹ�����
			new_poly.push_back (plane.to_2d(vconer[vp[sorted_index[j]]]));
	}
	return point_added;		//�������
}

int outer_bound::join_veb (const block_system &bs, const std::vector <int> &vebi, const bool& debug_info)
{
	//����fuck��������Ҫʵ��һ�����Ӽ����㷨�����Ǻ���ʵ�ֵ�ʱ������Ҫʵ������㷨����鷳��
	//���ұȽ����׳�bug������Ϊ��ʡ�£�����ֱ������Ա�һ��ķ������Ƚϰ�ȫʡ��һ��
	return join(outer_bound(bs,vebi),debug_info);
}

int outer_bound::join_cb (const block_system &bs, const int cbi, const bool& debug_info)
{
	return join_veb(bs,bs.vcb[cbi].vb,debug_info);
}

int outer_bound::join (const outer_bound &ob, const bool &debug_info)
{
	//���뺯���е�ob��*this���������ڲ�ͬ�����ϵͳ���������������vplane�е�ƽ����ܻ���ڷ�����������Ҫע��
	vector<vector<int>> common_pl;	//�ҳ�this��ob��vplane��һ���ж��ٲ��غϵ�ƽ�档common_pl[i]�е�ƽ�涼�����غϣ���Щƽ����������ʾ������ob�е�ƽ��������ԭ��ֵ�Ļ����ϼ���this->vplane.size().���ob�е�ĳ��ƽ���this�е�ĳ��ƽ���෴�������������Ϊ����
	common_pl.reserve (vplane.size () + ob.vplane.size ());		//common_pl[i][0]���µ�ƽ��vector
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
			if (common_pl[count][0]<vplane.size ())		//ע��common_pl[i][0]������С��0����this->vpolyΪ�յ�ʱ��ſ��ܵ���0
				p_temp = &vplane[common_pl[count][0]];	//���common_pl[count]����this�е�ƽ��
			else 
				p_temp = &ob.vplane [common_pl[count][0] - vplane.size ()];		//���common_pl��û��this�е�ƽ��
				
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
			common_pl.back ().push_back (i+(int)vplane.size ());		//���this->vplane��û�����ƽ��Ļ������һ�¡�
		}
	}

	vector<polygon_set_2> n_poly(common_pl.size ());		//��ʾ��common_pl[i]����ʾƽ�渺��Ķ����
	vector<polygon_set_2> p_poly(common_pl.size ());		//��ʾ��common_pl[i]����ʾƽ������Ķ���� ���ⷨ�߷������ƽ���෴��
	for (int i(0);i!=(int)common_pl.size ();++i)
	{
		for (int j(0);j!=(int)common_pl[i].size ();++j)
		{
			if (abs(common_pl[i][j]) < vplane.size ())			//��ƽ����this��
				if (outer_normal[common_pl[i][j]])
					n_poly[i].join(vpoly[common_pl[i][j]]);		//���ⷨ�߷�����common_pl[i]����ʾƽ��ͬ��Ķ������ӵ�n_poly��
				else 
					p_poly[i].join(vpoly[common_pl[i][j]]);
			else												//��ƽ����ob�С�
				if (common_pl[i][j] >= 0)						//��������±�ʾ���ƽ���common_pl[i][0]��ͬ��ֻ�е�this->vpoly.size() == 0 ʱ���п��ܵ���0
					if (ob.outer_normal[common_pl[i][j]-vplane.size ()])
						n_poly[i].join (ob.vpoly[common_pl[i][j]-vplane.size ()]);
					else
						p_poly[i].join (ob.vpoly[common_pl[i][j]-vplane.size ()]);
				else											//��������±�ʾ���ƽ���common_pl[i][0]�෴�����������common_pl[i][0]ֻ����this�е�ƽ�棬��������ob�е�ƽ��
					if (ob.outer_normal[-common_pl[i][j]-vplane.size ()])
						p_poly[i].join (change_proj_plane(ob.vplane[-common_pl[i][j]-vplane.size ()],vplane[common_pl[i][0]], ob.vpoly[-common_pl[i][j]-vplane.size ()]));		//��֮�������ָ����ob����Ӧ��ƽ���this�з����෴�����ʱ�򷴹������.ע�⣬���ʱ��Ӧ����������һ��Ŀ�����Σ������2015-10-28 11:59
					else
						n_poly[i].join (change_proj_plane(ob.vplane[-common_pl[i][j]-vplane.size ()],vplane[common_pl[i][0]], ob.vpoly[-common_pl[i][j]-vplane.size ()]));
		}
	}

	if (debug_info)
		cout<<"Generating polygon"<<endl;
	vector<plane_3> new_vplane;		//�µ������߽��Ӧ��ƽ��
	new_vplane.reserve (vplane.size ()+ob.vplane.size ());
	vector<polygon_with_holes_2> new_vpoly;		//�µ��������������
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
				if (split_polygon_2(v_n_poly[j],vpwh) == 0)		//�����Ž�ÿ��������ʷ�һ�£��ʷֳɼ򵥶����
				{
					new_vpoly.push_back (polygon_with_holes_2());
					simplify_polygon_2(v_n_poly[j],vpoly.back ());		//���ǰһ��Ҫ�ǵ�Ҫ��һ�±߽�
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
	//��vpoly�еĶ�������һ�������㣬�Ա�֤������ڶ����a�Ķ������ڶ����b���ϵ����
	point_3 **vpoint = new point_3*[vpoly.size()];		//��¼ÿ������ζ�������άƽ���ϵ�����
	int *p_count = new int[vpoly.size ()];				//��¼ÿ��������ӵ�еĶ��������
	//��¼�����߽���ÿ������ε�ÿ������
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
	mbox_3 *vbox = new mbox_3[vpoly.size ()];		//Ϊÿ������ζ���ʼ��һ����С��Χ��
	for (int i(0);i!=(int)vpoly.size ();++i)
		vbox[i] = build_mbox_3(vpoint[i],p_count[i]);

	//�����ⲿ�����ľ��Ǹ�vpoly�е�ÿ������ζ����һ�½ǵ�
	for (int i(0);i!=(int)vpoly.size ();++i)
		add_coner_point(i,vpoint,p_count,vbox);		//��vpoly�е�ÿ������ζ���ӽǵ�

	//�ͷ��ڴ�
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
			maxval = tar[i];		//��ʼ��һ��
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