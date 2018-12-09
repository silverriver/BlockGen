#include "entrance_block.h"
#include "geometry.h"
#include <vector>
#include "outer_poly.h"
#include <fstream>
#include <iostream>
using namespace std;

namespace BI
{

//�ж�һ�����һ����Ԫ�����λ�ù�ϵ,vfi��һ������ֵ����ʾ���������ڸõ�Ԫ�����ϣ����¼�����λ���ļ���ƽ����
Oriented_side point_eblock_location (const block_system &bs, const eblock &eb, const point_3&tar, vector<int>& vfi)
{
	int n_on(0), n_in(0), n_out(0);
	vfi.clear();
	for (int i(0);i!=eb.vf.size ();++i)
	{
		switch (bs.vpl[eb.vf[i].plane_id].oriented_side(tar))
		{
		case CGAL::ON_ORIENTED_BOUNDARY:
			++n_on; vfi.push_back (i);
			break;
		case CGAL::ON_POSITIVE_SIDE:
			if (eb.vf[i].outer_normal)
				++n_out;
			else
				++n_in;
			break;
		case CGAL::ON_NEGATIVE_SIDE:
			if (eb.vf[i].outer_normal )
				++n_in;
			else
				++n_out;
		}
	}
	if (n_out>0)
		return Oriented_side::OUT_;
	if (n_on > 0)
		return Oriented_side::ON_;
	else
		return Oriented_side::IN_;
}

bool angle_edge::sub_angle::is_internal_intersected (const vector_3 &v) const
{
	sub_angle sa_temp;
	sa_temp.vnormal.push_back (v);
	return is_internal_intersected(sa_temp);
}

bool angle_edge::sub_angle::is_internal_intersected (const sub_angle &sa) const
{
	switch (intersection(sa,point_3(0,0,0), point_3(0,0,0)))
	{
	case Intersection_Status::Inter_Boundary:
		return false;
	case Intersection_Status::Inter_Unbounded:
		return true;
	case Intersection_Status::Inter_Internal:
		throw logic_error ("bool angle_edge::sub_angle::intersection (const vector_3& v), intersection have boundary");
	case Intersection_Status::Separate:
		throw logic_error ("bool angle_edge::sub_angle::intersection (const vector_3& v), do not intersection");
	}
}

Intersection_Status angle_edge::sub_angle::intersection (const sub_angle&sa, const point_3 &p_this, const point_3 &p_sa) const
{
	vector<plane_3> vpl;
	vpl.reserve (vnormal.size ()+sa.vnormal.size ());
	for (int i(0);i!=(int)vnormal.size ();++i)
		vpl.push_back (plane_3(p_this, vnormal[i]));
	for (int i(0);i!=(int)sa.vnormal.size ();++i)
		vpl.push_back (plane_3(p_sa, sa.vnormal[i]));
	return BI::intersection(vpl);
}

bool angle_edge::sub_edge::is_internal_intersected (const sub_edge& se) const
{
	switch (intersection(se, point_3(0,0,0), point_3(0,0,0)))
	{
	case Intersection_Status::Inter_Boundary:
		return false;
	case Intersection_Status::Inter_Unbounded:
		return true;
	case Intersection_Status::Inter_Internal:
		throw logic_error ("bool angle_edge::sub_edge::intersection (const sub_edge& se), intersection have boundary");
	case Intersection_Status::Separate:
		throw logic_error ("bool angle_edge::sub_edge::intersection (const sub_edge& se), halp space do not intersect");
	}
}

Intersection_Status angle_edge::sub_edge::intersection (const sub_edge&se, const point_3 &p_this, const point_3 &p_se) const
{
	vector<plane_3> vpl;
	vpl.reserve (4);
	vpl.push_back (plane_3(p_this,n1));
	vpl.push_back (plane_3(p_this,n2));
	vpl.push_back (plane_3(p_se,se.n1));
	vpl.push_back (plane_3(p_se,se.n2));
	return BI::intersection (vpl);
}

bool angle_edge::angle::is_internal_intersected (const angle &a) const
{
	switch (intersection(a,point_3(0,0,0),point_3(0,0,0)))
	{
	case Intersection_Status::Inter_Boundary:
		return false;
	case Intersection_Status::Inter_Unbounded:
		return true;
	case Intersection_Status::Inter_Internal:
		throw logic_error ("bool angle_edge::angle::is_internal_intersected (const angle &a), intersection have boundary");
	case Intersection_Status::Separate:
		throw logic_error ("bool angle_edge::angle::is_internal_intersected (const angle &a), do not intersection");
	default:
		throw logic_error ("bool angle_edge::angle::is_internal_intersected (const angle &a), return type error");
	}
}

bool angle_edge::angle::is_internal_intersected (const vector_3 &v) const
{
	angle a_temp;
	a_temp.vsa.push_back (angle_edge::sub_angle());
	a_temp.vsa.back ().vnormal.push_back (v);
	return is_internal_intersected(a_temp);
}

Intersection_Status angle_edge::angle::intersection (const angle &a, const point_3&p_this, const point_3&p_a) const
{
	int n_sep(0), n_internal(0), n_bound(0);
	for (int i(0);i!= (int) vsa.size ();++i)
	{
		for (int j(0);j!= (int) a.vsa.size ();++j)
		{
			switch (vsa[i].intersection(a.vsa[j], p_this, p_a))
			{
			case Intersection_Status::Inter_Boundary:
				n_bound++; break;
			case Intersection_Status::Inter_Internal:
				n_internal++; break;
			case Intersection_Status::Inter_Unbounded:
				return Intersection_Status::Inter_Unbounded;
			case Intersection_Status::Separate:
				n_sep++; break;
			}
		}
	}
	if (n_internal > 0)
		return Intersection_Status::Inter_Internal;
	if (n_bound > 0)
		return Intersection_Status::Inter_Boundary;
	else 
		return Intersection_Status::Separate;
}

bool angle_edge::edge::is_internal_intersected (const edge &e) const
{
	switch (intersection(e,point_3(0,0,0),point_3(0,0,0)))
	{
	case Intersection_Status::Inter_Boundary:
		return false;
	case Intersection_Status::Inter_Unbounded:
		return true;
	case Intersection_Status::Inter_Internal:
		throw logic_error ("bool angle_edge::edge::is_internal_intersected (const edge &e), intersection have boundary");
	case Intersection_Status::Separate:
		throw logic_error ("bool angle_edge::edge::is_internal_intersected (const edge &e), do not intersection");
	default:
		throw logic_error ("bool angle_edge::edge::is_internal_intersected (const edge &e), return type error");
	}
}

Intersection_Status angle_edge::edge::intersection (const edge& e, const point_3& p_this, const point_3 &p_e) const
{
	int n_sep(0), n_internal(0), n_bound(0);
	for (int i(0);i!= (int)vse.size ();++i)
	{
		for (int j(0);j!= (int)e.vse.size ();++j)
		{
			switch (vse[i].intersection(e.vse[j],p_this,p_e))
			{
			case Intersection_Status::Inter_Boundary:
				n_bound++; break;
			case Intersection_Status::Inter_Internal:
				n_internal++; break;
			case Intersection_Status::Inter_Unbounded:
				return Intersection_Status::Inter_Unbounded;
			case Intersection_Status::Separate:
				n_sep++; break;
			}
		}
	}
	if (n_internal >0)
		return Intersection_Status::Inter_Internal;
	if (n_bound>0)
		return Intersection_Status::Inter_Boundary;
	else
		return Intersection_Status::Separate;
}

int entrance_block::gen_cover (const outer_bound &ob_a, const angle_edge &ae_a, const outer_bound &ob_b, const angle_edge &ae_b, const point_3 &a0)
{
	vcover.clear();
	vcover.reserve (ob_b.vpoly.size ());
	//��������������ɽ�����������渲��
	//E(A(0), B(2))
	for (int i(0);i!=(int)ob_b.vpoly.size ();++i)
	{
		vector_3 outer_n (ob_b.vplane[i].orthogonal_vector());	//����B�������ⷨ�߷���
		if (!ob_b.outer_normal[i])
			outer_n = -outer_n;
		for (int j(0);j!=(int)ae_a.va.size ();++j)		//���ڿ���A��ÿ����ά��
		{
			if (!ae_a.va[j].is_internal_intersected(outer_n))	//����A�������ά�ǺͿ���B������Ӧ�İ�ռ�û�н���
			{
				//�����ά���������ռ䲻�ཻ,��ʱE(ae_a.va[j].e, ob_b.vpoly[i])�п�����һ������
				vcover.push_back (poly_frac());
				const point_3 & a = ae_a.vp[ae_a.va[j].e];
				construct_poly_frac(vcover.back (),ob_b.vpoly[i],ob_b.vplane[i],vector_3(a, a0),false);
			}
		}
	}

	//E(A(2), B(0))
	for (int i(0);i!=(int)ob_a.vpoly.size ();++i)
	{
		vector_3 outer_n (ob_a.vplane[i].orthogonal_vector());		//�ҵ�����A����߽��Ӧ������淨��
		if (!ob_a.outer_normal[i])
			outer_n = -outer_n;
		for (int j(0);j!=(int)ae_b.va.size ();++j)		//���ڿ���B��ÿ����ά��
		{
			if (!ae_b.va[j].is_internal_intersected(outer_n))	//����B�������ά�ǺͿ���A������Ӧ�İ�ռ�û�н���
			{
				//�����ά�Ǻ������ռ䲻���ڲ��ཻ����ʱE(ob_a.vpoly[i], ae_b.va[j].e)�п�����һ������
				vcover.push_back (poly_frac());
				const point_3 & b = ae_b.vp[ae_b.va[j].e];
				construct_poly_frac(vcover.back (),ob_a.vpoly[i],ob_a.vplane[i],vector_3(b.x()+a0.x(), b.y()+a0.y(), b.z()+a0.z()),true);
			}
		}
	}

	//E(A(1), B(1))
	for (int i(0);i!=(int)ae_a.ve.size ();++i)
	{
		for (int j(0);j!=(int)ae_b.ve.size ();++j)
		{
			if(!ae_a.ve[i].is_internal_intersected (ae_b.ve[j]))	//��������û�н���
			{
				const point_3 &a1 = ae_a.vp[ae_a.ve[i].e1];
				const point_3 &a2 = ae_a.vp[ae_a.ve[i].e2];
				const point_3 &b1 = ae_b.vp[ae_b.ve[j].e1];
				const point_3 &b2 = ae_b.vp[ae_b.ve[j].e2];
				vcover.push_back (poly_frac());
				vcover.back ().vp_ob.push_back (point_3 (b1.x()-a1.x()+a0.x(), b1.y()-a1.y()+a0.y(), b1.z()-a1.z()+a0.z()));
				vcover.back ().vp_ob.push_back (point_3 (b2.x()-a1.x()+a0.x(), b2.y()-a1.y()+a0.y(), b2.z()-a1.z()+a0.z()));
				vcover.back ().vp_ob.push_back (point_3 (b2.x()-a2.x()+a0.x(), b2.y()-a2.y()+a0.y(), b2.z()-a2.z()+a0.z()));
				vcover.back ().vp_ob.push_back (point_3 (b1.x()-a2.x()+a0.x(), b1.y()-a2.y()+a0.y(), b1.z()-a2.z()+a0.z()));
			}
		}
	}
	return 0;
}

void update_mbbox (const poly_frac::v_loop &loop,FT &minx, FT &maxx, FT &miny, FT &maxy, FT&minz, FT&maxz)
{
	for (int i(0);i!=loop.size ();++i)
	{
		if (loop[i].x() > maxx) maxx = loop[i].x();
		if (loop[i].x() < minx) minx = loop[i].x();
		if (loop[i].y() > maxy) maxy = loop[i].y();
		if (loop[i].y() < miny) miny = loop[i].y();
		if (loop[i].z() > maxz) maxz = loop[i].z();
		if (loop[i].z() < minz) minz = loop[i].z();
	}
}

int entrance_block::gen_entrance_block (const block_system &bs_, const cblock& cba, const cblock& cbb, const point_3 &a0,
							const outer_bound &ob_a, const angle_edge &ae_a, const outer_bound &ob_b, const angle_edge &ae_b, bool debug_info)
{
	//debug_info
	if (debug_info)
		cout<<"Generating Covers..."<<endl;

	gen_cover(ob_a,ae_a,ob_b,ae_b,a0);
	if (vcover.empty())
		return 1;		//�����쳣

	//��ȷ�����������ķ�Χ��
	FT minx,maxx,miny,maxy,minz,maxz;
	minx = maxx = vcover[0].vp_ob[0].x();
	miny = maxy = vcover[0].vp_ob[0].y();
	minz = maxz = vcover[0].vp_ob[0].z();
	for (int i(0);i!=vcover.size ();++i)
	{
		update_mbbox(vcover[i].vp_ob,minx,maxx,miny,maxy,minz,maxz);
		for (int j(0);j!=vcover[i].vp_ib.size ();++j)
			update_mbbox (vcover[i].vp_ib[j],minx,maxx,miny,maxy,minz,maxz);
	}

	//����һ����Χ�ķ����о�����
	if (minx == maxx || miny == maxy || minz ==maxz) 
		return 2;		//�����쳣
	bs.add_rect_domain(minx-1,miny-1,minz-1,maxx+1,maxy+1,maxz+1);

	//�Ѹ��Ǽӵ�bs�У�ע�⣬һ��Ҫ�ȳ�ʼ���о������ٳ�ʼ����϶��
	for (int i(0);i!= (int)vcover.size ();++i)
		bs.add_poly_frac(vcover[i]);
	
	//debug_info
	if (debug_info)
		cout<<vcover.size ()<<" covers generated. "<<"Identifying entrance blocks"<<endl;

	//ʶ��֮
	bs.identify_block();

	//debug_info
	if (debug_info)
		cout<<"Determining entrance blocks"<<endl;

	//Ȼ���ҳ��������
	vector<int> vcb_in;		//��¼���ɵ�cblock����Щ�ڽ��������
	vector<int> vcb_out;
	vcb_in.reserve (bs.vcb.size ()); vcb_out.reserve (bs.vcb.size ());
	
	for (int i(0);i!=(int) bs.vcb.size ();++i)		//����ʶ�����ɵ�ÿ�����Ͽ���
	{
		point_3 p_temp (bs.vcb[i].get_rand_p(bs.veb,bs.vpo));		//�ж�������Ͽ����Ƿ���E(A,B)�С�
		switch (intersection(bs_.vpl,bs_.veb,cba,cbb,vector_3(a0,p_temp)))
		{
		case Intersection_Status::Separate:
			vcb_out.push_back (i); break;
		default:
			vcb_in.push_back (i);
			break;
		}
	}

	inner.vb.reserve (vcb_in.size ()*8);
	for (int i(0);i!= (int)vcb_in.size ();++i)
		for (int j(0);j!=(int)bs.vcb[vcb_in[i]].vb.size ();++j)
			inner.vb.push_back (bs.vcb[vcb_in[i]].vb[j]);

	outer_b.construct_outer_bound(bs,inner.vb);

	return 0;
}

int entrance_block::gen_entrance_block (const block_system &bs_, const cblock& cba, const cblock& cbb, const point_3 &a0, bool debug_info)
{
	//���ɸ���
	outer_bound ob_a, ob_b;
	ob_a.construct_outer_bound (bs_, cba.vb);
	ob_b.construct_outer_bound (bs_, cbb.vb);

	angle_edge ae_a, ae_b;
	outer_bound2angle_edge(bs_,cba,ob_a,ae_a);
	outer_bound2angle_edge(bs_,cbb,ob_b,ae_b);

	return gen_entrance_block(bs_,cba,cbb,a0,ob_a,ae_a,ob_b,ae_b,debug_info);
}

//�ж�p1p2����߶��Ƿ���vf��vp����ʾ�Ķ�����ڣ��������Ƿ���һ���㣩��Ҫ��ö������͹����Σ����е㶼��pl�ϣ�����p1�϶������������ϣ��ڻ�߽��ϣ�
bool segment_in_polygon (const vector<point_3> &vp, const vector<int> &vf, const plane_3 &pl, const point_3& e1, const point_3& e2)
{
	int ie1(-1), ie2(-1);		//��¼e1��͹����ε��ļ�������
	vector<plane_3> vpl;
	vector<CGAL::Oriented_side> vflag;
	vpl.reserve (vf.size ());
	vflag.reserve (vf.size ());
	//�ж�e1���������ε����λ��
	for (int i(0);i!=vf.size ();++i)
	{
		int i1(i+1), i2(i+2);
		if (i1 >= (int)vf.size ())
			i1 = i1 - (int)vf.size ();
		if (i2 >= (int)vf.size ())
			i2 = i2 - (int)vf.size ();
		//�ж�e1�ڱ�i,i1���ڲ࣬�ϣ��������
		vpl.push_back (plane_3 (vp[vf[i]], CGAL::cross_product(pl.orthogonal_vector(),vector_3(vp[vf[i]], vp[vf[i1]]))));
		vflag.push_back (vpl.back ().oriented_side(vp[vf[i2]]));
		switch (vpl.back ().oriented_side(e1))
		{
		case CGAL::Oriented_side::ON_ORIENTED_BOUNDARY:
			if (ie1<0) {ie1 = i; break;}
			if (ie2<0) {ie2 = i; break;}
			throw logic_error("segment_in_polygon, e1 on three edges");
		case CGAL::Oriented_side::ON_NEGATIVE_SIDE:
			if (vflag.back () == CGAL::Oriented_side::ON_POSITIVE_SIDE)
				throw logic_error ("segment_in_polygon, e1 is not in the polygon");
			break;
		case CGAL::Oriented_side::ON_POSITIVE_SIDE:
			if (vflag.back () == CGAL::Oriented_side::ON_NEGATIVE_SIDE)
				throw logic_error ("segment_in_polygon, e1 is not in the polygon");
			break;
		}
		if (vflag.back () == CGAL::Oriented_side::ON_ORIENTED_BOUNDARY)
			throw logic_error ("segment_in_polygon, illegal convex polygon");
	}
	if (ie1 < 0 && ie2 < 0)
		return true;		//e1�ڶ�����ڲ��������߶������εĽ���϶�����һ����
	else if (ie2 < 0)
	{
		//e1��ie1����ʾ���������ϣ��Ƕ��㴦��
		CGAL::Oriented_side tag_temp (vpl[ie1].oriented_side(e2));
		if (tag_temp == CGAL::Oriented_side::ON_ORIENTED_BOUNDARY)
			return true;		//e2Ҳ����������
		if (tag_temp == vflag[ie1])
			return true;
		else
			return false;
	}
	else
	{
		//e1ͬʱλ��ie1��ie2����������
		CGAL::Oriented_side tag_temp1 (vpl[ie1].oriented_side(e2));
		CGAL::Oriented_side tag_temp2 (vpl[ie2].oriented_side(e2));
		
		if (tag_temp1 == CGAL::Oriented_side::ON_ORIENTED_BOUNDARY && tag_temp2 == CGAL::Oriented_side::ON_ORIENTED_BOUNDARY)
			throw logic_error ("segment_in_polygon, e1 e2 are the same point");
		else if (tag_temp1 == CGAL::Oriented_side::ON_ORIENTED_BOUNDARY)
			if (tag_temp2 == vflag[ie2]) return true;
			else return false;
		else if (tag_temp2 == CGAL::Oriented_side::ON_ORIENTED_BOUNDARY)
			if (tag_temp1 == vflag[ie1]) return true;
			else return false;
		else
			if (tag_temp1 == vflag[ie1] && tag_temp2 == vflag[ie2]) return true;
			else return false;
	}
}

int outer_bound2angle_edge (const block_system & bs, const cblock & cb,const outer_bound &ob, angle_edge &ae)
{
	const vector<polygon_with_holes_2> &vpoly = ob.vpoly;	//���������Ǳ����ɴ�����������һ������
	const vector<plane_3> &vplane = ob.vplane;
	
	ae.va.clear(); ae.ve.clear (); ae.vp.clear ();
	typedef vector<int> loop;			//���������Ļ�·
	typedef vector<loop> polygon_with_holes_3;		//��һ��loop���ⷨ�߷���

	vector<polygon_with_holes_3> vpwh3;
	vpwh3.reserve (vpoly.size ());

	//�Ƚ�ÿ���㶼ת������ά�ռ���
	for (int i(0);i!=(int)vpoly.size ();++i)
	{
		vpwh3.push_back (polygon_with_holes_3());
		vpwh3.back ().reserve (vpoly[i].number_of_holes()+1);
		vpwh3.back ().push_back (loop());
		vpwh3.back ().back ().reserve (vpoly[i].outer_boundary ().size ());

		const polygon_2 &outer(vpoly[i].outer_boundary ());
		for (int j(0);j!=(int)outer.size ();++j)
			vpwh3.back ().back ().push_back (add_set(ae.vp,vplane[i].to_3d(outer[j])));

		for (polygon_with_holes_2::Hole_const_iterator hi = vpoly[i].holes_begin(); hi!=vpoly[i].holes_end(); ++hi)
		{
			vpwh3.back ().push_back (loop());
			vpwh3.back ().back ().reserve (hi->size ());
			for (int j(0);j!=(int)hi->size ();++j)
				vpwh3.back ().back ().push_back (add_set(ae.vp,vplane[i].to_3d((*hi)[j])));
		}
	}

	//����ľ������ɸ�����ά��
	ae.va.reserve (ae.vp.size ());		//��������һ�������Ӧһ����
	typedef struct{int ebi; vector<int> vfi;} v_eb_info;		//����һ���ṹ����¼ĳ���������ĸ���Ԫ����(ebi)����Щ��(vfi)��
	vector<vector<v_eb_info>> v_angle_suba_info;		//��¼ÿ����ά�ǵ�ÿ���ӽǵĹ�����Ϣ.��Щ��Ϣ����ֱ�Ӱ�����ÿ����ά���У���������������ᵼ��angle_edge�е����ݺ�bs�е������޷�����������������ʱ���ǽ��������������ݷֿ���������֮�����ϡ����ܺ����ᷢ��������Ŭ����ͽ�͵ģ����������Ժ���Ҫ���ǵ������ˡ�2016-2-19 20:35
	v_angle_suba_info.reserve (ae.va.size ());

	for (int i(0);i!=ae.vp.size ();++i)		//����ÿ������
	{
		ae.va.push_back (angle_edge::angle());
		ae.va.back ().e = i;
		ae.va.back ().vsa.reserve (4);	//���ҹ���ÿ����ά�ǰ���4���ӽ�
		v_angle_suba_info.push_back (vector<v_eb_info>());
		v_angle_suba_info.back ().reserve (4);

		for (int j(0);j!=(int)cb.vb.size ();++j)		//�ж���������ÿ����Ԫ�����λ�ù�ϵ
		{
			v_eb_info info_temp;
			info_temp.vfi.reserve (bs.veb[cb.vb[j]].vf.size ());
			info_temp.ebi = cb.vb[j];
			switch (point_eblock_location (bs,bs.veb[cb.vb[j]],ae.vp[i],info_temp.vfi))
			{
			case Oriented_side::IN_:
				throw logic_error ("outer_bound2angle_edge, point in eblock");		//������㲻�����ڵ�Ԫ������
				break;
			case Oriented_side::ON_:
				ae.va.back ().vsa.push_back (angle_edge::sub_angle());			//��������ڵ�Ԫ�����ϣ���vfi�ж�Ӧ���������ɸö����һ���ӽ�
				ae.va.back ().vsa.back ().vnormal.reserve (info_temp.vfi.size ());
				for (int k(0);k!=(int)info_temp.vfi.size ();++k)		//����ӽ�
				{
					if (bs.veb[cb.vb[j]].vf[info_temp.vfi[k]].outer_normal )
						ae.va.back ().vsa.back ().vnormal.push_back (bs.vpl[bs.veb[cb.vb[j]].vf [info_temp.vfi[k]].plane_id].orthogonal_vector());
					else
						ae.va.back ().vsa.back ().vnormal.push_back (-bs.vpl[bs.veb[cb.vb[j]].vf [info_temp.vfi[k]].plane_id].orthogonal_vector());
				}
				v_angle_suba_info.back ().push_back (info_temp);		//��¼����ӽǵ�������Ϣ
			default:
				break;
			}
		}
	}

	//����������ɸ�����
	ae.ve.reserve (ae.va.size () + ob.vpoly.size () - 2);	//�ȼ�����v+f-2����
	for (int i(0);i!=(int) vpwh3.size ();++i)		//���ҳ����еı�
	{
		for (int j(0);j!=(int) vpwh3[i].size ();++j)
		{
			for (int k(0);k!=(int) vpwh3[i][j].size ();++k)
			{
				int k1(k+1);
				if (k1 >= (int)vpwh3[i][j].size ())
					k1 = 0;
				add_set(ae.ve,angle_edge::edge (vpwh3[i][j][k], vpwh3[i][j][k1]));
			}
		}
	}
	for (int i(0);i!=(int) ae.ve.size ();++i)		//����ÿ����ά��
	{
		vector_3 edge_vector (ae.vp[ae.ve[i].e1], ae.vp[ae.ve[i].e2]);		//�������������ά�߶�Ӧ������
		ae.ve[i].vse.reserve (3);		//����ÿ�����������ӱߣ���Ԥ��λ��
		vector<v_eb_info> & angle_suba_info = v_angle_suba_info[ae.ve[i].e1];	//�ҵ�ÿ����ά�ߵĵ�һ�������������Ϣ
																				//ע�⣬����ֻ�����˵�һ�������������Ϣ����Ϊ�����������ɱߵ���������νṹob�ǹ淶��֮��ġ����ԣ�����õڶ����������Ϣ���ɸ����ӱߣ���ô�����ɵĶ�ά����һ���ģ�������������ά�ߵ��ӱ߿��ܲ�����ͬ
		for(int j(0);j!= (int)angle_suba_info.size ();++j)		//��������������ڵ�ÿ����Ԫ����
		{
			vector<int> vfi;	//��¼edge_vector�������Ԫ�������Щ����
			vfi.reserve (2);	//Ԥ��2������Ϊ�����ϲ����ܳ���3����
			int ebi = angle_suba_info[j].ebi;		//��ʱ������ָʾ�����ڴ����ĸ���Ԫ����

			for (int k(0);k!=(int) angle_suba_info[j].vfi.size ();++k)
			{
				int vfi_temp = angle_suba_info[j].vfi[k];		//��ʱ������ָʾ�����ڴ���ebi���ĸ���
				int plane_id_temp = bs.veb[ebi].vf[vfi_temp].plane_id;
				if ((bs.vpl[plane_id_temp].orthogonal_vector() * edge_vector == 0) &&
					segment_in_polygon(bs.vpo,bs.veb[ebi].vf[vfi_temp].np,bs.vpl[plane_id_temp],ae.vp[ae.ve[i].e1], ae.vp[ae.ve[i].e2]))
				{
					vfi.push_back (vfi_temp);//��ƽ����������ƽ�У�����������Ҳ�͸�ƽ���ϵ��������ཻ��
				}
				
			}
			if (vfi.size () >= 3)
				throw logic_error("outer_bound2angle_edge, edge in more than three planes");
			if (vfi.size () == 1)
			{
				ae.ve[i].vse.push_back (angle_edge::sub_edge());		//���������������Ԫ�����ĳ������
				if (bs.veb[ebi].vf[vfi[0]].outer_normal)		//sub_edge��n1��n2ָ���ӱߵ��ⲿ��ע�⣬���ⲿ����
					ae.ve[i].vse.back ().n1 = ae.ve[i].vse.back ().n2 = bs.vpl[bs.veb[ebi].vf[vfi[0]].plane_id].orthogonal_vector();
				else
					ae.ve[i].vse.back ().n1 = ae.ve[i].vse.back ().n2 = -bs.vpl[bs.veb[ebi].vf[vfi[0]].plane_id].orthogonal_vector();
			}
			if (vfi.size () == 2)
			{
				ae.ve[i].vse.push_back (angle_edge::sub_edge());		//���������������Ԫ�����ĳ������
				if (bs.veb[ebi].vf[vfi[0]].outer_normal )		//sub_edge��n1��n2ָ���ӱߵ��ⲿ��ע�⣬���ⲿ����
					ae.ve[i].vse.back ().n1 = bs.vpl[bs.veb[ebi].vf[vfi[0]].plane_id].orthogonal_vector();
				else
					ae.ve[i].vse.back ().n1 = -bs.vpl[bs.veb[ebi].vf[vfi[0]].plane_id].orthogonal_vector();

				if (bs.veb[ebi].vf[vfi[1]].outer_normal )
					ae.ve[i].vse.back ().n2 = bs.vpl[bs.veb[ebi].vf[vfi[1]].plane_id].orthogonal_vector();
				else
					ae.ve[i].vse.back ().n2 = -bs.vpl[bs.veb[ebi].vf[vfi[1]].plane_id].orthogonal_vector();
			}
		}
		if (ae.ve[i].vse.empty ())
			throw logic_error ("outer_bound2angle_edge, vsubedge is empty");	//��ά�ߵ��ӱ߼��ϲ���Ϊ��
	}
	return 0;
}


//int outer_bound2angle_edge (const outer_bound &ob, angle_edge &ae)
//{

//

//	
//	
//	typedef pair<int,int> angle_2d;		//��һ�����Ͷ�����ʾһ����ά��
//	vector<vector<angle_2d>> vbound(ae.vp.size ());	//ÿ����ά�ǵı߽�
//
//	for (int i(0);i!=(int)vpwh3.size ();++i)
//	{
//		loop &outer = vpwh3[i][0];		//Ϊouter_boundary��һ������
//		vector_3 outer_normal (vplane[i].orthogonal_vector());
//		if (!ob.outer_normal[i]) outer_normal = -outer_normal;		//ȷ�����ƽ����ⷨ��
//		vector_3 oriented_n (oriented_normal(ae.vp,outer));			//����������������ȷ��outer�ϵĵ����ת����
//		bool flag (false);		//T��ʾouter�еĵ���ת������ⷨ�߷���һ�£�F��ʾouter�еĵ���ת������ⷨ�߷����෴
//		if (outer_normal*oriented_n > 0)
//			flag = true;		
//
//		for (int j(0);j!=outer.size ();++j)
//		{
//			vector<angle_2d> & bound = vbound[outer[j]];		//����outer�е�ÿ�����㣬ȷ���ö�����������ά�Ǳ߽�
//			int j_1 (j-1), j1(j+1);
//			if (j_1 <0)	j_1 = (int)outer.size ()-1;
//			if (j1 >= (int)outer.size ())	j1 = 0;
//			vector_3 n_temp(CGAL::cross_product(vector_3(ae.vp[outer[j]],ae.vp[outer[j1]]), vector_3(ae.vp[outer[j]], ae.vp[outer[j_1]])));
//			
//		}
//	}
//
//
//}

}
