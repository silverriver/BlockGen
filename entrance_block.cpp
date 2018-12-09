#include "entrance_block.h"
#include "geometry.h"
#include <vector>
#include "outer_poly.h"
#include <fstream>
#include <iostream>
using namespace std;

namespace BI
{

//判断一个点和一个单元块体的位置关系,vfi是一个返回值，表示如果这个点在该单元块体上，则记录这个点位于哪几个平面上
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
	//这个函数用来生成进入块体的外表面覆盖
	//E(A(0), B(2))
	for (int i(0);i!=(int)ob_b.vpoly.size ();++i)
	{
		vector_3 outer_n (ob_b.vplane[i].orthogonal_vector());	//块体B这个面的外法线方向
		if (!ob_b.outer_normal[i])
			outer_n = -outer_n;
		for (int j(0);j!=(int)ae_a.va.size ();++j)		//对于块体A的每个三维角
		{
			if (!ae_a.va[j].is_internal_intersected(outer_n))	//块体A的这个三维角和块体B这个面对应的半空间没有交集
			{
				//这个三维角与这个半空间不相交,这时E(ae_a.va[j].e, ob_b.vpoly[i])有可能是一个覆盖
				vcover.push_back (poly_frac());
				const point_3 & a = ae_a.vp[ae_a.va[j].e];
				construct_poly_frac(vcover.back (),ob_b.vpoly[i],ob_b.vplane[i],vector_3(a, a0),false);
			}
		}
	}

	//E(A(2), B(0))
	for (int i(0);i!=(int)ob_a.vpoly.size ();++i)
	{
		vector_3 outer_n (ob_a.vplane[i].orthogonal_vector());		//找到块体A这个边界对应的外表面法向
		if (!ob_a.outer_normal[i])
			outer_n = -outer_n;
		for (int j(0);j!=(int)ae_b.va.size ();++j)		//对于块体B的每个三维角
		{
			if (!ae_b.va[j].is_internal_intersected(outer_n))	//块体B的这个三维角和块体A这个面对应的半空间没有交集
			{
				//这个三维角和这个半空间不在内部相交，这时E(ob_a.vpoly[i], ae_b.va[j].e)有可能是一个覆盖
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
			if(!ae_a.ve[i].is_internal_intersected (ae_b.ve[j]))	//这两个边没有交集
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
		return 1;		//返回异常

	//先确定覆盖所处的范围。
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

	//生成一个大范围的方块研究区域
	if (minx == maxx || miny == maxy || minz ==maxz) 
		return 2;		//返回异常
	bs.add_rect_domain(minx-1,miny-1,minz-1,maxx+1,maxy+1,maxz+1);

	//把覆盖加到bs中（注意，一定要先初始化研究区域再初始化裂隙）
	for (int i(0);i!= (int)vcover.size ();++i)
		bs.add_poly_frac(vcover[i]);
	
	//debug_info
	if (debug_info)
		cout<<vcover.size ()<<" covers generated. "<<"Identifying entrance blocks"<<endl;

	//识别之
	bs.identify_block();

	//debug_info
	if (debug_info)
		cout<<"Determining entrance blocks"<<endl;

	//然后找出进入块体
	vector<int> vcb_in;		//记录生成的cblock中哪些在进入块体内
	vector<int> vcb_out;
	vcb_in.reserve (bs.vcb.size ()); vcb_out.reserve (bs.vcb.size ());
	
	for (int i(0);i!=(int) bs.vcb.size ();++i)		//对于识别生成的每个复合块体
	{
		point_3 p_temp (bs.vcb[i].get_rand_p(bs.veb,bs.vpo));		//判断这个复合块体是否在E(A,B)中。
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
	//生成覆盖
	outer_bound ob_a, ob_b;
	ob_a.construct_outer_bound (bs_, cba.vb);
	ob_b.construct_outer_bound (bs_, cbb.vb);

	angle_edge ae_a, ae_b;
	outer_bound2angle_edge(bs_,cba,ob_a,ae_a);
	outer_bound2angle_edge(bs_,cbb,ob_b,ae_b);

	return gen_entrance_block(bs_,cba,cbb,a0,ob_a,ae_a,ob_b,ae_b,debug_info);
}

//判断p1p2这个线段是否在vf和vp所表示的多边形内（即交点是否是一个点），要求该多边形是凸多边形，所有点都在pl上，并且p1肯定在这个多边形上（内或边界上）
bool segment_in_polygon (const vector<point_3> &vp, const vector<int> &vf, const plane_3 &pl, const point_3& e1, const point_3& e2)
{
	int ie1(-1), ie2(-1);		//记录e1在凸多边形的哪几个边上
	vector<plane_3> vpl;
	vector<CGAL::Oriented_side> vflag;
	vpl.reserve (vf.size ());
	vflag.reserve (vf.size ());
	//判断e1和这个多边形的相对位置
	for (int i(0);i!=vf.size ();++i)
	{
		int i1(i+1), i2(i+2);
		if (i1 >= (int)vf.size ())
			i1 = i1 - (int)vf.size ();
		if (i2 >= (int)vf.size ())
			i2 = i2 - (int)vf.size ();
		//判断e1在边i,i1的内侧，上，还是外侧
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
		return true;		//e1在多边形内部，所有线段与多边形的交点肯定不是一个点
	else if (ie2 < 0)
	{
		//e1在ie1所表示的这条边上（非顶点处）
		CGAL::Oriented_side tag_temp (vpl[ie1].oriented_side(e2));
		if (tag_temp == CGAL::Oriented_side::ON_ORIENTED_BOUNDARY)
			return true;		//e2也在这条边上
		if (tag_temp == vflag[ie1])
			return true;
		else
			return false;
	}
	else
	{
		//e1同时位于ie1和ie2这两个边上
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
	const vector<polygon_with_holes_2> &vpoly = ob.vpoly;	//编译器老是报错，干脆给这个变量起一个别名
	const vector<plane_3> &vplane = ob.vplane;
	
	ae.va.clear(); ae.ve.clear (); ae.vp.clear ();
	typedef vector<int> loop;			//顶点索引的回路
	typedef vector<loop> polygon_with_holes_3;		//第一个loop是外法线方向

	vector<polygon_with_holes_3> vpwh3;
	vpwh3.reserve (vpoly.size ());

	//先将每个点都转化到三维空间中
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

	//下面的就是生成各个三维角
	ae.va.reserve (ae.vp.size ());		//理论上是一个顶点对应一个角
	typedef struct{int ebi; vector<int> vfi;} v_eb_info;		//定义一个结构来记录某个顶点在哪个单元块体(ebi)的哪些面(vfi)上
	vector<vector<v_eb_info>> v_angle_suba_info;		//记录每个三维角的每个子角的构成信息.这些信息可以直接包含在每个三维角中，但是如果这样做会导致angle_edge中的数据和bs中的数据无法做到独立。所以暂时还是将这两个部分数据分开，减少类之间的耦合。可能后来会发现这样的努力是徒劳的，不过那是以后需要考虑的事情了。2016-2-19 20:35
	v_angle_suba_info.reserve (ae.va.size ());

	for (int i(0);i!=ae.vp.size ();++i)		//对于每个顶点
	{
		ae.va.push_back (angle_edge::angle());
		ae.va.back ().e = i;
		ae.va.back ().vsa.reserve (4);	//暂且估计每个三维角包含4个子角
		v_angle_suba_info.push_back (vector<v_eb_info>());
		v_angle_suba_info.back ().reserve (4);

		for (int j(0);j!=(int)cb.vb.size ();++j)		//判断这个顶点和每个单元块体的位置关系
		{
			v_eb_info info_temp;
			info_temp.vfi.reserve (bs.veb[cb.vb[j]].vf.size ());
			info_temp.ebi = cb.vb[j];
			switch (point_eblock_location (bs,bs.veb[cb.vb[j]],ae.vp[i],info_temp.vfi))
			{
			case Oriented_side::IN_:
				throw logic_error ("outer_bound2angle_edge, point in eblock");		//这个顶点不可能在单元块体内
				break;
			case Oriented_side::ON_:
				ae.va.back ().vsa.push_back (angle_edge::sub_angle());			//如果顶点在单元块体上，则vfi中对应的面可以组成该顶点的一个子角
				ae.va.back ().vsa.back ().vnormal.reserve (info_temp.vfi.size ());
				for (int k(0);k!=(int)info_temp.vfi.size ();++k)		//添加子角
				{
					if (bs.veb[cb.vb[j]].vf[info_temp.vfi[k]].outer_normal )
						ae.va.back ().vsa.back ().vnormal.push_back (bs.vpl[bs.veb[cb.vb[j]].vf [info_temp.vfi[k]].plane_id].orthogonal_vector());
					else
						ae.va.back ().vsa.back ().vnormal.push_back (-bs.vpl[bs.veb[cb.vb[j]].vf [info_temp.vfi[k]].plane_id].orthogonal_vector());
				}
				v_angle_suba_info.back ().push_back (info_temp);		//记录这个子角的生成信息
			default:
				break;
			}
		}
	}

	//下面就是生成各个边
	ae.ve.reserve (ae.va.size () + ob.vpoly.size () - 2);	//先假设有v+f-2条棱
	for (int i(0);i!=(int) vpwh3.size ();++i)		//先找出所有的边
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
	for (int i(0);i!=(int) ae.ve.size ();++i)		//对于每条三维边
	{
		vector_3 edge_vector (ae.vp[ae.ve[i].e1], ae.vp[ae.ve[i].e2]);		//先生成下这个三维边对应的向量
		ae.ve[i].vse.reserve (3);		//假设每个边有三个子边，先预留位置
		vector<v_eb_info> & angle_suba_info = v_angle_suba_info[ae.ve[i].e1];	//找到每个三维边的第一个顶点的生成信息
																				//注意，这里只考虑了第一个顶点的生成信息，因为假设用来生成边的外表面多边形结构ob是规范化之后的。所以，如果用第二个顶点的信息生成各个子边，那么所生成的二维边是一样的，但是组成这个二维边的子边可能不尽相同
		for(int j(0);j!= (int)angle_suba_info.size ();++j)		//对于这个顶点所在的每个单元块体
		{
			vector<int> vfi;	//记录edge_vector在这个单元块体的哪些面上
			vfi.reserve (2);	//预留2个，因为理论上不可能出现3个。
			int ebi = angle_suba_info[j].ebi;		//临时变量，指示现在在处理哪个单元块体

			for (int k(0);k!=(int) angle_suba_info[j].vfi.size ();++k)
			{
				int vfi_temp = angle_suba_info[j].vfi[k];		//临时变量，指示现在在处理ebi的哪个面
				int plane_id_temp = bs.veb[ebi].vf[vfi_temp].plane_id;
				if ((bs.vpl[plane_id_temp].orthogonal_vector() * edge_vector == 0) &&
					segment_in_polygon(bs.vpo,bs.veb[ebi].vf[vfi_temp].np,bs.vpl[plane_id_temp],ae.vp[ae.ve[i].e1], ae.vp[ae.ve[i].e2]))
				{
					vfi.push_back (vfi_temp);//该平面与这条边平行，并且这条边也和该平面上的三角形相交。
				}
				
			}
			if (vfi.size () >= 3)
				throw logic_error("outer_bound2angle_edge, edge in more than three planes");
			if (vfi.size () == 1)
			{
				ae.ve[i].vse.push_back (angle_edge::sub_edge());		//这条边落在这个单元块体的某个面上
				if (bs.veb[ebi].vf[vfi[0]].outer_normal)		//sub_edge中n1和n2指向子边的外部，注意，是外部！！
					ae.ve[i].vse.back ().n1 = ae.ve[i].vse.back ().n2 = bs.vpl[bs.veb[ebi].vf[vfi[0]].plane_id].orthogonal_vector();
				else
					ae.ve[i].vse.back ().n1 = ae.ve[i].vse.back ().n2 = -bs.vpl[bs.veb[ebi].vf[vfi[0]].plane_id].orthogonal_vector();
			}
			if (vfi.size () == 2)
			{
				ae.ve[i].vse.push_back (angle_edge::sub_edge());		//这条边落在这个单元块体的某条棱上
				if (bs.veb[ebi].vf[vfi[0]].outer_normal )		//sub_edge中n1和n2指向子边的外部，注意，是外部！！
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
			throw logic_error ("outer_bound2angle_edge, vsubedge is empty");	//三维边的子边集合不能为空
	}
	return 0;
}


//int outer_bound2angle_edge (const outer_bound &ob, angle_edge &ae)
//{

//

//	
//	
//	typedef pair<int,int> angle_2d;		//用一个整型对来表示一个二维角
//	vector<vector<angle_2d>> vbound(ae.vp.size ());	//每个三维角的边界
//
//	for (int i(0);i!=(int)vpwh3.size ();++i)
//	{
//		loop &outer = vpwh3[i][0];		//为outer_boundary起一个别名
//		vector_3 outer_normal (vplane[i].orthogonal_vector());
//		if (!ob.outer_normal[i]) outer_normal = -outer_normal;		//确定这个平面的外法向
//		vector_3 oriented_n (oriented_normal(ae.vp,outer));			//按照右手螺旋法则确定outer上的点的旋转方向
//		bool flag (false);		//T表示outer中的点旋转方向和外法线方向一致，F表示outer中的点旋转方向和外法线方向相反
//		if (outer_normal*oriented_n > 0)
//			flag = true;		
//
//		for (int j(0);j!=outer.size ();++j)
//		{
//			vector<angle_2d> & bound = vbound[outer[j]];		//对于outer中的每个顶点，确定该顶点所属的三维角边界
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
