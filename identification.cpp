#define S_DEBUG
#include "geometry.h"
#include "GBUtility.h"
#include "identification.h"
#include "safety_factor.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <CGAL\squared_distance_3.h>
#include <CGAL\Timer.h>
#include <CGAL\Memory_sizer.h>
#include <CGAL\Boolean_set_operations_2.h>
#include <exception>
using namespace std;

namespace BI
{


//����-1��ʾ��������ĵ�������3��
//����-2��ʾĳ����ֻ������С���������ڣ�����������ȷ������һ����
//����-3��ʾĳ��������ڴ����������ڣ�����������Ҳ����ȷ��һ���㣬������������б���Ϊ�˽�ʡ�鷳����Ϊ��������������������Ӧ��exact���㣬�������ɿ����������
int block_system::init_domain(const GB_domain &d,const double err)
{
#ifdef S_DEBUG
	CGAL::Memory_sizer m_size;
	CGAL::Timer timer;
	timer.start();
#endif
	clear();
	vector<int> v_plane;		//ָʾd��ÿ��ƽ����vpl�е�����
	vector<bool> plane_flag;	//���d�е�ÿ��ƽ���Ƿ����������,��ʵ�������������Ӧ�ý�facet_flag


	//����d��ÿһ��face��Ӧ��ƽ�淽��
	for (int i(0);i!=d.GetFaceCount();++i)
	{
#ifdef S_DEBUG
		//cout<<i<<endl;
#endif
		if (d.GetFacePointCount(i)<3)
			return -1;	//��������ĵ�������3�����϶�����

		//ֻ��ǰ���������ƽ�淽��
		plane_3 plane_temp(d.GetFacePoint(i,0), d.GetFacePoint(i,1), d.GetFacePoint(i,2));
		FT res(0);	//���������㵽ƽ�淽�̵Ĳв�
		for (int j(0);j!=d.GetFacePointCount(i);++j)
			res += CGAL::squared_distance(d.GetFacePoint(i,j),plane_temp);
		
		if (res == 0)
		{
			//������ƽ����������ϣ���ôֱ��������ƽ��
			v_plane.push_back (add_set_pl(vpl,plane_temp));
			plane_flag.push_back (true);	//�������
			continue;
		}
		else
		{
#ifdef S_DEBUG
			cout<<res<<endl;	//DEBUG
#endif
			plane_flag.push_back (false);
		}

		//�����������ɵ�ƽ��
		vpl.push_back (plane_temp);
		v_plane.push_back ((int)vpl.size ()-1);

		if (i==0)
			continue;	//��һ��ƽ��ֱ�����
		//�����ж��Ƿ�ǰ���Ƿ��Ѿ�����˺�plane_temp���������ƽ��
		int j(0);
		for (;j<i;++j)
		{
			//����ÿһ���Ѿ�����˵�facet�����plane_temp�����facet���ڵ�ƽ���Ƿ��غ�
			int count(0);
			for (;count!=d.GetFacePointCount(j);++count)
				if (CGAL::squared_distance(d.GetFacePoint(j,count),plane_temp) > err)
					break;	//�����ĳ����ľ����������ޣ���ô���facet��plane_tempһ�����غ�
			if (count == d.GetFacePointCount(j))
				break;		//plane_temp��j���ƽ���غ�
		}
		
		if (j!=i)
		{
			//���������plane_temp��j���ƽ���غϣ���ô�Ͳ����plane_temp
			vpl.pop_back();
			v_plane.back () = v_plane[j];
		}
	}

#ifdef S_DEBUG
	timer.stop();
	cout<<"Plane ganerated complete. Time:"<<timer.time()<<endl;
	cout<<"Memory:"<<m_size.virtual_size()/1024.0/1024.0<<"Mb"<<endl;
	timer.reset();
	timer.start();
#endif

	//���ÿһ��fram
	vector<int> po_index(d.GetPointCount(),-1);	//��¼d�е�ÿ������vpo�е���������ʼֵ��-1,��ʾ��û����ӵ�vpo�У�
	for (int i(0);i!=d.GetFramCount();++i)
	{
		veb.push_back (eblock());	//ÿ��fram����һ��element_block
		eblock &eb = veb.back ();
		eb.cblock_index = i;	//ÿ�����嶼��һ���µ�����
		vector<int> po_in_fram;				//��¼���fram���ж��ٸ���,��¼��d.vp�е�����
		vector<vector<int>> po_in_face;		//��¼fram�е�ÿ��������Щface��,ÿ�����Ӧһ��face��������������Ӧ��d.fa�е�����

		//Ѱ�Ҹ�fram�е�ÿ���㣬����Ѱ��ÿ�������ڵ�ƽ��
		for (int j(0);j!=d.GetFramFaceCount(i);++j)
		{
			for (int k(0);k!=d.GetFramFace(i,j).pn.size ();++k)
			{
				//���ÿ����
				size_t ori_size = po_in_fram.size ();
				size_t position = add_set(po_in_fram,d.GetFramFace(i,j).pn[k]);
				if (po_in_fram.size () > ori_size)
					po_in_face.push_back (vector<int>());	//��po_in_fram��������˵�
				po_in_face[position].push_back (d.GetFram(i).nf[j]);
			}
		}

		for (int j(0);j!=po_in_face.size ();++j)
		{
			if (po_in_face[j].size () < 3)
				return -2;		//�������С��3��ƽ����
			if (po_in_face[j].size () > 3)
				return -3;		//������ڴ���3��ƽ����
		}

		for (int j(0);j!=po_in_fram.size ();++j)
		{
			if(po_index[po_in_fram[j]] != -1)
				continue;	//���������Ѿ���ӹ��ˣ���ôֱ������

			int facet0 = po_in_face[j][0];
			int facet1 = po_in_face[j][1];
			int facet2 = po_in_face[j][2];

			if( (plane_flag[facet0] || vpl[v_plane[facet0]].has_on(d.p[po_in_fram[j]]))
				&& (plane_flag[facet1] || vpl[v_plane[facet1]].has_on(d.p[po_in_fram[j]]))
				&& (plane_flag[facet2] || vpl[v_plane[facet2]].has_on(d.p[po_in_fram[j]])))
			{
				po_index[po_in_fram[j]] = (int)vpo.size ();	//��������������ȷ��������ƽ���ϣ���ôֱ�����
				vpo.push_back (d.p[po_in_fram[j]]);
			}
			else
			{
				//����������������꣬Ȼ�����֮
				po_index[po_in_fram[j]] = add_set(vpo,intersection_3planes(vpl[v_plane[facet0]],vpl[v_plane[facet1]],vpl[v_plane[facet2]]));
			}
		}

		for (int j(0);j!=d.GetFramFaceCount(i);++j)
		{
			//���ÿ��face����Ԫ����
			eb.vf.push_back (eblock::face());
			eblock::face & new_f = eb.vf.back();
			new_f.frac_id = -(1+d.GetFaceType(d.GetFram(i).nf[j]));	//�ø�����������ͱ����������϶���
			new_f.plane_id = v_plane[d.GetFram(i).nf[j]];			//��¼��������ڵ�ƽ����
			for (int k(0);k!=d.GetFramFacePointCount(i,j);++k)
				new_f.np.push_back (po_index[d.GetFramFace(i,j).pn[k]]);
		}
		eb.attribute = d.GetDensity(i);

		//����ÿ������ⷨ�߷��򣬲��ҿ��ܻᷴת������Ķ�������˳��
		for (int j(0);j!=eb.vf.size ();++j)
		{
			int p_temp = eb.p_inb_outf(j);	//��һ���������ƽ���ϵĵ�

			CGAL::Oriented_side side = vpl[eb.vf[j].plane_id].oriented_side(vpo[p_temp]);
			if (side == CGAL::ON_POSITIVE_SIDE)
				eb.vf[j].outer_normal = false;	//�������һ����p_temp��ƽ������࣬��ô���ⷨ�߷���Ӧ��ָ����
			else if (side == CGAL::ON_NEGATIVE_SIDE)
				eb.vf[j].outer_normal = true;
			else
				throw logic_error ("init_domain, point on plane");

			//�ж�һ�������Ķ��������Ƿ�����ʱ������,�������������˸�facet��͹���ԭ�򡣣�
			vector_3 v_temp (oriented_normal(vpo,eb.vf[j].np,vpl[eb.vf[j].plane_id]));

			FT dot_res (v_temp*vpl[eb.vf[j].plane_id].orthogonal_vector());
			if (eb.vf[j].outer_normal == false)
				dot_res = -dot_res;

			if (dot_res < 0)
				reverse(eb.vf[j].np.begin (),eb.vf[j].np.end ());
		}
	}
#ifdef S_DEBUG
	timer.stop();
	cout<<"fram added complete time:"<<timer.time()<<endl;
	cout<<"Memory:"<<m_size.virtual_size()/1024.0/1024.0<<"Mb"<<endl;
#endif
	return 0;
}

//������飬����0��ʾ���ͨ����������������Ϣ
int block_system::init_check () const
{
	//ÿ��������͹��
	for (int i(0);i!=veb.size ();++i)
	{
		if (veb[i].vf.size ()<4)
			return 3;	//����4����
		if (!eb_convex(i))
			return 1;	//���eblock����͹����
		if (!eb_sealed(i))
			return 2;	//���eblock�����
		for (int j(0);j!=vpo.size ();++j)
		{
			Oriented_side res = eb_p_location(i,vpo[j]);
			if (res == Oriented_side::IN_)
				return 4;		//�е�������
			if (res == Oriented_side::OUT_)
				continue;
			if (!veb[i].p_is_b_index(j))
				return 5;
		}
	}
	//������ô����ɣ���Ȼ������һ�������ʱ��ȽϷ�ʱ�䣬����Ҳ���в��˼���
	return 0;
}

bool block_system::eb_convex (const int &i) const
{
	for (int j(0);j!=veb[i].vf.size ();++j)
	{
		if(veb[i].vf[j].np.size () <3)
			return false;		//һ��������3����
		//�ж����ƽ���ϵĵ㶼�����ƽ���ϣ�����������㶼�����ƽ���һ��
		const plane_3 &p = vpl[veb[i].vf[j].plane_id];
		for (int k(0);k!=veb[i].vf[j].np.size ();++k)
			if (!p.has_on (vpo[veb[i].vf[j].np[k]]))
				return false;	//������ϵĵ�Ӧ�ö������ƽ����
		for (int k(0);k!=veb[i].vf.size ();++k)
		{
			if (j==k)
				continue;
			for (int l(0);l!=veb[i].vf[k].np.size ();++l)
			{
				if (find(veb[i].vf[j].np.begin (),veb[i].vf[j].np.end(),veb[i].vf[k].np[l]) == veb[i].vf[j].np.end())
					if (p.has_on_positive_side(vpo[veb[i].vf[k].np[l]]) == veb[i].vf[j].outer_normal)
						return false;		//�����������Ӧ�������ƽ���һ��
			}
		}
	}
	return true;
}

bool block_system::eb_sealed (const int &i) const
{
	vector<segment_3> vs;
	for (int j(0);j!=veb[i].vf.size ();++j)
	{

		for (int k(0);k!=(int)veb[i].vf[j].np.size ();++k)
		{
			int k1(k+1);
			if (k1 == (int)veb[i].vf[j].np.size ())
				k1=0;
			vs.push_back (segment_3(vpo[veb[i].vf[j].np[k]],vpo[veb[i].vf[j].np[k1]]));
		}
	}

	vector<int> seg_index1 (vs.size (),0);	//��һ��vs���ж���segment�͸���seg��ͬ����������
	vector<int> seg_index2 (vs.size (),0);	//��һ��vs���ж���segment�͸���seg�����෴��

	for (int j(0);j!=vs.size ();++j)
	{
		if (vs[j].is_degenerate())
			return false;	//�߶����˵��غ�
		for (int k(0);k!=vs.size ();++k)
		{
			if (j==k)
				continue;
			if (vs[j]==vs[k])
				seg_index1[j]++;
			if (vs[j] == vs[k].opposite())
				seg_index2[j]++;
		}
	}
	for (int j(0);j!=vs.size ();++j)
	{
		if (seg_index1[j] != 0)
			return false;	//������һ��seg�����seg�غ�
		if (seg_index2[j] != 1)
			return false;	//���ڲ�ֹһ��seg�����seg����
	}
	return true;
}

FT block_system::eb_volume (const int &ei) const
{
	FT res(0);
	for (int i(0);i!=veb[ei].vf.size ();++i)
	{
		for (int j(1);j!=veb[ei].vf[i].np.size ();++j)
		{
			int j1 = j+1;
			if (j1>=veb[ei].vf[i].np.size ())
				j1 = 0;
			const point_3 &p1=vpo[veb[ei].vf[i].np[0]];
			const point_3 &p2=vpo[veb[ei].vf[i].np[j]];
			const point_3 &p3=vpo[veb[ei].vf[i].np[j1]];
			res = res+(p1.x()*p2.y()*p3.z() + p1.y()*p2.z()*p3.x() + p1.z()*p2.x()*p3.y() - p1.z()*p2.y()*p3.x() - p1.y()*p2.x()*p3.z() - p1.x()*p2.z()*p3.y());
		}
	}
	return res/6;
}

FT block_system::cb_volume (const int &cbi) const
{
	FT res (0);
	for (int i(0);i!=vcb[cbi].vb.size ();++i)
		res = res + eb_volume(vcb[cbi].vb[i]);
	return res;
}

double block_system::eb_volume_est (const int &ei) const
{
	double res(0);
	for (int i(0);i!=veb[ei].vf.size ();++i)
	{
		for (int j(1);j!=veb[ei].vf[i].np.size ();++j)
		{
			int j1 = j+1;
			if (j1>=veb[ei].vf[i].np.size ())
				j1 = 0;
			const point_3 &p1=vpo[veb[ei].vf[i].np[0]];
			const point_3 &p2=vpo[veb[ei].vf[i].np[j]];
			const point_3 &p3=vpo[veb[ei].vf[i].np[j1]];
			double x1(CGAL::to_double(p1.x())), x2 (CGAL::to_double(p2.x())), x3(CGAL::to_double(p3.x()));
			double y1(CGAL::to_double(p1.y())), y2 (CGAL::to_double(p2.y())), y3(CGAL::to_double(p3.y()));
			double z1(CGAL::to_double(p1.z())), z2 (CGAL::to_double(p2.z())), z3(CGAL::to_double(p3.z()));

			res += (x1*y2*z3 + y1*z2*x3 + z1*x2*y3 - z1*y2*x3 - y1*x2*z3 - x1*z2*y3);
		}
	}
	return res/double (6);
}

double block_system::cb_volume_est (const int &cbi) const
{
	double res(0);
	for (int i(0);i!=vcb[cbi].vb.size ();++i)
		res += eb_volume_est(vcb[cbi].vb[i]);
	return res;
}

Oriented_side block_system::eb_p_location (const int& ei, const point_3 &p)const
{
	const eblock &eb = veb[ei];
	vector<Oriented_side> vside(eb.vf.size ());
	for (int i(0);i!=eb.vf.size ();++i)
	{
		CGAL::Oriented_side side= vpl[eb.vf[i].plane_id].oriented_side(p);
		if (side == CGAL::Oriented_side::ON_ORIENTED_BOUNDARY)
			vside[i] = Oriented_side::ON_;
		else if (side == CGAL::Oriented_side::ON_NEGATIVE_SIDE)
			if (eb.vf[i].outer_normal == true)
				vside[i] = Oriented_side::IN_;
			else
				vside[i] = Oriented_side::OUT_;
		else
			if (eb.vf[i].outer_normal == true)
				vside[i] = Oriented_side::OUT_;
			else
				vside[i] = Oriented_side::IN_;
	}
	for (int i(0);i!=vside.size ();++i)
	{
		if (vside[i] == Oriented_side::OUT_)
			return Oriented_side::OUT_;
	}
	for (int i(0);i!=vside.size ();++i)
	{
		if (vside[i] == Oriented_side::ON_)
			return Oriented_side::ON_;
	}
	return Oriented_side::IN_;
}

int block_system::common_plane (const eblock& eb1, const eblock&eb2) const
{
	for (int i(0);i!=eb1.vf.size ();++i)
	{
		int pl_index = eb1.vf[i].plane_id;
		for (int j(0);j!=eb2.vf.size ();++j)
			if (eb2.vf[j].plane_id == pl_index)
				return pl_index;
	}
	return -1;
}

int block_system::disc_frac_parser (std::istream &infile, const FT &dx)
{
	if (!infile)
		return 1;
	int count (0);
	infile>>count;
	for (int i(0);i!=count;++i)
	{
		vdisc.push_back (disc_frac());
		disc_frac &df = vdisc.back ();
		FT x,y,z,dip_dir,dip;
		infile>>x>>y>>z>>dip_dir>>dip>>df.r>>df.aperture>>df.cohesion>>df.f_angle;
		df.center = point_3(x,y,z);
		df.plane_id = add_set_pl(vpl,plane_3(df.center,dip2normal(dip_dir,dip,dx)));	//�ǳ���ʱ�Ĳ���
	}
	sort(vdisc.begin (), vdisc.end (),radius_bigger);	//����϶���ŵ�ǰ��
	for (int i(0);i!=vdisc.size ();++i)
		vdisc[i].frac_id = i;
	return 0;
}

int block_system::poly_frac_parser (std::istream &infile)
{
	/*******************************
	�������϶�ļ���ʽ��
	����϶��
		����ÿ����϶��
		��϶�ڵ��� ��϶��ճ��ϵ����Ħ����
			����ÿ���ڵ㣺
			�ڵ�x, y, z����
	�������϶���Բ���͹�����
	*******************************/
	//�����ļ���ʽ�иĶ����պ���ʵ�����������
	if (!infile)
		return 1;
	int frac_count(0);
	infile>>frac_count;
	for (int i(0);i!=frac_count;++i)
	{
		poly_frac pf;
		int fp_count(0);
		infile>>fp_count;
		pf.vp_ob.reserve (fp_count);

		infile>>pf.aperture>>pf.cohesion>>pf.f_angle;
		for (int j(0);j!=fp_count;++j)
		{
			FT x,y,z;
			infile>>x>>y>>z;
			pf.vp_ob.push_back (point_3(x,y,z));
		}
		pf.frac_id = (int)vdisc.size () + i;
		add_poly_frac(pf);
	}
	return 0;
}

//--------------Ϊ�������������stl_parser������������
struct point_t 
{double x,y,z;};

struct triangle_t
{
	int t1,t2,t3;
	triangle_t (int t1_, int t2_, int t3_): t1(t1_),t2(t2_), t3(t3_){}
};

class point_equal
{
	double err;		//Ӧ����һ������
public:
	point_equal(const double& err_):err(err_){}
	bool operator() (const point_t&p1, const point_t&p2) const
	{
		if (abs(p1.x-p2.x)<err && abs(p1.y-p2.y)<err && abs(p1.z-p2.z)<err)
			return true;
		else
			return false;
	}
};

class tri_equal
{
public:
	tri_equal () {}
	bool operator() (const triangle_t&t1, const triangle_t& t2) const
	{
		if ((t1.t1 == t2.t1 && t1.t2 == t2.t2 && t1.t3 == t2.t3) ||
			(t1.t1 == t2.t1 && t1.t2 == t2.t3 && t1.t3 == t2.t2) ||
			(t1.t1 == t2.t2 && t1.t2 == t2.t1 && t1.t3 == t2.t3) ||
			(t1.t1 == t2.t2 && t1.t2 == t2.t3 && t1.t3 == t2.t1) ||
			(t1.t1 == t2.t3 && t1.t2 == t2.t2 && t1.t3 == t2.t1) ||
			(t1.t1 == t2.t3 && t1.t2 == t2.t1 && t1.t3 == t2.t2))
			return true;
		else
			return false;
	}
};
//---------------------

int block_system::stl_parser (std::istream &infile, const double &err)
{
	if (!infile)
		return 1;	//ֱ�ӷ���
	string kw1,kw2;

	vector<point_t> vp;
	vector<triangle_t> vt;
	char temp('a');
	while (temp!='\n')
		infile.get (temp);

	point_equal pe(err);
	tri_equal te;
	while (infile)
	{
		point_t p[3];
		infile>>kw1>>kw2;
		if (kw1 == "endsolid")
			break;
		if (kw1 != "facet")
			return 2;
		getline (infile,kw1);
		infile>>kw1>>kw2;
		if (kw1 !="outer" || kw2 != "loop")
			return 2;
		for (int i(0);i!=3;++i)
		{
			infile>>kw1;
			if (kw1 != "vertex")
				return 2;
			infile>>p[i].x>>p[i].y>>p[i].z;
		}
		infile >> kw1>>kw2;
		if (kw1 != "endloop" || kw2!= "endfacet")
			return 2;
		int n1 = add_set (vp, p[0], pe);
		int n2 = add_set (vp, p[1], pe);
		int n3 = add_set (vp, p[2], pe);
		add_set (vt,triangle_t (n1,n2,n3),te);
	}
	for (int i(0);i!=vt.size ();++i)
	{
		poly_frac pf;
		pf.aperture = pf.cohesion = pf.f_angle = -1;
		pf.vp_ob.push_back (point_3 (vp[vt[i].t1].x, vp[vt[i].t1].y, vp[vt[i].t1].z));
		pf.vp_ob.push_back (point_3 (vp[vt[i].t2].x, vp[vt[i].t2].y, vp[vt[i].t2].z));
		pf.vp_ob.push_back (point_3 (vp[vt[i].t3].x, vp[vt[i].t3].y, vp[vt[i].t3].z));
		if (add_poly_frac(pf) !=0)
			return 3;
	}
	return 0;
}

int block_system::add_poly_frac (const poly_frac &pf)
{
	vector_3 normal(oriented_normal(pf.vp_ob));
	if (normal == vector_3(0,0,0))		//�����������ι��棬��ֱ���������������
		return 1;
	vpoly.push_back (pf);
	vpoly.back ().plane_id = add_set_pl(vpl,plane_3(pf.vp_ob[0],normal));
	vpoly.back ().frac_id = (int)vdisc.size () + (int)vpoly.size () - 1;
	vpoly.back().proj_plane(vpl[vpoly.back ().plane_id]);		//������������϶�ϵĵ㶼ͶӰ��ͬһ����
	return 0;
}

int block_system::add_disc_frac (const disc_frac &df, const vector_3 &normal)
{
	vdisc.push_back (df);
	vdisc.back ().plane_id = add_set_pl (vpl,plane_3(df.center,normal));
	vdisc.back ().frac_id = (int)vdisc.size ()-1;
	return 0;
}

int block_system::add_rect_domain (const FT &x1, const FT &y1, const FT &z1, 
	const FT &x2, const FT &y2, const FT &z2)
{
	clear();
	vpo.push_back (point_3(x1,y1,z1));
	vpo.push_back (point_3(x1,y2,z1));
	vpo.push_back (point_3(x2,y2,z1));
	vpo.push_back (point_3(x2,y1,z1));
	vpo.push_back (point_3(x1,y1,z2));
	vpo.push_back (point_3(x1,y2,z2));
	vpo.push_back (point_3(x2,y2,z2));
	vpo.push_back (point_3(x2,y1,z2));

	vpl.push_back (plane_3(vpo[2],vector_3(0,0,-1)));
	vpl.push_back (plane_3(vpo[4],vector_3(0,0,1)));
	vpl.push_back (plane_3(vpo[4],vector_3(0,-1,0)));
	vpl.push_back (plane_3(vpo[2],vector_3(0,1,0)));
	vpl.push_back (plane_3(vpo[2],vector_3(1,0,0)));
	vpl.push_back (plane_3(vpo[4],vector_3(-1,0,0)));

	veb.push_back (eblock());
	eblock &eb=veb.back ();
	eb.cblock_index = 0;		//�µĿ�����
	eb.vf.push_back (eblock::face ());
	eb.vf.back ().frac_id = -(1+2);		//����GB�еļǺţ��������������
	eb.vf.back ().plane_id = 0;
	eb.vf.back ().np .push_back (0);eb.vf.back ().np .push_back (1);eb.vf.back ().np .push_back (2);eb.vf.back ().np .push_back (3);
	eb.vf.back ().outer_normal = true;

	eb.vf.push_back (eblock::face());
	eb.vf.back ().frac_id = -(1+2);		
	eb.vf.back ().plane_id = 1;
	eb.vf.back ().np .push_back (4);eb.vf.back ().np .push_back (7);eb.vf.back ().np .push_back (6);eb.vf.back ().np .push_back (5);
	eb.vf.back ().outer_normal = true;

	eb.vf.push_back (eblock::face());
	eb.vf.back ().frac_id = -(1+2);		
	eb.vf.back ().plane_id = 2;
	eb.vf.back ().np .push_back (0);eb.vf.back ().np .push_back (3);eb.vf.back ().np .push_back (7);eb.vf.back ().np .push_back (4);
	eb.vf.back ().outer_normal = true;

	eb.vf.push_back (eblock::face());
	eb.vf.back ().frac_id = -(1+2);		
	eb.vf.back ().plane_id = 3;
	eb.vf.back ().np .push_back (1);eb.vf.back ().np .push_back (5);eb.vf.back ().np .push_back (6);eb.vf.back ().np .push_back (2);
	eb.vf.back ().outer_normal = true;

	eb.vf.push_back (eblock::face());
	eb.vf.back ().frac_id = -(1+2);		
	eb.vf.back ().plane_id = 4;
	eb.vf.back ().np .push_back (2);eb.vf.back ().np .push_back (6);eb.vf.back ().np .push_back (7);eb.vf.back ().np .push_back (3);
	eb.vf.back ().outer_normal = true;

	eb.vf.push_back (eblock::face());
	eb.vf.back ().frac_id = -(1+2);		
	eb.vf.back ().plane_id = 5;
	eb.vf.back ().np .push_back (4);eb.vf.back ().np .push_back (5);eb.vf.back ().np .push_back (1);eb.vf.back ().np .push_back (0);
	eb.vf.back ().outer_normal = true;

	eb.attribute = 2.7;		//�������ֵ
	return 0;
}

int block_system::plane_cut_eblock (const int &targeti, const int &pli, const int &frac_id, const int &cblock_index)
{
	eblock & eb = veb[targeti];

	for (int i(0);i!=eb.vf.size ();++i)
		if (eb.vf[i].plane_id == pli)
			return 1;	//pli��eb��һ������棬�����иֱ������

	const plane_3 &pl = vpl[pli];
	eblock res1, res2;
	vector<vector<CGAL::Oriented_side>> fp_flag(eb.vf.size ());	//��¼eb��ÿ������pl����һ��
	for (int i(0);i!=eb.vf.size ();++i)
		fp_flag[i].resize (eb.vf[i].np.size ());

	for (int i(0);i!=eb.vf.size ();++i)
	{
		for (int j(0);j!=eb.vf[i].np.size ();++j)
			fp_flag[i][j] = pl.oriented_side(vpo[eb.vf[i].np[j]]);
	}

	vector<int> f_flag (eb.vf.size ());		//��¼ÿ������pl��ʲôλ�ã�1��ʾ���࣬-1��ʾ���࣬0��ʾһ������һ�븺��
	for (int i(0);i!=eb.vf.size ();++i)
	{
		int p_count(0),n_count(0),o_count(0);	//�������ֱ��ʾ���࣬���࣬�����ϵĵ���
		for (int j(0);j!=eb.vf[i].np.size ();++j)
		{
			switch (fp_flag[i][j])
			{
			case CGAL::Oriented_side::ON_POSITIVE_SIDE:
				++p_count; break;
			case CGAL::Oriented_side::ON_NEGATIVE_SIDE:
				++n_count; break;
			case CGAL::Oriented_side::ON_ORIENTED_BOUNDARY:
				++o_count; break;
			}
		}
		if (p_count == 0 && n_count == 0)
			throw logic_error ("plane_cut_eblock, face on plane, eb.plane_id ERROR");	//������ϵĵ�ȫ����pl�ϣ����ǲ����ܵģ��ں���һ��ʼ���Ѿ��������������
		if (o_count >=3)
			throw logic_error ("plane_cut_eblock, face on plane, eb.plane_id ERROR");	//����������������ϵĵ���pl�ϣ����ǲ����ܵġ�
		if (p_count == 0)
		{f_flag[i] = -1; continue;}
		if (n_count == 0)
		{f_flag[i] = 1; continue;}
		f_flag[i] = 0;
	}
	
	int p_count(0),n_count(0),c_count(0);		//�������ֱ��ʾ��pl���࣬���࣬���������
	for (int i(0);i!=eb.vf.size ();++i)
	{
		switch (f_flag[i])
		{
		case 1:
			++p_count; break;
		case -1:
			++n_count; break;
		case 0:
			++c_count; break;
		}
	}
	if (p_count == eb.vf.size () || n_count == eb.vf.size ())
		return 1;		//��ȫ����������߸��࣬�����и�

	//��ʽ�и�eb
	res1.attribute = res2.attribute = eb.attribute;
	res1.cblock_index = cblock_index;
	res2.cblock_index = cblock_index+1;
	typedef pair<int,int> int_segment;
	vector<int_segment>	vseg;	//������¼�����ϣ�Ҳ���������ɵ����ϵı�
	vector<int_segment> intersected_seg;		//������¼�Ѿ������������߶Ρ�
	vector<int> calculated_pi;					//������¼�����Ѿ������������߶εĽ������ĸ�
	for (int i(0);i!=eb.vf.size ();++i)
	{
		vector<int> cutseg; cutseg.reserve(2);		//��¼����������������ڶ�����
		eblock * res = 0;
		switch (f_flag[i])
		{
		case 1:
			res = &res1;		//res1��ʾ���������
		case -1:
			if (res == 0)
				res = &res2;				//res2��ʾ�ڸ������

			for (int j(0);j!=eb.vf[i].np.size ();++j)
				if (fp_flag[i][j] == CGAL::Oriented_side::ON_ORIENTED_BOUNDARY)
					cutseg.push_back (eb.vf[i].np[j]);		//��¼һ�¿����ڶ����ϵ��߶�

			res->vf.push_back (eblock::face());
			res->vf.back().swap(eb.vf[i]);		//����ebҲҪ���������eb.vf[i]�Ѿ�û�������ü�ֵ��ֱ�ӽ�����Ϣswap��������ʡʱ�䡣res1��ʾ������Ĳ���

			break;
		case 0:
			//�ָ������
			res1.vf.push_back (eb.vf[i]); eblock::face &res1f(res1.vf.back ()); res1f.np.clear();		//����汻�г�������,���������������Ϣ��ԭ������ͬ��ֻ�Ƕ������в�ͬ
			res2.vf.push_back (eb.vf[i]); eblock::face &res2f(res2.vf.back ()); res2f.np.clear();
			CGAL::Oriented_side prev_p_flag (fp_flag[i].back ());	//��¼ǰһ�������ıߣ�Ҳ�����ж��Ƿ���Ҫ���㽻��
			for (int j(0);j!=eb.vf[i].np.size ();++j)
			{
				eblock::face *resf = 0;
				switch (fp_flag[i][j])
				{	
				case CGAL::Oriented_side::ON_POSITIVE_SIDE:
					resf = &res1f;
				case CGAL::Oriented_side::ON_NEGATIVE_SIDE:
					if (resf == 0)
						resf = &res2f;
					if ((prev_p_flag != CGAL::Oriented_side::ON_ORIENTED_BOUNDARY) && (prev_p_flag !=fp_flag[i][j]))
					{
						int j_1 (j-1);	//ǰһ�����Ŀǰ�ĵ㲻��ͬһ��ƽ���ϣ��н�����������㽻��,
						if (j_1 < 0) j_1 = (int)eb.vf[i].np.size ()-1;
						//����һ������߶εĽ�����ǰ�����û��
						int cal_pi = find_n(intersected_seg,int_segment(min(eb.vf[i].np[j], eb.vf[i].np[j_1]), max(eb.vf[i].np[j], eb.vf[i].np[j_1])));
						int ints_p_index;
						if (cal_pi == intersected_seg.size ())
						{
							point_3 ints_p;		//����߶εĽ�����ǰû�м��������Ҫ���¼���
							if (intersection_segp(vpo[eb.vf[i].np[j]], vpo[eb.vf[i].np[j_1]],pl,ints_p) != 0)
								throw logic_error ("plane_cut_eblock, intersection point calculation error");
							ints_p_index = add_set(vpo,ints_p);
							//����һ��������������
							intersected_seg.push_back (int_segment(min(eb.vf[i].np[j], eb.vf[i].np[j_1]), max(eb.vf[i].np[j], eb.vf[i].np[j_1])));
							calculated_pi.push_back (ints_p_index);
						}
						else
							ints_p_index = calculated_pi[cal_pi];
						
						res1f.np.push_back (ints_p_index);
						res2f.np.push_back (ints_p_index);
						cutseg.push_back (ints_p_index);	//ins_p_index������ڶ�����
					}
					resf->np.push_back (eb.vf[i].np[j]); break;

				case CGAL::Oriented_side::ON_ORIENTED_BOUNDARY:
					res1f.np.push_back (eb.vf[i].np[j]);
					res2f.np.push_back (eb.vf[i].np[j]); 
					cutseg.push_back (eb.vf[i].np[j]);	 break;//������ڶ�����
				}
				prev_p_flag = fp_flag[i][j];	//����һ���������
			}
			if (cutseg.size () != 2)
				throw logic_error ("plane_cut_eblock, cutting plane error");	//�������һ������ֻ���������ڶ�����

			break;
		}
#ifdef S_DEBUG
		if (cutseg.size () >=3)
			throw logic_error ("plane_cut_eblock, more than three points on the plane");	//ע�⣬����cutseg�Ĵ�С������>=3����Ϊ������ǰ���Ѿ����ų�����.�Է���һ���Ǽ��һ��
#endif
		if(cutseg.size () != 2)
			continue;		//������ϵ��߶�û���ڶ����ϵġ�
		if (cutseg[0] == cutseg[1])
			throw logic_error ("plane_cut_eblock, same point");
		add_set(vseg, int_segment(min(cutseg[0],cutseg[1]), max(cutseg[0],cutseg[1])));	//����߶��ڶ�����
	}

	//���ɶ���
	if (vseg.size ()<2)
		throw logic_error ("plane_cut_eblock, section segment less than 3");
	
	vector<int> section;
	if (construct_loop(vseg,section) != 0)		//���ɵ������
		throw logic_error ("plane_cut_eblock, section generation fail");
	res1.vf.push_back (eblock::face()); res2.vf.push_back (eblock::face());		//�������
	res1.vf.back().frac_id = res2.vf.back ().frac_id = frac_id;
	res1.vf.back().plane_id = res2.vf.back ().plane_id = pli;
	res1.vf.back ().np.resize (section.size ());
	res2.vf.back ().np.resize (section.size ());
	res1.vf.back ().outer_normal = false;		//res1��pl�����࣬��ô�����洦���ⷨ�߷����Ӧ�ú�pl�������෴
	res2.vf.back ().outer_normal = true;		//res2��pl�ĸ��ࡣͬ��

	vector_3 normal (oriented_normal(vpo,section));
	FT dot_p (normal*pl.orthogonal_vector());
	if (dot_p > 0)
	{
		copy(section.begin (), section.end (), res2.vf.back ().np.begin ());		//normal��pl����Ӧ�üӵ�res2��,��Ϊ��ʱsection�е�����б�ʾ���ⷨ�߷���Ӧ����res2������ⷨ�ߣ���Ϊres2��pl�ĸ���
		reverse_copy(section.begin (),section.end (),res1.vf.back().np.begin ());	//�ѷ���֮������мӵ�res1��
	}
	else if (dot_p<0)
	{
		copy(section.begin (), section.end (), res1.vf.back ().np.begin ());		//normal��pl�������ñ�ʾres1���ⷨ�߷���
		reverse_copy(section.begin (),section.end (),res2.vf.back ().np.begin ());	//ͬ��
	}
	else
		throw logic_error ("plane_cut_eblock, normal and pl parallel");	//���ʱ����ô���ܳ���ƽ�У�����Ц��������������

	veb[targeti].swap(res1);
	veb.push_back (eblock()); veb.back ().swap (res2);	//��������ɵ�����������
#ifdef S_DEBUG
	if (!eb_sealed(targeti) || !eb_convex(targeti) || !eb_sealed((int)veb.size ()-1) || !eb_convex((int)veb.size () - 1))
		throw logic_error ("eblock generation fail");

#endif
	return 0;	//�ɹ�ʶ�� oh yeah
}

int block_system::shrk_plane (const int &pli)
{
	const plane_3 &pl = vpl[pli];
	vector<int> p_eb;	//��pli����ĵ�Ԫ����
	vector<int> p_ebf;	//����veb[p_eb[i]]�ĸ����pli���غ�
	vector<int> n_eb;	//��pli����ĵ�Ԫ����
	vector<int> n_ebf;	//����veb[n_eb[i]]�ĸ����pli���غ�
	p_eb.reserve (10); n_eb.reserve(10); p_ebf.reserve (10); n_ebf.reserve(10);	//10�������������õ�
	for (int i(0);i!=veb.size ();++i)
	{
		int j(0);
		for (;j!=veb[i].vf.size ();++j)
			if (veb[i].vf[j].plane_id == pli && veb[i].vf[j].frac_id>-4)
				break;	//�����ĳ��������pli��������ɵģ�������
		if (j!=veb[i].vf.size ())
		{
			if (veb[i].vf[j].outer_normal)
			{n_eb.push_back (i); n_ebf.push_back (j);}	//��ʱ�������������ⷨ�߷�������pli�����࣬��˿�������pli�ĸ���
			else
			{p_eb.push_back (i); p_ebf.push_back (j);}	//��ʱ�������������ⷨ�߷�������pli�ĸ��࣬��˿�������pli������
		}
	}
	if(p_eb.size ()==0 || n_eb.size ()==0)
		return 0;		//����eblock��pli��ͬһ�࣬����Ҫ�ϲ�����
	
	//����ƽ���ϵ���϶
	vector<int> discf_i;		//Բ����϶����
	vector<int> polyf_i;		//�������϶����
	//����pli������Щ��϶
	for (int i(0);i!=vdisc.size ();++i)
		if (vdisc[i].plane_id == pli)
			discf_i.push_back (i);
	for (int i(0);i!=vpoly.size ();++i)
		if (vpoly[i].plane_id == pli)
			polyf_i.push_back (i);

	polygon_set_2 cutter;	//�������ƽ���ϵ���϶�Ĳ�������Ϊһ����Ƭ
	polygon_2 temp_d;
	for (int i(0);i!=discf_i.size ();++i)
	{
		if (proj_polygon(pl,vdisc[discf_i[i]],temp_d) != 0)
			throw logic_error ("shrk_plane, constructing disc cutter fail");
#ifdef S_DEBUG
		//ofstream outfile ("cutter.txt");//debug
		//for (polygon_2::Edge_const_iterator i = temp.edges_begin(); i!= temp.edges_end();++i)
		//	out_line_2(outfile,i->source(),i->end ());
		//outfile.close();
#endif
		cutter.join(temp_d);
	}
	polygon_with_holes_2 temp_p;
	for (int i(0);i!=polyf_i.size ();++i)
	{
		int return_res = proj_polygon(pl,vpoly[polyf_i[i]],temp_p);
		if ( return_res!= 0)
		{
			cout<<"error_res:"<<return_res<<endl;
			ofstream outfile ("debug_info.txt");
			out_line_3_rhino(outfile,vpoly[polyf_i[i]]);
			outfile.close();
			throw logic_error ("shrk_plane, constructing poly cutter fail");
		}
		cutter.join(temp_p);
	}
	
	vector<polygon_2> p_poly, n_poly;	//�ֱ��ʾ������ĺ��ڸ���Ŀ�����pl�ϵĶ����
	vector<mbox_2> p_box, n_box;		//�ֱ��ʾp_poly��n_poly�еĶ���ε�mbb
	p_poly.reserve(p_eb.size ()); p_box.reserve(p_eb.size ()); n_poly.reserve (n_eb.size ()); n_box.reserve(n_eb.size ());
	for (int i(0);i!=p_eb.size ();++i)
	{
		p_poly.push_back (polygon_2());		//����p_eb�еĿ�����pl�ϵĶ����
		proj_polygon(pl,vpo,veb[p_eb[i]].vf[p_ebf[i]].np,p_poly.back ());
#ifdef S_DEBUG
		if (!p_poly.back ().is_simple() || p_poly.back ().is_collinear_oriented() || !p_poly.back ().is_convex())
			throw logic_error ("shrk_plane, constructing eb polygon fail");

		//ofstream outfile ("p_eb.txt");	//debug
		//for (polygon_2::Edge_const_iterator i = p_poly.back ().edges_begin(); i!= p_poly.back ().edges_end();++i)
		//	out_line_2(outfile,i->source(),i->end ());
		//outfile.close();

#endif
		if (p_poly.back ().is_clockwise_oriented())
			p_poly.back ().reverse_orientation();
		p_box.push_back (build_mbox_2(p_poly.back ()));
	}
	for (int i(0);i!=n_eb.size();++i)
	{
		n_poly.push_back (polygon_2());		//����n_eb�еĿ�����pl�ϵĶ����
		proj_polygon(pl,vpo,veb[n_eb[i]].vf[n_ebf[i]].np,n_poly.back ());
#ifdef S_DEBUG
		if (!n_poly.back ().is_simple() || n_poly.back ().is_collinear_oriented() || !n_poly.back ().is_convex())
			throw logic_error ("shrk_plane, constructing eb polygon fail");

		//ofstream outfile ("n_eb.txt");	//debug
		//for (polygon_2::Edge_const_iterator i = n_poly.back ().edges_begin(); i!= n_poly.back ().edges_end();++i)
		//	out_line_2(outfile,i->source(),i->end ());
		//outfile.close();
#endif
		if (n_poly.back ().is_clockwise_oriented())
			n_poly.back ().reverse_orientation();
		n_box.push_back (build_mbox_2(n_poly.back ()));
	}

	int counter(0);		//������һ���ϲ��˼�������
	for (int i(0);i!=p_eb.size ();++i)
	{
		for (int j(0);j!=n_eb.size ();++j)
		{
			if (veb[p_eb[i]].cblock_index == veb[n_eb[j]].cblock_index) continue;	//������eblock�Ѿ�����ͬһ�����Ͽ����ˣ��Ͳ����ٺ�����
			if (!CGAL::do_overlap(p_box[i],n_box[j])) continue;		//������eblock��Ӧ�Ķ�������߽�ж�����ȣ������Ͳ������ཻ��ֱ������
			polygon_set_2 temp_ps;
			temp_ps.insert(p_poly[i]);
			temp_ps.intersection(n_poly[j]);
			temp_ps.difference(cutter);
			if (temp_ps.is_empty ()) continue;	//������eblock��cutter�п�
			++counter;
			int old_index (veb[n_eb[j]].cblock_index), new_index (veb[p_eb[i]].cblock_index);
			for (int k(0);k!=veb.size ();++k)
				if (veb[k].cblock_index == old_index)
					veb[k].cblock_index  = new_index;
		}
	}
	return counter;
}

int block_system::identify_block (std::ostream & os)
{
	CGAL::Timer timer;
#ifdef S_DEBUG
	cout<<"identifying block"<<endl;
#endif
	vector<mbox_3> veb_box, vdf_box, vpf_box;	//�ȹ�������eblock����϶�ı߽��
	veb_box.reserve (veb.size ()); vdf_box.reserve (vdisc.size ()); vpf_box.reserve (vpoly.size ());
	for (int i(0);i!=veb.size ();++i)
		veb_box.push_back (build_mbox_3(vpo,veb[i]));
	for (int i(0);i!=vdisc.size ();++i)
		vdf_box.push_back (build_mbox_3(vdisc[i].center,vpl[vdisc[i].plane_id],vdisc[i].r));
	for (int i(0);i!=vpoly.size ();++i)
		vpf_box.push_back (build_mbox_3(vpoly[i]));

	for (int i(0);i!=veb.size ();++i)
		veb[i].cblock_index = i;	//��ʼ��һ��cblock_index
#ifdef S_DEBUG
	cout<<"cutting blocks"<<endl;
#endif
	timer.start();
	int max_cblock_index((int)veb.size ());		//��¼����eblock�У�cblock_index�����ֵ.��˼�ǣ�����eblock��cblock_indexλ�ö�С��max_cblock_index
	for (int i(0);i!=vpoly.size ();++i)
	{
#ifdef S_DEBUG
		cout<<"p:"<<i<<"/"<<vpoly.size ()<<endl;
#endif
		int veb_count = (int) veb.size ();		//�ȼ�¼һ��eblock������,��Ϊveb�п��ܻ�����µĿ���
		for (int j=0; j!=veb_count; ++j)
			if (CGAL::do_overlap(vpf_box[i],veb_box[j]) &&
				(plane_cut_eblock(j,vpoly[i].plane_id,vpoly[i].frac_id,max_cblock_index) == 0))
			{
				veb_box[j] = build_mbox_3(vpo,veb[j]); 
				veb_box.push_back (build_mbox_3(vpo,veb.back ()));
				max_cblock_index += 2;
			}
	}
	for (int i(0);i!=vdisc.size ();++i)
	{
#ifdef S_DEBUG
		cout<<"d:"<<i<<"/"<<vdisc.size ()<<endl;
#endif
		int veb_count = (int) veb.size ();		//�ȼ�¼һ��eblock������,��Ϊveb�п��ܻ�����µĿ���
		for (int j=0; j!=veb_count; ++j)
			if (CGAL::do_overlap(vdf_box[i],veb_box[j]) &&
				(plane_cut_eblock(j,vdisc[i].plane_id,vdisc[i].frac_id,max_cblock_index) == 0))
			{
				veb_box[j] = build_mbox_3(vpo,veb[j]); 
				veb_box.push_back (build_mbox_3(vpo,veb.back ()));
				max_cblock_index += 2;
			}
	}
	timer.stop();
	cutting_time = timer.time ();
	
#ifdef S_DEBUG
	cout<<"eblock generating complete. Time:"<<timer.time()<<endl;
	cout<<"combining blocks"<<endl;
	//for (int i(0);i!=vpo.size ();++i)
	//	cout<<vpo[i]<<endl;
#endif
	timer.reset(); timer.start();
	for (int i(0);i!=vpl.size ();++i)
	{
		int count = shrk_plane (i);
#ifdef S_DEBUG
		cout<<i<<"/"<<vpl.size ()<<' '<<count<<endl;
		//cout<<vpl[i]<<endl;
#endif
	}

	vector<int> cb_n;	//��¼һ���ж��ٸ����Ͽ���
	for (int i(0);i!=veb.size ();++i)
		add_set(cb_n,veb[i].cblock_index);
	vcb.clear();
	vcb.resize (cb_n.size ());
	for (int i(0);i!=cb_n.size ();++i)
		for (int j(0);j!=veb.size ();++j)
			if (veb[j].cblock_index == cb_n[i])
				vcb[i].vb.push_back (j);
	for (int i(0);i!=veb.size ();++i)
	{
		int j(0);
		for (;j!=cb_n.size ();++j)
			if (cb_n[j] == veb[i].cblock_index)
				break;
		veb[i].cblock_index = j;
	}
	timer.stop();

	shrk_time = timer.time ();
#ifdef S_DEBUG
	cout<<"combine eblocks complete Time:"<<timer.time ()<<endl;
#endif

	return 0;
}

int construct_loop (const std::vector<std::pair<int,int>> &vis, std::vector<int> &res)
{
	if (vis.size () < 3)
		return 1;
	res.clear();
	vector<bool> s_flag(vis.size (),false);	//���vis�е�pair�Ƿ��Ѿ����õ�
	res.push_back (vis[0].first);	res.push_back (vis[0].second);
	s_flag[0] = true;	//���pair�Ѿ�������
	int flag_count(1);	//��¼�ж��ٸ�pair�Ѿ�������

	while (flag_count<vis.size ())
	{
		int i(1);
		for (;i!=vis.size ();++i)
		{
			if (s_flag[i])
				continue;
			if (vis[i].first == res.back())
			{
				if (vis[i].second == res.front())
					break;
				res.push_back (vis[i].second);
				s_flag[i] = true; ++flag_count;
				
			}
			else if (vis[i].second == res.back())
			{
				if (vis[i].first == res.front())
					break;
				res.push_back(vis[i].first);
				s_flag[i] = true; ++flag_count;
			}
		}
		if (i != vis.size ())
			break;
	}
	if (flag_count != vis.size ()-1)
		return 2;		//���е�pair��������˻���û���������loop ��flag_count == vis.size()��,�������е�pairû��������Ѿ����������loop
	return 0;
}

int block_system::output_eb_stl (const int &ebi,std::ostream &outfile, const FT &x_, const FT &y_, const FT &z_) const
{
	if (!outfile)
		return 1;
	const eblock &eb(veb[ebi]);
	outfile<<"solid eblock"<<'\n';
	for (int j(0);j!=eb.vf.size ();++j)
	{
		if (eb.vf[j].np.size () <3)
			continue;
		vector_3 n(vpl[eb.vf[j].plane_id].orthogonal_vector());
		if (!eb.vf[j].outer_normal)
			n = -n;
		for (int k(1);k < (int)eb.vf[j].np.size ()-1;++k)
		{
			outfile<<"  facet normal "<<n.x()<<' '<<n.y()<<' '<<n.z()<<'\n';
			outfile<<"    outer loop\n";
			outfile<<"      vertex "<<vpo[eb.vf[j].np[0  ]].x() + x_<<' '<<vpo[eb.vf[j].np[0  ]].y() + y_<<' '<<vpo[eb.vf[j].np[0  ]].z() + z_<<'\n';
			outfile<<"      vertex "<<vpo[eb.vf[j].np[k  ]].x() + x_<<' '<<vpo[eb.vf[j].np[k  ]].y() + y_<<' '<<vpo[eb.vf[j].np[k  ]].z() + z_<<'\n';
			outfile<<"      vertex "<<vpo[eb.vf[j].np[k+1]].x() + x_<<' '<<vpo[eb.vf[j].np[k+1]].y() + y_<<' '<<vpo[eb.vf[j].np[k+1]].z() + z_<<'\n';
			outfile<<"    endloop"<<endl;
			outfile<<"  endfacet"<<endl;
		}
	}
	outfile<<"endsolid eblock"<<endl;
	return 0;
}

int block_system::output_cb_stl (const int &cbi, std::ostream &outfile, const FT &x_, const FT &y_, const FT &z_, const string &cb_name) const 
{
	if (!outfile)
		return 1;
	if (cb_name.empty ())
		outfile<<"solid cblock"<<'\n';
	else
		outfile<<"solid "<<cb_name<<'\n';
	for (int i(0);i!=vcb[cbi].vb.size ();++i)
	{
		const eblock &eb(veb[vcb[cbi].vb[i]]);
		for (int j(0);j!=eb.vf.size ();++j)
		{
			if (eb.vf[j].np.size () <3)
				continue;
			vector_3 n(vpl[eb.vf[j].plane_id].orthogonal_vector());
			if (!eb.vf[j].outer_normal)
				n = -n;
			for (int k(1);k < (int)eb.vf[j].np.size ()-1;++k)
			{
				outfile<<"  facet normal "<<n.x()<<' '<<n.y()<<' '<<n.z()<<'\n';
				outfile<<"    outer loop\n";
				outfile<<"      vertex "<<vpo[eb.vf[j].np[0  ]].x() + x_<<' '<<vpo[eb.vf[j].np[0  ]].y() + y_<<' '<<vpo[eb.vf[j].np[0  ]].z() + z_<<'\n';
				outfile<<"      vertex "<<vpo[eb.vf[j].np[k  ]].x() + x_<<' '<<vpo[eb.vf[j].np[k  ]].y() + y_<<' '<<vpo[eb.vf[j].np[k  ]].z() + z_<<'\n';
				outfile<<"      vertex "<<vpo[eb.vf[j].np[k+1]].x() + x_<<' '<<vpo[eb.vf[j].np[k+1]].y() + y_<<' '<<vpo[eb.vf[j].np[k+1]].z() + z_<<'\n';
				outfile<<"    endloop"<<endl;
				outfile<<"  endfacet"<<endl;
			}
		}
	}
	if (cb_name.empty ())
		outfile<<"endsolid cblock"<<endl;
	else
		outfile<<"endsolid "<<cb_name<<endl;

	return 0;
}

int block_system::output_explode_cb_stl (const int &cbi, std::ostream &outfile, const double &ratio) const 
{
	if (!outfile)
		return 1;
	const cblock &cb = vcb[cbi];
	vector<point_3> vp_eb;
	vp_eb.reserve (cb.vb.size ());

	for (int i(0);i!=cb.vb.size ();++i)
		vp_eb.push_back (veb[cb.vb[i]].get_rand_point(vpo));

	double x(0),y(0),z(0);
	for (int i(0);i!=vp_eb.size ();++i)
	{
		x = x + CGAL::to_double (vp_eb[i].x());
		y = y + CGAL::to_double (vp_eb[i].y());
		z = z + CGAL::to_double (vp_eb[i].z());
	}
	x = x / double(vp_eb.size ());
	y = y / double(vp_eb.size ());
	z = z / double(vp_eb.size ());

	for (int i(0);i!=cb.vb.size ();++i)
	{

		output_eb_stl(cb.vb[i],outfile, (vp_eb[i].x() - x)*ratio, (vp_eb[i].y() - y)*ratio, (vp_eb[i].z() - z)*ratio);
	}
	
	return 0;
}

point_3 block_system::get_rand_p() const
{
	double x(0),y(0),z(0);
	for (int i(0);i!=vpo.size ();++i)
	{
		x += CGAL::to_double(vpo[i].x());
		y += CGAL::to_double(vpo[i].y());
		z += CGAL::to_double(vpo[i].z());
	}
	return point_3 (x/double(vpo.size ()), y/double (vpo.size ()), z/double (vpo.size ()));
}

void block_system::output_info (std::ostream &outfile) const
{
	if (!outfile)
		return;
	outfile<<"***Block system info.***"<<endl;
	outfile<<"Point count (vpo.size()):\t"<<vpo.size ()<<endl;
	outfile<<"Plane count (vpl.size()):\t"<<vpl.size ()<<endl;
	outfile<<"Disc frac count (vdisc.size()):\t"<<vdisc.size ()<<endl;
	outfile<<"Poly frac count (vpoly.size()):\t"<<vpoly.size ()<<endl;
	outfile<<"Element block count (veb.size()):\t"<<veb.size ()<<endl;
	outfile<<"Complex block count (vcb.siz ()):\t"<<vcb.size ()<<endl;
	outfile<<"Cutting time:\t"<<cutting_time<<endl;
	outfile<<"Shrinking time:\t"<<shrk_time<<endl;
	outfile<<"Bound surface count:\t"<<vbsf.size ()<<endl;
}


//------------code for stability analysis---------------
int block_system::init_stabl_alys (const bool& debug)
{
	vcb_bsf.resize (vcb.size ());
	vpl_bsf.resize (vpl.size ());
	vbsf.reserve (vcb.size()*10);

	//����˼·���ǰ�ÿ�����Ͽ����������ҳ���
	for (int i(0);i!=vcb.size ();++i)
	{
		if (debug)
			cout<<"cb:"<<i<<"/"<<vcb.size ()<<"\n";
		vector<int> plane_index;	//��¼������Ͽ����漰����Щ��
		vector<polygon_set_2> p_poly;	//����Ϊplane_index[i]��ƽ������Ķ����
		vector<polygon_set_2> n_poly;	//����Ϊplane_index[i]��ƽ�渺��Ķ����

		for (int j(0);j!=vcb[i].vb.size ();++j)
		{
			const eblock& eb = veb[vcb[i].vb[j]];
			for (int k(0);k!=eb.vf.size ();++k)
			{
				int index = add_set (plane_index, eb.vf[k].plane_id);
				if (plane_index.size () > p_poly.size ())
				{
					//������plane_index��������µ�ƽ��
					p_poly.push_back (polygon_set_2());
					n_poly.push_back (polygon_set_2());
				}
				polygon_2 poly_temp;
				if (proj_polygon(vpl[eb.vf[k].plane_id],vpo,eb.vf[k].np,poly_temp) != 0)
					throw logic_error ("init_stabl_alys, projection polygon fail");
				if (poly_temp.is_clockwise_oriented ())
					poly_temp.reverse_orientation();
				if (eb.vf [k].outer_normal)
					n_poly[index].join (poly_temp);
				else 
					p_poly[index].join (poly_temp);
			}
		}
		for (int j(0);j!=plane_index.size ();++j)
		{
			polygon_set_2 temp_set (n_poly[j]);
			n_poly[j].difference (p_poly[j]);
			p_poly[j].difference (temp_set);

			if ((!(n_poly[j].is_valid())) || (!(p_poly[j].is_valid())))
				throw logic_error ("init_stabl_alys, resulting polygon is not valid");
			if (n_poly[j].is_empty () && p_poly[j].is_empty ())
				continue;	//����ǿյľ�ֱ������
			else
			{
				vector<polygon_with_holes_2> v_n_poly, v_p_poly;
				n_poly[j].polygons_with_holes(std::back_inserter (v_n_poly));
				p_poly[j].polygons_with_holes(std::back_inserter (v_p_poly));
				//û��simplify��split��������Σ������мӹ���ֱ�����
				for (int k(0);k!=v_n_poly.size ();++k)
				{
					//��Щ����ⷨ�߷�����ƽ��vpl[plane_index[j]]��������ͬ
					vbsf.push_back (bound_surface());
					vbsf.back ().poly = v_n_poly[k];
					vbsf.back ().pli = plane_index[j];
					vbsf.back ().cbi = i;
					vbsf.back ().outer_normal = true;
					vcb_bsf[i].push_back ((int)vbsf.size ()-1);
					vpl_bsf[plane_index[j]].push_back ((int)vbsf.size ()-1);
				}
				for (int k(0);k!=v_p_poly.size ();++k)
				{
					//��Щ����ⷨ�߷�����ƽ��vpl[plane_index[j]]�ķ�������ͬ
					vbsf.push_back (bound_surface());
					vbsf.back ().poly = v_p_poly[k];
					vbsf.back ().pli = plane_index[j];
					vbsf.back ().cbi = i;
					vbsf.back ().outer_normal = false;
					vcb_bsf[i].push_back ((int)vbsf.size ()-1);
					vpl_bsf[plane_index[j]].push_back ((int)vbsf.size ()-1);
				}
			}
		}
	}
	if (debug)
		cout<<"Construcing Adj.  Relation"<<endl;
	for (int i(0);i!=(int)vbsf.size ();++i)
		for (int j(i+1);j<(int) vbsf.size ();++j)
		{
			if (vbsf[i].pli != vbsf[j].pli)
				continue;
			if (vbsf[i].outer_normal == vbsf[j].outer_normal)
				continue;
			if (CGAL::do_intersect(vbsf[i].poly, vbsf[j].poly))
			{
				vbsf[i].vadj_bsf_i.push_back (j);
				vbsf[j].vadj_bsf_i.push_back (i);
				
				vbsf[i].vadj_cb_i.push_back (vbsf[j].cbi);
				vbsf[j].vadj_cb_i.push_back (vbsf[i].cbi);
			}
		}

	for (int i(0);i!=(int) vbsf.size ();++i)
		if (vbsf[i].vadj_bsf_i.empty())
			vbsf[i].type = FACE_TYPE::EXCAV;
		else
			vbsf[i].type = FACE_TYPE::FRAC;
	if (debug)
		cout<<"init_stabl_alys complete"<<endl;
	return 0;
}

//bool block_system::is_interlock (const int &cbi, const int &cbj) const
//{
//	vector<vector_3> v_fdir;		//��¼cbi��������������������cbj����������ڵ������
//	v_fdir.reserve (vcb_bsf[cbi].size());
//	for (int i(0);i!=vcb_bsf[cbi].size ();++i)
//	{
//		//����cbi��ÿ�������
//		int count(0);
//		for (;count!=(int)vbsf[vcb_bsf[cbi][i]].vadj_cb_i.size ();++count)
//		{
//			if (vbsf[vcb_bsf[cbi][i]].vadj_cb_i[count] == cbj)
//				break;
//		}
//		if (count != vbsf[vcb_bsf[cbi][i]].vadj_cb_i.size ())	
//		{
//			//����������cbj�����������
//			if (vbsf[vcb_bsf[cbi][i]].outer_normal)
//				v_fdir.push_back (-(vpl[vbsf[vcb_bsf[cbi][i]].pli].orthogonal_vector()));
//			else
//				v_fdir.push_back ((vpl[vbsf[vcb_bsf[cbi][i]].pli].orthogonal_vector()));
//		}
//	}
//	if (removability(v_fdir) == Removability_Status::Nonremovable)
//		return true;		//��ʱ������������������һ���
//	else
//		return false;		//��ʱ�����������ǿ��Էֿ��ģ�û������
//}
//
//bool block_system::is_interlock (const std::vector<int> &vcbi) const
//{
//	if (vcbi.size ()<=1)
//		return true;		//��������϶��Ǻ��Լ�������
//	
//
//}

bool block_system::is_removable (const int &cbi) const
{
	vector<vector_3> v_fdir;
	v_fdir.reserve (vcb_bsf[cbi].size ());
	for (int i(0);i!=(int)vcb_bsf[cbi].size ();++i)
	{
		if (vbsf[vcb_bsf[cbi][i]].type == FACE_TYPE::FRAC)
		{
			if (vbsf[vcb_bsf[cbi][i]].outer_normal)
				v_fdir.push_back (-(vpl[vbsf[vcb_bsf[cbi][i]].pli].orthogonal_vector()));
			else
				v_fdir.push_back (vpl[vbsf[vcb_bsf[cbi][i]].pli].orthogonal_vector());
		}
	}
	if (removability(v_fdir) == Removability_Status::Nonremovable)
		return false;
	else
		return true;
}

bool block_system::is_removable (const std::vector<int> &vcbi) const
{
	vector<vector_3> v_fdir;
	v_fdir.reserve (vcbi.size ()*5);
	for (int i(0);i!=vcbi.size ();++i)
	{
		for (int j(0);j!=vcb_bsf[vcbi[i]].size ();++j)
		{
			const bound_surface & bsf = vbsf[vcb_bsf[vcbi[i]][j]];
			if (bsf.type == FACE_TYPE::FRAC)
			{
				if (!is_contained(bsf.vadj_cb_i,vcbi))	//���������һ���Ӧ�Ŀ��岢û����ȫ������vcbi�У�������һ���޶���
					if (bsf.outer_normal)
						v_fdir.push_back (-(vpl[bsf.pli].orthogonal_vector()));
					else 
						v_fdir.push_back (vpl[bsf.pli].orthogonal_vector());
			}
		}
	}
	if (removability(v_fdir) == Removability_Status::Nonremovable)
		return false;
	else 
		return true;
}

bool block_system::is_removable (const std::vector<int> &vcbi, const std::vector <int> &vfcbi) const
{
	vector<vector_3> v_fdir;
	v_fdir.reserve (vcbi.size ()*5);
	for (int i(0);i!=vcbi.size ();++i)
	{
		for (int j(0);j!=vcb_bsf[vcbi[i]].size ();++j)
		{
			const bound_surface & bsf = vbsf[vcb_bsf[vcbi[i]][j]];
			if (bsf.type == FACE_TYPE::FRAC)
			{
				if ((!is_contained (bsf.vadj_cb_i,vcbi)) && (has_commonelem(bsf.vadj_cb_i,vfcbi)))
					if (bsf.outer_normal)
						v_fdir.push_back (-(vpl[bsf.pli].orthogonal_vector()));
					else 
						v_fdir.push_back (vpl[bsf.pli].orthogonal_vector());
			}
		}
	}
	if (removability(v_fdir) == Removability_Status::Nonremovable)
		return false;
	else 
		return true;
}

int block_system::test_output_STL (std::ostream&outfile, const int&index, const int &flag,bool single,const std::string &stl_name) const
{
	using CGAL::to_double;
	if (!outfile)
		return 1;
	outfile.precision (20);
	if (index>=0)
	{
		if(single)
		{
			if (stl_name.empty ())
				outfile<<"solid CB "<<index<<endl;
			else
				outfile<<"solid "<<stl_name<<endl;
		}
		for (int i(0);i!=vcb_bsf[index].size ();++i)
			output_STL_poly(outfile,vbsf[vcb_bsf[index][i]].poly,vpl[vbsf[vcb_bsf[index][i]].pli],
				vbsf[vcb_bsf[index][i]].outer_normal,0,0,0,false);
		if (single)
		{
			if (stl_name.empty ())
				outfile<<"endsolid CB "<<index<<endl;
			else
				outfile<<"endsolid "<<stl_name<<endl;
		}
	}
	else
	{
		if (stl_name.empty ())
			outfile<<"solid plane"<<endl;
		else 
			outfile<<"solid "<<stl_name<<endl;
		for (int i(0);i!=vpl_bsf[-index+1].size ();++i)
		{
			if (flag==1 && (vbsf[vpl_bsf[-index+1][i]].outer_normal))
				continue;
			if (flag==2 && !(vbsf[vpl_bsf[-index+1][i]].outer_normal))
				continue;
			output_STL_poly(outfile,vbsf[vpl_bsf[-index+1][i]].poly, vpl[vbsf[vpl_bsf[-index+1][i]].pli],
				vbsf[vpl_bsf[-index+1][i]].outer_normal,0,0,0 ,false);
		}
		if (stl_name.empty ())
			outfile<<"endsolid plane"<<endl;
		else 
			outfile<<"endsolid "<<stl_name<<endl;
	}
	return 0;
}

}