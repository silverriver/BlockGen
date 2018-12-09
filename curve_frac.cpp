#include "curve_frac.h"
#include <string>
#include <algorithm>

using namespace std;

namespace BI
{

void curve_frac::mark_domains (CDT &cdt)
{
	std::list<CDT::Face_handle> queue;
	queue.push_back (cdt.infinite_face ());
	while (!queue.empty ())
	{
		CDT::Face_handle fh = queue.front();
		queue.pop_front();
		if (fh->info().nesting_level == 0)
		{
			fh->info().nesting_level = 1;		//������Ϊ��Χ������
			for (int i(0);i!=3;++i)
			{
				CDT::Edge e(fh,i);
				CDT::Face_handle n = fh->neighbor(i);
				if (n->info ().nesting_level == 0)
					if (!cdt.is_constrained(e))
						queue.push_back (n);
			}
		}
	}
}

void curve_frac::add_constraint (CDT &cdt, const point_3 &p1, const point_3 &p2, const point_3 &p3, const plane_3 &pl)
{
	cdt.insert_constraint(pl.to_2d(p1), pl.to_2d(p2));
	cdt.insert_constraint(pl.to_2d(p2), pl.to_2d(p3));
	cdt.insert_constraint(pl.to_2d(p1), pl.to_2d(p3));
}

void curve_frac::add_constraint (CDT & cdt, const triangle_3 &t, const plane_3 &pl)
{
	add_constraint(cdt,t.vertex (0), t.vertex (1), t.vertex (2),pl);
}

void curve_frac::add_constraint (CDT &cdt, const std::vector<point_3> &vp, const tri_facet &f, const plane_3 &pl)
{
	add_constraint(cdt,vp[f.pn[0]], vp[f.pn[1]], vp[f.pn[2]],pl);
}

int curve_frac::read_STL (std::istream &infile,const FT & tag)
{
	if (!infile)
		return 1;
	vp.clear();
	vtri.clear();
	string keyword1,keyword2;
	char temp('a');
	while (temp!='\n')
		infile.get (temp);

	while (infile)
	{
		point_3 p[3];
		infile>>keyword1>>keyword2;
		std::transform(keyword1.begin (),keyword1.end (),keyword1.begin (),std::tolower);
		std::transform(keyword2.begin (),keyword2.end (),keyword2.begin (),std::tolower);
		if (keyword1 == "endsolid")
			break;
		if (keyword1 != "facet")
			return 2;		//�ļ���ʽ���Ϸ�
		getline(infile,keyword1);		//���Ե���������Ϣ
		infile>>keyword1>>keyword2;
		std::transform(keyword1.begin (),keyword1.end (),keyword1.begin (),std::tolower);
		std::transform(keyword2.begin (),keyword2.end (),keyword2.begin (),std::tolower);
		if (keyword1 != "outer" || keyword2 !="loop")
			return 2;
		for (int i(0);i!=3;++i)
		{
			infile>>keyword1;
			std::transform(keyword1.begin (),keyword1.end (),keyword1.begin (),std::tolower);
			if (keyword1 != "vertex")
				return 2;
			FT x,y,z;
			infile>>x>>y>>z;
			p[i] = point_3(x,y,z);
		}
		infile>>keyword1>>keyword2;
		std::transform(keyword1.begin (),keyword1.end (),keyword1.begin (),std::tolower);
		std::transform(keyword2.begin (),keyword2.end (),keyword2.begin (),std::tolower);
		if (keyword1 != "endloop" || keyword2 != "endfacet")
			return 2;
		int n1 = add_set (vp,p[0]);
		int n2 = add_set (vp,p[1]);
		int n3 = add_set (vp,p[2]);
		vtri.push_back (tri_facet(n1,n2,n3,tag));
	}
	return 0;
}

int curve_frac::output_STL (std::ostream &outfile)
{
	if (!outfile)
		return 1;
	outfile.precision (25);		//��������ȵ���һ��
	outfile<<"solid CurveFrac"<<endl;
	for (int i(0);i!=vtri.size ();++i)
	{
		vector_3 normal(plane_3(vp[vtri[i].pn[0]], vp[vtri[i].pn[1]], vp[vtri[i].pn[2]]).orthogonal_vector());
		double a (to_double (normal.x()));
		double b (to_double (normal.y()));
		double c (to_double (normal.z()));
		double d (sqrt(a*a+b*b+c*c));
		a/=d; b/=d; c/=d;
		outfile<<"  facet normal "<<a<<' '<<b<<' '<<c<<endl;
		outfile<<"    outer loop"<<endl;
		outfile<<"      vertex "<<vp[vtri[i].pn[0]].x()<<' '<<vp[vtri[i].pn[0]].y()<<' '<<vp[vtri[i].pn[0]].z()<<endl;
		outfile<<"      vertex "<<vp[vtri[i].pn[1]].x()<<' '<<vp[vtri[i].pn[1]].y()<<' '<<vp[vtri[i].pn[1]].z()<<endl;
		outfile<<"      vertex "<<vp[vtri[i].pn[2]].x()<<' '<<vp[vtri[i].pn[2]].y()<<' '<<vp[vtri[i].pn[2]].z()<<endl;
		outfile<<"    endloop"<<endl;
		outfile<<"  endfacet"<<endl;
	}
	outfile<<"endsolid CurveFrac"<<endl;
	return 0;
}

void curve_frac::insert_tri (vector <tri_facet> &new_vtri, const triangle_3 &t, const FT & tag)
{
	insert_tri(new_vtri,t.vertex (0),t.vertex (1), t.vertex (2),tag);
}

void curve_frac::insert_tri (std::vector <tri_facet> &new_vtri, const point_3 &p1, const point_3 &p2, const point_3 &p3, const FT & tag)
{
	new_vtri.push_back (tri_facet(add_set(vp,p1), add_set(vp,p2), add_set(vp,p3), tag));
}

void curve_frac::insert_tri (std::vector <tri_facet> &new_vtri, const CDT &cdt, const FT& tag, const plane_3 &pl)
{
	for (CDT::Finite_faces_iterator fi = cdt.finite_faces_begin(); fi!=cdt.finite_faces_end();++fi)
		if (fi->info().in_domain ())
			insert_tri(new_vtri, pl.to_3d(fi->vertex (0)->point ()), pl.to_3d(fi->vertex (1)->point ()), pl.to_3d(fi->vertex (2)->point ()),tag);
}

void curve_frac::join_tri (std::vector <tri_facet> &new_vtri, const triangle_3 &t, const FT &tag)
{
	join_tri(new_vtri, t.vertex(0),t.vertex (1), t.vertex (2), tag);
}

void curve_frac::join_tri (std::vector <tri_facet> &new_vtri, const point_3 &p1, const point_3 &p2, const point_3 &p3, const FT & tag)
{
	add_set(new_vtri,tri_facet(add_set(vp,p1), add_set(vp,p2), add_set(vp,p3), tag));
}

void curve_frac::join_tri (std::vector <tri_facet> &new_vtri, const CDT &cdt, const FT& tag,const plane_3 &pl)
{
	for (CDT::Finite_faces_iterator fi = cdt.finite_faces_begin(); fi!=cdt.finite_faces_end();++fi)
		if (fi->info ().in_domain ())
			join_tri(new_vtri, pl.to_3d(fi->vertex (0)->point ()), pl.to_3d(fi->vertex (1)->point ()), pl.to_3d(fi->vertex (2)->point ()),tag);
}

int curve_frac::join (const curve_frac & c)
{
	if (c.vp.empty() || c.vtri.empty ())
		return 1;
	vector<mbox_3> vmb1;		//��ʼ��һ����С��Χ��
	vector<mbox_3> vmb2;
	vmb1.reserve (vtri.size ());
	vmb2.reserve (c.vtri.size ());
	vector<int> temp_vp(3);
	for (int i(0);i!=vtri.size ();++i)
	{
		temp_vp[0] = vtri[i].pn[0]; temp_vp[1] = vtri[i].pn[1]; temp_vp[2] = vtri[i].pn[2];
		vmb1.push_back (build_mbox_3(vp,temp_vp));
	}
	for (int i(0);i!=c.vtri.size ();++i)
	{
		temp_vp[0] = c.vtri[i].pn[0]; temp_vp[1] = c.vtri[i].pn[1]; temp_vp[2] = c.vtri [i].pn[2];
		vmb2.push_back (build_mbox_3(c.vp,temp_vp));
	}

	vector<intersect_type> vinter;	//��¼�ཻ��Ԫ��
	vinter.reserve (vtri.size ()+c.vtri.size ());
	vector<vector<int>> inter_index1 (vtri.size ());		//vtri��ÿ�������ζ�Ӧһ��int�����飬ָʾ����������м����ཻԪ��
	vector<vector<int>> inter_index2 (c.vtri.size ());		//c.vtri��ÿ�������ζ�Ӧһ��int�����飬ָʾ����������м����ཻԪ��

	vector<int> vtag(vtri.size () + c.vtri.size ());		//��vtri��c.vtri�е�ÿ�������ζ�����һ����ǩ�������Щ���������ص��Ĳ��֡�����������������ص��Ĳ��֣���ô�����������εı�ǩ��ͬ
	for (int i(0);i!=vtri.size () + c.vtri.size ();++i)		//��ʼ��
		vtag[i] = i;

	//�����ཻԪ��
	for (int i(0);i!=vtri.size ();++i)
	{
		for (int j(0); j!=c.vtri.size ();++j)
		{
			if (!CGAL::do_overlap(vmb1[i],vmb2[j]))
				continue;
			triangle_3 t1(vp[vtri[i].pn[0]], vp[vtri[i].pn[1]], vp[vtri[i].pn[2]]);
			triangle_3 t2(c.vp[c.vtri[j].pn[0]], c.vp[c.vtri[j].pn[1]], c.vp[c.vtri[j].pn[2]]);
			
			vinter.push_back (CGAL::intersection(t1,t2));
			if (vinter.back ())
			{
				inter_index1[i].push_back ((int)vinter.size ()-1);		//������������ȷʵ�ཻ��
				inter_index2[j].push_back ((int)vinter.size ()-1);
				if (is_overlap(vinter.back ()))		//�����������������غϵĲ���
				{
					if (vtag[i] != vtag[j+vtri.size ()])
					{
						//ֻ��������tag����ͬ��ʱ���ִ����һ��������ʲôҲ����ִ��
						int old_tag (vtag[i]), new_tag(vtag[j+vtri.size ()]);
						for (int k(0);k!=vtag.size ();++k)
							if (vtag[k] == old_tag) vtag[k] = new_tag;
					}

				}
			}
			else
				vinter.pop_back();		//������������û���ཻ
		}
	}

	//��������Щ�������ص�
	vector<int> vdif_tag;					//��¼һ��һ���ж��ٸ���ͬ��tag
	vector<vector<int>> vtri_group;			//��¼��ÿһ����ͬ��tag��Ӧ��Щ������vtri_group[i][j]��¼����vdif_tag[i]���tag��Ӧ�ĵ�j�������ε�����
	vdif_tag.reserve (vtag.size ());
	vtri_group.reserve (vtag.size ());
	for (int i(0);i!=vtag.size ();++i)
	{
		int index = add_set (vdif_tag,vtag[i]);
		if (vdif_tag.size () > vtri_group.size ())
			vtri_group.push_back (vector<int>());
		vtri_group[index].push_back (i);		//��this��c�е������ε�����ͳһ��ʾ��c�������ε���������this����
	}

	vector<tri_facet> new_vtri;			//�µ�����������
	new_vtri.reserve (vtri.size ()+c.vtri.size ()+vinter.size ()*2);

	for (int i(0);i!=vtri_group.size ();++i)
	{
		//�ȴ���һЩ�������
		if (vtri_group[i].size () == 1)		//������漰�������������
		{
			if(vtri_group[i][0] < vtri.size ())		//��������������vtri��
			{
				if (inter_index1[vtri_group[i][0]].empty())		//�����������β���c���κ��������ཻ����ֱ����ӡ���Ӧ������������
				{
					new_vtri.push_back (vtri[vtri_group[i][0]]);
					continue;
				}
			}
			else 
			{
				if (inter_index2[vtri_group[i][0]-vtri.size ()].empty ())	//�����������c�ϣ����Ҳ���this�е��κ��������ཻ����ֱ����ӡ�
				{
					int index = vtri_group[i][0] - (int)vtri.size ();
					insert_tri(new_vtri,c.vp[c.vtri[index].pn[0]], c.vp[c.vtri[index].pn[1]], c.vp[c.vtri[index].pn[2]], c.vtri[index].tag);
					continue;		
				}
			}
		}
		
		plane_3 pl;		//������Щ���������ڵ�һ��ƽ�淽��
		FT tag;			//��¼��һ�������ε�tag�����沽���������ɵ������ε�tag�����ó����ֵ
		if (vtri_group[i][0] < vtri.size ())
		{
			pl = plane_3(vp[vtri[vtri_group[i][0]].pn[0]], vp[vtri[vtri_group[i][0]].pn[1]], vp[vtri[vtri_group[i][0]].pn[2]]);
			tag = vtri[vtri_group[i][0]].tag ;
		}
		else
		{
			int index = vtri_group[i][0] - (int)vtri.size ();
			pl = plane_3(c.vp[c.vtri [index].pn[0]], c.vp[c.vtri [index].pn[1]], c.vp[c.vtri [index].pn[2]]);
			tag = c.vtri[index].tag;
		}
		CDT cdt;
		for (int j(0);j!=vtri_group[i].size ();++j)		//���ڹ����ÿ��������
		{
			if (vtri_group[i][j] < vtri.size ())		//��������������this��
			{
				int index = vtri_group[i][j];
				add_constraint(cdt, vp, vtri[index], pl);
				add_constraint (cdt,vinter,inter_index1[index],pl);
			}
			else										//��������������c��
			{
				int index = vtri_group[i][j] - (int)vtri.size ();
				add_constraint (cdt,c.vp, c.vtri[index], pl);
				add_constraint (cdt,vinter, inter_index2[index],pl);
			}
		}
		if (vtri_group[i].size () > 1)
			mark_domains(cdt);		//��������ɸ��������غϵ��������Ҫ���һ����Щ����������Χ������
		insert_tri(new_vtri,cdt,tag,pl);
	}
	//for (int i(0);i!=vtri.size ();++i)
	//{
	//	if (inter_index1[i].empty())	
	//	{
	//		new_vtri.push_back (vtri[i]);
	//		continue;				//��⣬������������û�к�c�е�����һ���������ཻ����ֱ����������Ϊ��������ξ������curve_frac�У����ԾͲ�������ˡ���Ӧ���Ǵ󲿷ֵ����
	//	}
	//	CDT cdt;
	//	plane_3 pl (vp[vtri[i].pn[0]], vp[vtri[i].pn[1]], vp[vtri[i].pn[2]]);
	//	add_constraint(cdt,vp[vtri[i].pn[0]], vp[vtri[i].pn[1]], vp[vtri[i].pn[2]],pl);
	//	if (add_constraint(cdt,vinter,inter_index1[i],pl))	
	//		join_tri(new_vtri,cdt,vtri[i].tag ,pl);			//�����������κ���������ι���	,��������ε�ʱ����Ҫ����һ�������
	//	else
	//		insert_tri (new_vtri,cdt,vtri[i].tag ,pl);		//û�����������κ���������ι��棬��������ε�ʱ��ֱ�����
	//}
	//
	//for (int i(0);i!=c.vtri.size ();++i)
	//{
	//	if (inter_index2[i].empty())
	//	{
	//		insert_tri(new_vtri, c.vp[c.vtri[i].pn[0]], c.vp[c.vtri[i].pn[1]], c.vp[c.vtri[i].pn[2]], c.vtri[i].tag);	//����һ��ѭ����ͬ���ǣ���Ȼû�������κ�����������ཻ��������Ȼ��Ҫ�������������ӵ���������У���Ϊ��ԭ���������������
	//		continue;
	//	}
	//	CDT cdt;
	//	plane_3 pl (c.vp[c.vtri [i].pn[0]],c.vp[c.vtri [i].pn[1]],c.vp[c.vtri [i].pn[2]]);
	//	add_constraint (cdt,c.vp[c.vtri [i].pn[0]],c.vp[c.vtri [i].pn[1]],c.vp[c.vtri [i].pn[2]],pl);
	//	if (add_constraint(cdt,vinter,inter_index2[i],pl))
	//		join_tri(new_vtri,cdt,c.vtri[i].tag,pl);		//�����������κ���������ι���	,��������ε�ʱ����Ҫ����һ�������
	//	else
	//		insert_tri(new_vtri,cdt, c.vtri[i].tag,pl);	//û�����������κ���������ι��棬��������ε�ʱ��ֱ�����
	//}
	vtri.swap (new_vtri);
	return 0;		
	//�����Ͼ����������ཻ��������֮ǰʵ�ֵ�ʱ��Ҫ����̫���ˣ���Ϊ�����˺ܶ�߼��ĺ���
	//������һ�����ǺܺõĹ��ܣ�����������ɸ���������ͬһƽ���ཻ��������Щ�����εĲ������пն�����ô�ó����޷�����
	//��ʵ���н�������ģ��������������εĲ�����Ȼ��������cdt�аѲ������䲢���Ĳ����ų���ȥ������������̫�鷳��
	//�����������ص�����������������˻�������Ƚ���������������ʱ�Ȳ�������������˰ɡ�
}

bool curve_frac::add_constraint (CDT &cdt, const intersect_type &inter,const plane_3 &pl)
{
	if (const point_3* p = boost::get<point_3> (&*inter))		//���������εĽ�����һ��������
	{
		cdt.insert(pl.to_2d(*p));
		return false;
	}
	if (const segment_3* s = boost::get<segment_3> (&*inter))	//���������εĽ���һ���߶�
	{
		cdt.insert_constraint(pl.to_2d(s->start()), pl.to_2d(s->end()));
		return false;
	}
	if (const triangle_3 *t = boost::get<triangle_3> (&*inter))
	{
		add_constraint(cdt,*t,pl);			//���������εĽ�����һ��������
		return true;
	}
	const vector<point_3> *poly = boost::get<vector<point_3>>(&*inter);	//�ų������п��ܣ����������εĽ���ֻ������һ�������
	for (int k(0);k!=(int)poly->size ();++k)
	{
		int k1(k+1);
		if (k1 >=poly->size ()) k1 = 0;
		cdt.insert_constraint(pl.to_2d((*poly)[k]), pl.to_2d((*poly)[k1]));
	}
	return true;
}

bool curve_frac::add_constraint (CDT &cdt, const std::vector <intersect_type> &vinter, const plane_3 &pl)
{
	bool res (false);
	for (int i(0);i!=vinter.size ();++i)
		res = (add_constraint(cdt,vinter[i],pl) || res);
	return res;
}

bool curve_frac::add_constraint (CDT &cdt, const std::vector <intersect_type> &vinter, const std::vector <int> &vindex, const plane_3 &pl)
{
	bool res(false);
	for (int i(0);i!=vindex.size ();++i)
		res = (add_constraint(cdt,vinter[vindex[i]],pl) || res);
	return res;
}

bool curve_frac::is_overlap (const intersect_type &it) const
{
	if (!it)
		return false;
	if (const point_3 *p = boost::get<point_3> (&*it))
		return false;
	if (const segment_3 *p = boost::get<segment_3> (&*it))
		return false;
	return true;
}
}