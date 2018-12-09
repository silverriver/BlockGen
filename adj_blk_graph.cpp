#include "adj_blk_graph.h"
#include "geometry.h"
#include "outer_poly.h"
#include "safety_factor.h"
#include "fixed_surfaces.h"
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <boost\archive\text_oarchive.hpp>
#include <boost\archive\text_iarchive.hpp>
#include <boost\graph\biconnected_components.hpp>
#include <boost\graph\adjacency_list.hpp>
using namespace std;

namespace BI
{

int adj_blk_graph::output_vertices_type_info (const block_system &bs, ofstream &outfile)
{
	if (!outfile)
		return 1;
	graph::vertex_iterator v, v_end;
	boost::tie (v, v_end) = boost::vertices(m_graph);
	for (;v!=v_end;++v)
	{

		outfile<<bs.vcb[m_graph[*v].vcb[0]].vb.size ()<<"\t"<<bs.cb_volume_est(m_graph[*v].vcb[0])<<"\t";
		switch (m_graph[*v].type)
		{
		case EXCAV:
		case FIXED_AND_EXCAV:
			outfile<<"Excaved"; break;
		case UNDEXPOSED:
		case FIXED_AND_UNEXPOSED:
			outfile<<"Unexposed"; break;

		default:
			break;
		}
		outfile<<endl;
	}
	return 0;
}

int adj_blk_graph::init_graph (const block_system &bs,const bool &debug)
{
	//ÿ��cb����ʼ��Ϊһ�����㡣
	vector<graph::vertex_descriptor> v_gvd;
	v_gvd.reserve (bs.vcb.size());

	for (int i(0);i!=bs.vcb.size ();++i)
	{
		v_gvd.push_back (boost::add_vertex(m_graph));
		m_graph[v_gvd.back ()].vcb.push_back (i);
		m_graph[v_gvd.back ()].id = i;
		m_graph[v_gvd.back ()].volume = bs.cb_volume_est(i);
		m_graph[v_gvd.back ()].type = GROUP_TYPE::UNDEXPOSED;
		for (int j(0);j!=bs.vcb_bsf[i].size ();++j)
		{
			if (bs.vbsf[bs.vcb_bsf[i][j]].type == block_system::FACE_TYPE::EXCAV)
			{
				m_graph[v_gvd.back ()].type = GROUP_TYPE::EXCAV;		//��һ�����ǿ�������ô���������ǿ��ڿ���
				break;
			}
			else
				continue;
		}
	}

	//���ÿ���ڽӹ�ϵ �����ַ�ʽÿ���ڽӹ�ϵ�ᱻ������Σ���������Խ������̫��Ӱ�죩
	for (int i(0);i!=bs.vbsf.size ();++i)
	{
		pair<graph::edge_descriptor,bool> res;
		res.second = false;
		for (int j(0);j!=bs.vbsf [i].vadj_bsf_i.size ();++j)
		{
			//��������
			res = boost::add_edge (v_gvd[bs.vbsf[i].cbi],v_gvd[bs.vbsf[bs.vbsf[i].vadj_bsf_i[j]].cbi],m_graph);
			add_set(m_graph[res.first].vbsf_pair,int_pair(i,bs.vbsf[i].vadj_bsf_i[j]),int_pair_comp());
		}
	}
	return 0;
}

int adj_blk_graph::init_vob (const block_system &bs, std::vector<outer_bound> &vob) const
{
	vob.resize (bs.vcb.size ());
	graph::vertex_iterator v, v_end;
	for (boost::tie (v,v_end) = boost::vertices (m_graph);v!=v_end;++v)
		vob[m_graph[*v].id].construct_outer_bound(bs,m_graph[*v].vcb);

	return 0;
}

int adj_blk_graph::update_fixed_vertices (const block_system &bs, const fixed_surfaces & fs, const bool &debug_flag)
{
	if (debug_flag)
		cout<<"init_fixed_vertices()"<<endl;
	
	graph::vertex_iterator v, v_end;
	for (boost::tie (v, v_end) = boost::vertices (m_graph); v!=v_end; ++v)
	{
		if (is_fixed(bs,fs,*v))		//����ǹ̶���
		{
			switch (m_graph[*v].type)
			{
			case GROUP_TYPE::EXCAV:
				m_graph[*v].type = GROUP_TYPE::FIXED_AND_EXCAV; break;
			case GROUP_TYPE::UNDEXPOSED:
				m_graph[*v].type = GROUP_TYPE::FIXED_AND_UNEXPOSED; break;
			default:
				break;
			}
		}
		else		//������ǹ̶���
		{
			switch (m_graph[*v].type)
			{
			case GROUP_TYPE::FIXED_AND_EXCAV:
				m_graph[*v].type = GROUP_TYPE::EXCAV; break;
			case GROUP_TYPE::FIXED_AND_UNEXPOSED:
				m_graph[*v].type = GROUP_TYPE::UNDEXPOSED; break;
			default:
				break;
			}
		}
	}
	return 0;
}

adj_blk_graph::graph::vertex_descriptor adj_blk_graph::merge_vertices (graph::vertex_descriptor &v1, graph::vertex_descriptor &v2)
{
	//v1��v2������Ҳ���Ժϲ�
	//�Ⱥϲ���������������ԡ�
	m_graph[v1].vcb.insert (m_graph[v1].vcb.end (),m_graph[v2].vcb.begin (),m_graph[v2].vcb.end ());
	if (m_graph[v1].type == GROUP_TYPE::EXCAV)
		if (m_graph[v2].type == GROUP_TYPE::FIXED_AND_EXCAV || m_graph[v2].type == GROUP_TYPE::FIXED_AND_UNEXPOSED)
			m_graph[v1].type = GROUP_TYPE::FIXED_AND_EXCAV;
	if (m_graph[v1].type == GROUP_TYPE::FIXED_AND_UNEXPOSED)
		if (m_graph[v2].type == GROUP_TYPE::EXCAV || m_graph[v2].type == GROUP_TYPE::FIXED_AND_EXCAV)
			m_graph[v1].type = GROUP_TYPE::FIXED_AND_EXCAV;
	if (m_graph[v1].type == GROUP_TYPE::UNDEXPOSED)
		m_graph[v1].type = m_graph[v2].type;

	m_graph[v1].volume += m_graph[v2].volume;

	//Ȼ��v2������������ڽӹ�ϵת�Ƶ�v1����
	graph::out_edge_iterator e, e_end;
	boost::tie(e,e_end) = boost::out_edges(v2,m_graph);
	for (;e!=e_end ;++e)
	{
		if (v1 == boost::target (*e,m_graph))
			continue;
		graph::edge_descriptor added_e; bool flag;
		boost::tie(added_e,flag) = boost::add_edge(v1,boost::target (*e,m_graph),m_graph);
		//��ԭ���ߵ����Ը����������
		m_graph[added_e].vbsf_pair.insert (m_graph[added_e].vbsf_pair.end(),
			m_graph[*e].vbsf_pair .begin (),m_graph[*e].vbsf_pair .end ());
	}
	boost::clear_vertex(v2,m_graph);
	boost::remove_vertex (v2,m_graph);
	return v1;
}

adj_blk_graph::graph::vertex_descriptor adj_blk_graph::merge_vertices (const std::vector<graph::vertex_descriptor> &vv)
{
	if (vv.empty ())
		return 0;
	vector<graph::vertex_descriptor>::const_iterator i(vv.begin ());
	graph::vertex_descriptor rv (*i);		//ֻ��������㱻����
	++i;
	GROUP_TYPE new_type (m_graph[rv].type);

	//����ʣ�µ����ж���
	for (;i!=vv.end ();++i)
	{
		graph::vertex_descriptor v = *i;
		m_graph[rv].vcb.insert (m_graph[rv].vcb.end(),m_graph[v].vcb.begin (),m_graph[v].vcb.end ());
		if (new_type == GROUP_TYPE::EXCAV)
			if (m_graph[v].type == GROUP_TYPE::FIXED_AND_EXCAV || m_graph[v].type == GROUP_TYPE::FIXED_AND_UNEXPOSED)
				new_type = GROUP_TYPE::FIXED_AND_EXCAV;
		if (new_type == GROUP_TYPE::FIXED_AND_UNEXPOSED)
			if (m_graph[v].type == GROUP_TYPE::EXCAV || m_graph[v].type == GROUP_TYPE::FIXED_AND_EXCAV)
				new_type = GROUP_TYPE::FIXED_AND_EXCAV;
		if (new_type == GROUP_TYPE::UNDEXPOSED)
			new_type = m_graph[v].type;
		m_graph[rv].volume += m_graph[v].volume;
		
		graph::out_edge_iterator e,e_end;
		boost::tie (e,e_end) = boost::out_edges(v,m_graph);
		//����������������ӵ�ÿ����
		for (;e!=e_end;++e)
		{
			if (find(vv.begin (),vv.end(),boost::target (*e,m_graph)) != vv.end ())
				continue;		//��������vv����������ıߣ���ֱ��������Ҳ����˵����߻ᱻɾ��
			graph::edge_descriptor added_e; bool flag;
			boost::tie (added_e,flag) = boost::add_edge (rv,boost::target (*e,m_graph),m_graph);
			//��ԭ���Ǹ��ߵ����Ը��ӵ��������ӵı���
			m_graph[added_e].vbsf_pair.insert (m_graph[added_e].vbsf_pair .end(),
				m_graph[*e].vbsf_pair.begin (),m_graph[*e].vbsf_pair .end ());
		}
		boost::clear_vertex(v,m_graph);
		boost::remove_vertex(v,m_graph);
	}
	m_graph[rv].type = new_type;
	return rv;
}

mbox_3 adj_blk_graph::build_mbox_3(const block_system &bs, const graph::vertex_descriptor &v) const
{
	mbox_3 res;
	if (m_graph[v].vcb.empty ())
		return res;
	res = BI::build_mbox_3(bs,bs.vcb [m_graph[v].vcb[0]]);
	for (int i(0);i!=m_graph[v].vcb.size ();++i)
	{
		res += BI::build_mbox_3(bs,bs.vcb[m_graph[v].vcb[i]]);
	}
	return res;
}

void adj_blk_graph::output_node_and_edge_info (std::ostream &outfile) const
{
	if (!outfile)
		return;
	vector<int> count(4,0);
	graph::vertex_iterator v, v_end;
	for (boost::tie (v, v_end) = boost::vertices (m_graph); v!=v_end; ++v)
	{
		switch (m_graph[*v].type )
		{
		case GROUP_TYPE::EXCAV:
			count[0]++; break;
		case GROUP_TYPE::FIXED_AND_EXCAV:
			count[1]++; break;
		case GROUP_TYPE::FIXED_AND_UNEXPOSED:
			count[2]++; break;
		case GROUP_TYPE::UNDEXPOSED:
			count[3]++; break;
		default:
			break;
		}
	}

	outfile<<"***Stability analysis info.***"<<endl;
	outfile<<"Valid block groups:\t"<<boost::num_vertices(m_graph)<<endl;
	outfile<<"Excaved blocks (unfixed):\t"<<count[0]<<endl;
	outfile<<"Excaved blocks (fixed):\t"<<count[1]<<endl;
	outfile<<"Unexposed blocks (unfixed):\t"<<count[3]<<endl;
	outfile<<"Unexposed blocks (fixed):\t"<<count[2]<<endl;

}

void adj_blk_graph::output_chain_falling_info_CSV (const block_system &bs, std::ostream &outfile) const
{
	if (!outfile)
		return;
	outfile.precision (15);
	outfile<<"No.,vertices_count,volume,blk_indices"<<endl;
	for (int i(0);i!=v_cfalling.size ();++i)
	{
		outfile<<i<<","<<v_cfalling[i].size ()<<",";
		double vol(0);
		for (int j(0);j!=v_cfalling[i].size ();++j)
			for (int k(0);k!=m_graph[v_cfalling[i][j]].vcb.size ();++k)
				vol+= bs.cb_volume_est(m_graph[v_cfalling[i][j]].vcb[k]);
		outfile<<vol<<",";
		for (int j(0);j!=v_cfalling[i].size ();++j)
		{
			outfile<<"(";
			outfile<<"["<<m_graph[v_cfalling[i][j]].id<<"] ";
			for (int k(0);k!=m_graph[v_cfalling[i][j]].vcb.size ();++k)
				outfile<<m_graph[v_cfalling[i][j]].vcb[k]<<";";
			outfile<<")";
		}
		outfile<<endl;
	}
}

void adj_blk_graph::ouput_removable_blks (const block_system &bs, const vector_3 &r,std::ostream &outfile) const
{
	if (!outfile)
		return;
	outfile.precision (15);
	graph::vertex_iterator v, v_end;
	for (boost::tie (v,v_end) = boost::vertices (m_graph);v!=v_end;++v)
	{
		safety_factor sf;
		switch (m_graph[*v].type)
		{
		case GROUP_TYPE::EXCAV:
			sf = cal_sf(bs,*v,r);
			if (sf.type == safety_factor::S_TYPE::DOU_F_SLID ||
				sf.type == safety_factor::S_TYPE::SIG_F_SLID ||
				sf.type == safety_factor::S_TYPE::FALLING)
			{
				output_vertex_blk(bs,*v,outfile);
			}
		default:
			break;
		}
	}
}



bool adj_blk_graph::is_nested (const std::vector<graph::vertex_descriptor> & vv, const std::vector <graph::vertex_descriptor> & vap ,const bool &debug_flag) const
{
	if (debug_flag)
		cout<<"is_nested() "<<endl;

	int n_excav(0), n_fixed_unexp(0), n_unexp(0), n_fixed_and_excav(0);
	if (vv.size ()<=1)		//ֻ��һ���ڵ��û�нڵ㣬����ʵʩ�ϲ�����
		return false;

	int count(0);
	for (int i(0);i!= (int) vv.size ();++i)
	{
		switch (m_graph[vv[i]].type)
		{
		case GROUP_TYPE::EXCAV:
			n_excav++; count = i;break;
		case GROUP_TYPE::FIXED_AND_EXCAV:
			n_fixed_and_excav++; count = i;break;
		case GROUP_TYPE::FIXED_AND_UNEXPOSED:
			n_fixed_unexp++; break;
		case GROUP_TYPE::UNDEXPOSED:
			n_unexp++; break;
		default:
			break;
		}
	}

	if(debug_flag)
		cout<<"n_excav:"<<n_excav<<" n_fixed_and_unexposed:"<<n_fixed_unexp<<" n_unexp:"<<n_unexp<<" n_fixed_and_excav:"<<n_fixed_and_excav<<endl;

	if ( (n_excav + n_fixed_and_excav) > 1)		//�ж��������Ŀ��ڿ��������
		return false;	
	else if ((n_excav + n_fixed_and_excav) == 1)
	{
		if (find(vap.begin (),vap.end (),vv[count]) == vap.end ())
			return false;		//������Ψһ�Ŀ��ڿ��岻�ǹؽڽڵ�
		else 
			return true;
	}
	else
	{
		for (int i(0);i!=vv.size ();++i)
			if (find(vap.begin (), vap.end (),vv[i]) != vap.end ())		//�����Щ�ڵ�����һ���ؽڽڵ㣬����Ժϲ�
				return true;
		return false;		//���û�п��ڿ��壬����û�йؽڽڵ㣬�򲻺ϲ���
	}
}

bool adj_blk_graph::is_removable (const block_system &bs, const std::vector<graph::edge_descriptor> &ve, const graph::vertex_descriptor &v, const vector_3 &r, const bool & debug_flag) const
{
	if (debug_flag)
		cout<<"is_removable()"<<endl;

	safety_factor sf(cal_sf(bs,ve,v,r,debug_flag));
	if (sf.type == safety_factor::SIG_F_SLID || sf.type == safety_factor::DOU_F_SLID || sf.type == safety_factor::FALLING)		//ֻ�����ʱ���������r�������¿��ƶ�
		return true;
	else
		return false;
}

bool adj_blk_graph::is_fixed (const block_system &bs, const fixed_surfaces &fs, const graph::vertex_descriptor &v) const
{
	for (int i(0);i!=m_graph[v].vcb.size ();++i)
		if (fs.do_intersect (bs,bs.vcb[m_graph[v].vcb[i]]))		//ֻҪ��һ��cb��fs�ཻ������ǹ̶��ڵ�
			return true;
	return false;
}

safety_factor adj_blk_graph::cal_sf (const block_system& bs, const graph::vertex_descriptor&v, const vector_3 &r, const bool& debug_flag) const
{
	if (debug_flag)
		cout<<"cal_sf()"<<endl;
	vector<graph::edge_descriptor> ve;
	graph::out_edge_iterator e, e_end;
	boost::tie(e, e_end) = boost::out_edges(v,m_graph);
	ve.insert (ve.end (),e, e_end);
	return cal_sf(bs,ve,v,r,debug_flag);
}

safety_factor adj_blk_graph::cal_sf (const block_system& bs, const std::vector<graph::edge_descriptor> &ve, const graph::vertex_descriptor &v, const vector_3 &r, const bool & debug_flag) const
{
	if (debug_flag)
		cout<<"cal_sf()"<<endl;
	vector<vector_3> vfixed;
	for (int i(0);i!=ve.size ();++i)		//����ÿһ����
	{
		for (int j(0);j!=m_graph[ve[i]].vbsf_pair.size ();++j)		//���������ߵ�ÿһ����������
		{
			int first (m_graph[ve[i]].vbsf_pair [j].first), second (m_graph[ve[i]].vbsf_pair [j].second );		//�﷨��
			if (find (m_graph[v].vcb.begin (), m_graph[v].vcb.end (), bs.vbsf[first].cbi) != m_graph[v].vcb.end ())
			{
				//���first��ʾ�ı߽�������v��ĳ�������
				if (bs.vbsf[first].outer_normal )
					add_set (vfixed,-bs.vpl[bs.vbsf[first].pli].orthogonal_vector());
				else
					add_set (vfixed,bs.vpl[bs.vbsf[first].pli].orthogonal_vector());
			}
			else if (find (m_graph[v].vcb.begin (), m_graph[v].vcb.end (), bs.vbsf[second].cbi) != m_graph[v].vcb.end ())
			{
				//���second��ʾ�ı߽�������v��ĳ�������
				if (bs.vbsf[second].outer_normal)
					add_set (vfixed, -bs.vpl[bs.vbsf[second].pli].orthogonal_vector());
				else
					add_set (vfixed, bs.vpl[bs.vbsf[second].pli].orthogonal_vector());
			}
			else
				throw logic_error ("adj_blk_graph::is_removable, vector v is not on edge");
		}
	}
	safety_factor sf;
	cal_safety_factor(sf,vfixed,r);
	return sf;
}

bool adj_blk_graph::is_interlocked (const block_system &bs, const graph::edge_descriptor  &e) const
{
	graph::vertex_descriptor v1 = boost::target (e,m_graph), v2 = boost::source(e,m_graph);

	vector<vector_3> fixed_vector;
	const vector<int> &vcb1 = m_graph[v1].vcb;		//�﷨��
	const vector<int> &vcb2 = m_graph[v2].vcb;

	for (int i(0);i!=m_graph[e].vbsf_pair.size ();++i)
	{
		int bsf1(-1), bsf2(-1);		//�ֱ��ʾ��vcb1�ϻ���vcb2�ϵ������
		if (find (vcb1.begin (),vcb1.end (), bs.vbsf[m_graph[e].vbsf_pair[i].first].cbi) != vcb1.end ())
			bsf1 = m_graph[e].vbsf_pair[i].first;
		else if (bsf1 == -1 && find (vcb1.begin (),vcb1.end (), bs.vbsf[m_graph[e].vbsf_pair[i].second ].cbi) != vcb1.end ())
			bsf1 = m_graph[e].vbsf_pair[i].second;
		else 
			throw logic_error ("adj_blk_graph::is_interlocked, edge property error1");

		if (find (vcb2.begin (),vcb2.end (), bs.vbsf[m_graph[e].vbsf_pair[i].first].cbi) != vcb2.end ())
			bsf2 = m_graph[e].vbsf_pair[i].first ;
		else if (bsf2 == -1 && find (vcb2.begin (),vcb2.end (), bs.vbsf[m_graph[e].vbsf_pair[i].second].cbi) != vcb2.end ())
			bsf2 = m_graph[e].vbsf_pair[i].second ;
		else 
			throw logic_error("adj_blk_graph::is_interlocked, edge property error2");

		if (bs.vbsf[bsf1].pli != bs.vbsf[bsf2].pli)
			throw logic_error("adj_blk_graph::is_interlocked, edge property error3");		//˳����һ�´�

		if (bs.vbsf[bsf1].outer_normal)
			fixed_vector.push_back (bs.vpl[bs.vbsf[bsf1].pli].orthogonal_vector());		//ʹ�����е�����ָ��vcb2
		else
			fixed_vector.push_back (-bs.vpl[bs.vbsf[bsf1].pli].orthogonal_vector());
	}
	if (removability(fixed_vector) == Removability_Status::Removable )
		return false;
	else
		return true;		//�����һ������������ߵ��޶��²������ƶ�����������������������
}

bool adj_blk_graph::is_interlocked (const block_system &bs, const graph::vertex_descriptor & v1, const graph::vertex_descriptor & v2) const
{
	graph::edge_descriptor e; bool flag;
	boost::tie(e,flag) = boost::edge(v1,v2,m_graph);
	if (!flag)
		return false;		//���������㶼�����������Կ϶�����������.
	return is_interlocked(bs,e);
}

int adj_blk_graph::eliminate_nested_blk (const bool &debug_flag)
{
	if (debug_flag)
		cout<<"eliminate_nested_blk"<<endl;
	//ָ��һ���ⲿ��property_map���궨�����ڵ㡣biconnected_components����㷨��Ҫ�õ�
	boost::property_map<graph, std::size_t blk_junction::*>::type component = boost::get(&blk_junction::component_id,m_graph);		//��ʾ����������˫��ͨ������ӳ��
	vector<graph::vertex_descriptor> vap;		//���йؽڽڵ�ļ���

	std::size_t num_comps; back_insert_iterator<vector<graph::vertex_descriptor>> inter_temp(vap);
	//�����ṩһ��property map����������ӳ�䵽ĳ�������ϲ��У�boost::vertex_index_map���������BGL���ֲ�
	boost::tie(num_comps,inter_temp)= boost::biconnected_components(m_graph,component,back_inserter(vap),boost::vertex_index_map(boost::get(&blk_group::id,m_graph)));

	vector<vector<int>> v_ap_comps(vap.size ());				//��¼ÿ���ؽڽڵ�����Щ˫��ͨ������
	vector<vector<int>> v_comps_ap(num_comps);					//��¼ÿ��˫��ͨ����������Щ�ؽڽڵ�
	for (int i(0);i!=vap.size ();++i)
	{
		graph::out_edge_iterator e,e_end;
		boost::tie (e, e_end) = boost::out_edges(vap[i],m_graph);
		for (;e!=e_end;++e)
			add_set (v_ap_comps[i],(int)component[*e]);
		for (int j(0);j!=v_ap_comps[i].size ();++j)
			add_set(v_comps_ap[v_ap_comps[i][j]],i);
	}
	
	//�ֱ��¼һ�¸���˫��ͨ�����а�����Щ�ڵ�
	vector<vector<graph::vertex_descriptor>> vvv(num_comps);

	graph::edge_iterator e, e_end;
	boost::tie (e,e_end) = boost::edges(m_graph);
	for (;e!=e_end;++e)
	{
		int n( (int)component[*e]);
		add_set(vvv[n], boost::target(*e, m_graph));
		add_set(vvv[n], boost::source(*e, m_graph));
	}

	vector<int> v_flag (num_comps);		//����˫��ͨ����Ӧ����ô�ϲ�
	for (int i(0);i!=num_comps;++i)
	{
		if (is_nested(vvv[i],vap))		
			v_flag[i] = i;			//���˫��ͨ������Ҫ�ϲ�
		else
			v_flag[i] = -1;			//���˫��ͨ��������Ҫ�ϲ�
	}

	//�������Ӧ�úϲ�������˫��ͨͼ����ڵ㣬��ô�ýڵ�϶��ǹؽڽڵ㣬����������˫��ͨͼӦ��һ��ϲ�
	for (int i(0);i!=num_comps;++i)
	{
		if (v_flag[i]<0) continue;		//���˫��ͨ�鲻�úϲ�
		for (int j(0);j!=v_comps_ap[i].size ();++j)		//�������˫��ͨͼ�е����йؽڽڵ�
		{
			for (int k(0);k!=v_ap_comps[v_comps_ap[i][j]].size ();++k)		//��������ؽڽڵ����ڵ�����˫��ͨ����
			{
				int temp_i (v_ap_comps[v_comps_ap[i][j]][k]);		//���˫��ͨ����������
				if (temp_i == i)
					continue;
				if (v_flag[temp_i] < 0)
					continue;		//���˫��ͨͼ���úϲ���ֱ������
				//������v_flag�е�������v_flag[i]��ȫ���ĳ�v_flag[temp_i]
				int tar (v_flag[i]), src(v_flag[temp_i]);
				if (tar == src) continue;
				for (int l(0);l!=num_comps;++l)
					if (v_flag[l] == tar) v_flag[l] = src;
			}
		}
	}

	vector<int> v_group;		//��¼һ���м����ϲ���,ÿ���ϲ���Ӧ�������ɸ�˫��ͨ������� 
	for (int i(0);i!=num_comps;++i)
	{
		if (v_flag[i] < 0)
			continue;
		else
			add_set(v_group,v_flag[i]);
	}

	vector<vector<graph::vertex_descriptor>> vv_nested_group(v_group.size ());		//��¼��Щ�ڵ�Ӧ�ñ��ϲ�
	for (int i(0);i!=v_group.size ();++i)
	{
		for (int j(0);j!=num_comps;++j)
		{
			if (v_flag[j] != v_group[i]) continue;		//�ҵ�flag����v_group[i]����һ��
			for (int k(0);k!=vvv[j].size ();++k)
				add_set(vv_nested_group[i],vvv[j][k]);
		}
	}

	//����������vv_nested_group�еĸ��������е�Ԫ�ز�Ӧ�����غϵ�
	for (int i(0);i!=vv_nested_group.size ();++i)
		merge_vertices(vv_nested_group[i]);

	return (int)vv_nested_group.size ();
}

int adj_blk_graph::eliminate_interlocked_blk (const block_system& bs, const int& search_depth, const bool &debug_flag)
{
	if (debug_flag)
		cout<<"eliminate_interlocked_blk()"<<endl;

	if (search_depth != 1)	//Ŀǰֻ֧��һ���ߵ��������
		throw logic_error("adj_blk_graph::eliminate_interlocked_blk (const int& search_depth), only support 1 edge by now");
	vector<vector<graph::vertex_descriptor>> vvv;		//��¼��Щ����Ӧ�ñ��ϲ���vvv�а�����ÿ��������������Ӧ�úϲ��Ķ�������

	//�������ǰ����m_graph�в������ظ��ı�,����һ���ߵ���β���㲻����ͬ
	graph::edge_iterator e,e_end;
	for (boost::tie(e,e_end) = boost::edges(m_graph);e!=e_end;++e)
	{
		if (!is_interlocked(bs,*e))
			continue;
		graph::vertex_descriptor tar = boost::target (*e, m_graph), src = boost::source (*e, m_graph);

		int tari(-1), srci(-1);		//��¼tar��src�ֱ����ĸ������б���
		for (int i(0); i!=vvv.size ();++i)
		{
			if(tari == -1)
			{
				tari = find_n(vvv[i],tar);
				if (tari == vvv[i].size ())
					tari = -1;
				else
					tari = i;
			}
			if (srci == -1)
			{
				srci = find_n (vvv[i], src);
				if (srci == vvv[i].size ())
					srci = -1;
				else 
					srci = i;
			}
			if (tari != -1 && srci != -1)
				break;
		}
		if (tari == -1 && srci == -1)	//��������tar��src��û�б���ӹ�
		{
			vvv.push_back (vector<graph::vertex_descriptor>());
			vvv.back ().push_back (tar);
			vvv.back ().push_back (src);
			continue;
		}
		else if (tari == -1)		//��һ���Ѿ�����ӵ�vvv����
		{
			vvv[srci].push_back (tar); continue;
		}
		else if (srci == -1)
		{
			vvv[tari].push_back (src); continue;
		}
		else
		{
			if (srci == tari)		//�������㱻��ӵ�ͬһ�����У���ֱ������
				continue;
			//�������Ѿ�������ˣ��������ڲ�ͬ�����С���Ҫ����������ϲ�����
			int maxi,mini;
			maxi = max(tari,srci); mini = min(tari,srci);
			vvv[mini].insert (vvv[mini].end (),vvv[maxi].begin (),vvv[maxi].end ());
			vector<vector<graph::vertex_descriptor>>::iterator iter(vvv.begin ());
			for (int i(0);i!=vvv.size ();++i,++iter)
				if (i == maxi)
					break;
			vvv.erase (iter);
		}
	}
	size_t ori_vcount = boost::num_vertices(m_graph);

	for (int i(0);i!=vvv.size ();++i)
		merge_vertices(vvv[i]);
	return (int) (ori_vcount - boost::num_vertices(m_graph));
}

int adj_blk_graph::determine_chain_falling_blocks (const block_system &bs,const vector_3 &r, const bool &debug_flag)
{
	if (debug_flag)
		cout<<"determine_chain_falling_blocks()"<<endl;

	v_cfalling.clear ();
	//��¼ÿ��������ȶ�״̬�����������屻�̶���״̬��
	graph::vertex_iterator v, v_end;
	vector<safety_factor> vsf(bs.vcb.size ());		//m_graph���������ô������㡣������������Ϣ������Ϊ�˷�������ô��
	for (boost::tie(v, v_end) = boost::vertices (m_graph); v!=v_end;++v)
	{
		switch (m_graph[*v].type )
		{
		case GROUP_TYPE::FIXED_AND_EXCAV:
		case GROUP_TYPE::FIXED_AND_UNEXPOSED:
			vsf[m_graph[*v].id].type = safety_factor::S_TYPE::FIXED; break;		//�����������������Ϣ�������о��аɣ����ں���һ���ڴ� 
		case GROUP_TYPE::UNDEXPOSED:
			vsf[m_graph[*v].id].type = safety_factor::S_TYPE::UNEXPOSED; break;
		case GROUP_TYPE::EXCAV:
			vsf[m_graph[*v].id] = cal_sf(bs,*v,r); break;
		default:
			throw logic_error ("adj_blk_graph::determine_chain_falling_blocks, wrong group type");
		}
	}

	size_t nv = boost::num_vertices(m_graph);
	size_t count(0);
	for (boost::tie (v, v_end) = boost::vertices (m_graph);v!=v_end;++v)
	{
		if (debug_flag)
			cout<<count++<<"/"<<nv<<endl;
		if (m_graph[*v].type != GROUP_TYPE::EXCAV)
		{
			//if (debug_flag)
			//	cout<<"Not Excav "<<endl;
			continue;		//�����ĸ����Ͷ��ǲ����ƶ��ġ����ٲ���triggerһ����������
		}
		
		//if (debug_flag)
		//	cout<<"EXCAV "<<endl;
	
		vector<graph::edge_descriptor> ve;
		graph::out_edge_iterator e,e_end;
		for (boost::tie(e,e_end) = boost::out_edges(*v, m_graph);e!=e_end;++e)
			ve.push_back (*e);
	
		//if (debug_flag)
		//	cout<<"ve.size():"<<ve.size ()<<endl;
	
		if (!(vsf[m_graph[*v].id].type == safety_factor::S_TYPE::DOU_F_SLID ||
			vsf[m_graph[*v].id].type == safety_factor::S_TYPE::SIG_F_SLID ||
			vsf[m_graph[*v].id].type == safety_factor::S_TYPE::FALLING))
		{
			//if (debug_flag)
			//	cout<<"Not keyblock "<<endl;
			continue;		//���������屾���ǹؼ����壬��϶���������һ����������
		}
	
		//if (debug_flag)
		//	cout<<"Keyblock "<<endl;
	
		//*v�������ᶨ��һ���µ������������,����֮ǰ��û�б�����������������������С�
		v_cfalling.push_back (vector<graph::vertex_descriptor>());
		v_cfalling.back ().push_back (*v);
	
		list<graph::vertex_descriptor> v_cand;		//�ܵ�Ӱ��Ķ���
		graph::adjacency_iterator adj_v, adj_v_end;
		boost::tie (adj_v, adj_v_end) = boost::adjacent_vertices(*v,m_graph);
		v_cand.insert (v_cand.end (),adj_v,adj_v_end);
	
		//if (debug_flag)
		//	cout<<"v_cand.size ():"<<v_cand.size ()<<endl;
	
		while (!v_cand.empty ())
		{
			graph::vertex_descriptor v_temp(v_cand.front ());
			v_cand.pop_front();
	
			if (m_graph[v_temp].type == GROUP_TYPE::FIXED_AND_EXCAV || m_graph[v_temp].type == GROUP_TYPE::FIXED_AND_UNEXPOSED)
				continue;				//�̶��Ľڵ�Ͳ��������
			if (m_graph[v_temp].type == GROUP_TYPE::EXCAV)
				if (vsf[m_graph[v_temp].id].type == safety_factor::S_TYPE::SIG_F_SLID ||
					vsf[m_graph[v_temp].id].type == safety_factor::S_TYPE::DOU_F_SLID ||
					vsf[m_graph[v_temp].id].type == safety_factor::S_TYPE::FALLING)
					continue;			//����������ǹؼ����壬���������
			if (find(v_cfalling.back ().begin (), v_cfalling.back ().end (), v_temp) != v_cfalling.back ().end ())
				continue;				//������Ѿ��������������������ˡ��Ͳ��������
	
			vector<graph::edge_descriptor> ve;
			graph::out_edge_iterator e,e_end;
			for (boost::tie (e,e_end) = boost::out_edges(v_temp,m_graph);e!=e_end;++e)
				if (find (v_cfalling.back ().begin (), v_cfalling.back ().end (), boost::target(*e, m_graph)) == v_cfalling.back ().end ())
					ve.push_back (*e);		//�������ߵ���һ��û����v_cfalling.back()�С�����һ���޶���
			if (is_removable(bs,ve,v_temp,r))
			{
				v_cfalling.back ().push_back (v_temp);		//����������r���������ƶ�
				//������Ӱ�쵽�Ŀ��嶼��ӵ�v_cand��
				graph::adjacency_iterator adj_v, adj_v_end;
				boost::tie (adj_v, adj_v_end) = boost::adjacent_vertices(v_temp,m_graph);
				v_cand.insert (v_cand.end (),adj_v,adj_v_end);
			}
		}
	}
	return 0;
}

int adj_blk_graph::output_vertex_blk (const block_system &bs, const graph::vertex_descriptor &v,std::ostream &outfile) const
{
	if (!outfile)
		return 1;
	outer_bound ob;
	ob.construct_outer_bound_using_vcb(bs,m_graph[v].vcb);
	ob.output_STL(outfile);
	return 0;
}

int adj_blk_graph::output_cfalling_blks_quick (const block_system &bs, std::ostream &os, const std::vector<graph::vertex_descriptor> &vv) const
{
	if (!os)
		return 1;
	for (int i(0);i!=vv.size ();++i)
	{
		os<<"solid "<<to_string (m_graph[vv[i]].id)<<endl;
		for (int j(0);j!=m_graph[vv[i]].vcb.size ();++j)
			bs.test_output_STL(os,m_graph[vv[i]].vcb[j],0,false);
		os<<"endsolid "<<to_string (m_graph[vv[i]].id)<<endl;
	}
	return 0;
}

int adj_blk_graph::output_cfalling_blks (const block_system &bs, std::ostream  &os, const std::vector<graph::vertex_descriptor> &vv) const
{
	if (!os)
		return 1;
	for (int i(0);i!=vv.size ();++i)
		output_vertex_blk(bs,vv[i],os);
	return 0;
}

int adj_blk_graph::output_CSV (std::ostream &outfile, const ELE_TYPE& et) const
{
	if (!outfile)
		return 1;

	if (et == ELE_TYPE::Node)
	{
		outfile<<"id,vcb,type, volume"<<endl;
		graph::vertex_iterator v,v_end;
		boost::tie(v,v_end) = boost::vertices (m_graph);
		for (;v!=v_end;++v)
		{
			outfile<<m_graph[*v].id<<",[";
			for (int i(0);i!=m_graph[*v].vcb.size ();++i)
				outfile<<m_graph[*v].vcb[i]<<" ";
			outfile<<"],";
			switch (m_graph[*v].type)
			{
			case GROUP_TYPE::EXCAV:
				outfile<<"EXCAV,"; break;
			case GROUP_TYPE::FIXED_AND_EXCAV:
				outfile<<"FIXED_AND_EXCAV,"; break;
			case GROUP_TYPE::FIXED_AND_UNEXPOSED:
				outfile<<"FIXED_AND_UNEXPOSED,"; break;
			case GROUP_TYPE::UNDEXPOSED:
				outfile<<"UNEXPOSED,"; break;
			default:
				break;
			}
			outfile<<m_graph[*v].volume<<endl;
		}
		return 0;
	}
	else
	{
		outfile<<"Source,Target,Type,bsf_pairs"<<endl;
		graph::edge_iterator e, e_end;
		boost::tie (e, e_end) = boost::edges (m_graph);
		for (;e!=e_end;++e)
		{
			outfile<<m_graph[boost::source (*e,m_graph)].id<<","
				<<m_graph[boost::target (*e,m_graph)].id<<",Undirected,[";
			for (vector<int_pair>::const_iterator i = m_graph[*e].vbsf_pair.begin ();
				i!=m_graph[*e].vbsf_pair.end ();++i)
			{
				outfile<<"("<<i->first<<" "<<i->second<<") ";
			}
			outfile<<"]"<<endl;
		}
		return 0;
	}
}

int adj_blk_graph::output_safety_factor (std::ostream &outfile, const block_system &bs) const
{
	if (!outfile)
		return 1;

	outfile<<"id,volume,weight(kg),type,sliding_f,n1,n2,f1,f2,safety_factor"<<endl;
	graph::vertex_iterator v, v_end;
	double density (2700);
	double angle (30.0*M_PI/180.0);
	boost::tie (v,v_end) = boost::vertices (m_graph);
	for (;v!=v_end;++v)
	{
		outfile<<m_graph[*v].id<<","<<m_graph[*v].volume<<","<<m_graph[*v].volume*density<<",";
		safety_factor sf;
		vector_3 r (0,0,-m_graph[*v].volume*density*9.8);
		switch (m_graph[*v].type)
		{
		case GROUP_TYPE::EXCAV :
			sf = cal_sf (bs,*v,r);
			switch (sf.type )
			{
			case safety_factor::S_TYPE::DOU_F_SLID:
				outfile<<"Double_face,"<<sf.sliding_force<<","<<sf.normal_forcei<<","<<sf.normal_forcej<<","
					<<sf.normal_forcei*tan (angle)<<","<<sf.normal_forcej*tan (angle)<<","
					<<(sf.normal_forcei*tan (angle)+ sf.normal_forcej*tan (angle))/sf.sliding_force<<endl;
				break;
			case safety_factor::S_TYPE::SIG_F_SLID:
				outfile<<"Single_face,"<<sf.sliding_force<<","<<sf.normal_forcei<<","<<0<<","
					<<sf.normal_forcei*tan (angle)<<","<<0<<",-"<<endl;
				break;
			case safety_factor::S_TYPE::FALLING:
				outfile<<"Falling,"<<0<<","<<0<<","<<0<<","<<0<<","<<0<<",-"<<endl;
				break;
			case safety_factor::S_TYPE::NON_REMOVABLE:
				outfile<<"Non_removable,"<<0<<","<<0<<","<<0<<","<<0<<","<<0<<",-"<<endl;
				break;
			case safety_factor::S_TYPE::REMOVABLE_BUT_SAFE:
				outfile<<"Removable_but_safe,"<<0<<","<<0<<","<<0<<","<<0<<","<<0<<",-"<<endl;
				break;
			default:
				break;
			}
			break;
		case GROUP_TYPE::FIXED_AND_EXCAV:
			outfile<<"FIXED_AND_EXCAV,"<<0<<","<<0<<","<<0<<","<<0<<","<<0<<",-"<<endl; break;
		case GROUP_TYPE::FIXED_AND_UNEXPOSED:
			outfile<<"FIXED_AND_UNEXPOSED,"<<0<<","<<0<<","<<0<<","<<0<<","<<0<<",-"<<endl; break;
		case GROUP_TYPE::UNDEXPOSED:
			outfile<<"UNEXPOSED,"<<0<<","<<0<<","<<0<<","<<0<<","<<0<<",-"<<endl; break;
		default:
			break;
		}
	}
	return 0;
}



int adj_blk_graph::output_blk_group_quick (std::ostream & os, const block_system &bs, const GROUP_TYPE &type) const
{
	graph::vertex_iterator v, v_end;
	boost::tie (v, v_end) = boost::vertices(m_graph);
	for (;v!=v_end;++v)
		if (m_graph[*v].type == type)
		{
			os<<"solid "<<to_string(m_graph[*v].id)<<endl;
			for (int i(0);i!=m_graph[*v].vcb.size ();++i)
				bs.test_output_STL(os,m_graph[*v].vcb[i],0,false);
			os<<"endsolid "<<to_string(m_graph[*v].id)<<endl;
		}
	return 0;
}

safety_factor adj_blk_graph::test_sf_of_vertex (const block_system &bs, const int &id, const vector_3 &r) const
{
	safety_factor sf;
	vector<vector_3> vfixed;
	graph::vertex_iterator v, v_end;
	for (boost::tie (v,v_end) = boost::vertices (m_graph); v!=v_end;++v)
		if (m_graph[*v].id == id)
		{
			graph::out_edge_iterator e,e_end;
			for (boost::tie (e,e_end) = boost::out_edges(*v,m_graph); e!=e_end;++e)
			{
				for (int i(0);i!= m_graph[*e].vbsf_pair .size ();++i)
				{
					int first (m_graph[*e].vbsf_pair [i].first), second (m_graph[*e].vbsf_pair [i].second );
					if (find (m_graph[*v].vcb.begin (), m_graph[*v].vcb.end(), bs.vbsf[first].cbi) != m_graph[*v].vcb.end())
					{
						if (bs.vbsf [first].outer_normal )
							add_set (vfixed,-bs.vpl[bs.vbsf[first].pli ].orthogonal_vector());
						else
							add_set (vfixed,bs.vpl[bs.vbsf[first].pli ].orthogonal_vector());
					}
					else if (find (m_graph[*v].vcb.begin (), m_graph[*v].vcb.end(), bs.vbsf[second].cbi) != m_graph[*v].vcb.end())
					{
						if (bs.vbsf [second].outer_normal)
							add_set (vfixed,-bs.vpl[bs.vbsf[second].pli].orthogonal_vector());
						else
							add_set (vfixed,bs.vpl[bs.vbsf[second].pli].orthogonal_vector());
					}
					else 
						throw logic_error ("test fun.");
				}
			}
			cal_safety_factor(sf,vfixed,r);
			break;
		}
	return sf;
}

}
