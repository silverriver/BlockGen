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
	//每个cb被初始化为一个顶点。
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
				m_graph[v_gvd.back ()].type = GROUP_TYPE::EXCAV;		//有一个面是开挖面那么这个块体就是开挖块体
				break;
			}
			else
				continue;
		}
	}

	//添加每个邻接关系 （这种方式每个邻接关系会被添加两次，不过不会对结果产生太大影响）
	for (int i(0);i!=bs.vbsf.size ();++i)
	{
		pair<graph::edge_descriptor,bool> res;
		res.second = false;
		for (int j(0);j!=bs.vbsf [i].vadj_bsf_i.size ();++j)
		{
			//添加这个边
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
		if (is_fixed(bs,fs,*v))		//如果是固定的
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
		else		//如果不是固定的
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
	//v1和v2不相邻也可以合并
	//先合并这两个顶点的属性。
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

	//然后将v2与其他顶点的邻接关系转移到v1上来
	graph::out_edge_iterator e, e_end;
	boost::tie(e,e_end) = boost::out_edges(v2,m_graph);
	for (;e!=e_end ;++e)
	{
		if (v1 == boost::target (*e,m_graph))
			continue;
		graph::edge_descriptor added_e; bool flag;
		boost::tie(added_e,flag) = boost::add_edge(v1,boost::target (*e,m_graph),m_graph);
		//将原来边的属性附到这个边上
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
	graph::vertex_descriptor rv (*i);		//只有这个顶点被保留
	++i;
	GROUP_TYPE new_type (m_graph[rv].type);

	//对于剩下的所有顶点
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
		//对于这个顶点所连接的每个边
		for (;e!=e_end;++e)
		{
			if (find(vv.begin (),vv.end(),boost::target (*e,m_graph)) != vv.end ())
				continue;		//如果这边是vv中两个顶点的边，则直接跳过，也就是说这个边会被删除
			graph::edge_descriptor added_e; bool flag;
			boost::tie (added_e,flag) = boost::add_edge (rv,boost::target (*e,m_graph),m_graph);
			//将原来那个边的属性附加到这个新添加的边上
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
	if (vv.size ()<=1)		//只有一个节点或没有节点，不用实施合并操作
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

	if ( (n_excav + n_fixed_and_excav) > 1)		//有多于两个的开挖块体就跳过
		return false;	
	else if ((n_excav + n_fixed_and_excav) == 1)
	{
		if (find(vap.begin (),vap.end (),vv[count]) == vap.end ())
			return false;		//如果这个唯一的开挖块体不是关节节点
		else 
			return true;
	}
	else
	{
		for (int i(0);i!=vv.size ();++i)
			if (find(vap.begin (), vap.end (),vv[i]) != vap.end ())		//如果这些节点中有一个关节节点，则可以合并
				return true;
		return false;		//如果没有开挖块体，并且没有关节节点，则不合并。
	}
}

bool adj_blk_graph::is_removable (const block_system &bs, const std::vector<graph::edge_descriptor> &ve, const graph::vertex_descriptor &v, const vector_3 &r, const bool & debug_flag) const
{
	if (debug_flag)
		cout<<"is_removable()"<<endl;

	safety_factor sf(cal_sf(bs,ve,v,r,debug_flag));
	if (sf.type == safety_factor::SIG_F_SLID || sf.type == safety_factor::DOU_F_SLID || sf.type == safety_factor::FALLING)		//只有这个时候才算是在r的驱动下可移动
		return true;
	else
		return false;
}

bool adj_blk_graph::is_fixed (const block_system &bs, const fixed_surfaces &fs, const graph::vertex_descriptor &v) const
{
	for (int i(0);i!=m_graph[v].vcb.size ();++i)
		if (fs.do_intersect (bs,bs.vcb[m_graph[v].vcb[i]]))		//只要有一个cb和fs相交则其就是固定节点
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
	for (int i(0);i!=ve.size ();++i)		//对于每一条边
	{
		for (int j(0);j!=m_graph[ve[i]].vbsf_pair.size ();++j)		//对于这条边的每一个多边形组合
		{
			int first (m_graph[ve[i]].vbsf_pair [j].first), second (m_graph[ve[i]].vbsf_pair [j].second );		//语法糖
			if (find (m_graph[v].vcb.begin (), m_graph[v].vcb.end (), bs.vbsf[first].cbi) != m_graph[v].vcb.end ())
			{
				//如果first表示的边界多边形是v中某个块体的
				if (bs.vbsf[first].outer_normal )
					add_set (vfixed,-bs.vpl[bs.vbsf[first].pli].orthogonal_vector());
				else
					add_set (vfixed,bs.vpl[bs.vbsf[first].pli].orthogonal_vector());
			}
			else if (find (m_graph[v].vcb.begin (), m_graph[v].vcb.end (), bs.vbsf[second].cbi) != m_graph[v].vcb.end ())
			{
				//如果second表示的边界多边形是v中某个块体的
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
	const vector<int> &vcb1 = m_graph[v1].vcb;		//语法糖
	const vector<int> &vcb2 = m_graph[v2].vcb;

	for (int i(0);i!=m_graph[e].vbsf_pair.size ();++i)
	{
		int bsf1(-1), bsf2(-1);		//分别表示在vcb1上或者vcb2上的外表面
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
			throw logic_error("adj_blk_graph::is_interlocked, edge property error3");		//顺带查一下错

		if (bs.vbsf[bsf1].outer_normal)
			fixed_vector.push_back (bs.vpl[bs.vbsf[bsf1].pli].orthogonal_vector());		//使得所有的向量指向vcb2
		else
			fixed_vector.push_back (-bs.vpl[bs.vbsf[bsf1].pli].orthogonal_vector());
	}
	if (removability(fixed_vector) == Removability_Status::Removable )
		return false;
	else
		return true;		//如果另一个块体在这个边的限定下不可以移动，则这两个块体是锁定的
}

bool adj_blk_graph::is_interlocked (const block_system &bs, const graph::vertex_descriptor & v1, const graph::vertex_descriptor & v2) const
{
	graph::edge_descriptor e; bool flag;
	boost::tie(e,flag) = boost::edge(v1,v2,m_graph);
	if (!flag)
		return false;		//这两个顶点都不相连，所以肯定不是锁定的.
	return is_interlocked(bs,e);
}

int adj_blk_graph::eliminate_nested_blk (const bool &debug_flag)
{
	if (debug_flag)
		cout<<"eliminate_nested_blk"<<endl;
	//指定一个外部的property_map来标定各个节点。biconnected_components这个算法中要用到
	boost::property_map<graph, std::size_t blk_junction::*>::type component = boost::get(&blk_junction::component_id,m_graph);		//表示各个边所在双连通分量的映射
	vector<graph::vertex_descriptor> vap;		//所有关节节点的集合

	std::size_t num_comps; back_insert_iterator<vector<graph::vertex_descriptor>> inter_temp(vap);
	//必须提供一个property map将各个顶点映射到某个整数上才行（boost::vertex_index_map）。具体见BGL的手册
	boost::tie(num_comps,inter_temp)= boost::biconnected_components(m_graph,component,back_inserter(vap),boost::vertex_index_map(boost::get(&blk_group::id,m_graph)));

	vector<vector<int>> v_ap_comps(vap.size ());				//记录每个关节节点在哪些双连通分量中
	vector<vector<int>> v_comps_ap(num_comps);					//记录每个双连通分量包含哪些关节节点
	for (int i(0);i!=vap.size ();++i)
	{
		graph::out_edge_iterator e,e_end;
		boost::tie (e, e_end) = boost::out_edges(vap[i],m_graph);
		for (;e!=e_end;++e)
			add_set (v_ap_comps[i],(int)component[*e]);
		for (int j(0);j!=v_ap_comps[i].size ();++j)
			add_set(v_comps_ap[v_ap_comps[i][j]],i);
	}
	
	//分别记录一下各个双连通分量中包含哪些节点
	vector<vector<graph::vertex_descriptor>> vvv(num_comps);

	graph::edge_iterator e, e_end;
	boost::tie (e,e_end) = boost::edges(m_graph);
	for (;e!=e_end;++e)
	{
		int n( (int)component[*e]);
		add_set(vvv[n], boost::target(*e, m_graph));
		add_set(vvv[n], boost::source(*e, m_graph));
	}

	vector<int> v_flag (num_comps);		//各个双连通分量应该怎么合并
	for (int i(0);i!=num_comps;++i)
	{
		if (is_nested(vvv[i],vap))		
			v_flag[i] = i;			//这个双连通分量需要合并
		else
			v_flag[i] = -1;			//这个双连通分量不需要合并
	}

	//如果两个应该合并的两个双连通图共享节点，那么该节点肯定是关节节点，并且这两个双连通图应该一起合并
	for (int i(0);i!=num_comps;++i)
	{
		if (v_flag[i]<0) continue;		//这个双连通组不用合并
		for (int j(0);j!=v_comps_ap[i].size ();++j)		//对于这个双连通图中的所有关节节点
		{
			for (int k(0);k!=v_ap_comps[v_comps_ap[i][j]].size ();++k)		//对于这个关节节点所在的所有双连通分量
			{
				int temp_i (v_ap_comps[v_comps_ap[i][j]][k]);		//这个双连通分量的索引
				if (temp_i == i)
					continue;
				if (v_flag[temp_i] < 0)
					continue;		//这个双连通图不用合并，直接跳过
				//将数组v_flag中的数等于v_flag[i]的全部改成v_flag[temp_i]
				int tar (v_flag[i]), src(v_flag[temp_i]);
				if (tar == src) continue;
				for (int l(0);l!=num_comps;++l)
					if (v_flag[l] == tar) v_flag[l] = src;
			}
		}
	}

	vector<int> v_group;		//记录一共有几个合并组,每个合并组应该有若干个双连通分量组成 
	for (int i(0);i!=num_comps;++i)
	{
		if (v_flag[i] < 0)
			continue;
		else
			add_set(v_group,v_flag[i]);
	}

	vector<vector<graph::vertex_descriptor>> vv_nested_group(v_group.size ());		//记录哪些节点应该被合并
	for (int i(0);i!=v_group.size ();++i)
	{
		for (int j(0);j!=num_comps;++j)
		{
			if (v_flag[j] != v_group[i]) continue;		//找到flag等于v_group[i]的那一组
			for (int k(0);k!=vvv[j].size ();++k)
				add_set(vv_nested_group[i],vvv[j][k]);
		}
	}

	//理论上这里vv_nested_group中的各个向量中的元素不应该有重合的
	for (int i(0);i!=vv_nested_group.size ();++i)
		merge_vertices(vv_nested_group[i]);

	return (int)vv_nested_group.size ();
}

int adj_blk_graph::eliminate_interlocked_blk (const block_system& bs, const int& search_depth, const bool &debug_flag)
{
	if (debug_flag)
		cout<<"eliminate_interlocked_blk()"<<endl;

	if (search_depth != 1)	//目前只支持一个边的深度搜索
		throw logic_error("adj_blk_graph::eliminate_interlocked_blk (const int& search_depth), only support 1 edge by now");
	vector<vector<graph::vertex_descriptor>> vvv;		//记录哪些顶点应该被合并，vvv中包含的每个顶点向量都是应该合并的顶点向量

	//下面代码前提是m_graph中不能有重复的边,并且一个边的首尾顶点不能相同
	graph::edge_iterator e,e_end;
	for (boost::tie(e,e_end) = boost::edges(m_graph);e!=e_end;++e)
	{
		if (!is_interlocked(bs,*e))
			continue;
		graph::vertex_descriptor tar = boost::target (*e, m_graph), src = boost::source (*e, m_graph);

		int tari(-1), srci(-1);		//记录tar和src分别在哪个顶点列表中
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
		if (tari == -1 && srci == -1)	//两个顶点tar和src都没有被添加过
		{
			vvv.push_back (vector<graph::vertex_descriptor>());
			vvv.back ().push_back (tar);
			vvv.back ().push_back (src);
			continue;
		}
		else if (tari == -1)		//有一个已经被添加到vvv中了
		{
			vvv[srci].push_back (tar); continue;
		}
		else if (srci == -1)
		{
			vvv[tari].push_back (src); continue;
		}
		else
		{
			if (srci == tari)		//两个顶点被添加到同一个组中，则直接跳过
				continue;
			//两个都已经被添加了，但是是在不同的组中。需要把这两个组合并起来
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
	//记录每个顶点的稳定状态。在其他块体被固定的状态下
	graph::vertex_iterator v, v_end;
	vector<safety_factor> vsf(bs.vcb.size ());		//m_graph中最多有这么多个顶点。其中有冗余信息，但是为了方便先这么着
	for (boost::tie(v, v_end) = boost::vertices (m_graph); v!=v_end;++v)
	{
		switch (m_graph[*v].type )
		{
		case GROUP_TYPE::FIXED_AND_EXCAV:
		case GROUP_TYPE::FIXED_AND_UNEXPOSED:
			vsf[m_graph[*v].id].type = safety_factor::S_TYPE::FIXED; break;		//这种情况下有冗余信息，不过有就有吧，不在乎那一点内存 
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
			continue;		//其他四个类型都是不可移动的。至少不能trigger一个连锁垮落
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
			continue;		//如果这个块体本身不是关键块体，则肯定不能引发一个连锁垮落
		}
	
		//if (debug_flag)
		//	cout<<"Keyblock "<<endl;
	
		//*v这个顶点会定义一个新的连锁垮落组合,不管之前有没有被包含在其他连锁垮落组合中。
		v_cfalling.push_back (vector<graph::vertex_descriptor>());
		v_cfalling.back ().push_back (*v);
	
		list<graph::vertex_descriptor> v_cand;		//受到影响的顶点
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
				continue;				//固定的节点就不用添加了
			if (m_graph[v_temp].type == GROUP_TYPE::EXCAV)
				if (vsf[m_graph[v_temp].id].type == safety_factor::S_TYPE::SIG_F_SLID ||
					vsf[m_graph[v_temp].id].type == safety_factor::S_TYPE::DOU_F_SLID ||
					vsf[m_graph[v_temp].id].type == safety_factor::S_TYPE::FALLING)
					continue;			//如果这个结点是关键块体，则不用添加了
			if (find(v_cfalling.back ().begin (), v_cfalling.back ().end (), v_temp) != v_cfalling.back ().end ())
				continue;				//这个点已经在这个连锁垮落组合中了。就不再添加了
	
			vector<graph::edge_descriptor> ve;
			graph::out_edge_iterator e,e_end;
			for (boost::tie (e,e_end) = boost::out_edges(v_temp,m_graph);e!=e_end;++e)
				if (find (v_cfalling.back ().begin (), v_cfalling.back ().end (), boost::target(*e, m_graph)) == v_cfalling.back ().end ())
					ve.push_back (*e);		//如果这个边的另一端没有在v_cfalling.back()中。则是一个限定边
			if (is_removable(bs,ve,v_temp,r))
			{
				v_cfalling.back ().push_back (v_temp);		//这个块体会在r的驱动下移动
				//将其所影响到的块体都添加到v_cand中
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
