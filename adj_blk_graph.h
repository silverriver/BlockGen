#pragma once
#include "boost\graph\adjacency_list.hpp"
#include "boost\graph\graph_traits.hpp"
#include "identification.h"
#include "safety_factor.h"
#include "geometry.h"
#include "fixed_surfaces.h"
#include "outer_poly.h"
#include <vector>
#include <set>
#include <utility>
#include <ostream>
#include <fstream>

//����ļ��н����˿�����ڽӱ�
//�ײ�������GBLʵ�֣�Ŀ����Ϊ�˱�ʾ������ڽӹ�ϵ

namespace BI
{

class adj_blk_graph
{
public:
	enum GROUP_TYPE {EXCAV, UNDEXPOSED, FIXED_AND_EXCAV, FIXED_AND_UNEXPOSED};	//FIXED_AND_EXCAV, FIXED_AND_UNEXPOSED����ʾ�̶��ڵ�
	enum ELE_TYPE {Edge, Node};
private:
	typedef std::pair<int,int> int_pair;
	

	class int_pair_comp
	{
	public:
		bool operator () (const int_pair &a, const int_pair &b) const
		{
			if (a.first == b.first && a.second == b.second)
				return true;
			if (a.second == b.first && a.first == b.second )
				return true;
			return false;
		}
	};

	struct blk_group
	{
		std::vector<int> vcb;		//���������ϰ�����Щ����
		GROUP_TYPE type;					//���������ϵ����͡�
		int id;						//��ʾ���������ϵ�ID
		double volume;
	};

	struct blk_junction
	{
		std::vector<int_pair> vbsf_pair;	//�����߶�Ӧ�����߽�����,vbsf_pair[i].first����߽����ζ�Ӧsource����˼������������source����Ӧ���Ͽ���ı߽磩��
		std::size_t component_id;			//��ʾ��������ڵ�˫��ͨ����
		int id;						//��ʾ����ߵ�ID
	};

public:
	typedef boost::adjacency_list<boost::setS,boost::listS,
		boost::undirectedS,blk_group,blk_junction> graph;
private:

	graph m_graph;

	//�ϲ����ͼ����������,��v2���ڽӹ�ϵת�Ƶ�v1�ϣ�Ȼ��v2ɾ��.���غϲ�֮��Ķ���
	graph::vertex_descriptor merge_vertices (graph::vertex_descriptor &v1, graph::vertex_descriptor &v2);
	graph::vertex_descriptor merge_vertices (const std::vector<graph::vertex_descriptor> &vv);
	std::vector<std::vector<graph::vertex_descriptor> > v_cfalling;		//��������ÿ��Ԫ�ر�ʾһ�������������

	//-------------------��ʼ��������һЩ��������------------------------------------
public:		
	
	adj_blk_graph (const block_system &bs,const bool &debug_flag = false) {init_graph(bs, debug_flag);}
	adj_blk_graph (const block_system &bs, const fixed_surfaces &fs, const bool &debug_flag = false) {init_graph(bs,debug_flag); update_fixed_vertices (bs,fs,debug_flag);}

	//Ϊÿ��������ϳ�ʼ��һ�������ṹ������������������ÿ�������id����Ϊ�ö�����vob�е�����,������ʱû���õ�2016-8-18 17:35
	int init_vob (const block_system &bs, std::vector<outer_bound> &vob) const;

	//���ݿ���ļ������Գ�ʼ��ͼ
	int init_graph (const block_system &bs,const bool &debug = false);	
	//���ݹ̶����λ�ø��¹̶��ڵ㣬��ԭ�����ǹ̶��Ľڵ���ܱ�ɹ̶��ģ�ԭ���̶��Ľڵ���ܱ�ɲ��̶���
	int update_fixed_vertices (const block_system &bs, const fixed_surfaces & fs, const bool &debug_flag = false);

	//ͨ���ϲ�����ʹ�ð�����ϵ�Ŀ���ϲ�Ϊͬһ������
	//������������˫��ͨ�����͹ؽڽڡ�Ȼ��ϲ���Ӧ�������ɸ����������
	//����ֵΪ�ϲ��˼���˫��ͨ����
	//flag��ʾ�Ƿ���ʾ������Ϣ
	int eliminate_nested_blk (const bool &flag = false);

	//search_depthָ����������ȡ���Ϊ����û�г�����㷨���Կ����ҵ�interlock��ֻ��������Эһ��
	//����֮��֤�����������С�ڵ���search_depth����ͨ������ϣ�����֤�䲻��interlock�ġ�
	//����ֵ�Ǽ����˼����ڵ�
	//Ŀǰֻ֧�������������㣬Ҳ����search_depth=1�����
	int eliminate_interlocked_blk (const block_system &bs, const int& search_depth, const bool &debug_flag = false);

	//Ӧ���ǱȽ���Ҫ�ĺ����������жϸ����������������
	//ÿ���ؼ������Ӧһ�������������������������������в����������Ĺؼ����塣����ν�������������˼�ǵ��������廬��ʱ��ԭ�����ܿ���Ŀ���Ҳ�������ƶ���
	//��ʶ��õ���������ϴ�����v_cfalling�У�ÿ��������������һ�����ڿ��塣����ÿ���������п�������˳��Ϳ���˳��һ��
	int determine_chain_falling_blocks (const block_system &bs, const vector_3 &r,  const bool &debug_flag = false) ;

	//������,�������id�ڵ���ȶ�ϵ��
	BI::safety_factor test_sf_of_vertex (const block_system &bs, const int &id, const vector_3 &r) const;

	//----------------------�������---------------------------------
public:			
	int output_vertices_type_info (const block_system &bs, std::ofstream &outfile);		//���������������

	//���ĳ�������ʾ�Ŀ�����ϣ����ʱ������Щ������
	int output_vertex_blk (const block_system &bs, const graph::vertex_descriptor &v,std::ostream &outfile) const;

	//����ú���,���vv�б�ʾ��һϵ�к�������������Ŀ��岢�������������ṹ������ֱ��������������
	//ע�⣬�������ֻ���������ٲ鿴������״�����ɵ�stl�ļ���������¶�߽�ͷ����б߽�
	int output_cfalling_blks_quick (const block_system &bs, std::ostream &os, const std::vector<graph::vertex_descriptor> &vv) const;
	
	//�������type���͵Ŀ���.����Ŀ��岢�������������ṹ������ֱ�������������Σ�����ɿ�������¶�߽�ͷ����ͱ߽�
	//ע�⣬�������ֻ���������ٲ鿴������״�����ɵ�stl�ļ���������¶�߽�ͷ����б߽�
	int output_blk_group_quick (std::ostream & os, const block_system &bs, const GROUP_TYPE &type) const;

	//���ָ������������壬���ʱ����ÿ������������ṹ
	int output_cfalling_blks (const block_system &bs, std::ostream &os, const std::vector<graph::vertex_descriptor> &vv) const;
	int output_cfalling_blks (const block_system &bs, std::ostream &os, const int &count) const {return output_cfalling_blks(bs,os,v_cfalling[count]);}
	int output_cfalling_blks (const block_system &bs, std::ostream &os) const
	{
		for (int i(0);i!=v_cfalling.size ();++i)
			output_cfalling_blks(bs,os,v_cfalling[i]);
		return 0;
	}


	//���CSV�ļ���ʾ���ͼ������etָʾ��������ݻ��Ǳ�����
	int output_CSV (std::ostream &outfile, const ELE_TYPE& et) const;

	//������������������ȶ�ϵ��
	int output_safety_factor (std::ostream &outfile, const block_system &bs) const;

	//��������ͨͼ������ͳ����Ϣ��������ٽڵ㣬���ٹ̶��ڵ���ڽڵ��
	void output_node_and_edge_info (std::ostream &outfile) const;

	//�����������ͳ����Ϣ��ÿ�������������Ϣ�����,����bs���������������������
	void output_chain_falling_info_CSV (const block_system &bs, std::ostream &outfile) const;

	////���ÿ����������ʾ���弯�ϵ���Ϣ������id�����ͣ����������壬������ȶ��������
	////��output_CSV�������Ϣ���ظ�
	//void output_blks_info (const block_system &bs, std::ostream &outfile) const;

	//������п��ƶ�����
	void ouput_removable_blks (const block_system &bs, const vector_3 &r,std::ostream &os) const;
	//----------------��ȡ����������Ϣ�Ĳ���----------------
public:			

	size_t get_cfalling_count () const {return v_cfalling.size ();};
	std::vector<graph::vertex_descriptor>& get_cfalling_group (int index) {return v_cfalling[index];}
	const std::vector<graph::vertex_descriptor>& get_cfalling_group (int index) const {return v_cfalling[index];}

private:

	//�ж�vv�м�¼�Ľڵ��Ƿ�����Ƕ�ģ�����ǣ��򷵻�true
	//vap��¼�������еĹؽڽڵ�
	//precondition: vv����Ľڵ㲻��ͬ����������ͬһ��˫��ͨ����
	//flag��ʾ�Ƿ���ʾ������Ϣ
	bool is_nested (const std::vector<graph::vertex_descriptor> & vv, const std::vector <graph::vertex_descriptor> & vap, const bool &debug_flag = false) const;

	//�жϽڵ�v��ve���޶����ܷ���r����������ƶ�,Ҫ��v������ve��Щ����
	//ע�⣬��block_system��Ҳ�����Ƶĺ���������û�и���������������û����block_system�е�ʵ�֡�
	bool is_removable (const block_system &bs, const std::vector<graph::edge_descriptor> &ve, const graph::vertex_descriptor &v, const vector_3 &r, const bool & debug_flag = false) const;

	//�жϸ����Ľڵ��Ƿ��ǹ̶��ڵ�
	bool is_fixed (const block_system &bs, const fixed_surfaces &fs, const graph::vertex_descriptor &v) const;

	//�����ȶ�ϵ��.ve��ʾvֻ��ve��Щ�ߵ��޶��¡����û�и���ve�����ʾv�������бߵ��޶���
	safety_factor cal_sf (const block_system& bs, const graph::vertex_descriptor&v, const vector_3 &r, const bool& debug_flag = false) const;
	safety_factor cal_sf (const block_system& bs, const std::vector<graph::edge_descriptor> &ve, const graph::vertex_descriptor &v, const vector_3 &r, const bool & debug_flag = false) const;

	//�ж�v1��v2�Ƿ��������ģ�Ŀǰֻ���ж���������
	//ע�⣬��block_system��Ҳ�����Ƶĺ���������û��ʹ��block_system�е�ʵ�֡���ΪĿǰΪֹ����ʵ�ֶ���δ��ɽ׶� 2016-6-16 10:49
	bool is_interlocked (const block_system &bs, const graph::vertex_descriptor & v1, const graph::vertex_descriptor & v2) const;
	bool is_interlocked (const block_system &bs, const graph::edge_descriptor  &e) const;

	mbox_3 build_mbox_3(const block_system &bs, const graph::vertex_descriptor &v) const;
};


}