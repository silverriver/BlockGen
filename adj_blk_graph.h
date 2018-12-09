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

//这个文件中建立了块体的邻接表。
//底层运算用GBL实现，目的是为了表示块体的邻接关系

namespace BI
{

class adj_blk_graph
{
public:
	enum GROUP_TYPE {EXCAV, UNDEXPOSED, FIXED_AND_EXCAV, FIXED_AND_UNEXPOSED};	//FIXED_AND_EXCAV, FIXED_AND_UNEXPOSED都表示固定节点
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
		std::vector<int> vcb;		//这个块体组合包含哪些块体
		GROUP_TYPE type;					//这个块体组合的类型。
		int id;						//表示这个块体组合的ID
		double volume;
	};

	struct blk_junction
	{
		std::vector<int_pair> vbsf_pair;	//这条边对应几个边界多边形,vbsf_pair[i].first这个边界多边形对应source（意思是这个多边形是source所对应复合块体的边界）。
		std::size_t component_id;			//表示这个边所在的双连通分量
		int id;						//表示这个边的ID
	};

public:
	typedef boost::adjacency_list<boost::setS,boost::listS,
		boost::undirectedS,blk_group,blk_junction> graph;
private:

	graph m_graph;

	//合并这个图的两个顶点,将v2的邻接关系转移到v1上，然后将v2删掉.返回合并之后的顶点
	graph::vertex_descriptor merge_vertices (graph::vertex_descriptor &v1, graph::vertex_descriptor &v2);
	graph::vertex_descriptor merge_vertices (const std::vector<graph::vertex_descriptor> &vv);
	std::vector<std::vector<graph::vertex_descriptor> > v_cfalling;		//该数组中每个元素表示一个连锁垮落组合

	//-------------------初始化函数和一些操作函数------------------------------------
public:		
	
	adj_blk_graph (const block_system &bs,const bool &debug_flag = false) {init_graph(bs, debug_flag);}
	adj_blk_graph (const block_system &bs, const fixed_surfaces &fs, const bool &debug_flag = false) {init_graph(bs,debug_flag); update_fixed_vertices (bs,fs,debug_flag);}

	//为每个块体组合初始化一个外表面结构，方便进行输出，其中每个顶点的id域被作为该顶点在vob中的索引,不过暂时没有用到2016-8-18 17:35
	int init_vob (const block_system &bs, std::vector<outer_bound> &vob) const;

	//根据块体的几何特性初始化图
	int init_graph (const block_system &bs,const bool &debug = false);	
	//根据固定面的位置更新固定节点，即原来不是固定的节点可能变成固定的，原来固定的节点可能变成不固定的
	int update_fixed_vertices (const block_system &bs, const fixed_surfaces & fs, const bool &debug_flag = false);

	//通过合并顶点使得包含关系的块体合并为同一个块体
	//首先生成所有双连通分量和关节节。然后合并相应顶点生成各个块体组合
	//返回值为合并了几组双连通分量
	//flag表示是否显示调试信息
	int eliminate_nested_blk (const bool &flag = false);

	//search_depth指的是搜索深度。因为现在没有成熟的算法可以快速找到interlock。只能这样妥协一下
	//调用之后保证对于任意个数小于等于search_depth的连通顶点组合，都保证其不是interlock的。
	//返回值是减少了几个节点
	//目前只支持搜索两个顶点，也就是search_depth=1的情况
	int eliminate_interlocked_blk (const block_system &bs, const int& search_depth, const bool &debug_flag = false);

	//应该是比较重要的函数，用来判断各个连锁垮落分量。
	//每个关键块体对应一个连锁垮落分量。并且连锁垮落分量中不包含其他的关键块体。（所谓连锁垮落分量意思是当给定块体滑落时，原本不能垮落的块体也发生了移动）
	//将识别得到的连锁组合储存在v_cfalling中，每个数组中至少有一个开挖块体。并且每个子数组中块体排列顺序和垮落顺序一致
	int determine_chain_falling_blocks (const block_system &bs, const vector_3 &r,  const bool &debug_flag = false) ;

	//测试用,计算给定id节点的稳定系数
	BI::safety_factor test_sf_of_vertex (const block_system &bs, const int &id, const vector_3 &r) const;

	//----------------------输出函数---------------------------------
public:			
	int output_vertices_type_info (const block_system &bs, std::ofstream &outfile);		//输出块体的类型体积

	//输出某个顶点表示的块体组合，输出时生成这些块体所
	int output_vertex_blk (const block_system &bs, const graph::vertex_descriptor &v,std::ostream &outfile) const;

	//输出用函数,输出vv中表示的一系列函数，但是输出的块体并不产生其外表面结构，而是直接输出外表面多边形
	//注意，这个函数只能用来快速查看块体形状，生成的stl文件可能有外露边界和非流行边界
	int output_cfalling_blks_quick (const block_system &bs, std::ostream &os, const std::vector<graph::vertex_descriptor> &vv) const;
	
	//输出所有type类型的块体.输出的块体并不生产其外表面结构，而是直接输出外表面多边形，因此由可能有外露边界和非流型边界
	//注意，这个函数只能用来快速查看块体形状，生成的stl文件可能有外露边界和非流行边界
	int output_blk_group_quick (std::ostream & os, const block_system &bs, const GROUP_TYPE &type) const;

	//输出指定连锁垮落块体，输出时构建每个块体的外表面结构
	int output_cfalling_blks (const block_system &bs, std::ostream &os, const std::vector<graph::vertex_descriptor> &vv) const;
	int output_cfalling_blks (const block_system &bs, std::ostream &os, const int &count) const {return output_cfalling_blks(bs,os,v_cfalling[count]);}
	int output_cfalling_blks (const block_system &bs, std::ostream &os) const
	{
		for (int i(0);i!=v_cfalling.size ();++i)
			output_cfalling_blks(bs,os,v_cfalling[i]);
		return 0;
	}


	//输出CSV文件表示这个图。其中et指示输出点数据还是边数据
	int output_CSV (std::ostream &outfile, const ELE_TYPE& et) const;

	//计算和输出各个块体的稳定系数
	int output_safety_factor (std::ostream &outfile, const block_system &bs) const;

	//输出这个连通图包含的统计信息，比如多少节点，多少固定节点或开挖节点等
	void output_node_and_edge_info (std::ostream &outfile) const;

	//输出连锁块体统计信息，每个连锁块体的信息都输出,其中bs是用来计算块体组合体积的
	void output_chain_falling_info_CSV (const block_system &bs, std::ostream &outfile) const;

	////输出每个顶点所表示块体集合的信息，包括id，类型，所包含块体，体积，稳定性情况等
	////和output_CSV的输出信息有重复
	//void output_blks_info (const block_system &bs, std::ostream &outfile) const;

	//输出所有可移动块体
	void ouput_removable_blks (const block_system &bs, const vector_3 &r,std::ostream &os) const;
	//----------------获取连锁块体信息的部分----------------
public:			

	size_t get_cfalling_count () const {return v_cfalling.size ();};
	std::vector<graph::vertex_descriptor>& get_cfalling_group (int index) {return v_cfalling[index];}
	const std::vector<graph::vertex_descriptor>& get_cfalling_group (int index) const {return v_cfalling[index];}

private:

	//判断vv中记录的节点是否是内嵌的，如果是，则返回true
	//vap记录的是所有的关节节点
	//precondition: vv里面的节点不相同，并且来自同一个双连通分量
	//flag表示是否显示调制信息
	bool is_nested (const std::vector<graph::vertex_descriptor> & vv, const std::vector <graph::vertex_descriptor> & vap, const bool &debug_flag = false) const;

	//判断节点v在ve的限定下能否沿r方向的外力移动,要求v必须在ve这些边上
	//注意，在block_system中也有类似的函数，但是没有给定外力方向。这里没有用block_system中的实现。
	bool is_removable (const block_system &bs, const std::vector<graph::edge_descriptor> &ve, const graph::vertex_descriptor &v, const vector_3 &r, const bool & debug_flag = false) const;

	//判断给定的节点是否是固定节点
	bool is_fixed (const block_system &bs, const fixed_surfaces &fs, const graph::vertex_descriptor &v) const;

	//计算稳定系数.ve表示v只在ve这些边的限定下。如果没有给定ve，则表示v在其所有边的限定下
	safety_factor cal_sf (const block_system& bs, const graph::vertex_descriptor&v, const vector_3 &r, const bool& debug_flag = false) const;
	safety_factor cal_sf (const block_system& bs, const std::vector<graph::edge_descriptor> &ve, const graph::vertex_descriptor &v, const vector_3 &r, const bool & debug_flag = false) const;

	//判断v1和v2是否是锁定的，目前只能判断两个顶点
	//注意，在block_system中也有类似的函数，这里没有使用block_system中的实现。因为目前为止两个实现都是未完成阶段 2016-6-16 10:49
	bool is_interlocked (const block_system &bs, const graph::vertex_descriptor & v1, const graph::vertex_descriptor & v2) const;
	bool is_interlocked (const block_system &bs, const graph::edge_descriptor  &e) const;

	mbox_3 build_mbox_3(const block_system &bs, const graph::vertex_descriptor &v) const;
};


}