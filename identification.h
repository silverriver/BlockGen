#pragma once
#include <vector>
#include "geometry.h"
#include "GBUtility.h"
#include <istream>
#include <fstream>
#include <utility>
#include <string>
namespace BI
{

//这个类主要处理块体的几何特性，不考虑诸如固定面之类的力学特性
class block_system
{
	block_system (const block_system & bs) {}
	block_system& operator= (const block_system& bs) {}
public:
	std::vector<point_3>	vpo;
	std::vector<plane_3>	vpl;
	std::vector<disc_frac>	vdisc;
	std::vector<poly_frac>	vpoly;		//多边形裂隙的ID，也就是frac_id排在圆盘裂隙的ID后面
	std::vector<eblock>		veb;
	std::vector<cblock>		vcb;
	double cutting_time;
	double shrk_time;

	block_system ():cutting_time(-1.0),shrk_time(-1.0){}

	void clear ()
	{vpo.clear(); vpl.clear(); vdisc.clear(); vpoly.clear(); veb.clear(); vcb.clear();}

	//清空这个blocksystem的内容，然后用d来初始化,返回1表示因为误差问题初始化失败
	int init_domain (const GB_domain &d, const double err=0.0000001);
	int init_check () const;

	//判断一个单元块体是否是凸体，并且点都exactly在指定平面上，并且外法线方向是否正确。true表示是凸体并且外法线方向正确
	bool eb_convex (const int &i) const;
	
	//判断一个块体是否封闭，也就是查看其所有边界是否有开放边界，返回true表示没有开放边界
	bool eb_sealed (const int &i) const;

	FT eb_volume (const int &ebi) const;

	FT cb_volume (const int &cbi) const;

	//返回块体体积的估计值，因为FT在连续累加的时候会出错。
	double eb_volume_est (const int &ebi) const;	
	double cb_volume_est (const int &cbi) const;

	//判断p在块体的什么位置，前提是保证ei是凸体的情况
	Oriented_side eb_p_location (const int& ei, const point_3 &p) const;

	//放回两个eblock的共用平面，意思是两个块体如果有一个共用平面，那么就返回平面的索引，如果没有，则返回一个小于0的数
	int common_plane (const eblock& eb1, const eblock&eb2) const;

	//处理圆盘裂隙，读取disc的信息，并初始化，返回0表示成功
	int disc_frac_parser (std::istream &infile, const FT &dx);

	//处理多边形裂隙，读取poly的信息，并初始化，返回0表示成功
	int poly_frac_parser (std::istream &infile);

	//处理一个stl文件，将其中的每一个三角形面片看做是一个多边形裂隙,err表示读取时的容许误差
	//新添加裂隙的frac_id排列在vpoly.back().frac_id的后面
	//返回0表示解析成功
	int stl_parser (std::istream &infile, const double &err = 0.00000001);

	//添加一个多边形裂隙。输入时需要用户保证pf不包含退化的外表面和内表面边界
	int add_poly_frac (const poly_frac &pf);

	//添加一个圆盘裂隙。注意，应该先把所有的圆盘裂隙添加完再添加多边形裂隙。不然frac_id就乱了.normal是表示df所在平面的法向向量
	//所谓的添加裂隙也就是把裂隙添加到相应的数组中，并且给frac_id和plane_id赋值
	int add_disc_frac (const disc_frac &df,const vector_3 &normal);

	//初始化一个多面体形状的研究区域,需要指定两个角点
	int add_rect_domain (const FT &x1, const FT &y1, const FT &z1, 
		const FT &x2, const FT &y2, const FT &z2);

	//试图用平面vpl[pli]去切割veb[targeti]，新生成的块体被添加到veb中，并且其中一个新生成的块体替换veb[targeti]的位置
	//有可能切割也有可能不切割。如果切割了，新生成的面裂隙id被赋值为frac_id，新生成的两个eblock中的复合块体索引被赋值为cblock_index和c_block_index+1
	//返回0表示切割成功。
	//返回1表示没法切割
	int plane_cut_eblock (const int &targeti, const int &pli, const int &frac_id, const int &cblock_index);

	//将平面缩小到原来的尺寸，pli表示需要缩小的平面的索引
	//会查找平面上的裂隙，并且查找平面两侧的eblock，并且计算各个块体的面的交集。
	//返回值表示合并的单元块体的数量。
	//如果eblock[i]的某个面的frac_id<=-4，则该面不会被合并，否则，该面则有可能被合并
	//假设各个veb[i].cblock_index都有各自的意义，已经被初始化
	int shrk_plane (const int &pli);

	//假设裂隙已经都设置好了，并且研究区域也已经读取完毕了，就剩下识别块体了
	//os是用来输出调试信息的，目前什么也没有输出
	int identify_block (std::ostream & os = std::ofstream());

	//输出stl文件，只是简单的stl文件，可能有非流行边界和外露边界，返回0表示输出成功,x_, y_, z_表明的是偏移量
	int output_eb_stl (const int &ebi, std::ostream &outfile, const FT &x_ = 0, const FT &y_ = 0, const FT &z_ = 0) const;
	int output_cb_stl (const int &cbi, std::ostream &outfile, const FT &x_ = 0, const FT &y_ = 0, const FT &z_ = 0, const std::string &cb_name = "") const;

	//输出stl文件，组成cbi的单元块体都以exploded view给出
	int output_explode_cb_stl (const int &cbi, std::ostream &outfile, const double &ratio = 0.1) const;

	//目的是得到这个块体系统的中点，其实就是把vpo中的点加起来平均一下。
	//没有用到严格运算，因为用不到。关键是FT的实现中有bug
	point_3 get_rand_p() const;

	//输出这个bs的一些基本信息
	void output_info (std::ostream &outfile) const;

//private:
	//这部分内容是为了分析各个块体的稳定性而准备的。
	//并不打算给出外部接口。2016-5-28 11:00
	
	enum FACE_TYPE {EXCAV, FRAC,};
	struct bound_surface 
	{
		//表示一个复合块体的一个外表面
		polygon_with_holes_2 poly;		//外表面多边形
		int pli;						//这个外表面多边形所在平面在vpl中的索引
		int cbi;						//这个外表面属于哪个复合块体
		bool outer_normal;				//true：这个外表面的外法线方向与pli这个平面的正向相同
		std::vector <int> vadj_bsf_i;	//所有与这个外表面多边形相邻的其他外表面多边形索引（这些多边形肯定有相同的pli，并且与该外表面多边形内部相交）。
		std::vector <int> vadj_cb_i;	//所有与这个外表面多边形相邻的复合块体
		FACE_TYPE type;						//记录这个面的类型，比如是否是固定面，或者是哪个裂隙面生成的。（先预留一个位置，暂时没有用到这个内容。2016-5-28 11:23）
	};

	std::vector<bound_surface> vbsf;				//这个块体系统中所有的外表面多边形
	std::vector<std::vector<int>> vcb_bsf;			//每个复合块体都对应一个外表面多边形的集合
	std::vector <std::vector <int>> vpl_bsf;		//每个平面对应一个外表面多边形的集合

public:
	//初始化，为分析块体系统稳定系数做准备。（为上述数据赋值）
	int init_stabl_alys (const bool& debug);
	
	//不完善，但是保留待用，没准以后会用的到2016-6-16 11:46
	//bool is_interlock (const int &cbi, const int &cbj) const;
	//bool is_interlock (const std::vector<int> &vcbi) const;		//判断一组块体是否是锁定的，这样的判断过程有指数复杂度

	//没有用到过这一部分来计算稳定性，因为在正常情况下，在一个bs中还没有合并nested block。所以导致直接用这些函数分析得到的稳定性不一定是准确的。2016-6-16 11:09
	bool is_removable (const int &cbi) const;
	bool is_removable (const std::vector<int> &vcbi) const;
	//判断被vfcbi限定住的块体集合vcbi是否可移动
	bool is_removable (const std::vector<int> &vcbi, const std::vector <int> &vfcbi) const;


	//index>=0 表示输出某个复合块体的所有外表面。 index<0则-index+1表示输出某个平面上所有外表面，此时flag才有意义（1:输出该平面正向块体外表面，2：输出该平面反向块体外表面， 3：都输出）
	//single表示是否生成独立的stl文件。只有index>=0时 这个值才会被用到，single=true表示每个块体成为一个独立的stl文件，否则只输出其表面三角形
	int test_output_STL (std::ostream&outfile, const int&index, const int &flag = 3, bool single = true,const std::string &stl_name = "") const;

private:
	//过滤不能生成块体的裂隙，其中vbsf表示model_domain的外表面多边形的集合，该函数会对vdisc和vpoly中的裂隙进行修改
	//该函数只能用在只有圆盘裂隙的情况下
	//以后待定的任务2016-8-20 13:59
	int filter_disc_fractures (std::vector<bound_surface> &vbsf);

public:
	friend std::ostream& operator<< (std::ostream& os, const block_system &bs);
	friend std::istream& operator>> (std::istream& is, block_system &bs);
	friend std::ostream& operator<< (std::ostream &os, const block_system::bound_surface& bsf);
	friend std::istream& operator>> (std::istream& is, block_system::bound_surface& bsf);
	friend std::ostream& operator<< (std::ostream& os, const block_system::FACE_TYPE &ft);
	friend std::istream& operator>> (std::istream& is, block_system::FACE_TYPE &ft);
	friend bool operator== (const block_system& bs1, const block_system& bs2);
	friend class adj_blk_graph;		//防火防盗防友元
};

//利用vis中的点的信息，识别一个首尾相接的序列。假设vis[i].first < vis[i].second
//返回0表示构建成功，返回其他值表示构建失败
//假设vis中的元素不重复
int construct_loop (const std::vector<std::pair<int,int>> &vis, std::vector<int> &res);

inline bool operator== (const block_system& bs1, const block_system& bs2)
{
	return (bs1.cutting_time == bs2.cutting_time &&
		bs1.shrk_time == bs2.shrk_time &&
		bs1.vcb == bs2.vcb &&
		bs1.vdisc == bs2.vdisc &&
		bs1.veb == bs2.veb &&
		bs1.vpl == bs2.vpl &&
		bs1.vpo == bs2.vpo &&
		bs1.vpoly == bs2.vpoly);
}


//判断va中的每个元素是否都包含在vb中
inline bool is_contained (std::vector<int> va, std::vector<int> vb)
{
	for (int i(0);i!= (int) va.size ();++i)
	{
		int count(0);
		for (; count!= (int) vb.size ();++count)
		{
			if (va[i] == vb[count])
				break;
		}
		if (count == (int) vb.size ())
			return false;		//这个va[i]没有被包含在vb中
	}
	return true;
}

//判断va中是否存在某个元素包含在vb中
inline bool has_commonelem (std::vector<int> va, std::vector<int> vb)
{
	for (int i(0);i!=va.size ();++i)
	{
		for (int j(0);j!=vb.size ();++j)
			if (va[i] == vb[j])
				return true;
	}
	return false;
}

}