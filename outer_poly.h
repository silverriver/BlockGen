#pragma once
#include "geometry.h"
#include "identification.h"
#include <vector>
#include <ostream>
#include <string>

namespace BI
{
struct angle_edge;
struct safety_factor;
class entrance_block;

class outer_bound
{
public:
	friend int cal_safety_factor (const outer_bound &ob, safety_factor & sf);
	friend int outer_bound2angle_edge (const block_system & bs, const cblock & cb,const outer_bound &ob, angle_edge &ae);
	friend class entrance_block;
//private:  DEBUG
	std::vector<polygon_with_holes_2> vpoly;
	std::vector<plane_3> vplane;		//要求对于vplane[i]. 任意其他平面都不等于vplane[i]或vplane[i].opposite()
	std::vector<bool> outer_normal;	//true：这个多边形的外法线方向与vplane[i]平面的正向相同
public:
	outer_bound () {}
	outer_bound (const block_system &bs,const int &cbi,const bool&debug_info=false)
	{construct_outer_bound(bs,cbi,debug_info);}
	outer_bound (const block_system &bs, const std::vector<int>& vebi,const bool &debug_info = false)
	{construct_outer_bound(bs,vebi,debug_info);}

	
	//清空原有信息，构建边界。返回0表示构建成功，事实上这个函数除了抛出异常就是返回0
	//输入函数中vebi指的是需要构建外表面边界的单元块体的索引向量
	int construct_outer_bound (const block_system &bs, const std::vector <int> &vebi, const bool& debug_info = false);
	//下面这个函数实现的时候是直接调用上面这个函数的
	int construct_outer_bound (const block_system &bs, const int &cbi, const bool& debug_info = false);
	//生成若干个cb的外表面
	int construct_outer_bound_using_vcb (const block_system &bs, const std::vector <int> &vcbi, const bool& debug_info = false);

	//输出stl文件，返回0表示输出成功。输出的文件只是近似而已
	int output_STL (std::ostream &outfile, const FT& tranx = 0.0, const FT& trany = 0.0, const FT& tranz = 0.0, const std::string& ob_name="") const;
	int output_exp_STL (std::ostream &outfile, const double &ratio = 0.2, const FT& tranx = 0.0, const FT& trany = 0.0, const FT& tranz = 0.0, const std::string& ob_name="") const;		//ob_name没有用到
	int output_CAD_script (std::ostream &outfile,const bool& split_ploy = false, const bool &output_normal=false) const;

	//将生成的外表面结构输出成BlockAni_1.3.3所显示的结构
	int output_BlockAni_format (std::ostream &outfile, const int &id) const;

	//在原来外表面边界的基础上添加一些单元块体。具体算法比join这个函数高效一点，节省了一些内存复制操作
	int join_veb (const block_system &bs, const std::vector <int> &vebi, const bool& debug_info = false);
	//这个函数的目的是在已经建立好的bs的基础上添加另外一个复合块体。实现方法是直接调用函数joint_veb
	int join_cb (const block_system &bs, const int cbi, const bool& debug_info = false);

	//求这个块体的最小包围盒的中心点
	point_3 get_center () const;

	int get_polycount () const
	{return (int)vpoly.size ();}

	//将给定外表面边界与该边界结合起来
	int join (const outer_bound &ob, const bool &debug_info = false);

private:
	//私有函数，用来给vpoly[poly_count]这个多边形加上其他多边形的角点
	//vpoint表示vpoly中各个多边形的角点
	//p_count表示vpoly中各个多边形焦点的数量
	//vbox表示vpoly中各个多边形的最小包围盒
	//上面这些信息都可以用的时候生成，但是为了省时间，需要在调用这个函数的时候手动cache这些信息
	//返回值是在这个多边形上添加了多少点
	int add_coner_point (const int &poly_n,  point_3 **vpoint,  int *p_count,  mbox_3 *vbox);

	//这个函数用来给给定多边形添加点，其中original中的点会再次转化到三维空间中，虽然这一步可以优化掉，但是会涉及更复杂的指针接口，为了图省事，暂时忽略2015-10-25 19:37
	//同样为了省事，表示潜在添加点时，这里的接口直接用的是点的坐标，而不是点的索引。
	//vconer表示需要添加的点，这些点应该都在plane这个平面上，并且vconer中的点应该相互不相同
	//res表示新生成的多边形。返回值表示在original的基础上新添加点的数量
	int add_point_on_poly (const polygon_2 &original, const plane_3& plane, const std::vector<point_3> & vconer, polygon_2 &res);

	//这个函数用来整合生成的外表面边界，也就是将各个边的角点添加在其他外表面多边形上
	//在调用这个函数时最好先简化一个各个边界多边形simplify_polyon_2
	int integrate_boundary ();

	//flag 表示用什么输出格式，false表示不把这个多边形输出成单个solid。true表示把这个多边形输出成单个solid
	int output_STL_poly (std::ostream &outfile, const int &pi, const bool& flag = false, const FT& tranx = 0.0, const FT& trany = 0.0, const FT& tranz = 0.0) const;

	//求某个多边形的中心点，也就是各个顶点的平均值
	point_3 get_poly_center (const int &pi) const;

public:
	friend std::ostream& operator<< (std::ostream &os, const outer_bound& ob);
	friend std::istream& operator>> (std::istream &is, outer_bound &ob);
	friend bool operator== (const outer_bound&ob1, const outer_bound& ob2);
};

//给定的tar是需要排序的值(所有值大约0)，index作为返回值给出，从最小元素索引一直到最大元素索引
void sort_index (const std::vector<FT> &tar, std::vector<int> &index);

inline bool operator== (const outer_bound&ob1, const outer_bound& ob2)
{
	return (ob1.outer_normal == ob2.outer_normal &&
		ob1.vplane == ob2.vplane &&
		ob1.vpoly == ob2.vpoly);
}

}



