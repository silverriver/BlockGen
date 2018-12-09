#pragma once
#include "geometry.h"
#include "identification.h"
#include <vector>
#include <ostream>
#include <string>

namespace BI
{

class inner_bound
{
public:
	std::vector<polygon_with_holes_2> vpoly;
	std::vector<plane_3> vplane;

	inner_bound (){}
	int construct_inner_bound (const block_system &bs, const std::vector<int> &vebi, const bool& debug_info = false);
	//下面这个函数实现的时候是直接调用上面这个函数的
	int construct_inner_bound (const block_system &bs, const int &cbi, const bool& debug_info = false);

	int output_STL ( std::ostream &outfile, const FT& tranx=0.0, const FT& trany=0.0, const FT& tranz=0.0, const std::string& ob_name="") const;
	int output_exp_STL (std::ostream &outfile, const double &ratio = 0.2, const FT& tranx = 0.0, const FT& trany = 0.0, const FT& tranz = 0.0, const std::string& ob_name="") const;		//ob_name没有用到

public:
	friend std::ostream& operator<< (std::ostream &os, const inner_bound& ib);
	friend std::istream& operator>> (std::istream &is, inner_bound &ib);
	friend bool operator== (const inner_bound&ib1, const inner_bound& ib2);
private:
	//flag 表示用什么输出格式，false表示不把这个多边形输出成单个solid。true表示把这个多边形输出成单个solid
	int output_STL_poly (std::ostream &outfile, const int &pi, const bool& flag = false, const FT& tranx = 0.0, const FT& trany = 0.0, const FT& tranz = 0.0) const;
	
	//求某个多边形的中心点，也就是各个顶点的平均值
	point_3 get_poly_center (const int &pi) const;
};

inline bool operator== (const inner_bound&ib1, const inner_bound& ib2)
{
	return (ib1.vplane == ib2.vplane &&
		ib1.vpoly == ib2.vpoly);
}

}