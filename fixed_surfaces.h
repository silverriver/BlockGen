#pragma once
#include "geometry.h"
#include "identification.h"
#include <istream>
#include <ostream>
#include <vector>

//这个头文件主要表示一些固定面。
//每个固定面是一个三维多边形，目前先用三角形表示。
//每个和该三角形相交的块体都应该被标记为固定块体
namespace BI
{
class fixed_surfaces
{
	std::vector<triangle_3> vtri;
public:

	fixed_surfaces () {}
	void add_tri (const triangle_3 &tri) {vtri.push_back(tri);}
	int STL_parser (std::istream& infile);		//语法解析器喽
	bool do_intersect (const std::vector<point_3> &vpo, const std::vector<plane_3> &vpl, const eblock &eb) const;
	bool do_intersect (const std::vector<point_3> &vpo, const std::vector<plane_3> &vpl, const std::vector<eblock> &veb, const cblock &cb) const;
	bool do_intersect (const block_system &bs, const eblock &eb) const;
	bool do_intersect (const block_system &bs, const cblock &cb) const;
	
	void output_triangles_line (std::ostream &outfile) const
	{
		for (int i(0);i!=vtri.size();++i)
		{
			out_line_3(outfile,vtri[i].vertex(0),vtri[i].vertex(1));
			out_line_3(outfile,vtri[i].vertex(2),vtri[i].vertex(1));
			out_line_3(outfile,vtri[i].vertex(0),vtri[i].vertex(2));
		}
	}
};

}