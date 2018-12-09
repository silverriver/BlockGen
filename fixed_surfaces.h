#pragma once
#include "geometry.h"
#include "identification.h"
#include <istream>
#include <ostream>
#include <vector>

//���ͷ�ļ���Ҫ��ʾһЩ�̶��档
//ÿ���̶�����һ����ά����Σ�Ŀǰ���������α�ʾ��
//ÿ���͸��������ཻ�Ŀ��嶼Ӧ�ñ����Ϊ�̶�����
namespace BI
{
class fixed_surfaces
{
	std::vector<triangle_3> vtri;
public:

	fixed_surfaces () {}
	void add_tri (const triangle_3 &tri) {vtri.push_back(tri);}
	int STL_parser (std::istream& infile);		//�﷨�������
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