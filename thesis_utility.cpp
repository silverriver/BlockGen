////used for thesis modification
////2017-5-27
//
//#include "identification.h"
//#include "geometry.h"
//#include "GBUtility.h"
//#include "outer_poly.h"
//#include "curve_frac.h"
//#include "entrance_block.h"
//#include "serialize.h"
//#include "safety_factor.h"
//#include "inner_poly.h"
//#include "adj_blk_graph.h"
//#include "fixed_surfaces.h"
//#include <iostream>
//#include <fstream>
//#include <string>
//#include <sstream>
//#include <utility>
//#include <algorithm>
//#include <CGAL\Timer.h>
//#include <CGAL\Memory_sizer.h>
//#include <cstdlib>
//#include <map>
//
//using namespace std;
//using BI::operator<<;
//using BI::operator>>;
//
////简单粗暴，没时间了
//bool eb_in_sphere (const BI::point_3& center, const BI::FT& radius, const BI::block_system& bs, const int& eb_index)
//{
//	BI::FT r_square(radius*radius);
//	for (int i(0);i!=bs.veb[eb_index].vf.size();++i)
//		for (int j(0);j!=bs.veb[eb_index].vf[i].np.size();++j)
//		{
//			if (CGAL::squared_distance(center, bs.vpo[bs.veb[eb_index].vf[i].np[j]]) < r_square)
//				return true;		//有一个点在球内这个块体就在球内
//		}
//	return false;		//没有点在球内
//}
//
//bool cb_in_sphere (const BI::point_3& center, const BI::FT& radius, const BI::block_system& bs, const int& cb_index)
//{
//	for (int i(0);i!=bs.vcb[cb_index].vb.size();++i)
//	{
//		if (eb_in_sphere (center, radius, bs, bs.vcb[cb_index].vb[i]))
//			return true;
//	}
//	return false;
//}
//
//bool cb_in_sphere (const BI::point_3& center, const BI::FT& radius, const BI::block_system& bs, const vector<int>& cb_point)
//{
//	BI::FT square_r (radius*radius);
//	for (int i(0);i!=cb_point.size();++i)
//	{
//		if (CGAL::squared_distance(center, bs.vpo[cb_point[i]]) < square_r)
//			return true;
//	}
//	return false;
//}
//
//void get_cb_points (const BI::block_system& bs, const int& cb_index, vector<int>&res)
//{
//	res.clear();
//	for (int i(0);i!=bs.vcb[cb_index].vb.size ();++i)
//		for (int j(0);j!=bs.veb[bs.vcb[cb_index].vb[i]].vf.size ();++j)
//			for (int k(0);k!=bs.veb[bs.vcb[cb_index].vb[i]].vf[j].np.size ();++k)
//				BI::add_set(res, bs.veb[bs.vcb[cb_index].vb[i]].vf[j].np[k]);
//}
//
//bool cb_in_box (const BI::point_3& center, const BI::FT& radius, const BI::block_system& bs, const vector<int>& cb_point)
//{
//	for (int i(0);i!=cb_point.size();++i)
//	{
//		const BI::point_3 &p (bs.vpo[cb_point[i]]);
//		if (CGAL::abs(center.x()-p.x())<radius && CGAL::abs(center.y()-p.y())<radius && CGAL::abs(center.z()-p.z())<radius)
//			return true;
//	}
//	return false;
//}
//
//
//int main()
//{
//	ifstream cmd("thesis_mod_input.txt");
//	string filename;
//	cmd>>filename;
//
//	int format(0);
//	cmd>>format;
//	cout<<"format:"<<format<<endl;
//
//	double x(0.0), y(0.0), z(0.0);	//中心点的坐标
//	cmd>>x>>y>>z;
//	cout<<"x,y,z:"<<x<<','<<y<<','<<z<<endl;
//
//	double len(0.0);		//最大单元体的边长
//	cmd>>len;
//	cout<<"len:"<<len<<endl;
//
//	int count(0);		//一共计算几份
//	cmd>>count;
//	cout<<"count:"<<count<<endl;
//
//	string outfile_name;
//	cmd>>outfile_name;
//
//	cmd.close();
//	ifstream in_bs;
//
//	if (format == 1)
//	{
//		in_bs.open (filename, ios::binary);
//		CGAL::set_binary_mode(in_bs);
//	}
//	else
//		in_bs.open(filename);
//
//	if (!in_bs)
//	{
//		cout<<"Open "<<filename<<" error"<<endl;
//		return 0;
//	}
//	else
//		cout<<"Open "<<filename<<" success"<<endl;
//
//	BI::block_system bs;
//	cout<<"reading block system"<<endl;
//	in_bs>>bs;
//	in_bs.close();
//	cout<<"read block system complete"<<endl;
//
//	cout<<"building bounding box"<<endl;
//	vector<BI::mbox_3> vbox;
//	vector<vector<int>> vcb_point(bs.vcb.size ());
//	vbox.reserve(bs.vcb.size ());
//	for (int i(0);i!=bs.vcb.size ();++i)
//	{
//		get_cb_points(bs, i, vcb_point[i]);
//		vbox.push_back(BI::build_mbox_3(bs.vpo, vcb_point[i], 0.000001));
//	}
//	cout<<"boudning box built complete"<<endl;
//
//	ofstream outfile (outfile_name);
//	BI::FT step(len/count);
//
//	outfile<<"center:"<<x<<','<<y<<','<<z<<endl;
//	outfile<<"len:"<<len<<endl;
//	outfile<<"count:"<<count<<endl;
//	outfile<<"step:"<<CGAL::to_double(step)<<endl;
//	outfile<<"count:"<<count<<endl;
//
//	BI::point_3 center (x,y,z);
//	outfile<<"r;eb_count;density"<<endl;
//	for (int i(0);i!=count;++i)
//	{
//		cout<<i<<"/"<<count<<endl;
//		BI::FT r(step*(i+1));
//		double r_ (CGAL::to_double(r)+0.0000001);
//		BI::mbox_3 center_box(x-r_, y-r_, z-r_, x+r_, y+r_, z+r_);
//		int eb_count(0);
//		for (int j(0);j!=bs.vcb.size ();++j)
//		{
//			if (CGAL::do_intersect(center_box, vbox[j]) && cb_in_box(center, r, bs, vcb_point[j]))
//				eb_count++;
//		}
//		outfile<<CGAL::to_double(r)<<';'<<eb_count<<';'<<CGAL::to_double(eb_count/(4.0/3.0*3.141592653*r*r*r))<<endl;
//	}
//	outfile.close();
//	cout<<"calculation complete"<<endl;
//}