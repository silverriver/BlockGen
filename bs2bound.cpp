////bs2bound.exe
////读取块体系统，生成并输出外表面边界和内表面边界，输出所有块体的表面边界
////bs2bound_input.txt格式：
////bs文件名
////bs文件格式 （1表示是binary，2表示是ascii）
////爆炸系数   
////输出文件格式(1:生成内表面，2：生成外表面，3：生成整体内表面和外表面， 4：生成的内表面中将各个多边形分开，5：生成的外表面中将各个多边形分开，6：生成的内表面和外表面中将各个多边形分开)
////输出文件名的前缀
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
//int main()
//{
//	ifstream cmd ("bs2bound_input.txt");
//
//	if(!cmd)
//	{
//		cout<<"can't open bs2bound_input.txt"<<endl;
//		return 0;
//	}
//
//	string filename;
//	cmd>>filename;
//
//	int format (0);		//bs文件的格式，1表示是binary，2表示是ascii
//	cmd>>format;
//
//	double exp_ratio(0);
//	cmd>>exp_ratio;		//爆炸系数
//
//	int out_flag(0);			//out_flag记录为每个cblock生成什么信息(1:生成内表面，2：生成外表面，3：生成整体内表面和外表面
//								//									 4：生成的内表面中将各个多边形分开，5：生成的外表面中将各个多边形分开，6：生成的内表面和外表面中将各个多边形分开)
//	cmd>>out_flag;
//
//	string out_perfix;
//	cmd>>out_perfix;		//输出文件的前缀
//
//	cmd.close();
//
//	ifstream in_bs;
//	if(format ==1)
//	{
//		in_bs.open (filename, ios::binary );
//		CGAL::set_binary_mode (in_bs);
//	}
//	else
//		in_bs.open (filename);
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
//	in_bs>>bs;
//	in_bs.close();
//	cout<<"Read block system complete"<<endl;
//
//	ofstream out_ob, out_ib;
//
//	vector<BI::point_3> cb_p;
//	BI::point_3 bs_center(bs.get_rand_p());
//
//	cb_p.reserve (bs.vcb.size ());
//	for (int i(0);i!=bs.vcb.size ();++i)
//		cb_p.push_back (bs.vcb[i].get_rand_p(bs.veb,bs.vpo));
//
//	switch (out_flag)
//	{
//	case 3:
//	case 1:
//		out_ib.open (out_perfix+"_ib.stl");				//生成内表面
//		for (int i(0);i!=bs.vcb.size ();++i)
//		{
//			stringstream ss;
//			ss<<"IB"<<i;
//			BI::inner_bound ib;
//			ib.construct_inner_bound(bs,i);
//			ib.output_STL(out_ib,(cb_p[i].x()-bs_center.x())*exp_ratio, (cb_p[i].y()-bs_center.y())*exp_ratio, (cb_p[i].z()-bs_center.z())*exp_ratio,ss.str ());
//		}
//		cout<<"IB generation complete"<<endl;
//		if (out_flag == 1)
//			break;
//	case 2:
//		out_ob.open (out_perfix+"_ob.stl");			//生成外表面
//		for (int i(0);i!=bs.vcb.size ();++i)
//		{
//			stringstream ss;
//			ss<<"OB"<<i;
//			BI::outer_bound ob;
//			ob.construct_outer_bound (bs,i);
//			ob.output_STL(out_ob, (cb_p[i].x()-bs_center.x())*exp_ratio, (cb_p[i].y()-bs_center.y())*exp_ratio, (cb_p[i].z()-bs_center.z())*exp_ratio,ss.str ());
//		}
//		cout<<"OB generation complete"<<endl;
//		break;
//	case 6:
//	case 4:
//		out_ib.open (out_perfix+"_ib_separated.stl");				//生成内表面
//		for (int i(0);i!=bs.vcb.size ();++i)
//		{
//			stringstream ss;
//			ss<<"IB"<<i;
//			BI::inner_bound ib;
//			ib.construct_inner_bound(bs,i);
//			ib.output_exp_STL(out_ib, 0.0, (cb_p[i].x()-bs_center.x())*exp_ratio, (cb_p[i].y()-bs_center.y())*exp_ratio, (cb_p[i].z()-bs_center.z())*exp_ratio,ss.str());
//		}
//		cout<<"IB generation complete poly separated"<<endl;
//		if (out_flag == 4)
//			break;
//	case 5:
//		out_ob.open (out_perfix+"_ob_separated.stl");			//生成外表面
//		for (int i(0);i!=bs.vcb.size ();++i)
//		{
//			stringstream ss;
//			ss<<"OB"<<i;
//			BI::outer_bound ob;
//			ob.construct_outer_bound (bs,i);
//			ob.output_exp_STL(out_ob, 0.0, (cb_p[i].x()-bs_center.x())*exp_ratio, (cb_p[i].y()-bs_center.y())*exp_ratio, (cb_p[i].z()-bs_center.z())*exp_ratio,ss.str());
//		}
//		cout<<"OB generation complete poly separated"<<endl;
//		break;
//	default:
//		cout<<"successfully end"<<endl;
//		break;
//	}
//	out_ib.close();
//	out_ob.close();
//	cout<<"boundary structure generation complete"<<endl;
//	return 0;
//}
