////bs2bound.exe
////��ȡ����ϵͳ�����ɲ���������߽���ڱ���߽磬������п���ı���߽�
////bs2bound_input.txt��ʽ��
////bs�ļ���
////bs�ļ���ʽ ��1��ʾ��binary��2��ʾ��ascii��
////��ըϵ��   
////����ļ���ʽ(1:�����ڱ��棬2����������棬3�����������ڱ��������棬 4�����ɵ��ڱ����н���������ηֿ���5�����ɵ�������н���������ηֿ���6�����ɵ��ڱ����������н���������ηֿ�)
////����ļ�����ǰ׺
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
//	int format (0);		//bs�ļ��ĸ�ʽ��1��ʾ��binary��2��ʾ��ascii
//	cmd>>format;
//
//	double exp_ratio(0);
//	cmd>>exp_ratio;		//��ըϵ��
//
//	int out_flag(0);			//out_flag��¼Ϊÿ��cblock����ʲô��Ϣ(1:�����ڱ��棬2����������棬3�����������ڱ���������
//								//									 4�����ɵ��ڱ����н���������ηֿ���5�����ɵ�������н���������ηֿ���6�����ɵ��ڱ����������н���������ηֿ�)
//	cmd>>out_flag;
//
//	string out_perfix;
//	cmd>>out_perfix;		//����ļ���ǰ׺
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
//		out_ib.open (out_perfix+"_ib.stl");				//�����ڱ���
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
//		out_ob.open (out_perfix+"_ob.stl");			//���������
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
//		out_ib.open (out_perfix+"_ib_separated.stl");				//�����ڱ���
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
//		out_ob.open (out_perfix+"_ob_separated.stl");			//���������
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
