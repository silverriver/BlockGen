////gen_bs.exe
////ʶ�����
////��stl����ʾ�������е�ÿ����������Ƭ����Ϊһ���������϶
////gen_bs_input.txt��ʽ��
////Ҫʶ������ģ����
//////����ÿ��ģ�ͣ�
////model_domain�ļ���   
////����Բ����϶�ļ�  
//////ÿ��Բ����϶�ļ�������...   
////�����������϶�ļ�  
//////ÿ���������϶�ļ�������   
////����stl�ļ�  
//////ÿ��stl�ļ�������   
////����ļ����ļ���   
////����ļ��ĸ�ʽ��1��binary��2��ascii��
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
//int main ()
//{
//	ifstream cmd ("gen_bs_input.txt");
//	if (!cmd)
//	{
//		cout<<"can not find gen_bs_input.txt"<<endl;
//		return 0;
//	}
//
//	int nmodel(0);
//	cmd>>nmodel;
//	for (int imodel(0);imodel!=nmodel;++imodel)
//	{
//		cout<<"-------"<<imodel<<"-------"<<endl;
//		int nstl(0), ndisc(0), npoly(0);
//
//		string domain_name;
//		cmd>>domain_name;		//model_domain�ļ����ļ���
//		cmd>>ndisc;				//����Ҫ��ȡ����Բ����϶�ļ�
//		vector<string> vdisc_fn;
//		vdisc_fn.resize (ndisc);		//ÿ��Բ����϶�ļ�������
//		for (int i(0);i!=ndisc;++i)
//			cmd>>vdisc_fn[i];
//
//		if (!cmd)
//		{
//			cout<<"wrong format"<<endl;
//			return 0;
//		}
//
//		cmd>>npoly;
//		vector<string> vpoly_fn;
//		vpoly_fn.resize(npoly);
//		for (int i(0);i!=npoly;++i)
//			cmd>>vpoly_fn[i];
//
//		if (!cmd)
//		{
//			cout<<"wrong format"<<endl;
//			return 0;
//		}
//
//		cmd>>nstl;		//����Ҫ��ȡ����stl�ļ�
//		vector<string> vstl_fn;
//		vstl_fn.resize (nstl);		//ÿ��stl�ļ�������
//		for (int i(0);i!=nstl;++i)
//			cmd>>vstl_fn[i];
//
//		if (!cmd)
//		{
//			cout<<"wrong format"<<endl;
//			return 0;
//		}
//		
//		string output_name;
//		cmd>>output_name;		//ָ������ļ����ļ���
//		int output_format;
//		cmd>>output_format;		//ָ������ļ��ĸ�ʽ��1��binary��2��ascii��
//
//		if (!cmd)
//		{
//			cout<<"wrong format"<<endl;
//			return 0;
//		}
//
//		ifstream infile (domain_name);
//		if (!infile)
//		{
//			cout<<"Open domain file "<<domain_name<<" error"<<endl;
//			return 0;
//		}
//		else 
//			cout<<"Open domain file "<<domain_name<<" success"<<endl;
//
//		GB_domain d;
//		d.Input(infile);
//		infile.close();
//		BI::block_system bs;
//		bs.init_domain(d);
//
//		for (int i(0);i!=ndisc;++i)
//		{
//			infile.open (vdisc_fn[i]);
//			if (!infile)
//			{
//				cout<<"Open disc file "<<vdisc_fn[i]<<" error"<<endl;
//				continue;
//			}
//			else 
//				cout<<"Open disc file "<<vdisc_fn[i]<<" success"<<endl;
//
//			if (bs.disc_frac_parser(infile,d.Get_dx_axis()) ==0)
//				cout<<"Parse "<<vdisc_fn[i]<<" success"<<endl;
//			else
//			{
//				cout<<"Parse "<<vdisc_fn[i]<<" fail"<<endl;
//				continue;
//			}
//			infile.close();
//		}
//
//		for (int i(0);i!=npoly;++i)
//		{
//			infile.open (vpoly_fn[i]);
//			if (!infile)
//			{
//				cout<<"Open poly file "<<vpoly_fn[i]<<" error"<<endl;
//				continue;
//			}
//			else
//				cout<<"Open poly file "<<vpoly_fn[i]<<" success"<<endl;
//			if (bs.poly_frac_parser(infile) == 0)
//				cout<<"Parse "<<vpoly_fn[i]<<" success"<<endl;
//			else
//			{
//				cout<<"Parse "<<vpoly_fn[i]<<" fail"<<endl;
//				continue;
//			}
//			infile.close();
//		}
//
//		for (int i(0);i!=nstl;++i)
//		{
//			infile.open (vstl_fn[i]);
//			if (!infile)
//			{
//				cout<<"Open stl file "<<vstl_fn[i]<<" error"<<endl;
//				continue;
//			}
//			else
//				cout<<"Open stl file "<<vstl_fn[i]<<" success"<<endl;
//
//			if (bs.stl_parser(infile) == 0)
//				cout<<"Parse "<<vstl_fn[i]<<" success"<<endl;
//			else
//			{
//				cout<<"Parse "<<vstl_fn[i]<<" fail"<<endl;
//				continue;
//			}
//			infile.close();
//		}
//
//		cout<<"Read frac file complete...identifying blocks"<<endl;
//		bs.identify_block (cout);
//
//		cout<<"Outputing"<<endl;
//
//		if (output_format == 1)
//		{
//			ofstream outfile (output_name, ios::binary );
//			CGAL::set_binary_mode(outfile);
//			outfile<<bs;
//			outfile.close();
//		}
//		else 
//		{
//			ofstream outfile (output_name);
//			outfile<<bs;
//			outfile.close();
//		}
//		ofstream outfile (output_name+"_info.txt");
//		bs.output_info(outfile);
//		outfile.close();
//		cout<<"Models output complete"<<endl;
//	}
//	return 0;
//}