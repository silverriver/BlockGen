/////////////////////////////////////////////////////////////////////////////////////////////
////                                                                                       //
//// BlockCut                                                                              //
////                                                                                       //
//// A block identification program.                                                       //
////                                                                                       //
//// Version 1.0                                                                           //
//// 2017-5-18                                                                             //
////                                                                                       //
//// BlockCut is freely available through http://www.github.com/silverriver/BlockCut       //
////   It may be copied, modified, and redistributed for non-commercial use.               //
////   Please consult the file LICENSE for the detailed copyright notices.                 //
////                                                                                       //
/////////////////////////////////////////////////////////////////////////////////////////////
//
////#define BLKGLIBRARY
//
//// Uncomment the following line to disable assert macros. These macros were intented to be 
////   used to catch bugs in the code.
//#define NDEBUG
//
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
////using BI::operator<<;
////using BI::operator>>;
////下面是原始代码。2015-09-02
////int main ()
////{
////	ifstream infile("model_domain.dat");
////	GB_domain d;
////	d.Input(infile);
////	BI::block_system bs;
////	bs.init_domain(d);
////	cout<<"initial element blocks:"<<bs.veb.size ()<<endl;
////
////	CGAL::Real_timer timer;
////	if (bs.init_check()==0)
////		cout<<"read modeling_domain.dat success"<<endl;
////	else
////	{
////		cout<<"read modeling_domain.dat fail"<<endl;
////		return 0;
////	}
////
////	infile.close();
////	infile.open ("contri_fracture_xyzabr.dat");
////	timer.reset();
////	timer.start();
////	if (bs.disc_frac_parser(infile,d.Get_dx_axis()) == 0)
////		cout<<"parser disc frac file success time:"<<timer.time()<<" disc count:"<<bs.vdisc.size ()<<endl;
////	else
////	{
////		cout<<"parser disc frac file fail"<<endl;
////		return 0;
////	}
////	timer.stop(); timer.reset(); timer.start();
////	bs.identify_block();
////	ofstream outfile;
////	timer.stop();
////	cout<<bs.vcb.size ()<<" complex blocks are generated time:"<<timer.time()<<endl;
////	cout<<bs.veb.size ()<<" element blocks are generated"<<endl;
////	
////	//核对一下单元块体的体积是否正确
////	BI::FT vol(0);
////	for (int i(0);i!=bs.veb.size ();++i)
////		vol = vol+bs.eb_volume(i);
////	cout<<"total volume of element blocks:"<<vol<<endl;
////
////	cout<<"how many cblock to generate? (enter -1 to output all)"<<endl;
////	int cbn;
////	cin>>cbn;
////	timer.reset(); timer.start();
////	outfile.open ("cblock.stl");
////	if (cbn<0)
////		for (int i(0);i!=bs.vcb.size ();++i)
////		{
////			BI::outer_bound ob;
////			ob.construct_outer_bound(bs,i,false);
////			ob.output_STL(outfile);
////		}
////	else
////	{
////		vector<int> vcbi;
////		int cbi;
////		cout<<"please enter "<<cbn<<" numbers"<<endl;
////		for (int i(0);i!=cbn;++i)
////		{
////			cin>>cbi;
////			vcbi.push_back (cbi);
////		}
////		for (int i(0);i!=cbn;++i)
////		{
////			BI::outer_bound ob;
////			ob.construct_outer_bound(bs,vcbi[i],false);
////			ob.output_STL(outfile);
////		}
////	}
////	outfile.close();
////	timer.stop();
////	cout<<"output STL file complete time:"<<timer.time ()<<endl;
////
////	return 0;
////}
//
////下面是为了测试单元函数与运算时间和内存的关系得到的代码
////int main ()
////{
////	int n_exp (0);
////	cout<<"how many DFN to model:";
////	cin>>n_exp;
////	string dfn_file_name;
////	cout<<"enter the name of DFN files:";
////	cin>>dfn_file_name;
////	double exp_ratio;
////	cout<<"explode ratio [-1 to use the defult (0.25)]:";
////	cin>>exp_ratio;
////	if (exp_ratio < 0)
////		exp_ratio = 0.25;
////	int vol_flag;
////	cout<<"required to analysis vol? (1-yes, -1-no):";
////	cin>>vol_flag;
////	ofstream outfile("time_memory.txt");
////	outfile.precision(13);
////	for (int curr(0);curr!=n_exp;++curr)
////	{
////		try
////		{
////			//debug
////			cout<<"curr:"<<endl;
////			cout<<curr<<endl;
////			cout<<"memory:"<<endl;
////			CGAL::Memory_sizer mem;
////			cout<<"timer"<<endl;
////			CGAL::Timer timer;
////			timer.stop(); timer.reset();
////
////			outfile<<"----------model"<<curr<<"--------"<<endl;
////			cout<<"initlizing..."<<endl;
////			//read model_domain
////			timer.start();
////			ifstream infile("model_domain.dat");
////			GB_domain d;
////			d.Input(infile);
////			BI::block_system bs;
////			bs.init_domain(d);
////			timer.stop();
////			outfile<<"fram_count:\t\t\t"<<bs.veb.size ()<<endl<<"init_fram_time:\t\t"<<timer.time()<<endl;
////
////			if (bs.init_check()==0)
////				cout<<"read modeling_domain.dat success"<<endl;
////			else
////			{
////				cout<<"read modeling_domain.dat fail"<<endl;
////				continue;
////			}
////			infile.close();
////
////			//计算一下读取的模型的中心点
////			BI::FT minx,maxx,miny,maxy,minz,maxz;
////			minx = maxx = bs.vpo[0].x();
////			miny = maxy = bs.vpo[0].y();
////			minz = maxz = bs.vpo[0].z();
////			for (int i(0);i!=bs.vpo.size ();++i)
////			{
////				if (bs.vpo[i].x() < minx)
////					minx = bs.vpo[i].x();
////				if (bs.vpo[i].x() > maxx)
////					maxx = bs.vpo[i].x();
////				if (bs.vpo[i].y() < miny)
////					miny = bs.vpo[i].y();
////				if (bs.vpo[i].y() > maxy)
////					maxy = bs.vpo[i].y();
////				if (bs.vpo[i].z() < minz)
////					minz = bs.vpo[i].z();
////				if (bs.vpo[i].z() > maxz)
////					maxz = bs.vpo[i].z();
////			}
////			BI::point_3 domain_center((minx+maxx)/2, (miny+maxy)/2, (minz+maxz)/2);
////
////			//read contri_fracture_xuzabr.dat
////			ostringstream ost;
////			ost<<dfn_file_name<<curr<<".dat";
////			infile.open (ost.str ());
////			timer.reset();
////			timer.start();
////			if (bs.disc_frac_parser(infile,d.Get_dx_axis()) == 0)
////			{
////				outfile<<"parser_disc_time:\t"<<timer.time()<<endl;
////				outfile<<"disc_count:\t\t\t"<<bs.vdisc.size ()<<endl;
////			}
////			else
////			{
////				outfile<<"parser disc frac file fail"<<endl;
////				continue;
////			}
////
////			//识别块体
////			cout<<"identifying blocks"<<endl;
////			timer.stop(); timer.reset(); timer.start();
////			bs.identify_block();
////			timer.stop();
////
////			outfile<<"Memory_used:\t\t"<<mem.resident_size()/1024.0/1024.0<<endl;
////			outfile<<"BI_Time:\t\t\t"<<timer.time()<<endl;
////			outfile<<"cutting_time:\t\t"<<bs.cutting_time<<endl;
////			outfile<<"shrk_time:\t\t\t"<<bs.shrk_time<<endl;
////			outfile<<"CB_count:\t\t\t"<<bs.vcb.size ()<<endl;
////			outfile<<"EB_count:\t\t\t"<<bs.veb.size ()<<endl;
////			outfile<<"Point_count:\t\t"<<bs.vpo.size ()<<endl;
////			outfile<<"Plane_count:\t\t"<<bs.vpl.size ()<<endl;
////			cout<<"Time:"<<timer.time()<<endl;
////
////			//之所以吧下面这一段注释掉是因为在运行的时候碰到了一个很奇怪的bug。当算到第6个DFN的时候会出错。也就是说算306个裂隙的时候会出错，是CGAL内核的异常。没有细研究到底是什么情况。
////			//然后我把第二次vol=0; 这个语句换成了BI::FT vol1(0); 并且在下面做了一下修改，然后就没有问题了。总之是很诡异的bug。有时间可以研究一下，没准还能给CGAL提一个feedback之类的。 2015-9-29 14:42
////			//BI::FT vol(0);
////			//for (int i(0);i!=eb_vol.size ();++i)
////			//	vol = vol + eb_vol[i];
////			//outfile<<"total volume of element blocks:"<<vol<<endl;
////			//vol = 0;
////			//for (int i(0);i!=cb_vol.size ();++i)
////			//	vol = vol + cb_vol[i];
////			//outfile<<"total volume of complex blocks:"<<vol<<endl;
////
////			//输出STL文件
////			cout<<"generating STL file"<<endl;
////			timer.reset(); timer.start();
////			ost.str ("");
////			ost<<curr<<"_cblock.stl";
////			ofstream stl_file (ost.str());
////
////			ost.str ("");
////			ost<<curr<<"_cblock_exploded.stl";
////			ofstream stl_file_exp(ost.str());		//用来输出exploded的stl文件
////
////			vector<int> cb_facen;					//用来记录每个cb有几个面
////			cb_facen.reserve (bs.vcb.size ());	
////
////			for (int i(0);i!=bs.vcb.size ();++i)
////			{
////				BI::outer_bound ob;
////				ob.construct_outer_bound(bs,i,false);
////				ob.output_STL(stl_file);
////
////				cb_facen.push_back (ob.get_polycount());
////
////				//exp_ratio表示在explode view中块体被移动多远
////				BI::point_3 p_temp(ob.get_center());
////				ob.output_STL(stl_file_exp,CGAL::to_double(p_temp.x()-domain_center.x())*exp_ratio,
////					CGAL::to_double(p_temp.y()-domain_center.y())*exp_ratio,
////					CGAL::to_double(p_temp.z()-domain_center.z())*exp_ratio);
////			}
////			stl_file.close();
////			stl_file_exp.close();
////			timer.stop();
////			outfile<<"STL_file_time:\t\t"<<timer.time ()<<endl;
////			cout<<"Time:"<<timer.time()<<endl;
////
////			//计算块体体积
////			if (vol_flag > 0)
////			{
////				//输出eb和cb的统计信息
////				ofstream eb_file,cb_file;
////				ost.str ("");
////				ost<<curr<<"_eb_stat.txt";
////				eb_file.open (ost.str ());
////				ost.str ("");
////				ost<<curr<<"_cb_stat.txt";
////				cb_file.open (ost.str ());
////				eb_file.precision(13); cb_file.precision(13);
////				eb_file<<"\"point_count\"\t \"face_count\"\t \"volume\""<<endl;
////				for (int i(0);i!=bs.veb.size ();++i)
////				{
////					eb_file<<bs.veb[i].get_point_count()<<"\t"
////						<<bs.veb[i].vf.size ()<<"\t"<<bs.eb_volume(i)<<endl;
////				}
////				eb_file.close ();
////
////				cb_file<<"\"face_count\"\t \"eb_count\"\t \"volume\""<<endl;
////				for (int i(0);i!=bs.vcb.size ();++i)
////				{
////					BI::FT vol(0);
////					for (int j(0);j!=bs.vcb[i].vb.size ();++j)
////						vol += bs.eb_volume(bs.vcb[i].vb[j]);
////					cb_file<<cb_facen[i]<<"\t"<<bs.vcb[i].vb.size ()<<"\t"<<vol<<endl;
////				}
////				cb_file.close();
////				cout<<"cb_file output complete"<<endl;
////			}
////		}
////		catch (exception &e)
////		{
////			cout<<"ERROR!"<<endl<<e.what ()<<endl;
////			continue;
////		}
////	}
////	outfile.close();
////	return 0;
////}
//
////下面的代码是为了生成一个例子，因为要涉及到较细的划分，所以不能通过读取有限精度文件的方法来识别块体
////int main ()
////{
////	int cut_n(0);
////	cout<<"input n:";
////	cin>>cut_n;
////	double exp_ratio;
////	cout<<"explode ratio [-1 to use the defult (0.35)]:";
////	cin>>exp_ratio;
////	if (exp_ratio < 0)
////		exp_ratio = 0.35;
////
////	ofstream outfile("time_memory.txt");
////	outfile.precision (13);
////	
////	CGAL::Memory_sizer mem;
////	CGAL::Timer timer;
////	timer.stop (); timer.reset ();
////
////	BI::FT minx(0),maxx(10),miny(0),maxy(10),minz(0),maxz(10);
////	BI::block_system bs;
////	bs.add_rect_domain(minx,miny,minz,maxx,maxy,maxz);
////
////	if (bs.init_check()==0)
////		cout<<"generate modeling_domain.dat success"<<endl;
////	else
////	{
////		cout<<"generate modeling_domain.dat fail"<<endl;
////		return 0;
////	}
////
////	BI::point_3 domain_center((minx+maxx)/2, (miny+maxy)/2, (minz+maxz)/2);
////	
////	//生成多边形裂隙
////	bs.vpoly.reserve (cut_n*3+1);
////	BI::FT offset(maxx-minx);
////	for (int i(1);i!=cut_n+1;++i)
////	{
////		offset = offset/2.0;
////		BI::poly_frac poly;
////		poly.aperture = 0.01;
////		poly.cohesion = 0.01;
////		poly.f_angle = 0.1;
////		poly.vp.reserve (4);
////
////		poly.vp.push_back (BI::point_3(minx-1,miny-1,minz+offset));
////		poly.vp.push_back (BI::point_3(maxx+1,miny-1,minz+offset));
////		poly.vp.push_back (BI::point_3(maxx+1,maxy+1,minz+offset));
////		poly.vp.push_back (BI::point_3(minx-1,maxy+1,minz+offset));
////		poly.frac_id = bs.vdisc .size ()+3*i;
////		poly.plane_id = BI::add_set_pl(bs.vpl,BI::plane_3 (poly.vp[0],BI::oriented_normal(poly.vp)));
////		bs.vpoly.push_back (poly);
////
////		poly.vp[0] = BI::point_3(minx-1, miny+offset, minz-1);
////		poly.vp[1] = BI::point_3(maxx+1, miny+offset, minz-1);
////		poly.vp[2] = BI::point_3(maxx+1, miny+offset, maxz+1);
////		poly.vp[3] = BI::point_3(minx-1, miny+offset, maxz+1);
////		poly.frac_id = bs.vdisc .size () + 3*i+1;
////		poly.plane_id = BI::add_set_pl(bs.vpl,BI::plane_3(poly.vp[0],BI::oriented_normal(poly.vp)));
////		bs.vpoly.push_back (poly);
////
////		poly.vp[0] = BI::point_3(minx+offset, miny-1, minz-1);
////		poly.vp[1] = BI::point_3(minx+offset, maxy+1, minz-1);
////		poly.vp[2] = BI::point_3(minx+offset, maxy+1, maxz+1);
////		poly.vp[3] = BI::point_3(minx+offset, miny-1, maxz+1);
////		poly.frac_id = bs.vdisc .size () + 3*i+2;
////		poly.plane_id = BI::add_set_pl(bs.vpl,BI::plane_3(poly.vp[0],BI::oriented_normal(poly.vp)));
////		bs.vpoly.push_back (poly);
////	}
////	cout<<"generate poly. frac. complete"<<endl;
////	outfile<<"poly_disc_count:\t\t"<<bs.vpoly.size ()<<endl;
////	timer.start ();
////	bs.identify_block ();
////	timer.stop();
////
////
////	outfile<<"Memory_used:\t\t"<<mem.resident_size()/1024.0/1024.0<<endl;
////	outfile<<"BI_Time:\t\t\t"<<timer.time()<<endl;
////	outfile<<"cutting_time:\t\t"<<bs.cutting_time<<endl;
////	outfile<<"shrk_time:\t\t\t"<<bs.shrk_time<<endl;
////	outfile<<"CB_count:\t\t\t"<<bs.vcb.size ()<<endl;
////	outfile<<"EB_count:\t\t\t"<<bs.veb.size ()<<endl;
////	outfile<<"Point_count:\t\t"<<bs.vpo.size ()<<endl;
////	outfile<<"Plane_count:\t\t"<<bs.vpl.size ()<<endl;
////	cout<<"identify blocks complete Time:"<<timer.time()<<endl;
////
////	//输出stl文件
////	cout<<"generationg STL file"<<endl;
////	timer.reset ();timer.start();
////	ofstream stl_file("cube.stl");
////	ofstream stl_file_exp("cube_exp.stl");
////
////	vector<int> cb_facen;
////	cb_facen.reserve (bs.vcb.size ());
////	for (int i(0);i!=bs.vcb.size ();++i)
////	{
////		BI::outer_bound ob;
////		ob.construct_outer_bound (bs,i,false);
////		ob.output_STL(stl_file);
////		cb_facen.push_back (ob.get_polycount());
////
////		BI::point_3 p_temp(ob.get_center());
////		ob.output_STL(stl_file_exp,(p_temp.x()-domain_center.x())*exp_ratio,
////					(p_temp.y()-domain_center.y())*exp_ratio,
////					(p_temp.z()-domain_center.z())*exp_ratio);
////	}
////	stl_file.close();
////	stl_file_exp.close ();
////	timer.stop();
////	outfile<<"STL_file_time:\t\t"<<timer.time ()<<endl;
////	cout<<"Time:"<<timer.time()<<endl;
////
////	//计算块体体积
////	ofstream eb_file("cube_eb_stat.txt"), cb_file("cube_cb_state.txt");
////	eb_file.precision (13); cb_file.precision (13);
////
////	eb_file<<"\"point_count\"\t \"face_count\"\t \"volume\""<<endl;
////	for (int i(0);i!=bs.veb.size ();++i)
////	{
////		eb_file<<bs.veb[i].get_point_count()<<"\t"
////			<<bs.veb[i].vf.size ()<<"\t"<<bs.eb_volume(i)<<endl;
////	}
////	eb_file.close ();
////	
////	cb_file<<"\"face_count\"\t \"eb_count\"\t \"volume\""<<endl;
////	for (int i(0);i!=bs.vcb.size ();++i)
////	{
////		BI::FT vol(0);
////		for (int j(0);j!=bs.vcb[i].vb.size ();++j)
////			vol += bs.eb_volume(bs.vcb[i].vb[j]);
////		cb_file<<cb_facen[i]<<"\t"<<bs.vcb[i].vb.size ()<<"\t"<<vol<<endl;
////	}
////	cb_file.close();
////	cout<<"cb_file output complete"<<endl;
////
////	outfile.close();
////
////	return 0;
////}
//
////作一些小的测试
////int main ()
////{
////	ifstream infile ("model_domain.dat");
////	GB_domain d;
////	d.Input (infile);
////	BI::block_system bs;
////	bs.init_domain(d);
////	if (bs.init_check()==0)
////		cout<<"read modeling_domain.dat success"<<endl;
////	else
////	{
////		cout<<"read modeling_domain.dat fail"<<endl;
////		return 0;
////	}
////	infile.close();
////
////	BI::poly_frac pf;
////
////	pf.vp_ob.push_back (BI::point_3(-2,2,2));
////	pf.vp_ob.push_back (BI::point_3(-2,8,2));	//front
////	pf.vp_ob.push_back (BI::point_3(-2,8,8));
////	pf.vp_ob.push_back (BI::point_3(-2,2,8));
////	bs.add_poly_frac(pf);
////
////	pf.vp_ob.clear();
////	pf.vp_ob.push_back (BI::point_3(-8,2,2));
////	pf.vp_ob.push_back (BI::point_3(-8,8,2));
////	pf.vp_ob.push_back (BI::point_3(-8,8,8));	//back
////	pf.vp_ob.push_back (BI::point_3(-8,2,8));
////	bs.add_poly_frac(pf);
////
////	pf.vp_ob.clear();
////	pf.vp_ob.push_back (BI::point_3(-2,2,2));
////	pf.vp_ob.push_back (BI::point_3(-8,2,2));	//left
////	pf.vp_ob.push_back (BI::point_3(-8,2,8));
////	pf.vp_ob.push_back (BI::point_3(-2,2,8));
////	bs.add_poly_frac(pf);
////
////	pf.vp_ob.clear();
////	pf.vp_ob.push_back (BI::point_3(-2,8,2));
////	pf.vp_ob.push_back (BI::point_3(-8,8,2));
////	pf.vp_ob.push_back (BI::point_3(-8,8,8));	//right
////	pf.vp_ob.push_back (BI::point_3(-2,8,8));
////	bs.add_poly_frac(pf);
////
////	pf.vp_ob.clear();
////	pf.vp_ob.push_back (BI::point_3(-2,2,2));
////	pf.vp_ob.push_back (BI::point_3(-2,8,2));	//bottom
////	pf.vp_ob.push_back (BI::point_3(-8,8,2));
////	pf.vp_ob.push_back (BI::point_3(-8,2,2));
////	pf.vp_ib.push_back (BI::poly_frac::v_loop ());
////	pf.vp_ib[0].push_back (BI::point_3(-4,4,2));
////	pf.vp_ib[0].push_back (BI::point_3(-4,6,2));
////	pf.vp_ib[0].push_back (BI::point_3(-6,6,2));
////	pf.vp_ib[0].push_back (BI::point_3(-6,4,2));
////	bs.add_poly_frac(pf);
////
////	pf.vp_ob.clear(); pf.vp_ib.clear();
////	pf.vp_ob.push_back (BI::point_3(-2,2,8));
////	pf.vp_ob.push_back (BI::point_3(-2,8,8));	//top
////	pf.vp_ob.push_back (BI::point_3(-8,8,8));
////	pf.vp_ob.push_back (BI::point_3(-8,2,8));
////	bs.add_poly_frac(pf);
////
////	pf.vp_ob.clear();
////	pf.vp_ob.push_back (BI::point_3(-4,4,2));
////	pf.vp_ob.push_back (BI::point_3(-4,6,2));	//front
////	pf.vp_ob.push_back (BI::point_3(-4,6,4));
////	pf.vp_ob.push_back (BI::point_3(-4,4,4));
////	bs.add_poly_frac(pf);
////
////	pf.vp_ob.clear();
////	pf.vp_ob.push_back (BI::point_3(-6,4,2));
////	pf.vp_ob.push_back (BI::point_3(-6,6,2));	//back
////	pf.vp_ob.push_back (BI::point_3(-6,6,4));
////	pf.vp_ob.push_back (BI::point_3(-6,4,4));
////	bs.add_poly_frac(pf);
////
////	pf.vp_ob.clear();
////	pf.vp_ob.push_back (BI::point_3(-4,4,2));
////	pf.vp_ob.push_back (BI::point_3(-4,4,4));	//left
////	pf.vp_ob.push_back (BI::point_3(-6,4,4));
////	pf.vp_ob.push_back (BI::point_3(-6,4,2));
////	bs.add_poly_frac(pf);
////
////	pf.vp_ob.clear();
////	pf.vp_ob.push_back (BI::point_3(-4,6,2));
////	pf.vp_ob.push_back (BI::point_3(-4,6,4));	//right
////	pf.vp_ob.push_back (BI::point_3(-6,6,4));
////	pf.vp_ob.push_back (BI::point_3(-6,6,2));
////	bs.add_poly_frac(pf);
////
////	pf.vp_ob.clear();
////	pf.vp_ob.push_back (BI::point_3(-4,4,4));
////	pf.vp_ob.push_back (BI::point_3(-4,6,4));	//top
////	pf.vp_ob.push_back (BI::point_3(-6,6,4));
////	pf.vp_ob.push_back (BI::point_3(-6,4,4));
////	bs.add_poly_frac(pf);
////
////	ofstream outfile ("poly.txt");
////	for (int i(0);i!= (int) bs.vpoly.size ();++i)
////	{
////		BI::out_line_3(outfile, bs.vpoly[i]);
////		outfile<<endl;
////	}
////	outfile.close();
////
////
////	cout<<"identifying blocks"<<endl;
////	bs.identify_block();
////
////	cout<<"eb count: "<<bs.veb.size ()<<endl;
////	cout<<"cb count: "<<bs.vcb.size ()<<endl;
////	cout<<"Generating STL file"<<endl;
////
////	outfile.open ("eb.stl");
////	for (int i(0);i!= (int) bs.veb.size ();++i)
////	{
////		bs.output_eb_stl (i,outfile);
////	}
////	outfile.close();
////
////	outfile.open ("cb.stl");
////	for (int i(0);i!= (int) bs.vcb.size ();++i)
////		bs.output_cb_stl(i,outfile);
////	outfile.close();
////
////	BI::entrance_block entrb;
////	cout<<"constructing entrance block info"<<endl;
////	cout<<entrb.gen_entrance_block(bs,bs.vcb[0],bs.vcb[1],bs.vcb[0].get_rand_p(bs.veb,bs.vpo))<<endl;
////	cout<<"entrb.bs.vpoly.size():"<<entrb.bs.vpoly.size()<<endl;
////	cout<<"entrb.vcover.size():"<<entrb.vcover.size()<<endl;
////	cout<<"entrb.bs.veb.size ()"<<entrb.bs.veb.size ()<<endl;
////	cout<<"entrb.bs.vcb.size ()"<<entrb.bs.vcb.size ()<<endl;
////	cout<<"entrb.inner.vb.size ()"<<entrb.inner.vb.size ()<<endl;
////	outfile.open ("cover.txt");
////	for (int i(0);i!=entrb.vcover.size ();++i)
////	{
////		BI::out_line_3_rhino(outfile,entrb.vcover[i]);
////		outfile<<endl;
////	}
////	outfile.close();
////	outfile.open ("vpoly.txt");
////	for (int i(0);i!=entrb.bs.vpoly.size ();++i)
////		BI::out_line_3_rhino(outfile,entrb.bs.vpoly[i]);
////	outfile.close();
////
////	outfile.open("entrence_blk.stl");
////	entrb.outer_b.output_STL(outfile);
////	outfile.close();
////
////	outfile.open ("entrb.bs.cb.stl");
////	for (int i(0);i!=entrb.bs.vcb.size ();++i)
////		entrb.bs.output_cb_stl(i,outfile);
////	outfile.close();
////
////	outfile.open ("entrb.bs.eb.stl");
////	for (int i(0);i!=entrb.bs.veb.size ();++i)
////		entrb.bs.output_eb_stl(i,outfile);
////	outfile.close();
////
////	return 0;
////}
//
//////再做一些小测试
////int main ()
////{
////	BI::curve_frac c1,c2;
////	ifstream infile1 ("curve1.stl"), infile2("curve2.stl");
////	c1.read_STL (infile1,1); c2.read_STL (infile2,2);
////	infile1.close(); infile2.close();
////	
////	c1.join (c2);
////	ofstream outfile ("result.stl");
////	c1.output_STL(outfile);
////	outfile.close();
////
////	BI::curve_frac c3;
////	ifstream infile3("result.stl");
////	cout<<c3.read_STL (infile3,3);
////
////	return 0;
////}
//
////demo：如何使用进入entrance_block这个类
////int main ()
////{
////	ifstream infile ("model_domain.dat");
////	ofstream outfile;
////	GB_domain d;
////	d.Input(infile);
////	infile.close();
////	BI::block_system bs;
////	
////	bs.init_domain(d);
////	if (bs.init_check()==0)
////		cout<<"read modeling_domain.dat success"<<endl;
////	else
////	{
////		cout<<"read modeling_domain.dat fail"<<endl;
////		return 0;
////	}
////
////	infile.open("contri_fracture_xyzabr.dat");
////	if (!infile)
////	{
////		cout<<"cann't read contri_fracture_xyzabr.dat"<<endl;
////		return 0;
////	}
////	
////	if (bs.disc_frac_parser (infile,d.Get_dx_axis()) == 0)
////		cout<<"parser disc frac file success "<<bs.vdisc.size ()<<endl;
////	else
////	{
////		cout<<"parser disc frac file fail"<<endl;
////		return 0;
////	}
////	infile.close();
////	bs.identify_block();
////	cout<<bs.vcb.size ()<<" cblocks obtained"<<endl;
////
////	outfile.open ("eblock.stl");
////	for (int i(0);i!=bs.veb.size ();++i)
////		bs.output_eb_stl(i,outfile);
////	outfile.close();
////
////	outfile.open ("cblock.stl");
////	for (int i(0);i!=bs.vcb.size ();++i)
////		bs.output_cb_stl(i,outfile);
////	outfile.close();
////
////	while (1)
////	{
////		cout<<"input indices for A and B"<<endl;
////		int a,b;
////		cin>>a>>b;
////		
////		BI::outer_bound ob_a, ob_b;
////		ob_a.construct_outer_bound(bs,a);
////		ob_b.construct_outer_bound(bs,b);
////
////		BI::angle_edge ae_a, ae_b;
////		BI::outer_bound2angle_edge(bs,bs.vcb[a],ob_a,ae_a);
////		BI::outer_bound2angle_edge(bs,bs.vcb[b],ob_b,ae_b);
////
////		//outfile.open("edge_a.txt");
////		//for (int i(0);i!=ae_a.ve.size ();++i)
////		//{
////		//	BI::out_line_3(outfile,ae_a.vp[ae_a.ve[i].e1],ae_a.vp[ae_a.ve[i].e2]);
////		//	outfile<<ae_a.ve[i].vse.size ()<<endl;
////		//	for (int j(0);j!=ae_a.ve[i].vse.size ();++j)
////		//		outfile<<"n1"<<ae_a.ve[i].vse[j].n1<<"| n2"<<ae_a.ve[i].vse[j].n2<<endl;
////		//}
////		//outfile.close();
////		//
////		//outfile.open ("edge_b.txt");
////		//for (int i(0);i!=ae_b.ve.size ();++i)
////		//{
////		//	BI::out_line_3(outfile,ae_b.vp[ae_b.ve[i].e1],ae_b.vp[ae_b.ve[i].e2]);
////		//	outfile<<ae_b.ve[i].vse.size ()<<endl;
////		//	for (int j(0);j!=ae_b.ve[i].vse.size ();++j)
////		//		outfile<<"n1"<<ae_b.ve[i].vse[j].n1<<"| n2"<<ae_b.ve[i].vse[j].n2<<endl;
////		//}
////		//outfile.close();
////
////		cout<<"ob_a.vpoly.size ()"<<ob_a.vpoly.size ()<<endl;
////		cout<<"ob_b.vpoly.size ()"<<ob_b.vpoly.size ()<<endl;
////		cout<<"ae_a.va.size ()"<<ae_a.va.size ()<<endl;
////		cout<<"ae_a.ve.size ()"<<ae_a.ve.size ()<<endl;
////		cout<<"ae_a.vp.size ()"<<ae_a.vp.size ()<<endl;
////		cout<<"ae_b.va.size ()"<<ae_b.va.size ()<<endl;
////		cout<<"ae_b.ve.size ()"<<ae_b.ve.size ()<<endl;
////		cout<<"ae_b.vp.size ()"<<ae_b.vp.size ()<<endl;
////
////
////		
////		outfile.open ("a.stl");
////		ob_a.output_STL(outfile);
////		outfile.close();
////
////		outfile.open ("b.stl");
////		ob_b.output_STL(outfile);
////		outfile.close();
////
////		BI::entrance_block eblk;
////		cout<<"generating entrance block"<<endl;
////		cout<<eblk.gen_entrance_block(bs,bs.vcb[a],bs.vcb[b],bs.vcb[a].get_rand_p(bs.veb,bs.vpo),ob_a,ae_a,ob_b,ae_b,true)<<endl;
////		
////		cout<<"eblk.bs.vpoly.size():"<<eblk.bs.vpoly.size()<<endl;
////		cout<<"eblk.vcover.size():"<<eblk.vcover.size()<<endl;
////		cout<<"eblk.bs.vcb.size():"<<eblk.bs.vcb.size ()<<endl;
////		cout<<"eblk.bs.veb.size():"<<eblk.bs.veb.size ()<<endl;
////		cout<<"eblk.inner.vb.size ()"<<eblk.inner.vb.size ()<<endl;
////		
////		outfile.open ("cover.txt");
////		for (int i(0);i!=eblk.bs.vpoly.size ();++i)
////			BI::out_line_3_rhino(outfile,eblk.bs.vpoly[i]);
////		outfile.close();
////
////		outfile.open("poly.txt");
////		for (int i(0);i!=eblk.bs.vpoly.size ();++i)
////			BI::out_line_3_rhino(outfile,eblk.bs.vpoly[i]);
////		outfile.close();
////		
////		outfile.open("entrblkcblock.stl");
////		for (int i(0);i!=eblk.bs.vcb.size ();++i)
////			eblk.bs.output_cb_stl(i,outfile);
////		outfile.close();
////
////		outfile.open ("entrance_blk.stl");
////		eblk.outer_b.output_STL(outfile);
////		outfile.close();
////
////		ofstream out_bin ("bin", ios::binary);
////		ofstream out_asc ("asc");
////		CGAL::set_binary_mode(out_bin);
////		out_asc.precision(20);
////		
////		out_bin<<eblk;
////		out_asc<<eblk;
////		
////		out_bin.close();
////		out_asc.close();
////		
////		BI::entrance_block eblk_bin,eblk_asc;
////		
////		ifstream in_bin ("bin", ios::binary );
////		ifstream in_asc ("asc");
////		CGAL::set_binary_mode(in_bin);
////		
////		in_bin>>eblk_bin;
////		in_asc>>eblk_asc;
////		
////		cout<<"(eblk_bin == eblk)"<<(eblk_bin == eblk)<<endl;
////		cout<<"(eblk_asc == eblk)"<<(eblk_asc == eblk)<<endl;
////		
////		in_bin.close();
////		in_asc.close();
////		
////		out_asc.open("ob_asc");
////		out_bin.open("ob_bin", ios::binary);
////		CGAL::set_binary_mode(out_bin);
////		
////		out_asc.precision(20);
////		out_asc<<ob_a;
////		out_bin<<ob_a;
////
////		out_asc.close();
////		out_bin.close();
////		
////		in_asc.open ("ob_asc");
////		in_bin.open ("ob_bin", ios::binary);
////		CGAL::set_binary_mode(in_bin);
////		BI::outer_bound ob_a_asc,ob_a_bin;
////		in_asc>>ob_a_asc;
////		in_bin>>ob_a_bin;
////		
////		in_asc.close();
////		in_bin.close();
////		cout<<"ob_a_asc==ob_a:"<<(ob_a_asc==ob_a)<<endl;
////		cout<<"ob_a_bin==ob_a:"<<(ob_a_bin==ob_a)<<endl;
////
////	}
////	
////}
//
////int main ()
////{
////	ifstream infile ("model_domain.dat");
////	ofstream outfile;
////	outfile.open ("time_memory.txt");
////	GB_domain d;
////	d.Input(infile);
////	infile.close();
////	BI::block_system bs;
////	CGAL::Timer timer;
////
////	double threshold(0),ratio_cb(0), ratio_ob(0);
////	cout<<"Enter volume thresholds for cblocks"<<endl;
////	cin>>threshold;
////	cout<<"Enter explode ratio for cblocks"<<endl;
////	cin>>ratio_cb;
////	cout<<"Enter explode ratio for outer boundaries"<<endl;
////	cin>>ratio_ob;
////
////	timer.reset ();
////
////	//reading model_domain
////	timer.start();
////	bs.init_domain(d);
////	if (bs.init_check()==0)
////		cout<<"read modeling_domain.dat success"<<endl;
////	else
////	{
////		cout<<"read modeling_domain.dat fail"<<endl;
////		return 0;
////	}
////	timer.stop();
////	outfile<<"Subdomain number: "<<bs.veb.size ()<<endl;
////	outfile<<"Initalization time: "<<timer.time ()<<endl;
////
////
////	infile.open("contri_fracture_xyzabr.dat");
////	if (!infile)
////	{
////		cout<<"cann't read contri_fracture_xyzabr.dat"<<endl;
////		return 0;
////	}
////
////	if (bs.disc_frac_parser (infile,d.Get_dx_axis()) == 0)
////		cout<<"parser disc frac file success "<<bs.vdisc.size ()<<endl;
////	else
////	{
////		cout<<"parser disc frac file fail"<<endl;
////		return 0;
////	}
////	infile.close();
////	outfile<<"Disc frac number: "<<bs.vdisc.size ()<<endl;
////
////	cout<<"Identifying blocks"<<endl;
////	timer.reset ();timer.start();
////	bs.identify_block();
////	timer.stop ();
////
////	outfile<<"Identification time:"<<timer.time ()<<endl;
////	outfile<<"Element-block count: "<<bs.veb.size ()<<endl;
////	outfile<<"Complex-block count: "<<bs.vcb.size ()<<endl;
////	vector<BI::FT> v_vol;
////	BI::FT total_vol(0);
////	v_vol.reserve (bs.vcb.size ());
////	for (int i(0);i!= (int) bs.vcb.size ();++i)
////	{
////		v_vol.push_back (bs.cb_volume(i));
////		total_vol = total_vol + v_vol.back ();
////	}
////	outfile<<"Model domain volume: "<<total_vol<<endl;
////	outfile<<"Exploded ratio for complex blocks: "<<ratio_cb<<endl;
////	outfile<<"Exploded ratio for outer boundaries: "<<ratio_ob<<endl;
////	for (int i(0);i!= (int) bs.vcb.size ();++i)
////	{
////		cout<<i<<"/"<<bs.vcb.size ()<<endl;;
////		if (v_vol[i]/total_vol < threshold)
////			continue;
////		ostringstream ost;
////		ost<<"cblock"<<i<<".stl";
////		ofstream out_stl(ost.str ());
////		if (bs.output_explode_cb_stl(i,out_stl,ratio_cb) != 0)
////			throw logic_error ("something wrong when outputing exploded cblocks");
////		out_stl.close();
////		ost.str("");
////		ost<<"ob_cblock"<<i<<".stl";
////		out_stl.open (ost.str ());
////		timer.reset (); timer.start();
////		BI::outer_bound ob(bs,i);
////		timer.stop();
////		ob.output_STL(out_stl);
////		out_stl.close ();
////
////		ost.str ("");
////		ost<<"ob_cblock_exp"<<i<<".stl";
////		out_stl.open (ost.str ());
////		ob.output_exp_STL(out_stl,ratio_ob);
////		out_stl.close();
////
////
////		outfile<<"------------Complex block:"<<i<<"------------"<<endl;
////		outfile<<"Element block count: "<<bs.vcb[i].vb.size ()<<endl;
////		outfile<<"OB construction time: "<<timer.time ()<<endl;
////		outfile<<"Volume: "<<v_vol[i]<<"; v_ratio: "<<v_vol[i]/total_vol<<endl;
////		outfile<<"Plane count: "<<ob.vplane.size ()<<endl;
////		outfile<<"Polygon count "<<ob.vpoly.size ()<<endl;
////
////
////	}
////	outfile.close();
////	return 0;
////}
//
////int main ()
////{
////	ifstream infile ("model_domain.dat");
////	GB_domain d;
////	d.Input(infile);
////	infile.close();
////	BI::block_system bs;
////
////	//reading model_domain
////	bs.init_domain(d);
////	if (bs.init_check()==0)
////		cout<<"read modeling_domain.dat success"<<endl;
////	else
////	{
////		cout<<"read modeling_domain.dat fail"<<endl;
////		return 0;
////	}
////
////	infile.open("contri_fracture_xyzabr.dat");
////	if (!infile)
////	{
////		cout<<"cann't read contri_fracture_xyzabr.dat"<<endl;
////		return 0;
////	}
////
////	if (bs.disc_frac_parser (infile,d.Get_dx_axis()) == 0)
////		cout<<"parser disc frac file success "<<bs.vdisc.size ()<<endl;
////	else
////	{
////		cout<<"parser disc frac file fail"<<endl;
////		return 0;
////	}
////	infile.close();
////
////	cout<<"Identifying blocks"<<endl;
////	bs.identify_block();
////
////	ofstream of_bin ("bs_bin", ios::binary);
////	ofstream of_asc ("bs_asc");
////	CGAL::set_binary_mode(of_bin);
////	
////	of_asc.precision(20);
////
////	of_bin<<bs;
////	of_asc<<bs;
////
////	cout<<"of_bin"<<bool(of_bin)<<of_bin.bad()<<of_bin.fail()<<of_bin.eof()<<of_bin.good()<<endl;
////	cout<<"of_asc"<<bool(of_asc)<<of_asc.bad()<<of_asc.fail()<<of_asc.eof()<<of_asc.good()<<endl;
////
////	of_bin.close();
////	of_asc.close();
////	cout<<"output success"<<endl;
////
////	ifstream if_bin ("bs_bin", ios::binary);
////	ifstream if_asc ("bs_asc");
////	CGAL::set_binary_mode(if_bin);
////
////	BI::block_system bs1,bs2;
////
////	cout<<"if_bin_before"<<bool(if_bin)<<(if_bin.eof())<<(if_bin.bad())<<(if_bin.fail())<<(if_bin.good())<<endl;
////	if_bin>>bs1;
////	cout<<"if_asc_before"<<bool(if_asc)<<(if_asc.eof())<<(if_asc.bad())<<(if_asc.fail())<<(if_asc.good())<<endl;
////	if_asc>>bs2;
////	
////	cout<<"read file success"<<endl;
////	cout<<"if_bin"<<bool(if_bin)<<(if_bin.eof())<<(if_bin.bad())<<(if_bin.fail())<<(if_bin.good())<<endl;
////	cout<<"if_asc"<<bool(if_asc)<<(if_asc.eof())<<(if_asc.bad())<<(if_asc.fail())<<(if_asc.good())<<endl;
////
////	cout<<"bs==bs1"<<(bs==bs1)<<endl;
////	cout<<"bs==bs2"<<(bs==bs2)<<endl;
////	if_bin.close();
////	if_asc.close();
////	return 0;
////}
//
////批量识别块体，并且获得各个块体的外表面结构和内表面结构
////int main ()
////{
////	int model_count(0);
////	cout<<"how many model to generate?"<<endl;
////	cin>>model_count;
////	string model_file_name;
////	cout<<"enter the name of modeling_domain"<<endl;
////	cin>>model_file_name;
////	string frac_file_name;
////	cout<<"enter the name of DFN files:"<<endl;
////	cin>>frac_file_name;
////
////	double exp_ratio;
////	cout<<"explode ratio (defult 0.25)"<<endl;
////	cin>>exp_ratio;
////
////	double filter_ratio;
////	cout<<"filter_ratio: (only filter_ratio numbers of blocks are generated. 0 means no block. 1.0 means all the blocks)"<<endl;
////	cin>>filter_ratio;
////
////	ofstream outfile ("time_memory.txt");
////	outfile.precision(20);
////
////	for (int curr(0); curr!=model_count; ++curr)
////	{
////		cout<<"current model: "<<curr<<endl;
////		CGAL::Memory_sizer mem;
////		CGAL::Timer timer;
////		timer.reset ();
////		//read model domain
////		outfile<<"Model "<<curr<<endl;
////		cout<<"initlizing..."<<endl;
////		timer.start();
////		ostringstream ost;
////		ost<<model_file_name<<curr<<".dat";
////		ifstream infile (ost.str ());
////		GB_domain d;
////		d.Input (infile);
////		BI::block_system bs;
////		bs.init_domain(d);
////		timer.stop();
////		outfile<<"fram_count:\t"<<bs.veb .size ()<<endl<<"init_fram_time:\t"<<timer.time ()<<endl;
////
////		if (bs.init_check() == 0)
////			cout<<"read domain file success"<<endl;
////		else
////		{
////			cout<<"read domain file fail"<<endl;
////			outfile<<"read domain file fail"<<endl;
////			continue;
////		}
////		infile.close();
////
////		//read frac file
////		ost.str ("");
////		ost<<frac_file_name<<curr<<".dat";
////		infile.open (ost.str ());
////		timer.reset (); timer.start ();
////		if (bs.disc_frac_parser (infile, d.Get_dx_axis()) == 0)
////		{
////			outfile<<"parser_disc_time:\t"<<timer.time ()<<endl;
////			outfile<<"disc_count:\t\t\t"<<bs.vdisc .size ()<<endl;
////		}
////		else
////		{
////			cout<<"parser disc frac file fail"<<endl;
////			outfile<<"parser disc frac file fail"<<endl;
////			continue;
////		}
////
////		//identify block
////		cout<<"identifying block... ";
////		timer.stop(); timer.reset (); timer.start();
////		bs.identify_block();
////		timer.stop();
////
////		outfile<<"Memory_used:\t\t"<<mem.resident_size()/1024.0/1024.0<<endl;
////		outfile<<"BI_Time:\t\t\t"<<timer.time()<<endl;
////		outfile<<"cutting_time:\t\t"<<bs.cutting_time<<endl;
////		outfile<<"shrk_time:\t\t\t"<<bs.shrk_time<<endl;
////		outfile<<"CB_count:\t\t\t"<<bs.vcb.size ()<<endl;
////		outfile<<"EB_count:\t\t\t"<<bs.veb.size ()<<endl;
////		outfile<<"Point_count:\t\t"<<bs.vpo.size ()<<endl;
////		outfile<<"Plane_count:\t\t"<<bs.vpl.size ()<<endl;
////		cout<<"Time:"<<timer.time()<<endl;
////
////		vector<double> v_vol;		//FT在进行重复性加减时会出错，会出很严重的错误，所以这里只用double
////		v_vol.reserve (bs.vcb.size ());
////		for (int i(0);i!= (int)bs.vcb.size ();++i)
////			v_vol.push_back (bs.cb_volume_est(i));
////
////		vector<double> v_vol_temp(v_vol);
////		sort(v_vol_temp.begin (),v_vol_temp.end());
////
////		BI::FT vol_th = v_vol_temp[int((v_vol.size ()-1) * (1-filter_ratio))];
////		ost.str ("");
////		ost<<model_file_name<<"_"<<frac_file_name<<curr;
////		string dir_name = ost.str();
////		ost.str("");
////		ost<<"mkdir \""<<dir_name<<"\"";
////		system (ost.str ().c_str());
////
////		ost.str ("");
////		ost<<dir_name<<"\\block_system.bin";
////		ofstream bs_bin (ost.str (), ios::binary);
////		CGAL::set_binary_mode(bs_bin);
////		bs_bin<<bs;
////		bs_bin.close();
////
////		ofstream ob_static (dir_name + "\\ob_static.txt");
////
////		ob_static<<"Exploded ratio:"<<exp_ratio<<endl;
////		ob_static<<"Filter_ratio:"<<filter_ratio<<endl;
////		cout<<"Generating Boundaries and outputting"<<endl;
////		for (int i(0);i!= (int)bs.vcb.size ();++i)
////		{
////			if (v_vol[i] < vol_th)
////				continue;
////			ob_static<<"cblock:"<<i<<endl;
////			ob_static<<"vb.size():"<<bs.vcb[i].vb.size ()<<endl;
////
////			ost.str ("");
////			ost<<dir_name<<"\\cblock"<<i<<".stl";
////			ofstream out_cblock (ost.str ());
////			bs.output_cb_stl(i,out_cblock);
////			out_cblock.close();
////
////			ost.str ("");
////			ost<<dir_name<<"\\cblock_exp"<<i<<".stl";
////			out_cblock.open (ost.str ());
////			bs.output_explode_cb_stl(i,out_cblock,exp_ratio);
////			out_cblock.close();
////
////			ost.str ("");
////			ost<<dir_name<<"\\outer_bound"<<i<<".stl";
////			out_cblock.open (ost.str());
////			timer.reset();timer.start ();
////			BI::outer_bound ob;
////			ob.construct_outer_bound(bs,i,false);
////			timer.stop ();
////			ob.output_STL(out_cblock);
////			out_cblock.close();
////			ob_static<<"ob generation time:"<<timer.time ()<<endl;
////			ob_static<<"ob polygon count:\t"<<ob.vpoly.size ()<<endl;
////
////			ost.str ("");
////			ost<<dir_name<<"\\outer_bound"<<i<<".bin";
////			out_cblock.open (ost.str(), ios::binary);
////			CGAL::set_binary_mode(out_cblock);
////			out_cblock<<ob;
////			out_cblock.close();
////
////			//test
////			//ifstream in_cblock(ost.str (), ios::binary);
////			//CGAL::set_binary_mode (in_cblock);
////			//BI::outer_bound ob_temp;
////			//in_cblock>>ob_temp;
////			//in_cblock.close();
////			//cout<<"ob==ob_temp:"<<(ob==ob_temp)<<endl;
////
////
////			ost.str("");
////			ost<<dir_name<<"\\inner_bound"<<i<<".stl";
////			out_cblock.open (ost.str());
////			CGAL::set_ascii_mode(out_cblock);
////			timer.reset ();timer.start();
////			BI::inner_bound ib;
////			ib.construct_inner_bound (bs,i,false);
////			timer.stop();
////			ib.output_STL(out_cblock);
////			out_cblock.close();
////			ob_static<<"ib generation time:"<<timer.time ()<<endl;
////			ob_static<<"ib polygon count: \t"<<ob.vpoly.size ()<<endl;
////
////			ost.str ("");
////			ost<<dir_name<<"\\inner_bound"<<i<<".bin";
////			out_cblock.open (ost.str(), ios::binary);
////			CGAL::set_binary_mode(out_cblock);
////			out_cblock<<ib;
////			out_cblock.close();
////
////			//test
////			//in_cblock.open (ost.str(), ios::binary);
////			//CGAL::set_binary_mode(in_cblock);
////			//BI::inner_bound ib_temp;
////			//in_cblock>>ib_temp;
////			//in_cblock.close();
////			//cout<<"ib==ib_temp:"<<(ib==ib_temp)<<endl;
////
////			ob_static<<endl;
////		}
////		ob_static.close();
////		outfile<<endl;
////	}
////	outfile.close();
////
////	return 0;
////}
//
//////bs2ob_ib.exe
////////block_system convert to outer_bound or inner_bound
//////生成指定块体的外表面或者内表面，只生成指定块体的外表面或者内表面
//////input.txt格式：
//////bs文件名
//////bs文件格式 （1表示是binary，2表示是ascii）
//////要输出几个cblock
//////对于每个cblock：cb的索引   cb输出文件的格式(1:生成内表面，2：生成外表面，3：都生成, 4:生成内表面，但是每个多边形分开输出，5：生成外表面，但是每个多边形分开输出，6：都生成，多边形分开输出)
////int main ()
////{
////	ifstream cmd("bs2ob_ib_input.txt");
////	string filename;	//bs文件名
////	cmd>>filename;
////
////	int format(0);
////	cmd>>format;		//bs文件的格式，1表示是binary，2表示是ascii
////
////	int ncb(0);
////	cmd>>ncb;			//要输出几个cblock
////
////	vector<int> cbi, vflag;		//cbi记录各个cblock的索引，vflag记录为每个cblock生成什么信息(1:生成内表面，2：生成外表面，3：都生成, 4:生成内表面，但是每个多边形分开输出，5：生成外表面，但是每个多边形分开输出，6：都生成，多边形分开输出)
////	cbi.resize (ncb);
////	vflag.resize (ncb);
////	for (int i(0);i!=ncb;++i)
////		cmd>>(cbi[i])>>(vflag[i]);
////
////	cmd.close();
////	ifstream in_bs;
////	if (format == 1)
////	{
////		in_bs.open(filename, ios::binary);
////		CGAL::set_binary_mode(in_bs);
////	}
////	else 
////		in_bs.open (filename);
////
////	if (!in_bs)
////	{
////		cout<<"Open "<<filename<<" error"<<endl;
////		return 0;
////	}
////	else
////		cout<<"Open "<<filename<<" success"<<endl;
////
////	BI::block_system bs;
////	in_bs>>bs;
////	in_bs.close();
////	cout<<"Read bs complete"<<endl;
////
////	for (int i(0);i!=ncb;++i)
////	{
////		cout<<"cb:"<<cbi[i];
////		if (cbi[i] <0 || cbi[i]>= (int)bs.vcb.size ())
////		{
////			cout<<" out of range"<<endl;
////			continue;
////		}
////		else
////			cout<<" outputing"<<endl;
////
////		ostringstream ost;
////		ofstream out_ib, out_ob;
////		BI::inner_bound ib;
////		BI::outer_bound ob;
////		switch (vflag[i])
////		{
////		case 3:
////		case 1:
////			ib.construct_inner_bound(bs,cbi[i]);		//输出内表面
////
////			ost.str("");
////			ost<<"inner_bound"<<cbi[i]<<".stl";
////			out_ib.open (ost.str());
////			CGAL::set_ascii_mode(out_ib);
////			ib.output_STL(out_ib);
////			out_ib.close();
////
////			cout<<"Inner_bound output complete"<<endl;
////			if (vflag[i] == 1)
////				break;
////		case 2:
////			ob.construct_outer_bound (bs,cbi[i]);		//输出外表面
////
////			ost.str ("");
////			ost<<"outer_bound"<<cbi[i]<<".stl";
////			out_ob.open (ost.str ());
////			CGAL::set_ascii_mode(out_ob);
////			ob.output_STL(out_ob);
////			out_ob.close();
////
////			cout<<"Outer_bound output complete"<<endl;
////			break;
////		case 6:
////		case 4:
////			ib.construct_inner_bound(bs,cbi[i]);		//输出内表面
////
////			ost.str("");
////			ost<<"inner_bound_sept"<<cbi[i]<<".stl";
////			out_ib.open (ost.str());
////			CGAL::set_ascii_mode(out_ib);
////			ib.output_exp_STL(out_ib, 0.0);
////			out_ib.close();
////
////			cout<<"Inner_bound output complete"<<endl;
////			if (vflag[i] == 4)
////				break;
////		case 5:
////			ob.construct_outer_bound (bs,cbi[i]);		//输出外表面
////
////			ost.str ("");
////			ost<<"outer_bound_sept"<<cbi[i]<<".stl";
////			out_ob.open (ost.str ());
////			CGAL::set_ascii_mode(out_ob);
////			ob.output_exp_STL(out_ob,0.0);
////			out_ob.close();
////
////			cout<<"Outer_bound output complete"<<endl;
////			break;
////		default:
////			cout<<"Wrong commond"<<endl;
////			break;
////		}
////	}
////
////	cout<<"Generate complete"<<endl;
////	return 0;
////}
//
//
//////bs2all_exp_ob_ib.exe
//////读取块体系统，生成并输出外表面边界和内表面边界，输出所有块体的表面边界
//////input.txt格式：
//////bs文件名
//////bs文件格式 （1表示是binary，2表示是ascii）
//////爆炸系数   输出文件格式(1:生成内表面，2：生成外表面，3：生成整体内表面和外表面， 4：生成的内表面中将各个多边形分开，5：生成的外表面中将各个多边形分开，6：生成的内表面和外表面中将各个多边形分开)
////int main()
////{
////	ifstream cmd ("bs2all_exp_ob_ib_input.txt");
////	string filename;
////	cmd>>filename;
////
////	int format (0);		//bs文件的格式，1表示是binary，2表示是ascii
////	cmd>>format;
////
////	double exp_ratio(0);
////	cmd>>exp_ratio;		//爆炸系数
////
////	int out_flag(0);			//out_flag记录为每个cblock生成什么信息(1:生成内表面，2：生成外表面，3：生成整体内表面和外表面
////								//									 4：生成的内表面中将各个多边形分开，5：生成的外表面中将各个多边形分开，6：生成的内表面和外表面中将各个多边形分开)
////	cmd>>out_flag;
////
////	cmd.close();
////	ifstream in_bs;
////	if(format ==1)
////	{
////		in_bs.open (filename, ios::binary );
////		CGAL::set_binary_mode (in_bs);
////	}
////	else
////		in_bs.open (filename);
////
////	if (!in_bs)
////	{
////		cout<<"Open "<<filename<<" error"<<endl;
////		return 0;
////	}
////	else
////		cout<<"Open "<<filename<<" success"<<endl;
////
////	BI::block_system bs;
////	in_bs>>bs;
////	in_bs.close();
////	cout<<"Read bs complete"<<endl;
////
////	ofstream out_ob, out_ib;
////
////	vector<BI::point_3> cb_p;
////	BI::point_3 bs_center(bs.get_rand_p());
////
////	cb_p.reserve (bs.vcb.size ());
////	for (int i(0);i!=bs.vcb.size ();++i)
////		cb_p.push_back (bs.vcb[i].get_rand_p(bs.veb,bs.vpo));
////
////	switch (out_flag)
////	{
////	case 3:
////	case 1:
////		out_ib.open ("ib_exp.stl");				//生成内表面
////		for (int i(0);i!=bs.vcb.size ();++i)
////		{
////			stringstream ss;
////			ss<<"IB"<<i;
////			BI::inner_bound ib;
////			ib.construct_inner_bound(bs,i);
////			ib.output_STL(out_ib,(cb_p[i].x()-bs_center.x())*exp_ratio, (cb_p[i].y()-bs_center.y())*exp_ratio, (cb_p[i].z()-bs_center.z())*exp_ratio,ss.str ());
////		}
////		cout<<"IB generation complete"<<endl;
////		if (out_flag == 1)
////			break;
////	case 2:
////		out_ob.open ("ob_exp.stl");			//生成外表面
////		for (int i(0);i!=bs.vcb.size ();++i)
////		{
////			stringstream ss;
////			ss<<"OB"<<i;
////			BI::outer_bound ob;
////			ob.construct_outer_bound (bs,i);
////			ob.output_STL(out_ob, (cb_p[i].x()-bs_center.x())*exp_ratio, (cb_p[i].y()-bs_center.y())*exp_ratio, (cb_p[i].z()-bs_center.z())*exp_ratio,ss.str ());
////		}
////		cout<<"OB generation complete"<<endl;
////		break;
////	case 6:
////	case 4:
////		out_ib.open ("ib_exp_polysepted.stl");				//生成内表面
////		for (int i(0);i!=bs.vcb.size ();++i)
////		{
////			stringstream ss;
////			ss<<"IB"<<i;
////			BI::inner_bound ib;
////			ib.construct_inner_bound(bs,i);
////			ib.output_exp_STL(out_ib, 0.0, (cb_p[i].x()-bs_center.x())*exp_ratio, (cb_p[i].y()-bs_center.y())*exp_ratio, (cb_p[i].z()-bs_center.z())*exp_ratio,ss.str());
////		}
////		cout<<"IB generation complete poly seperated"<<endl;
////		if (out_flag == 4)
////			break;
////	case 5:
////		out_ob.open ("ob_exp_polysepted.stl");			//生成外表面
////		for (int i(0);i!=bs.vcb.size ();++i)
////		{
////			stringstream ss;
////			ss<<"OB"<<i;
////			BI::outer_bound ob;
////			ob.construct_outer_bound (bs,i);
////			ob.output_exp_STL(out_ob, 0.0, (cb_p[i].x()-bs_center.x())*exp_ratio, (cb_p[i].y()-bs_center.y())*exp_ratio, (cb_p[i].z()-bs_center.z())*exp_ratio,ss.str());
////		}
////		cout<<"OB generation complete poly seperated"<<endl;
////		break;
////	default:
////		cout<<"successfully end"<<endl;
////		break;
////	}
////	out_ib.close();
////	out_ob.close();
////	cout<<"boundary structure generation complete"<<endl;
////	return 0;
////}
//
////bs2exp_cb.exe
////输出指定cblock的爆炸视图,也就是把指定的cblock中的单元块体以爆炸视图展示出来
////input.txt格式：
////bs文件名
////bs文件格式 （1表示是binary，2表示是ascii）
////要输出几个cblock
//////对于每个cblock： cb的索引   cb的爆炸系数
////int main ()
////{
////	ifstream cmd ("input.txt");
////	string filename;		//bs文件名
////	cmd>>filename;
////
////	int format (0);			//bs文件的格式，1表示是binary，2表示是ascii
////	cmd>>format;
////
////	int ncb(0);
////	cmd>>ncb;				//要输出几个cblock
////
////	vector<int> cbi;
////	vector<double> vratio;
////	cbi.reserve (ncb);
////	vratio.reserve (ncb);
////
////	for (int i(0);i!=ncb;++i)
////	{
////		if (!cmd)
////		{
////			cout<<"Input file has wrong format"<<endl;
////			return 0;
////		}
////		cmd>>(cbi[i]) >> (vratio[i]);
////	}
////	cmd.close();
////	ifstream in_bs;
////	if (format == 1)
////	{
////		in_bs.open (filename,ios::binary );
////		CGAL::set_binary_mode(in_bs);
////	}
////	else
////	{
////		in_bs.open (filename);
////	}
////
////	if (!in_bs)
////	{
////		cout<<"Open "<<filename<<" error"<<endl;
////		return 0;
////	}
////	else 
////		cout<<"Open "<<filename<<" success"<<endl;
////
////	BI::block_system bs;
////	in_bs>>bs;
////	in_bs.close();
////
////	for (int i(0);i!=ncb;++i)
////	{
////		cout<<"cb:"<<cbi[i];
////		if (cbi[i] <0 || cbi[i]>= (int)bs.vcb.size ())
////		{
////			cout<<" out of range"<<endl;
////			continue;
////		}
////		else
////			cout<<" outputing"<<endl;
////		if (vratio[i]<0)
////		{
////			cout<<"given ratio <0"<<endl;
////			continue;
////		}
////		ostringstream ost;
////		ost.str("");
////		ost<<"cblock_exp_"<<cbi[i]<<".stl";
////		ofstream out_cb(ost.str());
////		bs.output_explode_cb_stl(cbi[i],out_cb,vratio[i]);
////		out_cb.close();
////	}
////	cout<<"Output cblock complete"<<endl;
////	return 0;
////}
//
//////bs2entrblk.exe
//////输出指定cblock的entrance block。
//////input.txt格式：
//////bs文件名
//////bs文件格式 （1表示是binary，2表示是ascii）
//////要生成几个entrance block
//////对于每个进入块体：a的索引   b的索引   输出的stl文件格式（0:不输出，1:输出外表面，2:输出外表面，但是每个多边形都分开表示）   如何序列化输出生成的进入块体（0：不输出，1：输出ascii格式的文件，2：输出二进制格式的文件）
////int main()
////{
////	ifstream cmd ("input.txt");
////	string filename;		//bs文件名
////	cmd>>filename;
////
////	int format (0);			//bs文件的格式，1表示是binary，2表示是ascii
////	cmd>>format;
////
////	int neb(0);
////	cmd>>neb;				//要输出几个entrance block
////
////	vector<int> ai,bi,vflag, vflag2;
////	ai.resize (neb);		//ai是块体a的索引
////	bi.resize (neb);		//bi是块体b的索引
////	vflag.resize (neb);		//vflag表明如何输出entrance block的stl文件（0:不输出，1:输出外表面，2:输出外表面，但是每个多边形都分开表示）
////	vflag2.reserve (neb);	//vflag2表示是否序列化输出entrance block （0：不输出，1：输出ascii格式的文件，2：输出二进制格式的文件）
////
////	for (int i(0);i!=neb;++i)
////	{
////		if (!cmd)
////		{
////			cout<<"input file have wrong format"<<endl;
////			return 0;
////		}
////		cmd>>ai[i]>>bi[i]>>vflag[i]>>vflag2[i];
////	}
////	
////	cmd.close ();		//读取输出部分完毕
////
////	ifstream in_bs;
////	if (format == 1)
////	{
////		in_bs.open(filename, ios::binary);
////		CGAL::set_binary_mode(in_bs);
////	}
////	else 
////		in_bs.open (filename);
////
////	if (!in_bs)
////	{
////		cout<<"Open "<<filename<<" error"<<endl;
////		return 0;
////	}
////	else
////		cout<<"Open "<<filename<<" success"<<endl;
////
////	BI::block_system bs;
////	in_bs>>bs;
////	in_bs.close();
////	cout<<"Read bs complete"<<endl;
////
////	for (int i(0);i!=neb;++i)
////	{
////		cout<<"a:"<<ai[i]<<" b:"<<bi[i];
////		if (ai[i] < 0 || ai[i]>= (int)bs.vcb.size () || bi[i]<0 || bi[i]>=(int)bs.vcb.size ())
////		{
////			cout<<" index out of range"<<endl;
////			continue;
////		}
////		else
////			cout<<" outputing"<<endl;
////
////		ostringstream ost;
////		BI::entrance_block eb;
////		BI::point_3 a0 (bs.vcb[ai[i]].get_rand_p(bs.veb,bs.vpo));
////		eb.gen_entrance_block(bs,bs.vcb[ai[i]], bs.vcb[bi[i]],a0,true);
////		ost.str ("");
////		ost<<"etrblk_a"<<ai[i]<<"_b"<<bi[i]<<"_info.txt";
////		ofstream info(ost.str ());
////		
////		info<<"**********Entrance Block***********"<<endl;
////		info<<"a0: ("<<CGAL::to_double(a0.x())<<' '<<CGAL::to_double(a0.y())<<' '<<CGAL::to_double(a0.z())<<")"<<endl;
////		info<<"Cover count: "<<eb.vcover .size ()<<endl;
////		info.close ();
////
////		ost.str ("");
////		ost<<"etrblk_a"<<ai[i]<<"_b"<<bi[i]<<".stl";
////		ofstream out_stl;
////		ofstream out_bin;
////		switch (vflag[i])
////		{
////		case 0:
////			break;
////		case 1:
////			out_stl.open (ost.str());
////			eb.outer_b.output_STL(out_stl);
////			break;
////		case 2:
////			out_stl.open (ost.str());
////			eb.outer_b.output_exp_STL(out_stl,0);
////			break;
////		default:
////			cout<<"wrong command"<<endl;
////			break;
////		}
////		ost.str ("");
////		ost<<"etrblk_a"<<ai[i]<<"_b"<<bi[i]<<".bin";
////		switch (vflag2[i])
////		{
////		case 0:
////			break;
////		case 1:
////			out_bin.open (ost.str());
////			out_bin<<eb;
////			break;
////		case 2:
////			out_bin.open (ost.str(), ios::binary);
////			CGAL::set_binary_mode(out_bin);
////			out_bin<<eb;
////			break;
////		default:
////			cout<<"wrong command"<<endl;
////			break;
////		}
////		out_stl.close();
////		out_bin.close();
////	}
////	cout<<"Generate complete"<<endl;
////	return 0;
////}
//
//////entrblk2ob_stl.exe
//////读取一个进入块体，输出这个进入块体的外表面
//////input.txt格式：
//////要读取的进入块体文件数
////////对于每个进入块体：文件名  输入文件格式（1：bin格式，2：ascii格式）   输出文件格式（1：输出一个整体，2：每个面分别输出）
////int main ()
////{
////	ifstream cmd ("input.txt");
////	int neb (0);		//表明要读取几个进入块体文件
////	cmd>>neb;
////
////	vector<string> v_filename;
////	vector<int> vflag1,vflag2;
////	v_filename.resize (neb);		//记录每个进入块体文件的文件名
////	vflag1.resize (neb);			//这个块体文件是什么格式（1：bin格式，2：ascii格式）
////	vflag2.resize (neb);			//以什么方式输出stl文件（1：输出一个整体，2：每个面分别输出）
////
////	for (int i(0);i!=neb;++i)
////		cmd>>v_filename[i]>>vflag1[i]>>vflag2[i];
////
////	for (int i(0);i!=neb;++i)
////	{
////		ifstream in_eb;
////		if (vflag1[i] == 1)
////		{
////			in_eb.open(v_filename[i], ios::binary);
////			CGAL::set_binary_mode(in_eb);
////		}
////		else 
////			in_eb.open (v_filename[i]);
////
////		if (!in_eb)
////		{
////			cout<<"Open "<<v_filename[i]<<" error"<<endl;
////			continue;
////		}
////		else
////			cout<<"Open "<<v_filename[i]<<" success"<<endl;
////
////		BI::entrance_block eb;
////		in_eb>>eb;
////		in_eb.close();
////		ofstream out_stl (v_filename[i]+"_STL.stl");
////		switch (vflag2[i])
////		{
////		case 1:
////			eb.outer_b.output_STL(out_stl);
////			break;
////		case 2:
////			eb.outer_b.output_exp_STL(out_stl,0.0);
////			break;
////		default:
////			cout<<"wrong command"<<endl;
////			break;
////		}
////		out_stl.close();
////	}
////
////	return 0;
////}
//
//////gen_bs.exe
//////识别块体
//////将stl所表示的曲面中的每个三角形面片都视为一个多边形裂隙
//////input.txt格式：
//////要识别块体的模型数
////////对于每个模型：model_domain文件名   几个圆盘裂隙文件   每个圆盘裂隙文件的名字...   几个stl文件   每个stl文件的名字   输出文件的文件名   输出文件的格式（1：binary，2：ascii）
////int main ()
////{
////	ifstream cmd ("gen_bs_input.txt");
////	int nmodel(0);
////	cmd>>nmodel;
////	for (int imodel(0);imodel!=nmodel;++imodel)
////	{
////		cout<<"-------"<<imodel<<"-----------"<<endl;
////		int nstl(0), ndisc(0);;
////		if (!cmd)
////		{
////			cout<<"can not find gen_bs_input.txt"<<endl;
////			return 0;
////		}
////		string domain_name;
////		cmd>>domain_name;		//model_domain文件的文件名
////		cmd>>ndisc;			//表明要读取几个圆盘裂隙文件
////		vector<string> vdisc_fn;
////		vdisc_fn.resize (ndisc);		//每个圆盘裂隙文件的名字
////		for (int i(0);i!=ndisc;++i)
////			cmd>>vdisc_fn[i];
////
////
////		cmd>>nstl;		//表明要读取几个stl文件
////		vector<string> vpoly_fn;
////		vpoly_fn.resize (nstl);		//每个stl文件的名字
////		for (int i(0);i!=nstl;++i)
////			cmd>>vpoly_fn[i];
////
////		string output_name;
////		cmd>>output_name;		//指定输出文件的文件名
////		int output_format;
////		cmd>>output_format;		//指定输出文件的格式（1：binary，2：ascii）
////
////		cmd.close();
////
////		ifstream infile (domain_name);
////		if (!infile)
////		{
////			cout<<"Open domain file "<<domain_name<<" error"<<endl;
////			return 0;
////		}
////		else 
////			cout<<"Open domain file "<<domain_name<<" success"<<endl;
////
////		GB_domain d;
////		d.Input(infile);
////		infile.close();
////		BI::block_system bs;
////		bs.init_domain(d);
////
////		for (int i(0);i!=ndisc;++i)
////		{
////			infile.open (vdisc_fn[i]);
////			if (!infile)
////			{
////				cout<<"Open disc file "<<vdisc_fn[i]<<" error"<<endl;
////				continue;
////			}
////			else 
////				cout<<"Open disc file "<<vdisc_fn[i]<<" success"<<endl;
////
////			if (bs.disc_frac_parser(infile,d.Get_dx_axis()) ==0)
////				cout<<"Parse "<<vdisc_fn[i]<<" success"<<endl;
////			else
////			{
////				cout<<"Parse "<<vdisc_fn[i]<<" fail"<<endl;
////				continue;
////			}
////			infile.close();
////		}
////
////		//for (int i(0);i!=nstl;++i)
////		//{
////		//	infile.open (vpoly_fn[i]);
////		//	if (!infile)
////		//	{
////		//		cout<<"Open poly file "<<vpoly_fn[i]<<" error"<<endl;
////		//		continue;
////		//	}
////		//	else
////		//		cout<<"Open poly file "<<vpoly_fn[i]<<" success"<<endl;
////		//	if (bs.stl_parser(infile) == 0)
////		//		cout<<"Parse "<<vpoly_fn[i]<<" success"<<endl;
////		//	else
////		//	{
////		//		cout<<"Parse "<<vpoly_fn[i]<<" fail"<<endl;
////		//		continue;
////		//	}
////		//	infile.close();
////		//}
////
////		for (int i(0);i!=nstl;++i)
////		{
////			infile.open (vpoly_fn[i]);
////			if (!infile)
////			{
////				cout<<"Open poly file "<<vpoly_fn[i]<<" error"<<endl;
////				continue;
////			}
////			else
////				cout<<"Open poly file "<<vpoly_fn[i]<<" success"<<endl;
////
////			if (bs.poly_frac_parser(infile) == 0)
////				cout<<"Parse "<<vpoly_fn[i]<<" success"<<endl;
////			else
////			{
////				cout<<"Parse "<<vpoly_fn[i]<<" fail"<<endl;
////				continue;
////			}
////			infile.close();
////		}
////
////		cout<<"Read frac file complete...identifying"<<endl;
////		bs.identify_block (cout);
////
////		cout<<"Outputing"<<endl;
////
////		if (output_format == 1)
////		{
////			ofstream outfile (output_name, ios::binary );
////			CGAL::set_binary_mode(outfile);
////			outfile<<bs;
////			outfile.close();
////		}
////		else 
////		{
////			ofstream outfile (output_name);
////			outfile<<bs;
////			outfile.close();
////		}
////		ofstream outfile (output_name+"_info.txt");
////		bs.output_info(outfile);
////		outfile.close();
////		cout<<"Models output complete"<<endl;
////	}
////	return 0;
////}
//
//////cf_analy.exe
//////进行连锁垮落分析的部分，通过读取bs.bin进行分析，如果读取的块体系统文件没有初始化过稳定性分析部分的话会主动初始化并更新原来的bs.bin文件。
//////并且会读取"fixed_surfaces.stl"文件，并将所有与这个曲面相交的块体定义为固定块体
//////"cf_analy_input.txt"格式：
//////block_system文件名，bs文件格式，fixed_surface文件名
////int main ()
////{
////	ifstream cmd ("cf_analy_input.txt");
////	string bs_file_name;
////	cmd>>bs_file_name;
////
////	int format(0);
////	cmd>>format;
////
////	string fs_file_name;
////	cmd>>fs_file_name;
////
////	BI::block_system bs;
////	ifstream infile;
////	if (format == 1)
////	{
////		infile.open (bs_file_name,ios::binary);
////		CGAL::set_binary_mode(infile);
////	}
////	else
////		infile.open (bs_file_name);
////
////	if (!infile)
////	{
////		cout<<"Open "<<bs_file_name<<" error"<<endl;
////		return 0;
////	}
////	else
////		cout<<"Open "<<bs_file_name<<" success"<<endl;
////
////	infile>>bs;
////	if(infile || infile.eof())
////		cout<<"read bs complete"<<endl;
////	else 
////	{
////		cout<<"read bs error"<<endl;
////		return 0;
////	}
////	infile.close();
////
////	infile.open (fs_file_name);
////	BI::fixed_surfaces fs;
////	int flag = fs.STL_parser(infile);
////	cout<<flag<<endl;
////	if (flag ==0)
////		cout<<"Read fixed surfaces file complete"<<endl;
////	else
////		cout<<"Read fixed_surfaces.stl error"<<endl;
////	infile.close();
////
////	ofstream outfile;
////
////
////	if (bs.vbsf.empty())
////	{
////		cout<<"Initializing BS for stability analysis"<<endl;
////		cout<<bs.init_stabl_alys(true)<<endl;
////		cout<<"Initialize complete"<<endl;
////		ofstream outfile ("block_system.bin",ios::binary);
////		CGAL::set_binary_mode(outfile);
////		outfile<<bs;
////		outfile.close();
////	}
////	else
////	{
////		cout<<"Already initilaized"<<endl;
////	}
////
////	cout<<"constructing adjacent graph"<<endl;
////	BI::adj_blk_graph g(bs,fs);
////
////	
////
////	outfile.open ("bs_info.txt");
////	bs.output_info(outfile);
////	g.output_node_and_edge_info(outfile);
////	outfile.close ();
////
////	outfile.open ("node_ori.csv");
////	g.output_CSV(outfile,BI::adj_blk_graph::ELE_TYPE::Node);
////	outfile.close();
////	outfile.open ("edge_ori.csv");
////	g.output_CSV(outfile, BI::adj_blk_graph::ELE_TYPE::Edge );
////	outfile.close();
////	
////
////	cout<<"merging nested blocks"<<endl;
////	cout<<g.eliminate_nested_blk(true)<<endl;
////
////	cout<<"merging interlocked blocks"<<endl;
////	cout<<g.eliminate_interlocked_blk(bs,1,true)<<endl;
////	
////	ofstream out_excav ("excav_unfixed.stl"), out_fixed_and_unexp("unexp_fixed.stl"), out_unexp("unexp_unfixed.stl"), out_excav_and_fixed("excav_fixed.stl");
////	cout<<"outputing all kinds of blks"<<endl;
////	cout<<g.output_blk_group_quick(out_excav,bs,BI::adj_blk_graph::GROUP_TYPE::EXCAV);		//临时查看块体的形状，如果想要输出其外表面多边形还需要再进行计算
////	cout<<g.output_blk_group_quick(out_fixed_and_unexp,bs,BI::adj_blk_graph::GROUP_TYPE::FIXED_AND_UNEXPOSED);
////	cout<<g.output_blk_group_quick(out_unexp,bs,BI::adj_blk_graph::GROUP_TYPE::UNDEXPOSED);
////	cout<<g.output_blk_group_quick(out_excav_and_fixed,bs,BI::adj_blk_graph::GROUP_TYPE::FIXED_AND_EXCAV)<<endl;
////	out_excav.close ();out_fixed_and_unexp.close(); out_unexp.close (); out_excav_and_fixed.close ();
////	cout<<"output blks complete"<<endl;
////
////	outfile.open ("removable_blk.stl");
////	g.ouput_removable_blks(bs,BI::vector_3(0,0,-1),outfile);
////	outfile.close();
////
////	outfile.open ("safety_factor.csv");
////	g.output_safety_factor(outfile, bs);
////	outfile.close();
////
////	outfile.open ("node_merged.csv");
////	g.output_CSV(outfile,BI::adj_blk_graph::ELE_TYPE::Node);
////	outfile.close ();
////	outfile.open ("edge_merged.csv");
////	g.output_CSV(outfile, BI::adj_blk_graph::ELE_TYPE::Edge );
////	outfile.close();
////
////	cout<<"output CSV file complete"<<endl;
////	
////	cout<<"determining chain falling blocks"<<endl;
////	cout<<g.determine_chain_falling_blocks(bs,BI::vector_3(0,0,-1),true)<<endl;
////
////	outfile.open("chain_falling.csv");
////	g.output_chain_falling_info_CSV(bs,outfile);
////	outfile.close();
////
////	outfile.open ("chain_falling_blk.stl");
////	g.output_cfalling_blks(bs,outfile);
////	outfile.close();
////
////	cout<<"Outputing chain falling blocks"<<endl;
////
////	while (1)
////	{
////		cout<<"Enter id:";
////		int temp_id;
////		cin>>temp_id;
////		stringstream ss;
////		ss<<"Cfalling"<<temp_id<<".stl";
////		outfile.open (ss.str ());
////		cout<<g.output_cfalling_blks(bs,outfile,temp_id)<<endl;
////		outfile.close();
////	}
////	return 0;
////}
//
//////convert_bs_format.exe
////将块体系统直接输出，cb由eb表示
////int main ()
////{
////	ifstream cmd ("convert_bs_format_input.txt");
////	string bs_file_name;
////	cmd>>bs_file_name;
////
////	int format (0);
////	cmd>>format;
////
////	string outfile_name;
////	cmd>>outfile_name;
////
////	BI::block_system bs;
////	ifstream infile;
////	if (format == 1)
////	{
////		infile.open (bs_file_name,ios::binary );
////		CGAL::set_binary_mode(infile);
////	}
////	else
////		infile.open (bs_file_name);
////
////	if (!infile)
////	{
////		cout<<"Open "<<bs_file_name<<" error"<<endl;
////		return 0;
////	}
////	else
////		cout<<"Open "<<bs_file_name<<" success"<<endl;
////
////	infile>>bs;
////	if(infile || infile.eof())
////		cout<<"read bs complete"<<endl;
////	else 
////	{
////		cout<<"read bs error"<<endl;
////		return 0;
////	}
////	infile.close();
////
////	using CGAL::to_double;
////	ofstream outfile(outfile_name);
////	outfile.precision(15);
////
////
////
////	outfile<<bs.veb.size ()<<endl;
////	for (int i(0);i!=bs.veb.size ();++i)
////	{
////		map<int,int> m;
////		vector<int> temp_vp;
////		temp_vp.reserve (bs.veb[i].vf.size ()+20);
////		for (int j(0);j!=bs.veb[i].vf.size ();++j)
////			for (int k(0);k!=bs.veb[i].vf[j].np.size ();++k)
////				m[bs.veb[i].vf[j].np[k]] = BI::add_set (temp_vp, bs.veb[i].vf[j].np[k]);
////
////		outfile<<temp_vp.size ()<<endl;
////		for (int j(0);j!=temp_vp.size ();++j)
////			outfile<<to_double (bs.vpo[temp_vp[j]].x())<<" "<<to_double (bs.vpo[temp_vp[j]].y())<<" "<<to_double(bs.vpo[temp_vp[j]].z())<<endl;
////
////		outfile<<bs.veb[i].vf.size ()<<endl;
////		for (int j(0);j!=bs.veb[i].vf.size ();++j)
////		{
////			outfile<<bs.veb[i].vf[j].np.size ()<<" ";
////			for (int k(0);k!=bs.veb[i].vf[j].np.size ();++k)
////				outfile<<m.at (bs.veb[i].vf[j].np[k])<<" ";
////			outfile<<endl;
////		}
////	}
////
////	outfile<<endl<<bs.vcb.size ()<<endl;
////	for (int i(0);i!=bs.vcb.size ();++i)
////	{
////		outfile<<bs.vcb[i].vb.size ()<<" ";
////		for (int j(0);j!=bs.vcb[i].vb.size ();++j)
////			outfile<<bs.vcb[i].vb[j]<<" ";
////		outfile<<endl;
////	}
////
////	outfile.close();
////	cout<<"Output complete"<<endl;
////	return 0;
////}
//
////bs2outer_bound.exe
////将块体系统输出成外表面结构
////int main ()
////{
////	ifstream cmd ("bs2outer_bound_input.txt");
////	string bs_file_name;
////	cmd>>bs_file_name;
////
////	int format (0);
////	cmd>>format;
////
////	string outfile_name;
////	cmd>>outfile_name;
////
////	BI::block_system bs;
////	ifstream infile;
////	if (format == 1)
////	{
////		infile.open (bs_file_name,ios::binary );
////		CGAL::set_binary_mode(infile);
////	}
////	else
////		infile.open (bs_file_name);
////
////	if (!infile)
////	{
////		cout<<"Open "<<bs_file_name<<" error"<<endl;
////		return 0;
////	}
////	else
////		cout<<"Open "<<bs_file_name<<" success"<<endl;
////
////	infile>>bs;
////	if(infile || infile.eof())
////		cout<<"read bs complete"<<endl;
////	else 
////	{
////		cout<<"read bs error"<<endl;
////		return 0;
////	}
////	infile.close();
////
////	using CGAL::to_double;
////	ofstream outfile(outfile_name);
////	outfile.precision(15);
////
////	outfile<<bs.vcb.size ()<<endl;
////	for (int i(0);i!=bs.vcb.size ();++i)
////	{
////		BI::outer_bound ob (bs,i);
////		vector<BI::point_3> vp;
////		typedef vector<int> index_poly;
////		typedef vector<index_poly> index_pwh;
////		vector<index_pwh> vi(ob.vpoly.size ());
////
////		for (int j(0);j!=ob.vpoly.size ();++j)
////		{
////			const BI::polygon_2 & poly (ob.vpoly[j].outer_boundary ());
////			vi[j].push_back (index_poly());
////			for (BI::polygon_2::Vertex_const_iterator k = poly.vertices_begin();
////				k!=poly.vertices_end();++k)
////				vi[j].back ().push_back (BI::add_set(vp,ob.vplane[j].to_3d(*k)));
////
////			for (BI::polygon_with_holes_2::Hole_const_iterator hi = ob.vpoly[j].holes_begin();
////				hi!=ob.vpoly[j].holes_end();++hi)
////			{
////				vi[j].push_back (index_poly());
////				for (BI::polygon_2::Vertex_const_iterator k = hi->vertices_begin();
////					k!=hi->vertices_end();++k)
////					vi[j].back ().push_back (BI::add_set(vp,ob.vplane[j].to_3d(*k)));
////			}
////			if (!ob.outer_normal[j])
////				for (int k(0);k!=vi[j].size ();++k)
////					reverse(vi[j][k].begin (), vi[j][k].end ());
////		}
////
////		outfile<<vp.size ()<<endl;
////		for (int j(0);j!=vp.size ();++j)
////			outfile<<to_double (vp[j].x()) <<" "<<to_double (vp[j].y())<<" "<<to_double (vp[j].z())<<endl;
////		for (int j(0);j!=vp.size ();++j)
////		{
////			//输出每个角点由哪个单元块体的哪个平面组成
////			vector<pair<int,vector<int>>> angle_struct;
////			for (int k(0);k!=bs.vcb[i].vb.size ();++k)
////			{
////				pair<int,vector<int>> eb_index;		//寻找下vp[j]这个点在eb的哪些面上
////				eb_index.second .reserve (bs.veb[bs.vcb[i].vb[k]].vf.size ());
////				eb_index.first = bs.vcb[i].vb[k];
////				switch (BI::point_eblock_location (bs,bs.veb[bs.vcb[i].vb[k]],vp[j],eb_index.second))
////				{
////				case BI::Oriented_side::IN_:
////					throw logic_error ("point in eblock");
////					break;
////				case BI::Oriented_side::ON_:
////					angle_struct.push_back (eb_index);
////				default:
////					break;
////				}
////			}
////			outfile<<angle_struct.size ()<<endl;
////			for (int k(0);k!=angle_struct.size ();++k)
////			{
////				outfile<<(angle_struct[k].first) <<" "<<angle_struct[k].second.size ()<<" ";
////				for (int l(0);l!=angle_struct[k].second.size ();++l)
////					outfile<<angle_struct[k].second [l]<<" ";
////				outfile<<endl;
////			}
////		}
////		
////		outfile<<vi.size ()<<endl;
////		for (int j(0);j!=vi.size ();++j)
////		{
////			outfile<<vi[j].size ()<<endl;
////			for (int k(0);k!=vi[j].size ();++k)
////			{
////				outfile<<vi[j][k].size ()<<" ";
////				for (int l(0);l!=vi[j][k].size ();++l)
////					outfile<<vi[j][k][l]<<" ";
////				outfile<<endl;
////			}
////		}
////		outfile<<endl;
////	}
////
////	cout<<"Output complete"<<endl;
////	return 0;
////}
//
////bs2baniformat
////int main()
////{
////	ifstream cmd ("bs2baniformat.txt");
////	string bs_file_name;
////	cmd>>bs_file_name;
////
////	int format (0);
////	cmd>>format;
////
////	string outfile_name;
////	cmd>>outfile_name;
////
////	BI::block_system bs;
////	ifstream infile;
////	if (format == 1)
////	{
////		infile.open (bs_file_name,ios::binary );
////		CGAL::set_binary_mode(infile);
////	}
////	else
////		infile.open (bs_file_name);
////
////	if (!infile)
////	{
////		cout<<"Open "<<bs_file_name<<" error"<<endl;
////		return 0;
////	}
////	else
////		cout<<"Open "<<bs_file_name<<" success"<<endl;
////
////	infile>>bs;
////	if(infile || infile.eof())
////		cout<<"read bs complete"<<endl;
////	else 
////	{
////		cout<<"read bs error"<<endl;
////		return 0;
////	}
////	infile.close();
////
////	using CGAL::to_double;
////	ofstream outfile(outfile_name);
////	outfile.precision(15);
////
////	for (int i(0);i!=bs.vcb.size ();++i)
////	{
////		cout<<i<<"/"<<bs.vcb.size ()<<endl;
////		BI::outer_bound ob (bs,i);
////		ob.output_BlockAni_format(outfile,i);
////	}
////	outfile.close();
////
////	return 0;
////}
//
////给郑文师同学的版本
////int main ()
////{
////	ifstream infile ("model_domain.dat");
////	if (!infile)
////	{
////		cout<<"Open domain file error"<<endl;
////		return 0;
////	}
////	else 
////		cout<<"Open domain file success"<<endl;
////
////	GB_domain d;
////	d.Input(infile);
////	infile.close();
////	BI::block_system bs;
////	bs.init_domain(d);
////	
////
////	infile.open ("contri_fracture_xyzabr.dat");
////
////	if (bs.disc_frac_parser(infile,d.Get_dx_axis()) ==0)
////		cout<<"Read fracture file success"<<endl;
////	else
////	{
////		cout<<"Read fracture file fail"<<endl;
////		return 0;
////	}
////	infile.close();
////
////	
////	cout<<"Identifying blocks"<<endl;
////	bs.identify_block (cout);
////
////	bs.init_stabl_alys(true);
////
////	cout<<"Identify block complete"<<endl;
////	cout<<"Outputing"<<endl;
////	
////	ofstream outfile ("vol_info.txt");
////	outfile.precision (18);
////	outfile<<"Number_of_elementblock"<<"\t"<<"Volume"<<"\t"<<"Type"<<endl;
////	BI::adj_blk_graph g(bs);
////	
////	if (g.output_block_type(bs,outfile) == 0)
////		cout<<"Output complete"<<endl;
////	else
////		cout<<"Output fail"<<endl;
////
////	outfile.close ();
////
////	return 0;
////}
//
