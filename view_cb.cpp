//view_cb.exe
//输出指定cblock的爆炸视图,也就是把指定的cblock中的单元块体以爆炸视图展示出来
//view_cb_input.txt格式：
//bs文件名
//bs文件格式 （1表示是binary，2表示是ascii）
//输出文件的前缀
//要输出几个cblock
//对于每个cblock： cb的索引   cb的爆炸系数

#include "identification.h"
#include "geometry.h"
#include "GBUtility.h"
#include "outer_poly.h"
#include "curve_frac.h"
#include "entrance_block.h"
#include "serialize.h"
#include "safety_factor.h"
#include "inner_poly.h"
#include "adj_blk_graph.h"
#include "fixed_surfaces.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <utility>
#include <algorithm>
#include <CGAL\Timer.h>
#include <CGAL\Memory_sizer.h>
#include <cstdlib>
#include <map>

using namespace std;
using BI::operator<<;
using BI::operator>>;

int main ()
{
	ifstream cmd ("view_cb_input.txt");
	if(!cmd)
	{
		cout<<"can't open view_cb_input.txt"<<endl;
		return 0;
	}
	string filename;		//bs文件名
	cmd>>filename;

	int format (0);			//bs文件的格式，1表示是binary，2表示是ascii
	cmd>>format;

	string out_perfix;		//输出文件的前缀
	cmd>>out_perfix;

	int ncb(0);
	cmd>>ncb;				//要输出几个cblock

	vector<int> cbi;
	vector<double> vratio;
	cbi.reserve (ncb);
	vratio.reserve (ncb);

	for (int i(0);i!=ncb;++i)
	{
		if (!cmd)
		{
			cout<<"Input file has wrong format"<<endl;
			return 0;
		}
		cmd>>(cbi[i]) >> (vratio[i]);
	}
	cmd.close();

	ifstream in_bs;
	if (format == 1)
	{
		in_bs.open (filename,ios::binary );
		CGAL::set_binary_mode(in_bs);
	}
	else
	{
		in_bs.open (filename);
	}

	if (!in_bs)
	{
		cout<<"Open "<<filename<<" error"<<endl;
		return 0;
	}
	else 
		cout<<"Open "<<filename<<" success"<<endl;

	BI::block_system bs;
	in_bs>>bs;
	in_bs.close();

	for (int i(0);i!=ncb;++i)
	{
		cout<<"cb:"<<cbi[i];
		if (cbi[i] <0 || cbi[i]>= (int)bs.vcb.size ())
		{
			cout<<" out of range"<<endl;
			continue;
		}
		else
			cout<<" outputing"<<endl;
		if (vratio[i]<0)
		{
			cout<<"given ratio <0"<<endl;
			continue;
		}
		ostringstream ost;
		ost.str("");
		ost<<out_perfix<<"_cb"<<cbi[i]<<".stl";
		ofstream out_cb(ost.str());
		bs.output_explode_cb_stl(cbi[i],out_cb,vratio[i]);
		out_cb.close();
	}
	cout<<"Output cblock complete"<<endl;
	return 0;
}