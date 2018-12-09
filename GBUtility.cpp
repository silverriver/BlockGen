#include "GBUtility.h"
#include <string>
#include <fstream>
using namespace std;

int GB_domain::GetFramPointCount (const int& fram_index)const
{
	vector<int> temp;
	if(fr[fram_index].nf.empty ())
		return 0;
	for (int i(0);i!=fr[fram_index].nf.size ();++i)
		for (int j(0);j!=fa[fr[fram_index].nf [i]].pn .size ();++j)
		{
			int count(0);
			for (; count!=temp.size ();++count)
				if (temp[count] == fa[fr[fram_index].nf [i]].pn[j])
					break;
			if (count == temp.size ())
				temp.push_back (fa[fr[fram_index].nf [i]].pn[j]);
			else 
				continue;
		}
		return (int) temp.size ();
}
int GB_domain::Input(std::ifstream & infile)
{
	if (!infile)
		return 1;

	int itemp,itemp2;
	double dtemp;

	//倾向，开挖类型
	infile>>dx_axis>>exc_type;

	//总面数，总节点数
	int nf,np;
	infile>>nf>>np;

	//每个面的节点数
	vector<int> f_points;
	for (int i(0);i!=nf;++i)
	{
		if (!infile)
			return 1;
		infile>>itemp;
		f_points.push_back (itemp);
	}

	//每个面的唯一编号，其实可以忽略
	for (int i(0);i!=nf;++i)
		infile>>itemp;	

	//每个面的力学性质
	for (int i(0);i!=nf;++i)
	{
		if (!infile)
			return 1;
		infile>>itemp;
		face_type.push_back(itemp);
	}

	//每个面的顶点序号
	for (int i(0);i!=nf;++i)
	{
		fa.push_back(face());
		for (int j(0);j!=f_points[i];++j)
		{
			if (!infile)
				return 1;
			infile>>itemp2;
			fa[fa.size ()-1].pn.push_back(itemp2-1);
		}
	}

	//每个节点的三维坐标
	double x,y,z;
	for (int i(0);i!=np;++i)
	{
		if (!infile)
			return 1;
		infile>>x>>y>>z;
		p.push_back(BI::point_3(x,y,z));
	}

	//子区个数
	int nsubdomain;
	infile>>nsubdomain;

	//每个子区的岩体密度
	for (int i(0);i!=nsubdomain;++i)
	{
		if (!infile)
			return 1;
		infile>>dtemp;
		density.push_back(dtemp);
	}

	//每个子区的面数
	vector<int> fr_nface;
	for (int i(0);i!=nsubdomain;++i)
	{
		if (!infile)
			return 1;
		infile>>itemp;
		fr_nface.push_back(itemp);
	}

	//每个子区的面号序列
	for (int i(0);i!=nsubdomain;++i)
	{
		fr.push_back(fram());
		for (int j(0);j!=fr_nface[i];++j)
		{
			if (!infile)
				return 1;
			infile>>itemp;
			fr[fr.size ()-1].nf.push_back (itemp-1);
		}
	}
	infile.close();
	return 0;
}
