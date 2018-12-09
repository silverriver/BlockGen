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

	//���򣬿�������
	infile>>dx_axis>>exc_type;

	//���������ܽڵ���
	int nf,np;
	infile>>nf>>np;

	//ÿ����Ľڵ���
	vector<int> f_points;
	for (int i(0);i!=nf;++i)
	{
		if (!infile)
			return 1;
		infile>>itemp;
		f_points.push_back (itemp);
	}

	//ÿ�����Ψһ��ţ���ʵ���Ժ���
	for (int i(0);i!=nf;++i)
		infile>>itemp;	

	//ÿ�������ѧ����
	for (int i(0);i!=nf;++i)
	{
		if (!infile)
			return 1;
		infile>>itemp;
		face_type.push_back(itemp);
	}

	//ÿ����Ķ������
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

	//ÿ���ڵ����ά����
	double x,y,z;
	for (int i(0);i!=np;++i)
	{
		if (!infile)
			return 1;
		infile>>x>>y>>z;
		p.push_back(BI::point_3(x,y,z));
	}

	//��������
	int nsubdomain;
	infile>>nsubdomain;

	//ÿ�������������ܶ�
	for (int i(0);i!=nsubdomain;++i)
	{
		if (!infile)
			return 1;
		infile>>dtemp;
		density.push_back(dtemp);
	}

	//ÿ������������
	vector<int> fr_nface;
	for (int i(0);i!=nsubdomain;++i)
	{
		if (!infile)
			return 1;
		infile>>itemp;
		fr_nface.push_back(itemp);
	}

	//ÿ���������������
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
