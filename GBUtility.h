#pragma once
#include <vector>
#include "geometry.h"
namespace BI
{
	class block_system;
}

class GB_domain
{
public:
	friend class BI::block_system;
	struct fram			// 记录domain中fram的面编号，编号从0开始
	{
		std::vector<int> nf;
	};
	struct face
	{
		std::vector<int> pn;			//记录face的顶点编号，编号对应于block中的顶点。编号从0开始
	};
	GB_domain () :dx_axis (0.0){}
	int GetPointCount () const
	{return (int) p.size();}
	int GetFaceCount ()const 
	{return (int) fa.size ();}
	int GetFramCount ()const
	{return (int) fr.size ();}
	int GetFramPointCount (const int& fram_index)const;

	int GetFacePointCount (const int& face_index)const
	{return (int) fa[face_index].pn .size();}
	int GetFramFaceCount (const int&fram_index)const
	{return (int) fr[fram_index].nf.size ();}
	int GetFramFacePointCount (const int&fram_index,const int& face_index)const
	{return (int) fa[fr[fram_index].nf [face_index]].pn.size();}
	const BI::point_3 &GetPoint (const int&index)const
	{return p[index];}
	const BI::point_3 &GetFacePoint (const int &face_index,const int &point_index)const
	{return p[fa[face_index].pn [point_index]];}
	const BI::point_3 &GetFramFacePoint (const int &fram_index,const int&face_index,const int &point_index)const
	{return p[fa[fr[fram_index].nf [face_index]].pn [point_index]];}
	double Get_dx_axis ()const
	{return dx_axis;}
	int GetFaceType (const int &i)const
	{return face_type[i];}
	int GetFramFaceType (const int &fram_index,const int &face_index) const
	{
		return face_type[fr[fram_index].nf [face_index]];
	}
	double GetDensity (const int &fram_index)const
	{return density[fram_index];}
	int GetExcType () const
	{return exc_type;}
	const fram & GetFram (const int &fram_index) const {return fr[fram_index];}
	const face & GetFramFace (const int &fram_index, const int &face_index)const {return fa[fr[fram_index].nf[face_index]];}

	void AddPoint (const BI::point_3 &p_)
	{p.push_back (p_);}
	void AddFace (const face &f_)
	{fa.push_back (f_);}
	void AddFram (const fram &f_)
	{fr.push_back (f_);}
	void AddPointToFace (const int& point_index,const int&face_index)
	{fa[face_index].pn .push_back (point_index);}
	void AddEmptyFace ()
	{fa.push_back(face());}
	void AddEmptyFram ()
	{fr.push_back (fram());}
	void AddFaceToFram (const int&fram_index,const int&face_index)
	{fr[fram_index].nf.push_back (face_index);}
	void AddFaceType (const int &type)
	{face_type.push_back (type);}
	void AddDensity (const double &d_)
	{density.push_back (d_);}

	void Set_dx_axis(const double& d)
	{dx_axis =d;}
	void SetFaceType (const int&face_index,const int &type)
	{face_type[face_index]=type;}
	void SetDensity (const int &fram_index,const double & d_)
	{density[fram_index]=d_;}
	void SetExcType (const int &et)
	{exc_type = et;}

	void ClearAllItem ()
	{p.clear ();fa.clear ();fr.clear ();face_type.clear ();density .clear ();dx_axis = 0.0;}

	//返回0表明成功，其他值表明失败
	//不管核对数据文件的正确性
	int Output (std::string filename) const ;

	//返回0表明成功，其他值表明失败
	//读取的是model_domain.dat
	int Input (std::ifstream & infile);

	std::vector<BI::point_3> p;
	std::vector<face> fa;
	std::vector<fram> fr;
	std::vector<int> face_type;		//每个面的力学性质0=公用，1=固定，2=自由
	std::vector<double> density;			//每个子区的密度
	double dx_axis;
	int exc_type;				//开挖类型，边坡0，隧洞1，硐室2，复合3
};
