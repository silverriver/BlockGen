#pragma once
#include "geometry.h"
#include "identification.h"
#include <vector>
#include <ostream>
#include <string>

namespace BI
{
struct angle_edge;
struct safety_factor;
class entrance_block;

class outer_bound
{
public:
	friend int cal_safety_factor (const outer_bound &ob, safety_factor & sf);
	friend int outer_bound2angle_edge (const block_system & bs, const cblock & cb,const outer_bound &ob, angle_edge &ae);
	friend class entrance_block;
//private:  DEBUG
	std::vector<polygon_with_holes_2> vpoly;
	std::vector<plane_3> vplane;		//Ҫ�����vplane[i]. ��������ƽ�涼������vplane[i]��vplane[i].opposite()
	std::vector<bool> outer_normal;	//true���������ε��ⷨ�߷�����vplane[i]ƽ���������ͬ
public:
	outer_bound () {}
	outer_bound (const block_system &bs,const int &cbi,const bool&debug_info=false)
	{construct_outer_bound(bs,cbi,debug_info);}
	outer_bound (const block_system &bs, const std::vector<int>& vebi,const bool &debug_info = false)
	{construct_outer_bound(bs,vebi,debug_info);}

	
	//���ԭ����Ϣ�������߽硣����0��ʾ�����ɹ�����ʵ��������������׳��쳣���Ƿ���0
	//���뺯����vebiָ������Ҫ���������߽�ĵ�Ԫ�������������
	int construct_outer_bound (const block_system &bs, const std::vector <int> &vebi, const bool& debug_info = false);
	//�����������ʵ�ֵ�ʱ����ֱ�ӵ����������������
	int construct_outer_bound (const block_system &bs, const int &cbi, const bool& debug_info = false);
	//�������ɸ�cb�������
	int construct_outer_bound_using_vcb (const block_system &bs, const std::vector <int> &vcbi, const bool& debug_info = false);

	//���stl�ļ�������0��ʾ����ɹ���������ļ�ֻ�ǽ��ƶ���
	int output_STL (std::ostream &outfile, const FT& tranx = 0.0, const FT& trany = 0.0, const FT& tranz = 0.0, const std::string& ob_name="") const;
	int output_exp_STL (std::ostream &outfile, const double &ratio = 0.2, const FT& tranx = 0.0, const FT& trany = 0.0, const FT& tranz = 0.0, const std::string& ob_name="") const;		//ob_nameû���õ�
	int output_CAD_script (std::ostream &outfile,const bool& split_ploy = false, const bool &output_normal=false) const;

	//�����ɵ������ṹ�����BlockAni_1.3.3����ʾ�Ľṹ
	int output_BlockAni_format (std::ostream &outfile, const int &id) const;

	//��ԭ�������߽�Ļ��������һЩ��Ԫ���塣�����㷨��join���������Чһ�㣬��ʡ��һЩ�ڴ渴�Ʋ���
	int join_veb (const block_system &bs, const std::vector <int> &vebi, const bool& debug_info = false);
	//���������Ŀ�������Ѿ������õ�bs�Ļ������������һ�����Ͽ��塣ʵ�ַ�����ֱ�ӵ��ú���joint_veb
	int join_cb (const block_system &bs, const int cbi, const bool& debug_info = false);

	//������������С��Χ�е����ĵ�
	point_3 get_center () const;

	int get_polycount () const
	{return (int)vpoly.size ();}

	//�����������߽���ñ߽�������
	int join (const outer_bound &ob, const bool &debug_info = false);

private:
	//˽�к�����������vpoly[poly_count]�������μ�����������εĽǵ�
	//vpoint��ʾvpoly�и�������εĽǵ�
	//p_count��ʾvpoly�и�������ν��������
	//vbox��ʾvpoly�и�������ε���С��Χ��
	//������Щ��Ϣ�������õ�ʱ�����ɣ�����Ϊ��ʡʱ�䣬��Ҫ�ڵ������������ʱ���ֶ�cache��Щ��Ϣ
	//����ֵ������������������˶��ٵ�
	int add_coner_point (const int &poly_n,  point_3 **vpoint,  int *p_count,  mbox_3 *vbox);

	//������������������������ӵ㣬����original�еĵ���ٴ�ת������ά�ռ��У���Ȼ��һ�������Ż��������ǻ��漰�����ӵ�ָ��ӿڣ�Ϊ��ͼʡ�£���ʱ����2015-10-25 19:37
	//ͬ��Ϊ��ʡ�£���ʾǱ����ӵ�ʱ������Ľӿ�ֱ���õ��ǵ�����꣬�����ǵ��������
	//vconer��ʾ��Ҫ��ӵĵ㣬��Щ��Ӧ�ö���plane���ƽ���ϣ�����vconer�еĵ�Ӧ���໥����ͬ
	//res��ʾ�����ɵĶ���Ρ�����ֵ��ʾ��original�Ļ���������ӵ������
	int add_point_on_poly (const polygon_2 &original, const plane_3& plane, const std::vector<point_3> & vconer, polygon_2 &res);

	//������������������ɵ������߽磬Ҳ���ǽ������ߵĽǵ���������������������
	//�ڵ����������ʱ����ȼ�һ�������߽�����simplify_polyon_2
	int integrate_boundary ();

	//flag ��ʾ��ʲô�����ʽ��false��ʾ����������������ɵ���solid��true��ʾ��������������ɵ���solid
	int output_STL_poly (std::ostream &outfile, const int &pi, const bool& flag = false, const FT& tranx = 0.0, const FT& trany = 0.0, const FT& tranz = 0.0) const;

	//��ĳ������ε����ĵ㣬Ҳ���Ǹ��������ƽ��ֵ
	point_3 get_poly_center (const int &pi) const;

public:
	friend std::ostream& operator<< (std::ostream &os, const outer_bound& ob);
	friend std::istream& operator>> (std::istream &is, outer_bound &ob);
	friend bool operator== (const outer_bound&ob1, const outer_bound& ob2);
};

//������tar����Ҫ�����ֵ(����ֵ��Լ0)��index��Ϊ����ֵ����������СԪ������һֱ�����Ԫ������
void sort_index (const std::vector<FT> &tar, std::vector<int> &index);

inline bool operator== (const outer_bound&ob1, const outer_bound& ob2)
{
	return (ob1.outer_normal == ob2.outer_normal &&
		ob1.vplane == ob2.vplane &&
		ob1.vpoly == ob2.vpoly);
}

}



