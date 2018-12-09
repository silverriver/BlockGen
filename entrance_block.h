#pragma once
#include "geometry.h"
#include "identification.h"
#include "outer_poly.h"
#include <vector>


namespace BI
{

class outer_bound;

struct angle_edge
{
	struct sub_angle
	{
		//ÿ���ӽǱ�ʾ��һϵ�а�ռ�Ľ������ӽǶ���͹��
		std::vector <vector_3> vnormal;		//���е�ÿ��������ʾһ����ռ�,ÿ������vnormal[i]ָ���ӽǵ��ⲿ
		//����false����ӽǴ��ڣ��ڵ�ļ��ϲ�Ϊ�գ�������true��ʾ����ӽǲ����ڣ��ڵ㼯Ϊ�գ�
		bool is_degenerate () const 
		{
			std::vector <plane_3> vpl;
			vpl.reserve (vnormal.size ());
			for (int i(0);i!=(int)vnormal.size ();++i)
				vpl.push_back (plane_3(vnormal[i].x(), vnormal[i].y(), vnormal[i].z(),0));
			switch (BI::intersection(vpl))
			{
			case Intersection_Status::Inter_Boundary:
				return true;
			case Intersection_Status::Inter_Unbounded:
				return false;
			default:
				throw std::logic_error ("bool angle_edge::sub_angle::is_degenerate () const, LP results error");
			}
		}
		//��������ӽǺ�v����ʾ�İ�ռ��Ƿ��ཻ������vָ���ռ��ⲿ���ⲿ����
		//�ж�����ӽǺ͸����İ�ռ��Ƿ��ڲ��ཻ(������ͬʱ)
		bool is_internal_intersected (const vector_3 &v) const;
		//�ж�����ӽǺ͸������ӽ�sa�Ƿ��ڲ��ཻ(������ͬʱ)
		bool is_internal_intersected (const sub_angle &sa) const;
		
		//�ж�����Ǻ͸����Ľ�sa�Ƿ��ཻ
		//����p_this��ʾthis����ǵĶ��㣬p_sa��ʾ������sa�Ķ���
		Intersection_Status intersection (const sub_angle&sa, const point_3 &p_this, const point_3 &p_sa) const;
		
	};

	struct angle
	{
		int e;						//����ǵĶ��㣬���е���������ʾ
		std::vector <sub_angle> vsa;	//ÿ���Ǳ���ʾΪ�����ӽǵĽ���

		//�ж�����Ǻ͸����Ľ�a�Ƿ��ڲ��ཻ(������ͬʱ)
		bool is_internal_intersected (const angle &a) const;
		//�ж�����Ǻ͸����İ�ռ��Ƿ��ڲ��ཻ,vָ���ռ��ⲿ
		bool is_internal_intersected (const vector_3 &v) const;

		//�ж�����Ǻ͸����Ľ�a�Ƿ��ཻ
		//����p_this��ʾthis����ǵĶ��㣬p_a��ʾ������a�Ķ���
		Intersection_Status intersection (const angle &a, const point_3&p_this, const point_3&p_a) const;
	};

	struct sub_edge
	{
		vector_3 n1, n2;		//��ʾһ��͹�ӱߣ���n1��n2��ʾ�İ�ռ�Ľ�����ȷ��������n1��n2ָ���ӱߵ��ⲿ��ע�⣬���ⲿ����
								//���n1==n2�����ʾ����ӱ���һ����ռ�
		//�ж������ӱ��Ƿ��ڲ��ཻ�����������ӱ��еİ�ռ䶼��ԭ�㣬Ȼ����Щ��ռ��Ƿ����ڲ�������
		//�ж�����ӱߺ͸������ӱ�se�Ƿ��ڲ��ཻ(������ͬʱ)
		bool is_internal_intersected (const sub_edge& se) const;

		//�ж�����ߺ͸����ı�se�Ƿ��ཻ
		//����p_this��ʾthis��������ϵ�ĳ���㣬p_se��ʾ������se���ϵ�ĳ����
		Intersection_Status intersection (const sub_edge&se, const point_3 &p_this, const point_3 &p_se) const;
	};

	struct edge
	{
		int e1,e2;		//��ʾ�����ߵ��⣬�㶼��������ʾ������Ҫ��֤e1<e2
		std::vector<sub_edge> vse;		//��ʾ�����ߵ��ӱߣ�ע��vse�еķ�������һ����ֱ������e1e2
		edge (const int&e1_, const int&e2_)
		{
			if (e1_<e2_)
			{e1 = e1_; e2 = e2_;}
			else
			{e1 = e2_; e2 = e1_;}
		}

		//�ж�����ߺ͸����ı�e�Ƿ��ڲ��ཻ(������ͬʱ)
		bool is_internal_intersected (const edge &e) const;

		//�ж��������Ƿ��ڲ��ཻ��
		//�ж�����ߺ͸����ı�e�Ƿ��ཻ
		//����p_this��ʾthis��������ϵ�ĳ���㣬p_e��ʾ������e���ϵ�ĳ����
		Intersection_Status intersection (const edge& e, const point_3& p_this, const point_3 &p_e) const;
	};

	std::vector<point_3> vp;
	std::vector<angle> va;
	std::vector<edge> ve;
};

//struct cover
//{
//	enum TYPE {Edge2Edge, Vertex2Poly};
//	int ia,ib;
//
//};

inline bool operator== (const angle_edge::edge &edge1, const angle_edge::edge &edge2)
{
	if (edge1.e1 == edge2.e1 && edge1.e2 == edge2.e2)
		return true;
	else
		return false;
}

class entrance_block
{
public:
	std::vector<poly_frac> vcover;		//ÿ����϶��һ�����Ƕ��� vector of cover
	outer_bound outer_b;						//���������������Ϊ�������
	block_system bs;
	cblock inner;		//��������ڲ��ĵ㱻��ʾΪһ����Ԫ��������
private:
	int gen_cover (const block_system &bs, const int &cbi, const int &cbj);		//��ȡ���Ͽ������״����ʼ������vcov
	int gen_cover (const outer_bound &oba, const outer_bound &obb, const point_3 &a0);				//��outer_bound��ʾA��B,Ҫ��A��B�����Ѿ�regulated֮���
	int gen_cover (const outer_bound &ob_a, const angle_edge &ae_a, const outer_bound &ob_b, const angle_edge &ae_b, const point_3 &a0 = point_3 (0,0,0));
public:
	int gen_entrance_block (const block_system &bs_, const cblock& cba, const cblock& cbb, const point_3 &a0, bool debug_info = false);		//�ӿ���ƵĲ��Ǻܺã��д��Ľ�
	int gen_entrance_block (const block_system &bs_, const cblock& cba, const cblock& cbb, const point_3 &a0,
							const outer_bound &ob_a, const angle_edge &ae_a, const outer_bound &ob_b, const angle_edge &ae_b, bool debug_info = false);
	
	friend std::ostream& operator<< (std::ostream& os, const entrance_block& eb);
	friend std::istream& operator>> (std::istream& is, entrance_block& eb);
	friend bool operator== (const entrance_block& eb1, const entrance_block& eb2);
};

int outer_bound2angle_edge (const block_system & bs, const cblock & cb,const outer_bound &ob, angle_edge &ae);

inline bool operator== (const entrance_block& eb1, const entrance_block& eb2)
{
	return (eb1.vcover == eb2.vcover &&
		eb1.outer_b == eb2.outer_b &&
		eb1.bs == eb2.bs &&
		eb1.inner == eb2.inner );
}

//�ж�һ�����һ����Ԫ�����λ�ù�ϵ,vfi��һ������ֵ����ʾ���������ڸõ�Ԫ�����ϣ����¼�����λ���ļ���ƽ����
Oriented_side point_eblock_location (const block_system &bs, const eblock &eb, const point_3&tar, std::vector<int>& vfi);


}