#include "safety_factor.h"
#include <iostream>
#include <CGAL\QP_functions.h>

using namespace std;

namespace BI
{

Removability_Status removability (const std::vector<vector_3> &v_dir)
{
	//typedef CGAL::Quadratic_program<K::FT> Program;
	//typedef CGAL::Quadratic_program_solution<K::FT> Solution;
	//Program lp (CGAL::LARGER, true, -1, true, 1);
	//for (int i(0);i!=(int) v_dir.size ();++i)
	//{
	//	lp.set_a (0, i, v_dir[i].x());
	//	lp.set_a (1, i, v_dir[i].y());
	//	lp.set_a (2, i, v_dir[i].z());
	//	lp.set_b (i, 0);
	//}
	//lp.set_c(0,-1);
	//Solution s = CGAL::solve_quadratic_program(lp,K::FT());
	////debug
	//cout<<"fun val:"<<s.objective_value()<<endl;
	//cout<<"x:"<<(*s.variable_values_begin())<<endl;
	//cout<<"y:"<<(*(s.variable_values_begin()+1))<<endl;
	//cout<<"z:"<<(*(s.variable_values_begin()+2))<<endl;
	//if (s.is_infeasible())
	//	throw logic_error ("removability, LP infeasible");
	//if (s.is_unbounded())
	//	throw logic_error ("removability, LP unbounded");
	//if (s.objective_value() == 0)
	//	return Removability_Status::Nonremovable;		//��Ӧ�Ĳ���ʽֻ�����
	//if (s.objective_value() >0)
	//	return Removability_Status::Removable;			//��Ӧ�Ĳ���ʽ�з����
	//else 
	//	throw logic_error ("removability, LP has negative solution");		//�����ϲ������������

	if (v_dir.size () <=3)
		return Removability_Status::Removable;

	//���ڿ��ƶ����� <=> ĳ�����ǿ��ƶ�����
	for (int i(0);i!= (int) v_dir.size ();++i)
		for (int j(i+1);j < (int) v_dir.size ();++j)
		{
			vector_3 temp1 (CGAL::cross_product(v_dir[i], v_dir[j]));
			vector_3 temp2 (-temp1);

			if (temp1*temp1 == 0)	//������������
				continue;

			int count1 (0);
			for (;count1!=(int) v_dir.size ();++count1)
				if ((v_dir[count1]*temp1) <0)
					break;
			if (count1 == (int) v_dir.size ())
				return Removability_Status::Removable;		//��������һ�����ƶ�����

			int count2(0);
			for (;count2!=(int) v_dir.size ();++count2)
				if ((v_dir[count2]*temp2) <0)
					break;
			if (count2 == (int) v_dir.size ())
				return Removability_Status::Removable;
		}
	return Removability_Status::Nonremovable;
}

//vo��ʾ�������������淨�����е�����ָ������⣬���ǿ����ڡ�
int cal_safety_factor (safety_factor & sf, const std::vector<vector_3> &vo, const vector_3 &r, const double & tol, const bool &debug_flag)
{
	if (debug_flag)
		cout<<"cal_safety_factor()"<<endl;
	if (vo.empty ())
	{
		sf.facei = sf.facej = -1;		//���û���޶�����ô�϶�����r�ķ����ƶ���
		sf.moving_dir = r;
		sf.normal_forcei = sf.normal_forcej = -1;
		sf.resisting_forcei = sf.resisting_forcej = -1;
		sf.sf = -1;
		sf.type = safety_factor::FALLING;
		return 0;
	}
	//ͨ���������п��ܵ��ƶ������жϿ���Ŀ��ƶ���
	vector<vector_3> vpos;		//���п��ܵ��ƶ�����
	vpos.reserve (vo.size ()+vo.size ()*(vo.size ()+1)/2+1);
	vector<int> facei, facej;	//������¼vpos��ÿ����������Դ
	facei.reserve (vpos.size ()); facej.reserve (vpos.size());	

	vpos.push_back (r);
	facei.push_back (-1), facej.push_back (-1);

	for (int i(0);i!=vo.size ();++i)
	{
		vector_3 v_temp (CGAL::cross_product(CGAL::cross_product(vo[i],r), vo[i]));
		if (v_temp == -v_temp)
			continue;		//�����������
		vpos.push_back (v_temp);		//���������һ��ƽ�滬���ķ���
		vpos.push_back (-v_temp);		//���������п��ܳ�Ϊ���ƶ�����
		facei.push_back (i); facei.push_back (i); facej.push_back (-1); facej.push_back (-1);
	}
	for (int i(0);i!=vo.size ();++i)
	{
		for (int j(i+1); j<vo.size ();++j)
		{
			vector_3 v_temp (CGAL::cross_product(vo[i],vo[j]));
			if (v_temp == -v_temp)
				continue;		//�����������
			vpos.push_back (v_temp);		//�������������ƽ�滬���ķ���
			vpos.push_back (-v_temp);
			facei.push_back (i); facei.push_back (i); facej.push_back (j);facej.push_back (j);
		}
	}

	if (debug_flag)
		cout<<"vpos.size():"<<vpos.size ()<<endl;

	vector<int> vm_index;		//��vpos��ɸѡ���ƶ�����,ֻ��¼���ƶ�������vpos�е�������
	vm_index.reserve (vpos.size ());
	for (int i(0);i!=vpos.size ();++i)
	{
		//if (vpos[i].squared_length() <= tol)		
		//	continue;			//������������Ϊ���ƶ�����
		int count(0);
		for (;count!=vo.size ();++count)
			if (vo[count]*vpos[i]<0)
				break;
		if (count == vo.size ())
			vm_index.push_back (i);
	}
	
	if (debug_flag)
		cout<<"vm_index.size ()"<<vm_index.size ()<<endl;

	if (vm_index.empty ())
	{sf.type = safety_factor::S_TYPE::NON_REMOVABLE; return 0;}	//���û�п��ƶ��������ǲ����ƶ�����
	double max_dotres;		//�п��ƶ�������һ����Ҫ�ҳ����п��ܵ��ƶ�����
	int vec_index(0);		//��ʾ���п����ƶ�������vm_index�е�����
	max_dotres = CGAL::to_double(vpos[vm_index[vec_index]]*r/sqrt(CGAL::to_double(vpos[vm_index[vec_index]].squared_length())));
	for (int i(0);i!=vm_index.size ();++i)
	{
		double temp_res = CGAL::to_double(vpos[vm_index[i]]*r/sqrt(CGAL::to_double(vpos[vm_index[i]].squared_length())));
		if (temp_res > max_dotres)
		{ max_dotres=temp_res; vec_index = i;}
	}
	if (max_dotres <= 0)
	{ sf.type = safety_factor::S_TYPE::REMOVABLE_BUT_SAFE; return 0;}	//���û���κ�һ�������ƶ�������r��ͬ����Ȼ�ÿ����ǿ��ƶ��ģ����ǲ�����r��Ӱ�����ƶ�
	
	if (facei[vm_index[vec_index]] == -1)
	{
		sf.type = safety_factor::S_TYPE::FALLING;		//��r�ķ����ƶ�
		sf.facei = sf.facej = -1;
		sf.moving_dir = r;
		sf.sliding_force = sf.normal_forcei = sf.normal_forcej = -1;		//�ǻ������ϵ���Ϊ0
	}
	else if (facej[vm_index[vec_index]] == -1)
	{
		sf.type = safety_factor::S_TYPE::SIG_F_SLID;		//���滬��
		sf.facei = facei[vm_index[vec_index]];
		sf.facej = -1;
		sf.moving_dir = vpos[vm_index[vec_index]];
		sf.sliding_force = CGAL::to_double(r*sf.moving_dir/sqrt(CGAL::to_double (sf.moving_dir.squared_length())));
		sf.normal_forcei = CGAL::to_double (r*vo[sf.facei]/sqrt(CGAL::to_double (vo[sf.facei].squared_length())));
		sf.normal_forcej = 0;		//�ǻ������ϵ���Ϊ0
	}
	else
	{
		sf.type = safety_factor::S_TYPE::DOU_F_SLID;		//˫�滬��
		sf.facei = facei[vm_index[vec_index]];
		sf.facej = facej[vm_index[vec_index]];
		sf.moving_dir = vpos[vm_index[vec_index]];
		sf.sliding_force = CGAL::to_double(r*sf.moving_dir/sqrt(CGAL::to_double (sf.moving_dir.squared_length())));
		vector_3 ni = vo[sf.facei], nj = vo[sf.facej];
		ni = (1.0/sqrt(CGAL::to_double(ni.squared_length())))*ni;
		nj = (1.0/sqrt(CGAL::to_double(nj.squared_length())))*nj;
		sf.normal_forcei = CGAL::to_double(CGAL::cross_product(r,nj)*CGAL::cross_product(ni,nj)/CGAL::cross_product(ni,nj).squared_length ());
		sf.normal_forcej = CGAL::to_double(CGAL::cross_product(r,ni)*CGAL::cross_product(ni,nj)/CGAL::cross_product(ni,nj).squared_length ());
	}
	return 0;
	//�������2015-10-24.
	//��δ���ԣ��κ�bug���п��ܳ��֡�
}

};