#pragma once
#include "geometry.h"
#include "identification.h"
#include "outer_poly.h"
#include "excav_surface.h"
#include <vector>

namespace BI
{

enum Removability_Status {Removable, Nonremovable};
//�ж�һ�������Ŀ����Ƿ��ǿ��ƶ��ģ�����v_dir�еĸ���������ʾ������������Ӵ��Ĺ̶��棬����ָ������ڲ���
//ע�⣬ֻ��һ���ƶ���������Ҳ���ж�Ϊ���ƶ�����Ҳ����v_dir����������������ƽ�еģ�
//ע�⣬��intersection����������������鲢��ǡ������Ϊ��Щ��ռ������ཻ��0�㡣
//����������������ж��Ƿ�������Ч����΢��һ�㡣
Removability_Status removability (const std::vector<vector_3> &v_dir);		

class outer_bound;
struct safety_factor	//������¼������ȶ�����Ϣ
{
	enum S_TYPE {FIXED, SIG_F_SLID, DOU_F_SLID, NON_REMOVABLE, UNEXPOSED, REMOVABLE_BUT_SAFE, FALLING};
	double sf;			//��ȫϵ��
	INT facej, facei;	//���������϶���
	S_TYPE type;			//��������
	double normal_forcej, normal_forcei;
	double sliding_force;
	double resisting_forcei, resisting_forcej;
	vector_3 moving_dir;	//���ܵ��ƶ�����	������ֵ����Ϊ���������õ���ƽ��������
};

//����һ������ĺ�����vo��ʾ��������淨��������r��ʾ�������ܵĺ������ķ���
//vo��ʾ�������������淨�����е�����ָ������ڡ�
//Ҫ��vo�е�������ͬ���������������ȫ��ͬ����������ܻ��д���Ľ��
//����ֵsf��facei, facej��ʾ����淨����vo�еĵı�š�
//safety_factor�е�sf��resisting_forcei, resisting_forcej��û�м���
//���ص��ƶ�����ֻ����SIG_F_SLID, DOU_F_SLID, NON_REMOVABLE, REMOVABLE_BUT_SAFE, FALLING
//��������Ĺ��ܺ�removability�е����ƣ�Ӧ��Ҳ��ǰ�߶ࡣ������Ч�ʵ�����¿������������
//tol���ڷ�����ʱ�������С��
int cal_safety_factor (safety_factor & sf, const std::vector<vector_3> &vo, const vector_3 &r, const double & tol= 1E-10, const bool &debug_flag = false);



}