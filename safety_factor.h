#pragma once
#include "geometry.h"
#include "identification.h"
#include "outer_poly.h"
#include "excav_surface.h"
#include <vector>

namespace BI
{

enum Removability_Status {Removable, Nonremovable};
//判断一个给定的块体是否是可移动的，其中v_dir中的各个向量表示与这个这个块体接触的固定面，方向指向块体内部。
//注意，只有一个移动方向的情况也被判定为可移动。（也就是v_dir中两个向量是允许平行的）
//注意，用intersection这个函数做这样事情并不恰当。因为这些半空间总是相交于0点。
//这个函数单纯用来判断是否被锁定。效率稍微高一点。
Removability_Status removability (const std::vector<vector_3> &v_dir);		

class outer_bound;
struct safety_factor	//用来记录块体的稳定性信息
{
	enum S_TYPE {FIXED, SIG_F_SLID, DOU_F_SLID, NON_REMOVABLE, UNEXPOSED, REMOVABLE_BUT_SAFE, FALLING};
	double sf;			//安全系数
	INT facej, facei;	//滑动面的裂隙编号
	S_TYPE type;			//滑动类型
	double normal_forcej, normal_forcei;
	double sliding_force;
	double resisting_forcei, resisting_forcej;
	vector_3 moving_dir;	//可能的移动方向	。近似值，因为求解过程中用到了平方根运算
};

//重载一下上面的函数，vo表示所有外表面法向向量，r表示块体所受的合外力的方向。
//vo表示块体的所有外表面法向，所有的向量指向块体内。
//要求vo中的向量不同向。如果存在两个完全相同的向量则可能会有错误的结果
//返回值sf中facei, facej表示外表面法线在vo中的的编号。
//safety_factor中的sf，resisting_forcei, resisting_forcej域没有计算
//返回的移动类型只能是SIG_F_SLID, DOU_F_SLID, NON_REMOVABLE, REMOVABLE_BUT_SAFE, FALLING
//这个函数的功能和removability有点相似，应用也比前者多。不考虑效率的情况下可以用这个函数
//tol是在分析的时候容许的小误差。
int cal_safety_factor (safety_factor & sf, const std::vector<vector_3> &vo, const vector_3 &r, const double & tol= 1E-10, const bool &debug_flag = false);



}