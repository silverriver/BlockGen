#pragma once
#include "geometry.h"
#include "GBUtility.h"
#include "identification.h"
#include <vector>

namespace BI
{
//那么要开一个新坑了，这个类表示一个开挖面的类。意思是所有的开挖面都记录在这里
//其实还没有想好到底应该怎么利用这个类，初步想法是判断可移动性之前先用outer_bound和这个类比较一下，然后就可以找出固定面了（虽然目前想起来这个方法的效率实在是太低。但是先试试吧）
class excav_surface
{
public:
	std::vector <polygon_set_2> vpoly;	//表示每个平面上包含的开挖面边界
	std::vector <plane_3> vplane;		//这个vector中包含的平面都是不重合的
	//这个函数的实现先等一等，因为实现过程中还要涉及到处理输入数据中的误差
	excav_surface (const GB_domain &d);

	//将这个bs中所有单元块体组成的复合块体的外表面当成是开挖面。比较简单粗暴的方法
	excav_surface (const block_system &bs);
};

}