#pragma once
#include "geometry.h"
#include "GBUtility.h"
#include "identification.h"
#include <vector>

namespace BI
{
//��ôҪ��һ���¿��ˣ�������ʾһ����������ࡣ��˼�����еĿ����涼��¼������
//��ʵ��û����õ���Ӧ����ô��������࣬�����뷨���жϿ��ƶ���֮ǰ����outer_bound�������Ƚ�һ�£�Ȼ��Ϳ����ҳ��̶����ˣ���ȻĿǰ���������������Ч��ʵ����̫�͡����������԰ɣ�
class excav_surface
{
public:
	std::vector <polygon_set_2> vpoly;	//��ʾÿ��ƽ���ϰ����Ŀ�����߽�
	std::vector <plane_3> vplane;		//���vector�а�����ƽ�涼�ǲ��غϵ�
	//���������ʵ���ȵ�һ�ȣ���Ϊʵ�ֹ����л�Ҫ�漰���������������е����
	excav_surface (const GB_domain &d);

	//�����bs�����е�Ԫ������ɵĸ��Ͽ��������浱���ǿ����档�Ƚϼ򵥴ֱ��ķ���
	excav_surface (const block_system &bs);
};

}