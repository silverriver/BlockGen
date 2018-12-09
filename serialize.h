#pragma once
#include "geometry.h"
#include "identification.h"
#include "outer_poly.h"
#include "inner_poly.h"
#include <istream>
#include <ostream>
#include <vector>
#include <CGAL\Gmpq.h>
#include <CGAL\Vector_3.h>
#include <CGAL\gmp.h>
//����ļ�����һЩ�����л�����ĺ���
//�ѵ�����һ���������л�FT�����ڵ���Ҫ�����������ȼ���FT��Ӧ���ϸ�����
//��ΪĿǰ�׶�FT���õ������������lazy computation������ֱ����serialize��Щ���ѡ�
//��ȡʱ��Gmpq��Ϊ�м���������������ж�ȡGmpq��Ȼ��������ϸ����ͳ�ʼ��FT

//ʹ����֪��
//���л�T��
//���T���������ͣ����ֶ�����ostream::write��istream::in��ʵ�ֶ������������� ��ע�⣬��T��bool��������Ҫ�ֶ�ת��һ�£�һ����˵��true��ӦΪ'1'��false��Ӧ��'0'���ڶ�ȡʱ������任��
//���T�����������ͣ���ֱ����<<��>>��������ʵ��

//���л�vector<T>:
//���T���������ͣ������binary_in��binary_out��ʵ�ֶ������������� ����ʱ����T�Ƿ���bool������ֱ�ӵ�������������������<<��>>������ʵ��ASCII����������
//���T�����������ͣ���ֱ����<<��>>��������ʵ��

//���л�vector<vector<T>>
//���T���������ͣ������binary_in��binary_out��ʵ�ֶ�����������������<<��>>����������ASCII��ʽ���������� ��û��ʵ��bool���͵ĺ�������Ҫ��ʱ����д��
//���T�����������ͣ���ֱ����<<��>>��������ʵ��
//����ʱ��Ƶ�ʱ�����ӳ�磬��Ҫ�ò�������ʵ�ֶ����Ƶ������������ͼװ��һ�ѡ�û�����Ǻ����������е���ҡ�ʧ�ߡ������пմ���һ�¡�2016-6-9 17:42��
namespace BI
{
//���T���������ͣ����ò�����<<��>>ֻ�����ascii��ʽ���ļ����޷����ز�����ʹ�����binary��ʽ���ļ�����ʱ���os��binaryģʽ�������˲�Ҫ�����������ʹ��
//�������͵Ķ������������Ӧ����binary_in��binary_out��ʵ��

//���T�����������ͣ�������os��ģʽ�Զ�ѡ�������ʽ������ͬʱ����ASCII��ʽ�����л��Ͷ����Ƹ�ʽ�����л�
//���T���������ͣ���ֻ��������ASCII��ʽ�����л�����Ҫ�����������л�Ҫ�õ������binary_in��binary_out������
template <typename T>
std::ostream& operator<< (std::ostream& os, const std::vector <T>& vt)
{
	size_t size(vt.size ());
	switch (os.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		os.write( (char*)&size, sizeof (size));
		for (int i(0);i!=size;++i)
		{
			if (!os)
				return os;
			os<<vt[i];
		}
		return os;
	default:
		os<<vt.size ()<<std::endl;
		for (int i(0);i!=size;++i)
		{
			if (!os)
				return os;
			os<<vt[i]<<std::endl;
		}
		return os;
	}
}

template <typename T>
std::istream& operator>> (std::istream& is, std::vector<T>& vt)
{
	size_t size;
	switch (is.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		is.read ( (char*)&size, sizeof (size));
		break;
	case CGAL::IO::ASCII:
		is>>size;
		break;
	default:
		throw std::logic_error ("operator>> (std::istream&, std::vector<T>");
	}

	vt.resize(size);
	for (int i(0);i!=(int) size;++i)
		is>>vt[i];
	return is;
}

//���ڲ������ض�ȡbool���͵����ݣ���ר��Ϊ������ݶ���һ��������
//ע�⣬bool���������ͣ��������os��is��ģʽ��binary�Ļ���ֱ�ӱ���.
std::ostream& operator<< (std::ostream &os, const std::vector<bool>& vt)
{
	if (os.iword (CGAL::IO::mode) == CGAL::IO::BINARY)
		throw std::logic_error ("operator<< (std::ostream &os, const std::vector<bool>& vt), stream should be ascii");
	os<<vt.size ()<<std::endl;
	for (int i(0);i!=vt.size ();++i)
		if (vt[i])
			os<<1<<' ';
		else 
			os<<0<<' ';
	return os<<std::endl;
}

std::istream& operator>> (std::istream &is, std::vector <bool>&vt)
{
	if (is.iword(CGAL::IO::mode) == CGAL::IO::BINARY)
		throw std::logic_error ("operator>> (std::istream &is, const std::vector <bool>&vt), stream should be ascii");
	size_t size;
	is>>size;
	vt.resize(size);
	for (int i(0);i!=size;++i)
	{
		int temp;
		is>>temp;
		if (temp == 1)
			vt[i] = true;
		else
			vt[1] = false;
	}
	return is;
}

//���T�����õģ���������������������������л�����������ļ�����Ҫ�õ������������
template <typename T>
void binary_out (std::ostream&os, const std::vector <T>& vt)
{
	size_t size(vt.size ());
	os.write((char*) (&size), sizeof(size));
	for (int i(0);i!=size;++i)
		os.write( (const char*) (&(vt[i])), sizeof (vt[i]));
}

void binary_out (std::ostream&os, const std::vector<bool> &vt)
{
	size_t size (vt.size ());
	os.write((char*) (&size), sizeof (size));
	for (int i(0);i!=size;++i)
	{
		char temp;
		if (vt[i])
			temp = '1';
		else 
			temp = '0';
		os.write((char*) (&temp),sizeof (temp));
	}
}

template <typename T>
void binary_out (std::ostream&os, const std::vector<std::vector<T>> &vvt)	//������Ϊ��ʱ��Ҫ������һ���������ϸ���������˵����Ҫbool�������¶���һ�����������������������������Ϊ����û���õ���ʡ����2016-6-9 17:12
{
	size_t size (vvt.size ());
	os.write((char*) (&size), sizeof(size));;
	for (int i(0);i!=size;++i)
		binary_out(os,vvt[i]);
}

template <class T>
void binary_in (std::istream&is, std::vector<T> &vt)
{
	size_t size;
	is.read ( (char*) (&size), sizeof (size));
	vt.resize(size);
	for (int i(0);i!=size;++i)
		is.read ( (char*) (&(vt[i])), sizeof(vt[i]));
}

void binary_in (std::istream&is, std::vector<bool> &vt)
{
	size_t size;
	is.read ((char*) &size, sizeof (size));
	vt.resize(size);
	for (int i(0);i!=size;++i)
	{
		char temp;
		is.read ((char*) (&temp), sizeof (temp));
		switch (temp)
		{
		case '0':
			vt[i] = false; break;
		case '1':
			vt[i] = true; break;
		default:
			throw std::logic_error("binary_in (std::istream&is, std::vector<bool> &vt), file with wrong format");
		}	
	}
}

template <typename T>
void binary_in (std::istream&is, std::vector<std::vector<T>> &vvt)
{
	size_t size;
	is.read ( (char*) (&size), sizeof (size));
	vvt.resize(size);
	for (int i(0);i!=size;++i)
		binary_in(is,vvt[i]);
}

void write_gmpq_raw (std::ostream &os, const CGAL::Gmpq &value)
{
	void (*freefunc) (void *, size_t);
	mp_get_memory_functions(NULL, NULL, &freefunc);

	int8_t sign = mpq_sgn(value.mpq());
	os.write( (char*)&sign, sizeof (sign));

	size_t c;
	void * raw1 = mpz_export(NULL, &c, 1, 1, 0, 0, value.numerator().mpz());
	os.write( (char*)&c, sizeof(c) );
	os.write( (char*)raw1, c );
	freefunc(raw1,c);

	void * raw2 = mpz_export(NULL, &c, 1, 1, 0, 0, value.denominator().mpz());
	os.write( (char*)&c, sizeof(c) );
	os.write( (char*)raw2, c );
	freefunc(raw2,c);
}

void read_gmpq_raw (std::istream &in, CGAL::Gmpq &value)
{
	CGAL::Gmpz num, den;
	size_t size1,size2;

	int8_t sign;
	in.read ( (char*)&sign, sizeof (sign));

	in.read ( (char*)&size1, sizeof(size1));
	uint8_t * vbit1 = new uint8_t[size1];
	in.read ( (char*)vbit1, size1);
	mpz_import (num.mpz (),size1, 1, 1, 0, 0, vbit1);

	in.read ( (char*) &size2, sizeof (size2));
	uint8_t *vbit2 = new uint8_t[size2];
	in.read ( (char*)vbit2, size2);
	mpz_import (den.mpz (), size2, 1, 1, 0, 0, vbit2);

	value = CGAL::Gmpq(num,den);
	if (sign<0)
		value = -value;
	delete []vbit1;
	delete []vbit2;
}

//�õ���CGAL��Ϊio�������״̬
//��Ȼ���Ѿ�ȫ������CGAL�ˣ��Ͳ���Ҫ�Ʒ���
//ͨ���ж������������״̬���ж�������ȡ��ʽ
inline std::ostream& operator<<(std::ostream& os, const FT &x)
{
	switch (os.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		write_gmpq_raw(os,x.exact ());
		return os;
	default:
		return os<<x.exact();
	}
}

inline std::istream& operator>>(std::istream& is, FT &x)
{
	CGAL::Gmpq temp;
	switch (is.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		read_gmpq_raw(is,temp);
		break;
	case CGAL::IO::ASCII:
		is>>temp ;
		break;
	default:
		throw std::logic_error ("operator>> (istream& is, FT &x), stream must be ascii or binary mode");
		break;
	}
	if (is)
		x = FT(temp);
	return is;
}

inline std::ostream& operator<< (std::ostream& os,const point_3 &p)
{
	switch (os.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		return os<<p.x()<<p.y()<<p.z();
	default:
		return os<<p.x()<<' '<<p.y()<<' '<<p.z();
	}
}

inline std::istream& operator>> (std::istream& is, point_3 &p)
{
	FT x,y,z;
	is>>x>>y>>z;
	if (is)
		p = point_3(x,y,z);
	return is;
}

inline std::ostream& operator<< (std::ostream&os, const point_2 &p)
{
	switch (os.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		return os<<p.x()<<p.y();
	default:
		return os<<p.x()<<' '<<p.y();
		break;
	}
}

inline std::istream& operator>> (std::istream&is, point_2 &p)
{
	FT x,y;
	is>>x>>y;
	if (is)
		p = point_2(x,y);
	return is;
}

inline std::ostream& operator<< (std::ostream&os, const plane_3 &p)
{
	switch (os.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		return os<<p.a()<<p.b()<<p.c()<<p.d();
	default:
		return os<<p.a()<<' '<<p.b()<<' '<<p.c()<<' '<<p.d();
	}
}

inline std::istream& operator>> (std::istream&is, plane_3 &p)
{
	FT a,b,c,d;
	is>>a>>b>>c>>d;
	if (is)
		p = plane_3(a,b,c,d);
	return is;
}

inline std::ostream& operator<< (std::ostream&os, const vector_3 &v)
{
	switch (os.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		return os<<v.x()<<v.y()<<v.z();
	default:
		return os<<v.x()<<' '<<v.y()<<' '<<v.z();
	}
}

inline std::istream& operator>> (std::istream&is, vector_3 &v)
{
	FT x,y,z;
	is>>x>>y>>z;
	if (is)
		v = vector_3(x,y,z);
	return is;
}

inline std::ostream& operator<< (std::ostream&os, const polygon_2 &poly)
{
	size_t size(poly.size ());
	switch (os.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		os.write((char*) &size, sizeof (size));
		for (int i(0);i!=size;++i)
			os<<poly[i];
		break;
	default:
		os<<size<<std::endl;
		for (int i(0);i!=size;++i)
			os<<poly[i]<<std::endl;
		break;
	}
	return os;
}

inline std::istream& operator>> (std::istream& is, polygon_2 &poly)
{
	size_t size;
	switch (is.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		is.read((char*) &size, sizeof (size));
		break;
	case CGAL::IO::ASCII:
		is>>size;
		break;
	default:
		throw std::logic_error("operator>> (std::istream& is, const polygon_2 &poly), stream should be binary or ascii");
	}
	poly.clear();
	for (int i(0);i!=size;++i)
	{
		point_2 temp;
		is>>temp;
		poly.push_back (temp);
	}
	return is;
}

inline std::ostream& operator<< (std::ostream& os, const polygon_with_holes_2 &pwh)
{
	os<<pwh.outer_boundary();
	size_t count (pwh.number_of_holes());
	switch(os.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		os.write((char*) &(count), sizeof (count));
		break;
	default:
		os<<count<<std::endl;
		break;
	}
	for (polygon_with_holes_2::Hole_const_iterator hi = pwh.holes_begin();
		hi!=pwh.holes_end();++hi)
		os<<(*hi);
	return os;
}

inline std::istream& operator>> (std::istream& is, polygon_with_holes_2 &pwh)
{
	polygon_2 poly;
	is>>poly;
	pwh = polygon_with_holes_2(poly);
	size_t count;
	switch (is.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		is.read ((char*) &count, sizeof (count));
		break;
	case CGAL::IO::ASCII:
		is>>count;
		break;
	default:
		throw std::logic_error("operator>> (std::istream& is, const polygon_with_holes_2 &pwh), stream should be binary or ascii");
		break;
	}
	for (int i(0);i!=count;++i)
	{
		polygon_2 temp;
		is>>temp;
		pwh.add_hole (temp);
	}
	return is;
}

inline std::ostream& operator<< (std::ostream &os, const disc_frac &x)
{
	switch (os.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		os.write((char*) &(x.plane_id), sizeof (x.plane_id));
		os.write((char*) &(x.frac_id ), sizeof (x.frac_id));
		return os<<x.cohesion<<x.f_angle<<x.aperture<<x.r<<x.center;
	default:
		return os<<x.plane_id<<' '<<x.frac_id<<' '<<x.cohesion<<' '<<x.f_angle<<' '<<x.aperture<<' '<<x.r<<std::endl
			<<x.center;
	}
}

inline std::istream& operator>> (std::istream &is, disc_frac& x)
{
	switch (is.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		is.read ((char*) &(x.plane_id), sizeof (x.plane_id));
		is.read ((char*) &(x.frac_id), sizeof (x.frac_id));
		break;
	case CGAL::IO::ASCII:
		is>>x.plane_id>>x.frac_id;
		break;
	default:
		throw std::logic_error("operator>> (std::istream &is, disc_frac& x), stream should be binary or ascii");
	}
	is>>x.cohesion>>x.f_angle>>x.aperture>>x.r>>x.center;
	return is;
}

inline std::ostream& operator<< (std::ostream &os, const poly_frac &x)
{
	switch(os.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		os.write((char*) &(x.plane_id), sizeof (x.plane_id));
		os.write((char*) &(x.frac_id), sizeof (x.frac_id));
		return os<<x.cohesion<<x.f_angle<<x.aperture<<x.vp_ob<<x.vp_ib;
	default:
		return os<<x.plane_id<<' '<<x.frac_id<<' '<<x.cohesion<<' '<<x.f_angle<<' '<<x.aperture<<std::endl
			<<x.vp_ob<<x.vp_ib;
	}
}

inline std::istream& operator>> (std::istream &is, poly_frac& x)
{
	switch (is.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		is.read ((char*) &(x.plane_id), sizeof (x.plane_id));
		is.read ((char*) &(x.frac_id), sizeof (x.frac_id ));
		break;
	case CGAL::IO::ASCII:
		is>>x.plane_id>>x.frac_id;
		break;
	default:
		throw std::logic_error("operator>> (std::istream &is, poly_frac& x), stream should be binary or ascii");
	}
	is>>x.cohesion>>x.f_angle>>x.aperture>>x.vp_ob>>x.vp_ib;
	return is;
}

inline std::ostream& operator<< (std::ostream &os, const eblock::face &x)
{
	switch (os.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		binary_out(os,x.np);
		os.write ((char*) &(x.frac_id), sizeof (x.frac_id));
		os.write ((char*) &(x.plane_id), sizeof (x.plane_id));
		os.write ((char*) &(x.outer_normal), sizeof (x.outer_normal));
		return os;
	default:
		return os<<x.np<<' '<<x.frac_id<<' '<<x.plane_id<<' '<<x.outer_normal;
	}
}

inline std::istream& operator>> (std::istream &is, eblock::face& x)
{
	switch (is.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		binary_in(is,x.np);
		is.read ((char*) &(x.frac_id), sizeof (x.frac_id));
		is.read ((char*) &(x.plane_id), sizeof (x.plane_id));
		is.read ((char*) &(x.outer_normal), sizeof (x.outer_normal));
		break;
	case CGAL::IO::ASCII:
		is>>x.np>>x.frac_id>>x.plane_id>>x.outer_normal;
		break;
	default:
		throw std::logic_error("operator>> (std::istream &is, eblock::face& x), stream should be binary or ascii");
		break;
	}
	return is;
}

inline std::ostream& operator<< (std::ostream &os, const eblock &x)
{
	switch (os.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		os<<x.vf<<x.attribute;
		os.write ((char*) &(x.cblock_index), sizeof (x.cblock_index));
		return os;
	default:
		return os<<x.vf<<' '<<x.attribute<<' '<<x.cblock_index;
	}
}

inline std::istream& operator>> (std::istream &is, eblock& x)
{
	switch (is.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		is>>(x.vf)>>(x.attribute );
		is.read ((char*) &(x.cblock_index), sizeof (x.cblock_index));
		break;
	case CGAL::IO::ASCII:
		is>>(x.vf)>>(x.attribute)>>(x.cblock_index);
		break;
	default:
		throw std::logic_error("operator>> (std::istream &is, eblock& df), stream should be binary or ascii");
		break;
	}
	return is;
}

inline std::ostream& operator<< (std::ostream &os, const cblock &x)
{
	switch(os.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		binary_out(os,x.vb);
		return os;
	default:
		return os<<x.vb;
	}
}

inline std::istream& operator>> (std::istream &is, cblock& x)
{
	switch (is.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		binary_in(is,x.vb);
		break;
	case CGAL::IO::ASCII:
		is>>x.vb;
		break;
	default:
		throw std::logic_error("operator>> (std::istream &is, cblock& df), stream should be binary or ascii");
		break;
	}
	return is;
}

inline std::ostream& operator<< (std::ostream& os, const block_system &bs)
{
	switch (os.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		os<<bs.vpo<<bs.vpl<<bs.vdisc<<bs.vpoly<<bs.veb<<bs.vcb;
		os.write((char*) &(bs.cutting_time), sizeof(bs.cutting_time));
		os.write((char*) &(bs.shrk_time), sizeof (bs.shrk_time));
		os<<bs.vbsf;
		binary_out(os,bs.vcb_bsf);
		binary_out(os,bs.vpl_bsf);
		return os;
	default:
		return os<<bs.vpo<<bs.vpl<<bs.vdisc<<bs.vpoly<<bs.veb<<bs.vcb
			<<' '<<bs.cutting_time<<' '<<bs.shrk_time<<bs.vbsf<<bs.vcb_bsf<<bs.vpl_bsf;
	}
}

inline std::istream& operator>> (std::istream& is, block_system &bs)
{
	is>>bs.vpo>>bs.vpl>>bs.vdisc>>bs.vpoly>>bs.veb>>bs.vcb;
	switch (is.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		is.read ((char*) &(bs.cutting_time), sizeof (bs.cutting_time));
		is.read ((char*) &(bs.shrk_time), sizeof (bs.shrk_time));
		is>>bs.vbsf;
		if(!is)
			return is;
		binary_in(is,bs.vcb_bsf);
		binary_in(is,bs.vpl_bsf);
		break;
	case CGAL::IO::ASCII:
		is>>bs.cutting_time>>bs.shrk_time>>bs.vbsf>>bs.vcb_bsf>>bs.vpl_bsf;
		break;
	default:
		throw std::logic_error("operator>> (std::istream &is, block_system& df), stream should be binary or ascii");
		break;
	}
	return is;
}

inline std::ostream& operator<< (std::ostream &os, const outer_bound& ob)
{
	switch (os.iword (CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		os<<ob.vpoly<<ob.vplane;
		binary_out(os, ob.outer_normal);
		break;
	default:
		os<<ob.vpoly<<ob.vplane<<ob.outer_normal;
		break;
	}
	return os;
}

inline std::istream& operator>> (std::istream &is, outer_bound &ob)
{
	switch (is.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		is>>(ob.vpoly)>>(ob.vplane);
		binary_in(is,ob.outer_normal);
		break;
	case CGAL::IO::ASCII:
		is>>(ob.vpoly)>>(ob.vplane)>>(ob.outer_normal);
		break;
	default:
		throw std::logic_error("operator>> (std::istream &is, outer_bound& df), stream should be binary or ascii");
		break;
	}
	return is;
}

inline std::ostream& operator<< (std::ostream &os, const inner_bound& ib)
{
	return os<<ib.vpoly<<ib.vplane;
}

inline std::istream& operator>> (std::istream &is, inner_bound &ib)
{
	switch (is.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
	case CGAL::IO::ASCII:
		return is>>(ib.vpoly)>>(ib.vplane);
	default:
		throw std::logic_error("operator>> (std::istream &is, inner_bound& df), stream should be binary or ascii");
		break;
	}
	
}

inline std::ostream& operator<< (std::ostream& os, const block_system::FACE_TYPE &ft)
{
	char temp;
	switch (ft)
	{
	case block_system::FACE_TYPE::EXCAV:
		temp = '0'; break;
	case block_system::FACE_TYPE::FRAC:
		temp = '2'; break;
	default:
		throw std::logic_error ("operator<< (std::ostream& os, const block_system::FACE_TYPE &ft), unexpected type");
		break;
	}
	switch (os.iword (CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		os.write ((char*) (&temp), sizeof (temp));
		break;
	default:
		os<<temp<<std::endl;
		break;
	}
	return os;
}

inline std::istream& operator>> (std::istream& is, block_system::FACE_TYPE &ft)
{
	char temp;
	switch (is.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		is.read ((char*) &temp, sizeof (temp));
		break;
	case CGAL::IO::ASCII:
		is>>temp;
		break;
	default:
		throw std::logic_error ("operator>> (std::istream& is, block_system::FACE_TYPE &ft)");
	}
	switch (temp)
	{
	case '0':
		ft = block_system::FACE_TYPE::EXCAV; break;
	case '2':
		ft = block_system::FACE_TYPE::FRAC; break;
	default:
		throw std::logic_error ("operator>> (std::istream& is, block_system::FACE_TYPE &ft), wrong type");
		break;
	}
	return is;
}

inline std::ostream& operator<< (std::ostream &os, const block_system::bound_surface& bsf)
{
	switch (os.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		os<<bsf.poly;
		os.write((char*) &(bsf.pli), sizeof (bsf.pli));
		os.write((char*) &(bsf.cbi), sizeof (bsf.cbi));
		char temp;
		if (bsf.outer_normal)
			temp = '1';
		else
			temp = '0';
		os.write((char*) (&temp),sizeof (temp));
		binary_out(os,bsf.vadj_bsf_i);
		binary_out(os,bsf.vadj_cb_i);
		os<<bsf.type;
		break;
	default:
		os<<bsf.poly<<' '<<bsf.pli<<' '<<bsf.cbi<<' ';
		if (bsf.outer_normal)
			os<<'1';
		else
			os<<'0';
		os<<' '<<bsf.vadj_bsf_i<<bsf.vadj_cb_i<<' '<<bsf.type;
		break;
	}
	return os;
}

inline std::istream& operator>> (std::istream& is, block_system::bound_surface& bsf)
{
	is>>bsf.poly;
	char temp;
	switch (is.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		is.read ((char*) &(bsf.pli), sizeof (bsf.pli));
		is.read ((char*) &(bsf.cbi), sizeof (bsf.cbi));
		is.read ((char*) &(temp), sizeof (temp));
		if (temp == '1')
			bsf.outer_normal = true;
		else
			bsf.outer_normal = false;
		binary_in(is,bsf.vadj_bsf_i);
		binary_in(is,bsf.vadj_cb_i);
		is>>bsf.type;
		break;
	case CGAL::IO::ASCII:
		is>>bsf.pli>>bsf.cbi;
		is>>temp;
		if (temp == '1')
			bsf.outer_normal = true;
		else
			bsf.outer_normal = false;
		is>>bsf.vadj_bsf_i>>bsf.vadj_cb_i>>bsf.type;
		break;
	default:
		throw std::logic_error("operator>> (std::istream& is, block_system::bound_surface& bsf)");
	}
	return is;
}

inline std::ostream& operator<< (std::ostream& os, const entrance_block& eb)
{
	switch (os.iword(CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
		os<<eb.vcover<<eb.outer_b<<eb.bs<<eb.inner;
		break;
	default:
		os<<eb.vcover<<' '<<eb.outer_b<<' '<<eb.bs<<' '<<eb.inner;
		break;
	}
	return os;
}

inline std::istream& operator>> (std::istream& is, entrance_block& eb)
{
	switch (is.iword (CGAL::IO::mode))
	{
	case CGAL::IO::BINARY:
	case CGAL::IO::ASCII:
		is>>eb.vcover>>eb.outer_b>>eb.bs>>eb.inner;
		break;
	default:
		throw std::logic_error("operator>> (std::istream &is, entrance_block& df), stream should be binary or ascii");
		break;
	}
	return is;
}

}