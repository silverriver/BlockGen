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
//这个文件包含一些做序列化输出的函数
//难点是找一个方法序列化FT，现在的主要做法就是首先计算FT对应的严格类型
//因为目前阶段FT中用到了区间运算和lazy computation。所以直接做serialize有些困难。
//读取时用Gmpq作为中间变量，首先在流中读取Gmpq，然后用这个严格类型初始化FT

//使用须知：
//序列化T：
//如果T是内置类型，则手动调用ostream::write和istream::in来实现二进制输入和输出 （注意，当T是bool类型是需要手动转换一下，一般来说将true对应为'1'，false对应于'0'，在读取时进行逆变换）
//如果T不是内置类型，则直接用<<和>>操作符来实现

//序列化vector<T>:
//如果T是内置类型，则调用binary_in和binary_out来实现二进制输入和输出 （此时不管T是否是bool都可以直接调用这两个函数），用<<和>>操作符实现ASCII的输入和输出
//如果T不是内置类型，则直接用<<和>>操作符来实现

//序列化vector<vector<T>>
//如果T是内置类型，则调用binary_in和binary_out来实现二进制输入和输出。用<<和>>操作符进行ASCII格式的输入和输出 （没有实现bool类型的函数。需要的时候再写）
//如果T不是内置类型，则直接用<<和>>操作符来实现
//（当时设计的时候脑子抽风，想要用操作符来实现二进制的输入输出，试图装逼一把。没想后果是函数的作用有点混乱。失策。改天有空大修一下。2016-6-9 17:42）
namespace BI
{
//如果T是内置类型，则用操作符<<和>>只会输出ascii格式的文件，无法重载操作符使其输出binary格式的文件。此时如果os是binary模式会出错，因此不要在这种情况下使用
//内置类型的二进制输入输出应该用binary_in和binary_out来实现

//如果T不是内置类型，则会根据os的模式自动选择输出格式。可以同时进行ASCII格式的序列化和二进制格式的序列化
//如果T是内置类型，则只能用来做ASCII格式的序列化。想要做二进制序列化要用到下面的binary_in和binary_out函数。
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

//由于不能愉快地读取bool类型的数据，所专门为这个数据定义一个函数。
//注意，bool是内置类型，所以如果os和is的模式是binary的话则直接报错.
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

//如果T是内置的，则不能用上面的两个函数来做序列化输出二进制文件，需要用到下面这个函数
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
void binary_out (std::ostream&os, const std::vector<std::vector<T>> &vvt)	//这里因为临时需要加了这一个函数，严格意义上来说还需要bool类型重新定义一个函数做二进制输入输出，但是因为现在没有用到就省略了2016-6-9 17:12
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

//用到了CGAL中为io流定义的状态
//既然都已经全面依赖CGAL了，就不用要牌坊了
//通过判断输入输出流的状态来判断输出或读取格式
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