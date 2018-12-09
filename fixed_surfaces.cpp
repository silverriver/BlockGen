#include "fixed_surfaces.h"
#include "geometry.h"
#include "identification.h"
#include <string>
#include <vector>
#include <CGAL\intersections.h>

using namespace std;
namespace BI
{

int fixed_surfaces::STL_parser (istream& infile)
{
	if (!infile)
		return 1;
	string kw1, kw2;

	char temp('a');
	while (temp!='\n')
		infile.get (temp);

	while (infile)
	{
		double x[3],y[3],z[3];
		infile>>kw1>>kw2;
		if (kw1 == "endsolid")
			break;
		if (kw1 != "facet")
			return 2;
		getline (infile,kw1);
		infile>>kw1>>kw2;
		if (kw1 !="outer" || kw2 != "loop")
			return 2;
		for (int i(0);i!=3;++i)
		{
			infile>>kw1;
			if (kw1 != "vertex")
				return 2;
			infile>>x[i]>>y[i]>>z[i];
		}
		infile >> kw1>>kw2;
		if (kw1 != "endloop" || kw2!= "endfacet")
			return 2;
		vtri.push_back (triangle_3 (point_3 (x[0],y[0],z[0]),point_3(x[1],y[1],z[1]), point_3(x[2],y[2],z[2])));
		if (vtri.back ().is_degenerate())
			vtri.pop_back();
	}
	return 0;
}

bool fixed_surfaces::do_intersect (const std::vector<point_3> &vpo, const std::vector<plane_3> &vpl, const eblock &eb) const
{
	mbox_3 mbox(build_mbox_3(vpo,eb));
	int count(0);
	for (;count!=vtri.size ();++count)
		if (CGAL::do_intersect(vtri[count],mbox))
			break;
	if (count == vtri.size ())
		return false;		//和eb的mbox都不相交，直接跳出

	for (int i(0);i!=(int)eb.vf.size ();++i)
	{
		int t1 (eb.vf[i].np[0]);
		for (int j(1);j < (int) eb.vf[i].np.size ()-1;++j)
		{
			triangle_3 t_temp(vpo[t1], vpo[eb.vf[i].np[j]], vpo[eb.vf[i].np[j+1]]);
			for (int k(0);k!=(int) vtri.size ();++k)
				if (CGAL::do_intersect (t_temp,vtri[k]))
					return true;		//和eb的某个边界相交则算是相交
		}
	}

	for (int i(0);i!=(int) vtri.size ();++i)
	{
		int j(0);
		for (;j!=(int) eb.vf.size ();++j)
		{
			CGAL::Orientation ori = vpl[eb.vf[j].plane_id].oriented_side(vtri[i].vertex(0));
			if (ori == CGAL::Orientation::ON_NEGATIVE_SIDE && (!eb.vf[j].outer_normal))		//如果这个点在eb外面，则提前跳出
				break;
			else if (ori == CGAL::Orientation::ON_POSITIVE_SIDE && eb.vf[j].outer_normal )	//如果这个点在eb外面，则提前跳出
				break;
		}
		if (j == eb.vf.size ())		//这个三角形在包含在eb内。所以肯定相交
			return true;
	}

	return false;		//eb的边界不与vtri不相交，并且也不包含vtri中任何一个三角形的顶点,所以是不相交
}

bool fixed_surfaces::do_intersect (const std::vector<point_3> &vpo, const std::vector<plane_3> &vpl, const std::vector<eblock> &veb, const cblock &cb) const
{
	for (int i(0);i!=cb.vb.size ();++i)
		if (do_intersect (vpo,vpl,veb[cb.vb[i]]))
			return true;
	return false;
}

bool fixed_surfaces::do_intersect (const block_system &bs, const eblock &eb) const
{
	return do_intersect(bs.vpo ,bs.vpl,eb);
}
bool fixed_surfaces::do_intersect (const block_system &bs, const cblock &cb) const
{
	return do_intersect (bs.vpo, bs.vpl, bs.veb, cb);
}


}
