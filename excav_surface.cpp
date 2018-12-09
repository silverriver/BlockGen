#include "excav_surface.h"
#include "GBUtility.h"
#include "outer_poly.h"
#include <vector>
using namespace std;

namespace BI
{

excav_surface::excav_surface (const block_system &bs)
{
	vpoly.clear();
	vplane.clear();

	vector<int> veb;
	veb.reserve (bs.veb.size ());
	for (int i(0);i!=bs.veb.size ();++i)
		veb.push_back (i);
	outer_bound ob(bs,veb);

	vpoly.reserve (ob.vplane.size ());
	vplane.reserve (ob.vplane.size ());

	for (int i(0);i!=ob.vplane.size ();++i)
	{
		int plane_index (add_set_pl(vplane,ob.vplane[i]));
		if (vplane.size () > vpoly.size ())
			vpoly.push_back (polygon_set_2());
		vpoly[plane_index].join(ob.vpoly[i]);
	}
}

}