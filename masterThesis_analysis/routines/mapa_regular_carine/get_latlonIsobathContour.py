from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

def get(m, Xgrid, Ygrid, depgrid, isobath):

	"""
	extract lat and lon coordinates of given isobath depth
	Carine G. R. Costa, 2017
	"""

	plt.close('all')

	cs = m.contour(Xgrid, Ygrid, depgrid, isobath, latlon=True, colors='darkgray', linewidths=.5);

	p0 = cs.collections[0].get_paths()[0]
	v = p0.vertices
	x = v[:, 0]
	y = v[:, 1]

	loniso, latiso = m(x, y, inverse=True)

	plt.close('all')

	return loniso, latiso

# def polygon(x, y):

# 	# create path of closed polygon
# 	pol, patches = [], []
# 	for i in range(len(x)):
# 		pol.append((x[i], y[i]))
# 	p = path.Path(pol)
# 	polygon = Polygon(pol, closed=True)
# 	patches.append(polygon)
