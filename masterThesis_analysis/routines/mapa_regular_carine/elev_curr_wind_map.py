import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import cmocean
from scipy import interpolate
import get_latlonIsobathContour as latlonIsobath
import plot_basicMap as basicMap
import fnearest

def map(II, JJ, Xgrid, Ygrid, depgrid, elev, u, v, wu, wv, time, trun, x_psi, y_psi):

	"""

	----- from model_grid file: II, JJ, Xgrid, Ygrid, depgrid
	II and JJ are the I and J numbers printed on model_grid
	Xgrid.shape = (JJ-2, II-2)
	same for Ygrid, depgrid

	----- from gcmplt.cdf file: elev, u, v, wu, wv, time
	elev.shape = (NSTEPS, sigma, JJ, II)
	same for u, v, wu, wv

	----- user defined:
	trun: time step to plot figure
	x_psi, y_psi: corners of each grid cell
	x_psi.shape = (JJ-1, II-1)

	Carine G. R. Costa, 2017

	"""

	# positions to draw map scale
	lon1,lat1 = -48, -22.6
	lengthscale = 100

	# position to write max and mean currents
	xcur, ycur = 0.52, 0.22

	lai2 = -28
	lae2 = -22.3
	loi2 = -48.8
	loe2 = -41.2
	m200 = Basemap(projection='merc',llcrnrlat=lai2,urcrnrlat=lae2,llcrnrlon=loi2,urcrnrlon=loe2,resolution='i');
	meridians = np.arange(-48, -38, 2)
	parallels = np.arange(-28, -20, 2)

	# to make pcolor
	xplot2, yplot2 = m200(x_psi, y_psi)

	# remove everything deeper than 200 m
	elev200 = np.copy(elev[:, 1:-1, 1:-1])
	lon200, lat200 = np.copy(Xgrid), np.copy(Ygrid)
	u200, v200 = np.copy(u[:, :, 1:-1, 1:-1]), np.copy(v[:, :, 1:-1, 1:-1])
	wu200, wv200 = np.copy(wu[:, 1:-1, 1:-1]), np.copy(wv[:, 1:-1, 1:-1])
	for j in range(JJ-2):
		for i in range(II-2):
			if depgrid[j, i] > 200:
				elev200[:, j, i] = np.nan
				lon200[j, i], lat200[j, i] = np.nan, np.nan
				u200[:, :, j, i], v200[:, :, j, i] = np.nan, np.nan
				wu200[:, j, i], wv200[:, j, i] = np.nan, np.nan
	# remove 6 first and last rows
	elev200[:, :6, :] = np.nan
	elev200[:, -6:, :] = np.nan
	u200[:, :, :6, :] = np.nan
	u200[:, :, -6:, :] = np.nan
	v200[:, :, :6, :] = np.nan
	v200[:, :, -6:, :] = np.nan
	wu200[:, :6, :] = np.nan
	wu200[:, -6:, :] = np.nan
	wv200[:, :6, :] = np.nan
	wv200[:, -6:, :] = np.nan
	dep = np.copy(depgrid)
	dep[:6, :] = np.nan
	dep[-6:, :] = np.nan

	# interpolate on regular grid
	# desired number of lines (J) and  columns (I)
	jr, ir = 40, 40
	lons, lats, x, y = m200.makegrid(jr, ir, returnxy=True)
	ut, vt = np.ravel(u200[trun, 0, :, :]), np.ravel(v200[trun, 0, :, :])
	X1, Y1 = np.ravel(Xgrid), np.ravel(Ygrid)
	X1, Y1 = X1[~np.isnan(ut)], Y1[~np.isnan(ut)]
	ut, vt = ut[~np.isnan(ut)], vt[~np.isnan(vt)]
	wut, wvt = np.ravel(wu200[trun, :, :]), np.ravel(wv200[trun, :, :])
	wut, wvt = wut[~np.isnan(wut)], wvt[~np.isnan(wvt)]
	points = np.array([X1, Y1]).T
	uI = interpolate.griddata(points, ut, (lons, lats), method='linear')
	vI = interpolate.griddata(points, vt, (lons, lats), method='linear')
	wuI = interpolate.griddata(points, wut, (lons, lats), method='linear')
	wvI = interpolate.griddata(points, wvt, (lons, lats), method='linear')
	maxvI = np.nanmax(np.sqrt(uI**2 + vI**2))
	meanvI = np.nanmean(np.sqrt(uI**2 + vI**2))
	maxwI = np.nanmax(np.sqrt(wuI**2 + wvI**2))

	# remove regular data located on grid points deeper than 200 m
	loniso, latiso = latlonIsobath.get(m200, Xgrid, Ygrid, depgrid, [200])
	A = np.array([loniso, latiso]).T
	for j in range(jr):
		for i in range(ir):
			point = [lons[j, i], lats[j, i]]
			distance, index = fnearest.find(A, point)
			pnear = [A[index, 0], A[index, 1] ]
			# if point is east of isobath, interp values receive NAN
			if point[0] >= pnear[0] and point[1] < pnear[1]:
				wuI[j, i], wvI[j, i] = np.NAN, np.NAN
				uI[j, i], vI[j, i] = np.NAN, np.NAN
			elif point[1] < pnear[1]:
				wuI[j, i], wvI[j, i] = np.NAN, np.NAN
				uI[j, i], vI[j, i] = np.NAN, np.NAN

	# max elev over time
	me = np.nanmax(elev200)

	fig = basicMap.map(m200, Xgrid, Ygrid, meridians, parallels, lon1, lat1, lengthscale, gridname, ' ', ' ', ' ');
	basicMap.contourf(m200, Xgrid, Ygrid, elev200[trun, :, :], 'elev (m)', -me, me, cmocean.cm.balance);
	basicMap.contour(m200, Xgrid, Ygrid, dep);
	basicMap.wind(m200, lons, lats, wuI, wvI, maxwI)
	basicMap.quiver(fig, m200, lons, lats, uI, vI, maxvI, meanvI, xcur, ycur, 2, 'black')
	plt.title(str(time[trun]) + ' days');
	plt.savefig('elev_curr_wind_Map_' + str(trun) + '_.png', bbox_inches='tight', dpi=120);
	plt.show()
