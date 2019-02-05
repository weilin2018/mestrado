# routine created by MSc Carine Godoi

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.cm as cm
from mpl_toolkits.axes_grid.inset_locator import inset_axes

def sigma2stdl(variable, sigma, nstdl, depth, ang, h1, name, lon, lat, figpath):

	print('\ninterpolating sigma to standard levels: ' + name)

	# lat and lon limits to draw grid
	lai = np.floor(np.nanmin(lat))
	lae = np.ceil(np.nanmax(lat))
	loi = np.floor(np.nanmin(lon))
	loe = np.ceil(np.nanmax(lon))
	m = Basemap(projection='merc',llcrnrlat=lai,urcrnrlat=lae,llcrnrlon=loi,urcrnrlon=loe,resolution='i');
	X, Y = m(lon,lat)

	# nsteps, sigmalevels, lines, columns = variable.shape
	nsteps, sigmalevels, lines, columns = variable.shape

	variableI = np.zeros((nsteps, nstdl, lines, columns))
	variableI[:] = np.NAN

	if nstdl == 37:
		stdl = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 1000, 1200, 1500, 1800, 2000]
	elif nstdl == 23:
		stdl = [0, 10, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 1000, 1200, 1500, 1800, 2000]


	### --- vectorize data to decrease number of for loops --- ###

	vecvar = np.zeros((sigmalevels, nsteps*lines*columns))
	vecvarI = np.zeros((nstdl, nsteps*lines*columns))
	localdep = np.zeros((nsteps*lines*columns))
	cont = 0

	for n in np.arange(0, nsteps):
		# print('nstep ' + str(n) + ' of ' + str(nsteps-1))
		for l in np.arange(0, lines):
			for c in np.arange(0, columns):
				vecvar[:, cont] = variable[n, :, l, c]
				localdep[cont] = depth[l, c]
				cont += 1


	### --- interpolate sigma to standard levels --- ###

	# plot counter
	kplot = int(nsteps*lines*columns/10)
	k = 0

	for i in np.arange(0, nsteps*lines*columns):
		k += 1
		if k == kplot:
			print(str(np.round(i/(nsteps*lines*columns)*100)) + '%')
			k = 0

		if ~np.isnan(localdep[i]):
			# print('yes 1/2')

			# levels of data (m)
			depsigma = -localdep[i]*sigma

			# include surface with same value of first sigma level m to interpolate
			D = list(depsigma)
			D.insert(0, 0)

			# select profile and include surface
			profile = np.zeros(sigmalevels+1)
			profile[1:] = vecvar[:, i]
			profile[0] = profile[1]

			# watercolumn positions only
			watercolumn = stdl <= localdep[i]
			stdl2interp = np.array(stdl)[watercolumn]

			# interpolate to the same standard levels
			fsigma2stdl = interpolate.interp1d(D, profile)
			profileI = fsigma2stdl(stdl2interp)

			# stores at vectorized variable
			vecvarI[watercolumn, i] = profileI
			vecvarI[~watercolumn, i] = np.NAN

			# print('yes 2/2')
			# if k == kplot:
			# 	print('printing figure')
			# 	k = 0

			# 	plt.figure();
			# 	plt.plot(profile, -np.array(D), marker='o', linewidth=2.5);
			# 	plt.plot(profileI, -stdl2interp, marker='o', color='red');
			# 	plt.legend(['sigma','interp stdl'], loc='best');
			# 	x1, x2 = plt.xlim()
			# 	plt.plot([x1, x2], [-localdep[i], -localdep[i]], color='brown', linewidth=3);
			# 	plt.ylabel('depth (m)');
			# 	plt.xlabel(name)
			# 	plt.savefig('../figures/' + path + '/' + run + '/sigma2stdl_ex' + str(i) + name + '.png', bbox_inches='tight', dpi = 130);
			# 	plt.close();


	### --- back to original shape --- ###

	cont = 0

	for n in np.arange(0, nsteps):
		for l in np.arange(0, lines):
			for c in np.arange(0, columns):
				variableI[n, :, l, c] = vecvarI[:, cont]
				cont += 1

	# corners to NaN
	variableI[:, :, 1, -2] = np.NAN
	variableI[:, :, -2, -2] = np.NAN


	### --- cross-shore sections --- ###

	if name != 'temperature' and name != 'salinity':

		# for lin in [1, 12, 25, 37, 50, 62, 75, 87, 100, 112, 120, 125, 128, 130, 135]:
		for lin in np.arange(1, lines-1):
		# for lin in [129]:
			tim = nsteps-1
			# tim = 0
			istart = np.min(np.where(variableI[tim, 0, lin, :] != 0))
			i200 = np.min(np.where(depth[lin, :] >= 200))
			i100 = np.min(np.where(depth[lin, :] >= 100))

			# cols = h1[lin,istart:-1]
			cols = h1[lin, i100:-1]
			cols = np.cumsum(cols)
			# cols = (cols - cols[0])/1000
			cols = cols/1000
			# iend = np.nanmax(np.where(cols <= 140))
			# cols = cols[0:iend]

			if nstdl == 37:
				p = 34
			elif nstdl == 23:
				p = 20
			# var2plot = variableI[tim, 0:p, lin, istart:-1]
			var2plot = variableI[tim, 0:p, lin, i100:-1]
			# var2plot = var2plot[:, 0:iend]

			#c2 = (np.nanmax(np.abs(variableI[variableI!=0])))
			#if lin == 129:
			#	c2 = (np.nanmax(np.abs(var2plot[var2plot!=0])))
			#passo = 3*c2/(nstdl-1)
			#clevs = np.arange(-c2, c2+(passo/2), passo)
			c2 = 0.5#float('%.2f' %(np.abs(np.nanmax(var2plot))))
			passo = 0.05
			clevs = np.arange(-c2, c2+passo/2, passo)

			cmap2 = cm.bwr
			cmap2.set_bad('black','NaN')

			fig = plt.figure();
			ax = plt.subplot(111);
			plt.gca().patch.set_color('.25')
			plt.contourf(cols, -np.array(stdl[0:p]), var2plot, clevs, cmap=cmap2);
			cbar = plt.colorbar(format='%.2f');
			cbar.set_label('m/s', labelpad=-38, y=1.05, rotation=0);
			cs = plt.contour(cols, -np.array(stdl[0:p]), var2plot, clevs, colors='black');
			plt.clabel(cs, fontsize=9, inline=1, fmt='%.2f');
			plt.ylabel('depth (m)');
			plt.xlabel('distance from the 100 m isobath (km)');
			# plt.xlabel('distance from the coast (km)');
			plt.title(name[0:-9] + ' ' + name[-8:] + ', {:.1f}'.format(abs(lat[lin, i200])) + '$^\circ$S');
			if cols[-1] > 350:
				plt.xlim(0, 350);
				iend = istart + np.min(np.where(cols > 350))
			else:
				#iend = -1
				iend= np.min(np.where(depth[lin,:]>=1200))


			from mpl_toolkits.axes_grid.inset_locator import inset_axes
			inset_axes = inset_axes(ax, width="30%", height=1.1, loc=3);
			m.drawcoastlines();
			m.fillcontinents(color='wheat');
			lon1, lat1 = lon[lin, i100:iend], lat[lin, i100:iend]
			x, y = m(lon1, lat1)
			m.plot(X[1,:], Y[1,:], color='gray', linewidth=1);
			m.plot(X[-2,:], Y[-2,:], color='gray', linewidth=1);
			m.plot(X[:,-2], Y[:,-2], color='gray', linewidth=1);
			m.plot(x, y, color='red', linewidth=2);

			plt.savefig(figpath + '/sigma2stdl_' + name + '_J' + str(lin+1) + '.png', bbox_inches='tight', dpi = 150);
			plt.close();

	return variableI
