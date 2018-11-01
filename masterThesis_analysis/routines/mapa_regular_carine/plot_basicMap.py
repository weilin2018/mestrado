# maps...
# Carine G. R. Costa, 2017

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import cmocean

def map(m, lon, lat, meridians, parallels, lon1, lat1, lengthscale, gridname, lonsantos, latsantos, patches):

	fig,ax = plt.subplots(1);
	m.drawmeridians(meridians, labels=[0,0,0,1], color='0.2', dashes=[2,6], linewidth=0.3, zorder=95);
	m.drawparallels(parallels, labels=[1,0,0,0], color='0.2', dashes=[2,6], linewidth=0.3, zorder=96);
	m.drawmapscale(lon1, lat1, np.nanmean(lon), np.nanmean(lat), lengthscale, barstyle='fancy', zorder=100);

	if gridname == 'SSVBES':

		m.plot(lonsantos, latsantos, color='black', linewidth=.2, zorder=92);
		collection = PatchCollection(patches);
		ax.add_collection(collection);
		collection.set_zorder(91);
		collection.set_color([.9, .9, .9, 1]);

	elif gridname == 'SBB':

		m.fillcontinents(color=[.9, .9, .9, 1], zorder=92);
		m.drawcoastlines(color='black', linewidth=.2, zorder=93);
		m.drawstates(color='black', linewidth=.1, zorder=94);

	return fig


def pcolor(m, xplot, yplot, varback, name, vmin, vmax, cmap):

	cmap.set_under([1,1,1])
	p = m.pcolor(xplot, yplot, varback, vmin=vmin, vmax=vmax, cmap=cmap);
	cbar = m.colorbar(p);
	cbar.set_label(name, labelpad=-40, y=1.05, rotation=0);

def contourf(m, xplot, yplot, varback, name, vmin, vmax, cmap):

	passo = (vmax-vmin)/30
	clevs = np.arange(vmin, vmax+passo, passo)
	p = m.contourf(xplot, yplot, varback, clevs, latlon=True, cmap=cmap);
	cbar = m.colorbar(p, format='%.1f');
	cbar.set_label(name, labelpad=-40, y=1.05, rotation=0);

def contour(m, lon, lat, depgrid):

	cs = m.contour(lon, lat, depgrid, [50, 100, 200], latlon=True, colors='darkgray', linewidths=.5);
	plt.clabel(cs, inline=2, fontsize=12, fmt='%1i');

def scatter(m, xplot, yplot, j, i, cor):

	m.scatter(xplot[j, i], yplot[j, i], color=cor, s=70, marker='o', zorder=200);

def wind(m, lon, lat, wuplot, wvplot, maxw):

	q = m.quiver(lon[::4, ::4], lat[::4, ::4], wuplot[::4, ::4], wvplot[::4, ::4], latlon=True, headwidth=4, headlength=5, width=0.01, scale=maxw*8.5, alpha=.5, color='none', edgecolors=('k'), linewidths=(1));
	plt.quiverkey(q, 0.77, 0.23, maxw, r'max wind speed = $%.1f$ m/s' %maxw, labelpos='W');

def quiver(fig, m, lon, lat, uplot, vplot, maxv, meanv, xcur, ycur, p, cor):

	q = m.quiver(lon[::p, ::p], lat[::p, ::p], uplot[::p, ::p], vplot[::p, ::p], latlon=True, headwidth=4, headlength=5, width=0.002, scale=maxv*8.5, color=cor);
	plt.quiverkey(q, 0.77, 0.05, maxv, r'max surf current = $%.1f$ m/s' %maxv, labelpos='W');
	# fig.text(xcur, ycur,'max = ' + '{:.1f}'.format(maxv) + ' m/s', fontsize=10);
	fig.text(xcur, ycur-.04,'mean surf current = ' + '{:.1f}'.format(meanv) + ' m/s', fontsize=10);

def windOriGrid(m, lon, lat, wuplot, wvplot, maxw):

	q = m.quiver(lon[1:22:4, 40:-1:20], lat[1:22:4, 40:-1:20], wuplot[1:22:4, 40:-1:20], wvplot[1:22:4, 40:-1:20], latlon=True, headwidth=4, headlength=5, width=0.01, scale=maxw*8.5, alpha=.5, color='none', edgecolors=('k'), linewidths=(1));
	q = m.quiver(lon[27:40:10, 35:-1:20], lat[27:40:10, 35:-1:20], wuplot[27:40:10, 35:-1:20], wvplot[27:40:10, 35:-1:20], latlon=True, headwidth=4, headlength=5, width=0.01, scale=maxw*8.5, alpha=.5, color='none', edgecolors=('b'), linewidths=(1));
	q = m.quiver(lon[90:-1:20, 44:-1:22], lat[90:-1:20, 44:-1:22], wuplot[90:-1:20, 44:-1:22], wvplot[90:-1:20, 44:-1:22], latlon=True, headwidth=4, headlength=5, width=0.01, scale=maxw*8.5, alpha=.5, color='none', edgecolors=('r'), linewidths=(1));
	q = m.quiver(lon[120:125:4, ::10], lat[120:125:4, ::10], wuplot[120:125:4, ::10], wvplot[120:125:4, ::10], latlon=True, headwidth=4, headlength=5, width=0.01, scale=maxw*8.5, alpha=.5, color='none', edgecolors=('g'), linewidths=(1));
	plt.quiverkey(q, 0.77, 0.23, maxw, r'max wind speed = $%.1f$ m/s' %maxw, labelpos='W');

def quiverOriGrid(fig, m, lon, lat, uplot, vplot, maxv, meanv, xcur, ycur):

	# q = m.quiver(lon[::4,::4], lat[::4,::4], uplot[::4, ::4], vplot[::4, ::4], latlon=True, headwidth=4, headlength=5, scale=maxv*8.5, alpha=0.5);
	q = m.quiver(lon[0:22:2, ::8], lat[0:22:2, ::8], uplot[0:22:2, ::8], vplot[0:22:2, ::8], latlon=True, headwidth=4, headlength=5, width=0.002, scale=maxv*8.5);#, color='green');
	q = m.quiver(lon[22:40:5, ::8], lat[22:40:5, ::8], uplot[22:40:5, ::8], vplot[22:40:5, ::8], latlon=True, headwidth=4, headlength=5, width=0.002, scale=maxv*8.5);#, color='blue');
	q = m.quiver(lon[50:100:40, 40:-1:4], lat[50:100:40, 40:-1:4], uplot[50:100:40, 40:-1:4], vplot[50:100:40, 40:-1:4], latlon=True, headwidth=4, headlength=5, width=0.002, scale=maxv*8.5);#, color='red');

	q = m.quiver(lon[105:120:6, 50:-1:4], lat[105:120:6, 50:-1:4], uplot[105:120:6, 50:-1:4], vplot[105:120:6, 50:-1:4], latlon=True, headwidth=4, headlength=5, width=0.002, scale=maxv*8.5);#, color='blue');

	q = m.quiver(lon[120:-1:2, ::4], lat[120:-1:2, ::4], uplot[120:-1:2, ::4], vplot[120:-1:2, ::4], latlon=True, headwidth=4, headlength=5, width=0.002, scale=maxv*8.5);#, color='magenta');

	plt.quiverkey(q, 0.77, 0.05, maxv, r'$%.1f$ m/s' %maxv, labelpos='W');
	fig.text(xcur, ycur,'max = ' + '{:.1f}'.format(maxv) + ' m/s', fontsize=10);
	fig.text(xcur, ycur-.04,'mean = ' + '{:.1f}'.format(meanv) + ' m/s', fontsize=10);

