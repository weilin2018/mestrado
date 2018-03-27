'''
        INTERPOLACAO DOS DADOS DE GRADE REFINADA PARA MENOS REFINADA

Funcao 1 - check
- usar CDO para recortar nova regiao de dados e converter arquivos
grb2 para nc
    . somente arquivos com grade refinada !!!

Funcao 2:
- carregar grade menos e mais refinada

- griddata

-


após testada, adicionar funcao ao pacote!







# novos codigos:

. realizar download dos dados do banco de dados CFSv2

. recortar para um domínio específico

. converter de grb2 para .nc

. tomar as medias

. plotarem_tes

'''

def interpolar_grade(coarsedGrid,refinedGrid,variable):
    '''
        Funcao para interpolar da grade do CFSv2 (mais refinada) para a grade
        do CSFR (menos refinada).

        input:
            coarseGrid = lista com lons e lats da grade menos refinada
            refineGrid = lista com lons e lats da grade mais  refinada
            variable   = variavel a ser interpolada para nova grade

        output:
            interpolated_variable

        Importante: a coarseGrid deve ser maior que a refineGrid, ou seja,
        a grade menos refinada deve englobar totalmente a grade mais refinada.
        Assim a interpolação não irá gerar valores NaN.
    '''

    try:
        from scipy.interpolate import griddata
    except:
        print('Need to install scipy package to import griddata')



#############################################################################


def make_map(ax):

	m = Basemap(projection='merc', llcrnrlat=-35, urcrnrlat=-15, llcrnrlon=-56, urcrnrlon=-34, resolution='l')

	# m = pickle.load(open("/media/danilo/Danilo/mestrado/ventopcse/rotinas/sudesteBR.pkl", "r"))
	m.ax = ax

	m.drawcoastlines(linewidth=.8)
	m.drawmapboundary()

	# definir meridianos e paralelos para plotar no mapa
	meridians=np.arange(-50,-40,3)
	parallels=np.arange(-30,-20,2)
	# desenhar meridianos e paralelos conforme definido acima
	m.drawparallels(parallels,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')
	m.drawmeridians(meridians,labels=[True,False,False,True],fontsize=13,fontweight='bold',color='gray')


	return m


coarsed     = xr.open_dataset('/media/danilo/Danilo/mestrado/ventopcse/data/CFSR/1992_2011/wnd10m.gdas.200810.grb2.nc')
lon_coarsed = coarsed['lon'].values - 360
lat_coarsed = coarsed['lat'].values

refined = xr.open_dataset('/media/danilo/Danilo/mestrado/ventopcse/data/CFSv2/verao2014/cdas1.20140316.sfluxgrbf.grb2.nc')
lon_refined = refined['lon'].values - 360
lat_refined = refined['lat'].values

# plotar as duas grades para verificar os tamanhos
fig, ax = plt.subplots()

m = make_map(ax)

# coarsed
x,y = np.meshgrid(lon_coarsed,lat_coarsed)
m.plot(x,y,'k',latlon=True)
m.plot(x.T,y.T,'k',latlon=True)

# refined
x,y = np.meshgrid(lon_refined,lat_refined)
m.plot(x,y,'r',latlon=True,alpha=.2)
m.plot(x.T,y.T,'r',latlon=True,alpha=.2)

plt.show()

# devemos remover a primeira e ultima linha da grade refinada, pois elas estão
# para fora da grade grosseira e isso pode gerar problemas na interpolação
lon_refined = lon_refined[1:-1]
lat_refined = lat_refined[1:-1]

# extraindo as variáveis que serão interpoladas
wu = np.squeeze(ncdata['U_GRD_L103'].values)
wu = wu.mean(axis=0)[1:-1,1:-1]
wv = np.squeeze(ncdata['V_GRD_L103'].values)
wv = wv.mean(axis=0)[1:-1,1:-1]

# preparar os dados para interpolacao
lon_r,lat_r = np.meshgrid(lon_refined, lat_refined)
points = np.array([lon_r, lat_r])               # shape da variavel
lon_c,lat_c = np.meshgrid(lon_coarsed, lat_coarsed)
xi     = np.array([lon_c, lat_c])               # shape de interpolacao

# interpolar
wu_interp = griddata((lon_r.flatten(),lat_r.flatten()), wu.flatten(), (lon_c,lat_c), method='cubic')
wv_interp = griddata((lon_r.flatten(),lat_r.flatten()), wv.flatten(), (lon_c,lat_c), method='cubic')

# plotando
spd_interp = np.sqrt(wu_interp**2 + wv_interp**2)

plt.contourf(lon_c,lat_c,spd_interp)
plt.quiver(lon_c,lat_c,wu_interp,wv_interp)

plt.show()
