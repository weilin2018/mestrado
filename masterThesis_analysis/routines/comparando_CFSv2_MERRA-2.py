
# funcao para facilitar a vida na visualizacao dos dados
def plotar(lajeCut,merrCut):
    ## comparing two datasets
    fig,ax = plt.subplots(nrows=2,ncols=1,sharex=True)

    ax[0].plot(lajeCut.index,lajeCut.wu.values,label='CFSv2')
    ax[0].plot(merrCut.index,merrCut.wu.values,label='MERRA')
    ax[0].set_title('Eastward Component')

    ax[1].plot(lajeCut.index,lajeCut.wv.values,label='CFSv2')
    ax[1].plot(merrCut.index,merrCut.wv.values,label='MERRA')
    ax[1].set_title('Northward Component')

    ax[0].legend(loc='best')
    ax[1].legend(loc='best')

    # plt.suptitle('Laje dos Santos and MERRA nearest location data for DJFM/2014',fontsize=24)
    plt.suptitle('Comparison between CFSv2 and MERRA-2 \n at the nearest point between the two grids',fontsize=24)

    plt.show()

cfsv2 = pickle.load(open('/media/danilo/Danilo/mestrado/github/masterThesis_analysis/routines/pickles/cfsv2_laje.pickle','r'))

merra = pickle.load(open('/media/danilo/Danilo/mestrado/github/masterThesis_analysis/routines/pickles/merraLaje.pickle','r'))
