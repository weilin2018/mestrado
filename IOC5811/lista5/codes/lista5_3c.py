

import scipy.io as sio

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')

modosF = sio.loadmat('/home/tparente/danilo/mestrado/github/IOC5811/lista5/outputs/modos_fratantoni.mat')
modosB = sio.loadmat('/home/tparente/danilo/mestrado/github/IOC5811/lista5/outputs/modos_bruBrown.mat')
# recuperando os dados da solução numérica do matlab
Fp = modosF['Fp']
z = modosF['z']*(-1)
Ht = modosF['Ht']
nl = modosF['nl']

# plotando os dados
fig, ax = plt.subplots(nrows=1,ncols=2)

for i in range(0,6):
    ax[0].plot(Fp[:,i], z[0], label=str(i))

ax[0].set_xlim([-5, 5])
ax[0].set_xlabel(r'Modos ($F_j$)', fontsize=15)

ax[0].set_ylim([Ht, 0])
ax[0].set_ylabel(r'Profundidade ($H_i$) [m]', fontsize=15)
ax[0].set_title('Modelo de Fratantoni et al. (1995)',fontsize=15)

# recuperando os dados da solução numérica do matlab
Fp = modosB['Fp']
z = modosB['z']*(-1)
Ht = modosB['Ht']
nl = modosB['nl']


for i in range(0,6):
    ax[1].plot(Fp[:,i], z[0], label=str(i))

ax[1].set_xlim([-5, 5])
ax[1].set_xlabel(r'Modos ($F_j$)', fontsize=15)

ax[1].set_ylim([Ht, 0])
ax[1].set_ylabel(r'Profundidade ($H_i$) [m]', fontsize=15)
ax[1].set_title('Modelo de Bub & Brown (1996)',fontsize=15)


ax[1].legend(loc='best')

plt.suptitle(u'Ex 3.c - Modos Dinâmicos num Oceano de 6 Camadas', fontsize=20)


plt.show()
