#-*-coding:utf-8-*-
import numpy as np
import matplotlib.pyplot as plt
import seawater as sw

l = -2500. # metros
y = np.arange(1,l,-1)
H = 14. # metros

v = 2e-2 # corresponde a um vento de 8m/s
c = np.sqrt(9.8*H)

coef = -v**2/c

decai = []
for i in y:
    decai.append(2*i/np.sqrt(np.pi))

decai = np.asarray(decai)

eta = []

for val in decai:
    eta.append(coef - coef*val)

eta = np.asarray(eta)

fig, ax = plt.subplots(figsize=(8,10))
plt.plot(eta,y,label='modelo')
plt.xlabel(u'Elevação $\eta$ em metros', fontsize=20)
plt.ylabel(u'Distância $l$ em metros', fontsize=20)
plt.title(u'Variação do Nível do Mar', fontsize=22)
plt.ylim([-2500, 0])
#plt.xlim([-0.11, 0])
plt.show()



#### exemplo caderno
import numpy as np
import matplotlib.pyplot as plt
import seawater as sw

l = -2500. # metros, considerando minha origem na costa
y = np.arange(1,l,-1)
H = 14. # metros

v = 2e-2 # corresponde a um vento de 8m/s
c = np.sqrt(9.8*H)

f=np.abs(sw.f(-23.))

coef = v**2/(2*f*c)
R = c/np.abs(f)

etaNO = []

for i in y:
    etaNO.append(coef * np.e**(-i/R))

etaNO = np.asarray(etaNO)

plt.plot(etaNO,y,label='ex')
#plt.xlim([0,0.3])
plt.title('Sobre elevacao do nivel do mar na costa')
plt.xlabel(r'Elevacao $\eta$ em metros')
plt.ylabel(r'Distancia $l$ em metros')

plt.ylim([-2500, 0])

plt.legend()
plt.show()
