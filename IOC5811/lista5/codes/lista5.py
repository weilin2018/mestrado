"""

Objetivo:

    Calcular os modos dinâmicos num oceano de 6 camadas

    item A:
        calculo dos raios de deformação baroclínico pelo modelo
        de Fratantoni et al. (1995) e Bub & Brown (1996)

    item B:
        encontrar os modos discretos

    item C:
        descritivo

Original matlab's code:

clear all; format shortEng
% constantes 'ambiente'
f0 = 2*(2*pi/86400)*sind(5);
g = 9.8;

% Modelo de Fratantoni et al. (1995)
Hf = [80, 170, 175, 250, 325, 3000];
Sf = [24.97, 26.3, 26.83, 27.12, 27.32, 27.77];
epsf = [Inf, (Sf(2:end)-Sf(1:end-1))/1020];
Hmf = (Hf(1:end-1).*Hf(2:end))./(Hf(1:end-1)+Hf(2:end));
glf = g.*epsf;
glf(1) = g;
Rdf = (1/f0 .* [sqrt(glf(1).*sum(Hf)), sqrt(glf(2:end).*Hmf)]).*(10^-3)

% Modelo de Bub & Brown, (1996)
Hb = [150, 440, 240, 445, 225, 2500];
Sb = [24.13, 26.97, 27.28, 27.48, 27.74, 27.87];
epsb = [Inf, (Sb(2:end)-Sb(1:end-1))/1020];
Hmb = (Hb(1:end-1).*Hb(2:end))./(Hb(1:end-1)+Hb(2:end));
glb = g.*epsb;
glb(1) = g;
Rdb = (1/f0 .* [sqrt(glb(1).*sum(Hb)), sqrt(glb(2:end).*Hmb)]).*(10^-3)

format short

"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')

def calcula_epsion(Si, rho_0=1020):
    """
        calculo de epsilon, como eps = (sigma[i+1] - sigma[i])/rho_0
    """
    epsf = [(Si[1:] - Si[:-1])/1020]

    eps = np.zeros(6)
    eps[0]  = np.inf

    for i in range(1,len(Si)):
        eps[i] = epsf[0][i-1]

    return eps

def calcula_ortonormalidade(Hi):

    Hm = np.zeros(len(Hi)-1)                # i - 1

    Hm = (Hi[:-1] * Hi[1:]) / (Hi[:-1] + Hi[1:])

    return Hm

def calcula_Rdi(Si,Hi, lat=5):
    """
    Calcular radio de deformação baroclinico
    """

    f0 = 2*(2*np.pi/86400)*np.sin(5)            # plano Beta centrado em 5N
    g  = 9.8                                    # aceleração da gravidade
    # calcular epsilon
    eps = calcula_epsion(Si)
    # calcular condição de ortonormalidade
    Hm = calcula_ortonormalidade(Hi)

    # criando vetor de raios de deformação
    Rdi = np.zeros(6)
    # calculo do raio da camada 1
    Rdi[0] = 1/f0 * np.sqrt(g * Hm.sum())
    # calculando o Rdi da camada 2-6
    R = [ np.sqrt(gHi[1:] * Hm) ]

    for i in range(1,len(Si)):
        Rdi[i] = R[0][i-1]
    #
    # # calcular termo g*He
    # gHi = np.zeros(6)
    # gHi[0] = g
    # gHi[1:] = g * eps[1:]
    #
    # # calcular para a primeira camada (epsilon = inf)
    # Rdi[0] = (1/f0) * np.sqrt(gHi[0] * Hm.sum())
    # # calcular para camadas 1:i
    # for i in range(1,len(Hi)):
    #     Rdi[i] = (1/f0) * (np.sqrt( gHi[i] * Hm.sum() ))

    return Rdi






#####################################################################
#               MODELO DE FRATANTONI ET AL (1995)                   #
#####################################################################
Hi = np.asarray([80, 170, 175, 250, 325, 3000])                 # profundidade equivalente
Si = np.asarray([24.97, 26.30, 26.83, 27.12, 27.32, 27.77])     # densidade

Rdi = calcula_Rdi(Si, Hi)

#####################################################################
#                 MODELO DE BUB & BROWN (1996)                      #
#####################################################################
