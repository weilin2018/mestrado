#-*-coding: utf-8-*-

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import os
import string
import xray as xr
from scipy.spatial import cKDTree

import matplotlib
matplotlib.style.use('ggplot')

import sys
sys.path.append('../rotinas/artigoTGpack/')

import artigoTGpack as oceano

# definir o lag, em horas, q pode ser baseado em oceano.max_correlacao()
lag = 3 # em horas a partir do run00

# criar BASE_DIR
BASE_DIR = oceano.make_dir()
DATA_DIR = BASE_DIR + 'artigoTG/rotinas/tides/'
FILE_OUT = DATA_DIR+'lag_'+str(lag)+'h_fromOriginal'
FILE_IN  = DATA_DIR+'tide_complement.txt'

# período das componentes de maré aceitas no sECOM
Tcomponentes = [12., 12.42, 12.66, 23.94, 24.06, 25.86]

# leitura do arquivo com dados de maré (do run_data)
f = open(FILE_IN)
f = f.read()
f = f.split('\n')
data = f[2:]

# inicializar novo arquivo de texto
out = open(FILE_OUT,'w')
out.write(f[0]+'\n')
out.write(f[1]+'\n')

# loop 
for i in np.arange(0,len(data)-1,3):
    loc = data[i]
    com = data[i+1]
    deg = data[i+2]
    
    nova_fase  = []
    componente = 0
    
    for fase in string.split(deg,sep=' '):
        if len(fase)!=0:
            nfase = (lag/Tcomponentes[componente]) * 360
            nova_fase.append(float(nfase))
            componente += 1
            
    # imprimir em novo arquivo os dados modificados
    out.write(loc+'\n')
    out.write(com+'\n')
    out.write('%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n'%(nova_fase[0],nova_fase[1],nova_fase[2],nova_fase[3],nova_fase[4],nova_fase[5]))
    
# fechar arquivo
out.close()

###############################################################################
# troca de arquivos de entrada e saida para facilitar o entendimento
FILE_IN  = FILE_OUT
FILE_RUN = DATA_DIR+'run_data'

# definir a componente
M2 = True
S2 = False

if M2:
	FILE_OUT2 = DATA_DIR+'mare_M2_'+str(lag)+'h'
if S2:
	FILE_OUT2 = DATA_DIR+'mare_S2_'+str(lag)+'h'

out      = open(FILE_OUT2, 'w') # criando arquivo de saída

# ler run_data pulando um determinado numero de linhas até certo ponto
rowStart = 37
rowFinal = 1081

f1       = open(FILE_RUN)
f1       = f1.read()
f1       = f1.split('\n')
data     = f1[rowStart:rowFinal] # somente as componentes de maré

f 		 = open(FILE_IN)
f        = f.read()
f        = f.split('\n')
newData  = f[2:] # somente as componentes de maré

# imprimir as componentes de maré
for i in np.arange(0,len(newData)-1,3):
    loc  = newData[i]
    com  = newData[i+1]
    deg  = newData[i+2]
    
    nova_fase  = []
    componentes= string.split(com, sep=' ') # separando as amplitudes
    
    out.write(loc+'\n') # imprimindo os pontos de contorno
    
    if M2:
        out.write('%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n'%
              (0.0,float(componentes[6]),0.0,0.0,0.0,0.0))
    if S2:
        out.write('%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n'%
              (float(componentes[3]),0.0,0.0,0.0,0.0,0.0)) # imprimindo as amplitudes
     
    out.write(deg+'\n') # imprimindo as fases

# fechar arquivo    
out.close()

##############################################################################
# gerar o run_data final
f1       = open(FILE_RUN)
f1       = f1.read()
f1       = f1.split('\n')

header   = f1[:rowStart]
footer   = f1[rowFinal:]

f2       = open(FILE_OUT2)
f2       = f2.read()
f2       = f2.split('\n')

out 	 = open(DATA_DIR+'run_data'+str(lag)+'hLag','w')

# imprimir o cabeçalho do run_data (copiando do proprio run_data lido)
for line in header:
    out.write(line)
    out.write('\n')
   
# imprimir as componentes de maré
for line in f2[:-1]:
	out.write(line)
	out.write('\n')

out.write(f1[-1])


# imprimir o restante do run_data
for line in footer:
    out.write(line)
    out.write('\n')

# fechar arquivo    
out.close()

# limpando os arquivos gerados
os.system('rm %s %s' % (FILE_OUT2, FILE_OUT))