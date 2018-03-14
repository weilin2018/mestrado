
# ler model_grid sem terra
# pegar o II, JJ

import numpy as np


# diretorios de trabalho
BASE_DIR = '/media/danilo/Danilo/mestrado/github/artigoTG/'
SAVE_DIR = BASE_DIR+'grade/model_grid_comterra'
GRID_DIR = BASE_DIR+'grade/model_grid_semterra'

# importar model_grid
grade = np.loadtxt(GRID_DIR, skiprows=20)

II = 94 # tamanho real da grade
JJ = 297# tamanho da real da grade


Is = grade[:,0]                 # Is existentes na grade
Js = grade[:,1]                 # Js existentes na grade

cont = 0                        # contador de posicao

for I in np.arange(2,II-1):
    for J in np.arange(2,JJ-1):
