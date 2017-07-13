#!/usr/bin/env python2.7
#-*-coding:utf-8-*-

""" 
Código para criar um novo esquema de diretórios de uma lista:

					    lista_n/
					_______|_______
				   |       |       |
			     codes  outputs  data
			     	    ___|___
					   |       |
					pickles  *.png

Como usar:
	$ chmod +x criar_diretorio.py

	Selecione o diretório que deseja gerar a nova estrutura (ex: IOCXXX)

	Digite no terminal o nome do diretório base da estrutura (ex: listaX)


"""

import os
from Tkinter import Tk
from tkFileDialog import askopenfilename,askdirectory

Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
ROOT_DIR = os.getcwd()
BASE_DIR = askdirectory(initialdir=ROOT_DIR,title='Selecione o diretório para gerar a estrutura de diretórios:')

NAME_DIR = str(raw_input('Digite o nome do novo diretório a ser criado: '))

LIST_DIR = [NAME_DIR, NAME_DIR+'/codes', NAME_DIR+'/data', NAME_DIR+'/outputs', NAME_DIR+'/outputs/pickles']

os.system('clear')

for folder in LIST_DIR:
	print("Diretório %s criado!\n" % (folder))
	os.mkdir(os.path.join(BASE_DIR, folder))

