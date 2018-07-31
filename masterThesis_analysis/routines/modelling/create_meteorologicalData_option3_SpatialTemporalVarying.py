#!/usr/bin/env python
# -*- coding: utf-8 -*-

#IMPORTANTE LEMBRAR DE ABRIR O PYTHON E COLOCAR NA PASTA CERTA,
#PARA QUE OS ARQUIVOS SEJAM SALVOS NO LOCAL CORRETO
#obs.: pode ser antes de abrir o python

#cria os arquivos de vento para rodar no sECOM

import os

i = 110
j = 137
t = 253 # numero de timesteps do seu arquivo wind.f!!!(ju.py eh o --> t0/6)

g ="""
       INTEGER K, KI, KJ, NI, NJ, NT
       dimension I(%d,%d), J(%d,%d), TIME(%d)
       dimension U(%d,%d), V(%d,%d), BP(%d,%d)
       dimension SW(%d,%d), AT(%d,%d), RH(%d,%d)
       dimension CC(%d,%d), EX(%d,%d), QP(%d,%d), QE(%d,%d)

       NI = %d
       NJ = %d
       NT = %d

       open(20,file='synop_met',form='unformatted')
       open(50,file='metData')

       do 175 K=1,NT
         read(50,*)TIME(K)
         read(50,1)((I(KI,KJ),J(KI,KJ),U(KI,KJ),
     .      V(KI,KJ),SW(KI,KJ),AT(KI,KJ),RH(KI,KJ),BP(KI,KJ),
     .      CC(KI,KJ),EX(KI,KJ),QP(KI,KJ),QE(KI,KJ),KI=1,NI),KJ=1,NJ)
    1    format (2I5,10F10.3)

      write(*,*)'TIME=',TIME(K)

         write(20)TIME(K)
         write(20)((I(KI,KJ),J(KI,KJ),U(KI,KJ),
     .      V(KI,KJ),SW(KI,KJ),AT(KI,KJ),RH(KI,KJ),BP(KI,KJ),
     .      CC(KI,KJ),EX(KI,KJ),QP(KI,KJ),QE(KI,KJ),KI=1,NI),KJ=1,NJ)

  175  continue

       STOP

       END

"""

f = open('file.f','w+')
f.write(g % (i,j,
             i,j,
             t,
             i,j,
             i,j,
             i,j,
             i,j,
             i,j,
             i,j,
             i,j,
             i,j,
             i,j,
             i,j,
             i,j,t))
f.close()

#os.system('gfortran file.f')
#os.system('./a.out') #./a.out gera os arquivos
                     #./ executa qualqeur executavel
                     #a.out executavel do vento do modelo
