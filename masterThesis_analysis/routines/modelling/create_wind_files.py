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
! original script adapted from test case by Rafaela F. Nascimento (LHiCo-IOUSP june/2015) and Carine G. R. Costa (LHiCo-IOUSP july/2015)
! updated version by Carine G. R. Costa (Davidson Lab, november/2016)
! create "synop_wind_sbb" and "synop_hflx_sbb" files from "wind_syn_cptec2sbb" file

      INTEGER K, KI, KJ, NI, NJ, NT
      dimension I(%d,%d), J(%d,%d), TIME(%d)
      dimension U(%d,%d), V(%d,%d), P(%d,%d)
      dimension H(%d,%d)

      NI = %d
      NJ = %d
      NT = %d

      open(20,file='synop_wind',form='unformatted')
      open(30,file='synop_hflx',form='unformatted')
      open(50,file='vento')

      do 175 K=1,NT

!        read wind_syn_cptec2sbb file, time by time
         read(50,*)TIME(K)
         read(50,1)((I(KI,KJ),J(KI,KJ),U(KI,KJ),
     .      V(KI,KJ),P(KI,KJ),KI=1,NI),KJ=1,NJ)
    1    format (2I5,3F10.3)

         do KI=1,NI
           do KJ=1,NJ
             H(KI,KJ)=0
           end do
         end do

!        print data from wind_syn_cptec2sbb on screen
         write(*,*)'TIME =',TIME(K)
!         write(*,2)I(80,17),J(80,17),U(80,17)
!     .             ,V(80,17),P(80,17),H(80,17)

    2    format (2I5,4F10.3)

!        print on files synop_wind and synop_hflx
         write(20)TIME(K)
         write(20)((U(KI,KJ),V(KI,KJ),P(KI,KJ),KI=1,NI),KJ=1,NJ)
         write(30)TIME(K)
         write(30)((H(KI,KJ),KI=1,NI),KJ=1,NJ)

  175 continue

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
             i,j,t))
f.close()

#os.system('gfortran file.f')
#os.system('./a.out') #./a.out gera os arquivos
                     #./ executa qualqeur executavel
                     #a.out executavel do vento do modelo
