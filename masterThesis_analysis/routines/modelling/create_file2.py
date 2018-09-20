    #!/usr/bin/env python
    # -*- coding: utf-8 -*-

    #IMPORTANTE LEMBRAR DE ABRIR O PYTHON E COLOCAR NA PASTA CERTA,
    #PARA QUE OS ARQUIVOS SEJAM SALVOS NO LOCAL CORRETO
    #obs.: pode ser antes de abrir o python

    #cria os arquivos de vento para rodar no sECOM

    import os

    i = 110
    j = 137
    """
    t é a quantidade de timesteps do arquivo file.f

    Atenção aqui:
        . se você estiver para rodar uma simulação HOT START, então:
            t0 em create_wind_files.py deve conter a quantidade de horas
            já simuladas no aquecimento e, para calcular t desta rotina,
            deve-se subtrair essa quantidade de horas e então dividir pela
            frequência de input dos dados (ex: no caso de vento são 6 horas)

        . se você estiver para rodar uma simulação COLD START, então:
            t0 = 0
            t = t0/frequencia dos dados (vento = 6 horas)
    """
    t = 245 # numero de timesteps do seu arquivo wind.f!!!

    g ="""
    ! original script adapted from test case by Rafaela F. Nascimento (LHiCo-IOUSP june/2015) and Carine G. R. Costa (LHiCo-IOUSP july/2015)
    ! updated version by Carine G. R. Costa (Davidson Lab, november/2016)
    ! create "synop_wind_sbb" and "synop_hflx_sbb" files from "wind_syn_cptec2sbb" file

          INTEGER K, KI, KJ, NI, NJ, NT
          dimension I(%d,%d), J(%d,%d), TIME(%d)
          dimension I2(%d,%d), J2(%d,%d), TIME2(%d)
          dimension U(%d,%d), V(%d,%d), P(%d,%d)
          dimension H(%d,%d)

          NI = %d
          NJ = %d
          NT = %d

          open(20,file='synop_wind',form='unformatted')
          open(30,file='synop_hflx',form='unformatted')
          open(50,file='vento')
          open(60,file='calor')

          do 175 K=1,NT

    !        read wind_syn_cptec2sbb file, time by time
             read(50,*)TIME(K)
             read(50,1)((I(KI,KJ),J(KI,KJ),U(KI,KJ),
         .      V(KI,KJ),P(KI,KJ),KI=1,NI),KJ=1,NJ)
        1    format (2I5,3F10.3)

    !        read heat_syn_cptec2sbb file, time by time
             read(60,*)TIME2(K)
             read(60,2)((I2(KI,KJ),J2(KI,KJ),H(KI,KJ),KI=1,NI),KJ=1,NJ)
        2    format (2I5,1F10.3)

    !        do KI=1,NI
    !          do KJ=1,NJ
    !            H(KI,KJ)=0
    !          end do
    !        end do

    !        print data from wind_syn_cptec2sbb on screen
             write(*,*)'TIME =',TIME(K)
    !         write(*,2)I(80,17),J(80,17),U(80,17)
    !     .             ,V(80,17),P(80,17),H(80,17)

    !    2    format (2I5,4F10.3)

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
