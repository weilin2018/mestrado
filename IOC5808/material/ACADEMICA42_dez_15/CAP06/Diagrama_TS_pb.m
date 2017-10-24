%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA Diagrama_TS.m                           %
% PLOTA DADOS DE TEMPERATURA E DE SALINIDADE       % 
% NO CL�SSICO DIAGRAMA TS DE HELLAND HANSEN        %      
% GERANDO O DIAGRAMA DE ESTADO NA FORMA DE IMAGENS % 
% DE PARES DE PONTOS (T,S) ESPALHADOS              %
% PREPARADO POR FERNANDO PINHEIRO ANDUTTA EM 2010  %
% MAIS DETALHES SOBRE APLICA��O PARA MASSAS        % 
% DE �GUA ESTUARINAS EM Miranda et al. (2002, 2012)% 
% p. 233 e p. 237-238                              %
% ATEN��O: O PROGRAMA USA A SUBROTINA "SW_DENS0"   %
% DO MORGAN (1994)- CSIRO MARINE LABORATORIES,     %
% AUSTR�LIA, 222, 28 P. PARA CALCULAR A DENSIDADE  %
% DA �GUA ESTUARINA � PRESS�O HIDROST�TICA ZERO    % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % DEVIDO A R�PIDA VARIA��O NO TEMPO DAS MASSAS DE �GUA ESTUARINAS  %
 % O PROGRAMA FOI PREPARADO PARA PLOTAR O DIAGRAMA COM              %
 % DADOS AMOSTRADOS DURANTE MAR�S DE QUADRATURA E DE SIZ�GIA        %
 % OS ARQUIVOS DESSES DADOS S�O OS SEGUINTES:                       %
 % "quadratura.dat" e "siz�gia.dat"                                 % 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  LEITURA DOS DADOS PARA QUADRATURA E SIZ�GIA

   load quadratura.dat; 
   TQ=quadratura(:,5); 
   SQ=quadratura(:,4);

   load sizigia.dat;
   TS=sizigia(:,5); 
   SS=sizigia(:,4);

% Prepara��o do Diagrama T-S com as linhas isopicnais
% Ou seja, o Diagrama de Estado da massa de �gua estuarina

    Se=17:1:37;
    Te=26:1:30;
        
    [Sg,Tg]=meshgrid(Se,Te,10);
    dens=sw_dens0(Sg,Tg)-1000;
     
     figure
     c = contour(Se,Te,dens,0:2:32);
     colormap('bone');
     %shading flat

     clabel(c,0:1:28);
     
 % Plota os pares de pontos (SQ,TQ) da quadratura e superp�e
 % os pares de pontos (SS,TS) da siz�gia no Diagrama de Estado     
     
     
     hold
     plot(SQ,TQ,'k.','linewidth',1);
     
     plot(SS,TS,'ko','linewidth',1);
     
     hold
     
     title('Diagrama T-S')
     xlabel('Salinidade (psu)','fontsize',14)
     ylabel('Temperatura (^{o}C)','fontsize',14)
     
     
     print -dbitmap f:\academica42\cap06\diagrama_pb


