%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA vel_eq_9_75.m                                    %
% DIST�NCIA RELATIVA DA PENETRA��O DA CUNHA SALINA (x/Xc)   %
% EQUA��O DE FARMER & MORGAN (1953)                         %
% TEORIA Miranda et al. (2002, 2012 2nd)                    %
% (eq. 9.75 p. 333 ou 342)                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Simula��o da dist�ncia relativa de penetra��o da cunha salina 
%  em um estu�rio tipo cunha salina (x/Xc) 
%  Foi utilizada a equa��o de Farmer & Morgan (1953)
%  Programa��o de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Equil�brio dos componentes barotr�pico, descarga fluvial, tens�o
%  interfacial de atrito e o atrito interno na cunha salina e atrito m�ximo no fundo"
%  Tens�o de cisalhamento do vento foi considerada desprez�vel

clear
clc
close all

nr=1; % input('Numero da simulacao:');
Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

%COEFICIENTES, PAR�METROS 

Ho=12.0; % espessura da  camada na cabeceira do estu�rio (m) - Ho:');
hm=12.0; % Profundida da cunha salina na boca do estu�rio

% C�lculo da profundidade adimensional Hm

Hm=hm/Ho

% C�lculo da espessura relativa da cunha salina

H=0:0.1:1.0;

Xrel=((H./Hm).^2).*(3.0-2.0*(H./Hm))

figure

plot(Xrel,H,'k')

xlabel('Dist�ncia relativa (x/Xc)','fontsize',14);
ylabel('Espessura relativa (H/Hm)','fontsize',14);

gtext('De acordo com Farmer & Morgan (1953)','fontsize',14)

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

 
diary off



