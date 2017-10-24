%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA vel_eq_9_75.m                                    %
% DISTÂNCIA RELATIVA DA PENETRAÇÃO DA CUNHA SALINA (x/Xc)   %
% EQUAÇÃO DE FARMER & MORGAN (1953)                         %
% TEORIA Miranda et al. (2002, 2012 2nd)                    %
% (eq. 9.75 p. 333 ou 342)                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Simulação da distância relativa de penetração da cunha salina 
%  em um estuário tipo cunha salina (x/Xc) 
%  Foi utilizada a equação de Farmer & Morgan (1953)
%  Programação de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Equilíbrio dos componentes barotrópico, descarga fluvial, tensão
%  interfacial de atrito e o atrito interno na cunha salina e atrito máximo no fundo"
%  Tensão de cisalhamento do vento foi considerada desprezível

clear
clc
close all

nr=1; % input('Numero da simulacao:');
Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

%COEFICIENTES, PARÂMETROS 

Ho=12.0; % espessura da  camada na cabeceira do estuário (m) - Ho:');
hm=12.0; % Profundida da cunha salina na boca do estuário

% Cálculo da profundidade adimensional Hm

Hm=hm/Ho

% Cálculo da espessura relativa da cunha salina

H=0:0.1:1.0;

Xrel=((H./Hm).^2).*(3.0-2.0*(H./Hm))

figure

plot(Xrel,H,'k')

xlabel('Distância relativa (x/Xc)','fontsize',14);
ylabel('Espessura relativa (H/Hm)','fontsize',14);

gtext('De acordo com Farmer & Morgan (1953)','fontsize',14)

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

 
diary off



