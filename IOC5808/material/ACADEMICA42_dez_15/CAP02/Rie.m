%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA rie.m                                      %
% CÁLCULO DO NÚMERO DE RICHARDSON ESTUARINO           %
% ALESSANDRO LUVIZON BÉRGAMO - JULHO DE 1998/2001     %
% E LUIZ BRUNER DE MIRANDA                            %
% TEORIA Miranda et al. (2002, 2012 2nd edition)      %
% Equação 2.36 p. 85                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

clear

% O Programa calcula o número Rie para condições de   %
% maré de quadratura e de sizígia de uma estação fixa %
% utilizando dados da descarga fluvial (Qf)  e da     %
% geometria do estuário: largura B e profundidade h   %

diary e:\Academica42\cap02\DiarioRichardson.txt

%OBSERVAÇÃO: A VELOCIDADE GERADA PELA DESCARGA FLUVIAL (uf) 
%FOI SUBSTITUIDA PELA vELOCIDADE RESIDUAL (ua). %
%QUANDO Qf FOR UM DADO DO PROBLEMA, ESSA VELOCIDADE É CALCULADA
%PELA RAZÃO uf=Qf/B.h, com B e h DENOTANDO A LARGURA E A 
%PROFUNDIDADE DA SEÇÃO TRANSVERSAL 

%COEFICIENTES, PARÂMETROS (geometria, descarga fluvial)

g=9.8; % Aceleração da gravidade em (m/s2)
hq=6.49; % Profundidade na quadratura em (m)
hs=7.50; % Profundidade na sizígia em (m)
Dro=24.0; % Densidade média da coluna de água
%B=660; Largura da secção transversal
%Qf= 2.0+003; Descarga Fluvial em (m3/s)


% Carrega a densidade na estação C
% Quadratura

load roC_quad.dat
ros=roC_quad;

% Carrega a variação horária da velocidade (uC_quad) na estação C
% e o valor da velocidade residual (uaC_quad) na Quadratura

load uC_quad.dat
load uaC_quad.dat

% Cálculo de Rie na quadratura

uq=(uC_quad)-(uaC_quad);

uf=uaC_quad;

rq=(((ros(1)+ros(end)))/2+sum(ros(2:end-1)))/26;

u2q=uq.^2;

Mu2q=(((u2q(1)+u2q(end)))/2+sum(u2q(2:end-1)))/26;

Urmq_q=Mu2q.^0.5;

Rie_quad=(g*(Dro/rq)*uf*hq)/(Urmq_q^3)

% Mostra e grava o Rie da quadratura 

save e:\Academica42\cap02\Rie_quad.dat Rie_quad -ascii

% Carrega a densidade na estação C
% Sizígia

load roC_siz.dat
ros=roC_siz;

% Carrega a variação horária da velocidade (uC_siz) na estação C
% e o valor da velocidade residual (uaC_siz) na Sizígia

load uC_siz.dat
load uaC_siz.dat

% Cálculo de Rie na sizígia

us=uC_siz-uaC_siz;

uf=uaC_siz;

rs=(((ros(1)+ros(end)))/2+sum(ros(2:end-1)))/26;

u2s=us.^2;

Mu2s=(((u2s(1)+u2s(end)))/2+sum(u2s(2:end-1)))/26;

Urmq_s=Mu2s.^0.5;

Rie_siz=(g*(Dro/rs)*uf*hs)/(Urmq_s^3)

% Mostra e grava o Rie da quadratura 

save e:\Academica42\cap02\Rie_siz.dat Rie_siz -ascii

