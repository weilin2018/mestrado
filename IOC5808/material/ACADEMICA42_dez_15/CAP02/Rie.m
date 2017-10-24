%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA rie.m                                      %
% C�LCULO DO N�MERO DE RICHARDSON ESTUARINO           %
% ALESSANDRO LUVIZON B�RGAMO - JULHO DE 1998/2001     %
% E LUIZ BRUNER DE MIRANDA                            %
% TEORIA Miranda et al. (2002, 2012 2nd edition)      %
% Equa��o 2.36 p. 85                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

clear

% O Programa calcula o n�mero Rie para condi��es de   %
% mar� de quadratura e de siz�gia de uma esta��o fixa %
% utilizando dados da descarga fluvial (Qf)  e da     %
% geometria do estu�rio: largura B e profundidade h   %

diary e:\Academica42\cap02\DiarioRichardson.txt

%OBSERVA��O: A VELOCIDADE GERADA PELA DESCARGA FLUVIAL (uf) 
%FOI SUBSTITUIDA PELA vELOCIDADE RESIDUAL (ua). %
%QUANDO Qf FOR UM DADO DO PROBLEMA, ESSA VELOCIDADE � CALCULADA
%PELA RAZ�O uf=Qf/B.h, com B e h DENOTANDO A LARGURA E A 
%PROFUNDIDADE DA SE��O TRANSVERSAL 

%COEFICIENTES, PAR�METROS (geometria, descarga fluvial)

g=9.8; % Acelera��o da gravidade em (m/s2)
hq=6.49; % Profundidade na quadratura em (m)
hs=7.50; % Profundidade na siz�gia em (m)
Dro=24.0; % Densidade m�dia da coluna de �gua
%B=660; Largura da sec��o transversal
%Qf= 2.0+003; Descarga Fluvial em (m3/s)


% Carrega a densidade na esta��o C
% Quadratura

load roC_quad.dat
ros=roC_quad;

% Carrega a varia��o hor�ria da velocidade (uC_quad) na esta��o C
% e o valor da velocidade residual (uaC_quad) na Quadratura

load uC_quad.dat
load uaC_quad.dat

% C�lculo de Rie na quadratura

uq=(uC_quad)-(uaC_quad);

uf=uaC_quad;

rq=(((ros(1)+ros(end)))/2+sum(ros(2:end-1)))/26;

u2q=uq.^2;

Mu2q=(((u2q(1)+u2q(end)))/2+sum(u2q(2:end-1)))/26;

Urmq_q=Mu2q.^0.5;

Rie_quad=(g*(Dro/rq)*uf*hq)/(Urmq_q^3)

% Mostra e grava o Rie da quadratura 

save e:\Academica42\cap02\Rie_quad.dat Rie_quad -ascii

% Carrega a densidade na esta��o C
% Siz�gia

load roC_siz.dat
ros=roC_siz;

% Carrega a varia��o hor�ria da velocidade (uC_siz) na esta��o C
% e o valor da velocidade residual (uaC_siz) na Siz�gia

load uC_siz.dat
load uaC_siz.dat

% C�lculo de Rie na siz�gia

us=uC_siz-uaC_siz;

uf=uaC_siz;

rs=(((ros(1)+ros(end)))/2+sum(ros(2:end-1)))/26;

u2s=us.^2;

Mu2s=(((u2s(1)+u2s(end)))/2+sum(u2s(2:end-1)))/26;

Urmq_s=Mu2s.^0.5;

Rie_siz=(g*(Dro/rs)*uf*hs)/(Urmq_s^3)

% Mostra e grava o Rie da quadratura 

save e:\Academica42\cap02\Rie_siz.dat Rie_siz -ascii

