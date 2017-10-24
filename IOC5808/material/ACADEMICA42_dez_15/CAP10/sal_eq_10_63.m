%  Simulacao de perfis unidimensionais de salinidade [S=S(x)] 
%  Estuários parcialmente misturados estacionários 
%  Foi utilizada a equação de Arons & Stommel(1951) 
%  Detalhes adicionias sobre a teoria em Miranda et al. (2002, 2012 2nd ed.)
%  (eq. 10.63 p. 357 ou p. 365) ou no artigo original
%  Essa solução tem como condição de contorno na superfície fluxo difusivo 
%  de sal igual a "zero" no limite das zonas de mistura (ZM) e do rio (ZR)
%  localizada em (x=0) 
%  Programação de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Movimento gerado pelo efeito barotrópico da maré e sem atrito e
%  uniforme na profundidade. Foi utilizada a equação de conservação 
%  de sal simplificada (forma unidimensional), Teoria apresentada
%  em Miranda et al. (2002, 2011 2nd ed.), p. 274 (eq. 7.137). 

clear
clc
close all

%COEFICIENTES, PARÂMETROS (geometria, descarga fluvial)

nr=1; %input('Numero da simulacao:');
uf=0.2; % velocidade gerada pela descarga fluvial
So=36.0; % salinidade na boca do estuário
h=-10.0; %input('Prof. local (m) - H:');
T=43200.0; % período da maré em segundos;
L=-1.0+004;% distância longitudinal da ZM;
etao=0.5; % amplitude da velocidade gerada pela descarga fluvial;
B=3.0e+005; % coeficiente adimensional 

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)


% Cálculo do coeficiente F
% O resultado final da variação longitudinal da salinidade
% S=S(x) tem esse coeficiente (F) como parâmetro

F=((uf*h.^2*T)/(4*pi*B*etao.^2*L))


% Formulação matemática - salinidade média  - [S=(x)];
% Equação 10.63 p. 357 

X=0:0.02:1.0

% Cálculo do coeficiente cinemático de difusão longitudinal Kx
% equação 10.57 p. 356 (Miranda et al. 2011)

ome=2*pi/T;

Kx=((2*B*etao.^2*2*ome*X.^2)/h.^2)*L.^2

XX=(1-1./X)

S=So*exp(F*XX)

%(S/So)=SR

SR=exp(F*XX)


% Resultados

figure

%f1 = plot(SR,X,'k');
f1 = plot(X,SR,'k');
set(f1,'linewidth',3);

%legend('Baroclínico','Atrito máximo no fundo',4)
xlabel('Distância adimensional X','fontsize',14);
ylabel('Salinidade relativa, S/So','fontsize',14);

gtext('F=','fontsize',14')

I=num2str(F);
texto=['gtext(','I',')'];
eval(texto)

%Salvando o perfil vertical

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   


diary off


