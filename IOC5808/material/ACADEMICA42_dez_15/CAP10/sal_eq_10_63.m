%  Simulacao de perfis unidimensionais de salinidade [S=S(x)] 
%  Estu�rios parcialmente misturados estacion�rios 
%  Foi utilizada a equa��o de Arons & Stommel(1951) 
%  Detalhes adicionias sobre a teoria em Miranda et al. (2002, 2012 2nd ed.)
%  (eq. 10.63 p. 357 ou p. 365) ou no artigo original
%  Essa solu��o tem como condi��o de contorno na superf�cie fluxo difusivo 
%  de sal igual a "zero" no limite das zonas de mistura (ZM) e do rio (ZR)
%  localizada em (x=0) 
%  Programa��o de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Movimento gerado pelo efeito barotr�pico da mar� e sem atrito e
%  uniforme na profundidade. Foi utilizada a equa��o de conserva��o 
%  de sal simplificada (forma unidimensional), Teoria apresentada
%  em Miranda et al. (2002, 2011 2nd ed.), p. 274 (eq. 7.137). 

clear
clc
close all

%COEFICIENTES, PAR�METROS (geometria, descarga fluvial)

nr=1; %input('Numero da simulacao:');
uf=0.2; % velocidade gerada pela descarga fluvial
So=36.0; % salinidade na boca do estu�rio
h=-10.0; %input('Prof. local (m) - H:');
T=43200.0; % per�odo da mar� em segundos;
L=-1.0+004;% dist�ncia longitudinal da ZM;
etao=0.5; % amplitude da velocidade gerada pela descarga fluvial;
B=3.0e+005; % coeficiente adimensional 

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)


% C�lculo do coeficiente F
% O resultado final da varia��o longitudinal da salinidade
% S=S(x) tem esse coeficiente (F) como par�metro

F=((uf*h.^2*T)/(4*pi*B*etao.^2*L))


% Formula��o matem�tica - salinidade m�dia  - [S=(x)];
% Equa��o 10.63 p. 357 

X=0:0.02:1.0

% C�lculo do coeficiente cinem�tico de difus�o longitudinal Kx
% equa��o 10.57 p. 356 (Miranda et al. 2011)

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

%legend('Barocl�nico','Atrito m�ximo no fundo',4)
xlabel('Dist�ncia adimensional X','fontsize',14);
ylabel('Salinidade relativa, S/So','fontsize',14);

gtext('F=','fontsize',14')

I=num2str(F);
texto=['gtext(','I',')'];
eval(texto)

%Salvando o perfil vertical

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   


diary off


