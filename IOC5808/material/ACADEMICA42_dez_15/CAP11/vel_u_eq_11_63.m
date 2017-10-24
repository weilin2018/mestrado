%  Simulacao de perfis de velocidade estacionários u=u(x,Z)
%  Foi utilizada a equação de Miranda et al. (2011) (eq. 11.63 p. 382)
%  Programação de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Estuários parcialmente misturados ou bem misturados 
%  Equilíbrio dos componentes barotrópico e baroclínico com o atrito interno
%  com as forçantes vento e descarga fluvial e atrito de fundo moderado
%  Forçantes e parcelas que dependem da profundidade (Z) tiveram o sinal trocado
%  em relação à orientação utilizada no livro (Miranda et al., 2011), logo,
%  a profundidade adimensional Z está orientada positivamente para cima com origem na superfície -1<Z<0
%  Tensão de cisalhamento do vento tau>0 e tau<0 estuário abaixo e acima

clear
clc
close all

%COEFICIENTES, PARÂMETROS (geometria, descarga fluvial)

nr=1; %input('Número da simulacao:');
nz=1.5e-003; % 1.5e-002 input('Coef. de viscosidade turbulenta vertical m.m/s - Nz:');
tau=-0.10; %tensão de cisalhamento do vento em Pa, 
k=2.5e-003; %Coeficiente k de Prandle ou Rossiter; 
uo=0.8; %Amplitude da velocidade gerada pela maré;
uf=0.2; % 0.1 velocidade gerada pela descarga fluvial;
h=-10.0; %input('Prof. local (m) - H:');
rob=1025.0; %input a densidade na boca do estuário;
roc=1000.0; %input a densidade na cabeceira do estuário;
deltax=1.5e+003; % 1.5e+004 distância longitudinal;
g=9.8;% aceleracao da gravidade (m/s2)
ro=1020; %densidade da água do mar em kg/m3

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

% Gradiente longitudinal de densidade 

delro=(rob-roc);

rox=(delro./deltax)/ro

% Formulação matemática - Inclinação da superfície livre (etax)
% Equação 11.66 p. 382


cb1=-0.592*(h*rox);
cb2=2.212*(k*uo/g*h)*uf;
cb3=0.631*tau/(ro*g*h);

etax=cb1+cb2+cb3

% Formulação matemática - perfil longitudinal de velocidade - (u=u(x,Z));
% Equação 11.63 - p. 382.

Z=0.0:-0.1:-1.0

coef1=[(g*h.^2)/(k*uo)]*(rox)

coef2=[tau/(ro*k*uo)]

% Cálculo da velocidade na superfície [u(x,o)]
% Equação 11.64 - p. 382

uxo=0.058*coef1+0.63*uf+0.395*coef2


% Sem a tensão de cisalhamento do vento (tau=0)

uzv=(coef1)*(-0.167*Z.^3-0.296*Z.^2+0.058)+uf*(1.106*Z.^2+0.630)

vmedia=mean(uzv)

% Com a tensão de cisalhamento do vento (tau)

uzvv=(coef1)*(-0.167*Z.^3-0.296*Z.^2+0.058)+uf*(1.106*Z.^2+0.630)+(coef2)*(0.316*Z.^2+1.0*Z+0.395)

vmedia=mean(uzvv)


% Resultados

x=[0,0]
X=[0,-1]

figure

f1 = plot(uzv,Z,'k');
set(f1,'linewidth',3);
hold
f2 = plot(uzvv,Z,'k');
set(f2,'linewidth',2);

line(x,X)

legend('Barotrópico+Baroclínico','Com vento estuário acima',4)
xlabel('Componente, u(m s^{-1})','fontsize',14);
ylabel('Profundidade, Z','fontsize',14);

%Salvando os perfis verticais com e sem vento

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

diary off


 


