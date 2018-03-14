%  Simulacao de perfis de velocidade estacion�rios u=u(x,Z)
%  Foi utilizada a equa��o de Miranda et al. (2011) (eq. 11.63 p. 382)
%  Programa��o de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Estu�rios parcialmente misturados ou bem misturados 
%  Equil�brio dos componentes barotr�pico e barocl�nico com o atrito interno
%  com as for�antes vento e descarga fluvial e atrito de fundo moderado
%  For�antes e parcelas que dependem da profundidade (Z) tiveram o sinal trocado
%  em rela��o � orienta��o utilizada no livro (Miranda et al., 2011), logo,
%  a profundidade adimensional Z est� orientada positivamente para cima com origem na superf�cie -1<Z<0
%  Tens�o de cisalhamento do vento tau>0 e tau<0 estu�rio abaixo e acima

clear
clc
close all

%COEFICIENTES, PAR�METROS (geometria, descarga fluvial)

nr=1; %input('N�mero da simulacao:');
nz=1.5e-003; % 1.5e-002 input('Coef. de viscosidade turbulenta vertical m.m/s - Nz:');
tau=-0.10; %tens�o de cisalhamento do vento em Pa, 
k=2.5e-003; %Coeficiente k de Prandle ou Rossiter; 
uo=0.8; %Amplitude da velocidade gerada pela mar�;
uf=0.2; % 0.1 velocidade gerada pela descarga fluvial;
h=-10.0; %input('Prof. local (m) - H:');
rob=1025.0; %input a densidade na boca do estu�rio;
roc=1000.0; %input a densidade na cabeceira do estu�rio;
deltax=1.5e+003; % 1.5e+004 dist�ncia longitudinal;
g=9.8;% aceleracao da gravidade (m/s2)
ro=1020; %densidade da �gua do mar em kg/m3

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

% Gradiente longitudinal de densidade 

delro=(rob-roc);

rox=(delro./deltax)/ro

% Formula��o matem�tica - Inclina��o da superf�cie livre (etax)
% Equa��o 11.66 p. 382


cb1=-0.592*(h*rox);
cb2=2.212*(k*uo/g*h)*uf;
cb3=0.631*tau/(ro*g*h);

etax=cb1+cb2+cb3

% Formula��o matem�tica - perfil longitudinal de velocidade - (u=u(x,Z));
% Equa��o 11.63 - p. 382.

Z=0.0:-0.1:-1.0

coef1=[(g*h.^2)/(k*uo)]*(rox)

coef2=[tau/(ro*k*uo)]

% C�lculo da velocidade na superf�cie [u(x,o)]
% Equa��o 11.64 - p. 382

uxo=0.058*coef1+0.63*uf+0.395*coef2


% Sem a tens�o de cisalhamento do vento (tau=0)

uzv=(coef1)*(-0.167*Z.^3-0.296*Z.^2+0.058)+uf*(1.106*Z.^2+0.630)

vmedia=mean(uzv)

% Com a tens�o de cisalhamento do vento (tau)

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

legend('Barotr�pico+Barocl�nico','Com vento estu�rio acima',4)
xlabel('Componente, u(m s^{-1})','fontsize',14);
ylabel('Profundidade, Z','fontsize',14);

%Salvando os perfis verticais com e sem vento

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

diary off


 


