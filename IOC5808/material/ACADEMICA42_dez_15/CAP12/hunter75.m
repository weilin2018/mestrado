%  Simulacao de perfis de velocidade estacion�rios u=u(x,Z)
%  Programa��o de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Estu�rios parcialmente misturados 
%  Equil�brio dos componentes barocl�nico com o atrito interno...
%  com as for�antes vento e descarga fluvial e atrito de fundo maximo
%  A profundidade adimensional Z est� orientada positivamente para cima com origem no
%  fundo -1<Z<0
%  Tens�o de cisalhamento do vento tau>0 e tau<0 estu�rio abaixo e acima

clear
clc
close all

%COEFICIENTES, PAR�METROS (geometria, descarga fluvial)

nr=1; %input('N�mero da simulacao:');
Nz=1.60e-02; % 3.6e-03 input('Coef. de viscosidade turbulenta vertical m.m/s - Nz:');
tau=0.5; % 0.5 tens�o de cisalhamento do vento em Pa, 
%k=2.5e-003; %Coeficiente k de Prandle ou Rossiter; 
%uo=0.8; %Amplitude da velocidade gerada pela mar�;
uf=0.02; % 0.02 velocidade gerada pela descarga fluvial;
h=-10.5; %input('Prof. local (m) - H:');
rro=1000.0; %input a densidade na boca do estu�rio;
%roc=1000.0; %input a densidade na cabeceira do estu�rio;
delsal=2.3e-02; % 2.2e-002 dist�ncia longitudinal;
g=9.8;% aceleracao da gravidade (m/s2)
beta=7.0e-004; %densidade da �gua do mar em kg/m3

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

% Gradiente longitudinal de densidade 

%delro=(rob-roc);

%rox=(delro./deltax)/ro

% Formula��o matem�tica - Inclina��o da superf�cie livre (etax)
% Equa��o 11.66 p. 382


%cb1=-0.592*(h*rox);
%cb2=2.212*(k*uo/g*h)*uf;
%cb3=0.631*tau/(ro*g*h);

%etax=cb1+cb2+cb3

% Formula��o matem�tica - perfil longitudinal de velocidade - (u=u(x,Z));
% Equa��o 11.63 - p. 382.

Z=0.0:-0.1:-1.0

coef1=(-1/(48*Nz))*[(g*beta*delsal*h.^3)]

coef2=-[(1.5*uf)]

coef3=(1/4)*[(tau*h)/(rro*Nz)]

% C�lculo da velocidade na superf�cie [u(x,o)]
% Equa��o 11.64 - p. 382

%uxo=0.058*coef1+0.63*uf+0.395*coef2


% Sem a tens�o de cisalhamento do vento (tau=0)

%uzv=(coef1)*(-0.167*Z.^3-0.296*Z.^2+0.058)+uf*(1.106*Z.^2+0.630)

%vmedia=mean(uzv)

% Cisalhamento do vento, descarga fluvial+gradiente barocl�nico

uzvv=(coef1)*(-8.0*Z.^3-9.0*Z.^2+1.0)+coef2*(1.0*Z.^2-1.0)+coef3*(3*Z.^2+4.0*Z+1.0)

%vmedia=mean(uzvv)


% Resultados

x=[0,0]
X=[0,-1]

figure

%f1 = plot(uzv,Z,'k');
%set(f1,'linewidth',3);
%hold
f1 = plot(uzvv,Z,'k');
set(f1,'linewidth',2);

line(x,X)

legend('Barocl�nic+River+Wind',4)
xlabel('Component, u(m s^{-1})','fontsize',14);
ylabel('Depth, Z','fontsize',14);

%Salvando os perfis verticais com e sem vento

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

diary off


 


