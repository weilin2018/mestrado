%  SIMULAÇÃ0 DE PERFIS VERTICAIS DA VELOCIDADE LONGITUDINAL
%  Equilíbrio do componente barotrópico estacionário com o atrito interno
%  descarga fluvial e tensão do vento (estuário acima (<0) e estuário abaixo (>0).   
%  (Andutta & Miranda, 2009)


clear
clc
close all

% ENTRADA DOS PARÂMETROS

nr=1; %input('Numero da simulacao:');
nz=0.001; %input('Coef. de viscosidade turbulenta vertical m.m/s - Nz:');
deta=-0.002; %input('Inclinação da superfície livre - eta:');
dx=10000; %input('Delta x (m) - D:');
tau=0.005;%tensão de cisalhamento do vento em Pa;
h=10; %input('Prof. local (m) - H:');

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

%Z Orientado positivamente para cima com origem na superfície -1<Z<0

Z=0:.1:1;

g=9.8;                      % aceleracao da gravidade (m/s2)

ro=1020.          %densidade da água do mar em kg/m3


% FORMULAÇÃO MATEMÁTICA - PERFIL DE VELOCIDADE LONGITUDINAL - (U)

Detax=deta./(dx);

coef1=[(h^(2)*g)/(2*nz)]*Detax

coef2=[(h^(2)*g)/(2*nz)]*Detax

coef3=[(tau*h)/(ro*nz)]


uzv=coef1*(Z.^(2)-1.)

uzvv=coef1*(Z.^(2)-1.)+coef3*(-Z+1)


x=[0,0]
X=[0,-1]

% RESULTADOS PERFIL DE VELOCIDADE LONGITUDINAL

figure


f1 = plot(uzv,-Z,'k');
set(f1,'linewidth',2);
hold
f2 = plot(uzvv,-Z,'k')
set(f2,'linewidth',3);
line(x,X)
legend('Barotropic','Gradient Barotrop+wind',3)
xlabel('Component, u(m s^{-1})','fontsize',14);
ylabel('Deoth, Z','fontsize',14);


Y=num2str(nr);
salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

 
%diary off
 diary


