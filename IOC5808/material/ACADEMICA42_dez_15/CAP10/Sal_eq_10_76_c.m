%  ESTUÁRIO BEM MISTURADO
%  Simulacao de perfil vertical estacionário de Salinidade: u=u(x,Z)
%  Nz - Coeficiente cinemático de viscosidade turbulento
%  Kz -     "           "         difusão          " 
%  deta="eta" - inclinação estacionária da superfície livre (maré)
%  delS = Gradiente longitudinal médio de salinidade 
%  Etax = Gradiente barotrópico de pressão
%  Com ou sem a influência da tensão de cisalhamento do vento

clear
clc
close all

nr=1; %input('Numero da simulacao:');
nz=0.001; %input('Coef. Cin. de Visc. turbulento m.m/s - Nz:');
kz=0.01; %input('Coef. Cin. de Dif. turbulento m.m/s - Kz:');
Sb=36.0; % input (salinidade na boca do estuário);
Sc=1.0; % input (salinidade na cabeceira do estuário);
S0=30; % Salinidade na superfície
eta=1.0 %input('Inclinação da superfície livre - eta>0 ou eta<0');
Tau=0.0 % input a tensão de cisalhamento do vento; 
dx=10000; %input('Delta x (m) - D:');
ro=1010; % densidade ;

h=10; %input('Prof. local (m) - H:');

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

%Z Orientado positivamente para baixo com origem na superfície -1<Z<0

Z=0:.1:1

g=9.8; % aceleracao da gravidade (m/s2)


% Formulação matemática - Perfil Estacionário de Salinidade.


DelSx=(Sb-Sc)*(.001)/dx;

NK=nz*kz;

etax=eta./(dx);

coef=(DelSx)*(etax);

coef1=[(g*h^(4)*coef)/(2*NK)];

coef2=[(g*h^(4)*coef)/(3*NK)];

coef3=[(DelSx*h^(3)*Tau)/(2*ro*NK)]

Sz=-coef1*(1./12)*(Z.^(4)-(1/2)*Z.^(2))-(coef2+coef3)*Z;

%S=S0+Sz

S=S0-Sz

SM=mean(S)

x=[SM,SM]
X=[0,-1]

figure

f1 = plot(S,-Z,'k');
set(f1,'linewidth',2);

line(x,X)
 
legend('Vertical Salinity Profile',1)
xlabel('Salinity (^o/_o_o)','fontsize',14);
ylabel('Depth Z','fontsize',14);

Y=num2str(nr);
salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

 
diary off


