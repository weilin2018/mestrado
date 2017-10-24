%  Simulação de perfis de salinidade estacionários # %  (Andutta & Miranda, 2012)
%  Estuário parcialmente misturado 
%  Foram utilizadas as equações de Miranda et al. (2002) (eqs. 11.88, p. 389)
%  Forçantes e parcelas que dependem da profundidade (Z) tiveram o sinal trocado
%  em relação à orientação utilizada no livro, logo,
%  com Z orientado positivamente para cima com origem na superfície -1<Z<0


clear
clc
close all

%COEFICIENTES, PARÂMETROS (geometria, descarga fluvial)

nr=1; %input('Numero da simulacao:');
nz=1.0e-003; %Coeficiente cinemático de viscosidade;
kz=1.5e+001; %'Coeficiente cinemático de difusão turbulenta m.m/s);
tau=0.02; %Tensão de cisalhamento do vento;
uf=0.1;%velocidade gerada pela descarga fluvial;
Sm=20.0; % Salinidade média;
h=-8.0; %input('Prof. local (m) - H:');
sb=36.0; %Salinidade na boca do estuário;
sc=1.0; %Salinidade na cabeceira do estuário;
rob=1025.0;%input a densidade na boca do estuário;
roc=1000.0;%input a densidade na cabeceira do estuário;
deltax=5.0e+002;%distância longitudinal;
uo=0.5;% Velocidade gerada pela maré;
k=2.5e-003; %Coeficiente de Prandle
g=9.8;% aceleracao da gravidade (m/s2)
ro=1020;% densidade da água do mar em kg/m3

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    


% Gradientes de densidade e de salinidade 

delro=(rob-roc);

rox=(delro./deltax)/ro

delsal=(sb-sc);

salx=(delsal./deltax);

% Formulação matemática - Cálculo de u(x,0) 
% Equação 11.47 p. 378


cb1=-((0.0208*g*h.^3)/(nz))*rox;
cb2=(0.25*tau*h)/(ro*nz);

uxo=cb1+1.0*uf+0.25*cb2

% Formulação matemática - perfil vertical de Salinidade - (S=S(x,Z));
% Com escorregamento de fundo - Equação 11.89 p. 389.

Z=0.0:-0.1:-1.0

coef1=[(g*h.^2)/(k*uo)*(rox)]

coef2=[1./(ro*k*uo)]

sztaub=Sm+(coef1)*uxo*(-0.144*Z.^5-0.425*Z.^4+0.5*Z.^2-0.105)+uf*(0.091*Z.^5+0.36*Z.^4-0.057)

sm=mean(sztaub)

x1=[sm,sm]
x2=[-1,0]

% Resultados

figure

f1 = plot(sztaub,Z,'k');
set(f1,'linewidth',3);
legend('Baroclínico/u-atrito máximo',3)
xlabel('Salinidade - S(^{o}/{o})','fontsize',14);
ylabel('Profundidade, Z','fontsize',14);

line(x1,x2)

%Salvando os perfis verticais

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

 
diary off


 


