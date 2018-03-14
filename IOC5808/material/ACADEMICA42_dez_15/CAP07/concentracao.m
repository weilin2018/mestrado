%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA concentracao.m                                      %
% LUIZ bRUNER DE MIRANDA                                       %
% Simulação da concentração de substâncias não conservativas   %
% Teoria de Stommel (1953) e Fischer et al. (1979)             %
% Maiz detaalhes em Miranda et al. (2002, 2012 2nd edição      %
% (eqs. 7.146, 7.147 e 7.148 p. 282-283)                       %
% ESTUÁRIO UNIDIMENSIONAL E TIPO BEM MISTURADO                 % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%COEFICIENTES, PARÂMETROS (geometria, descarga fluvial)

nr=1; % input('Numero da simulacao:');
co=0.13; % Concentração inicial em kg/m3 
Kxc=60.0e+0002; % Coef. de difusão turbulenta horizotal m.m/s - Kxc');
Qf=100.0; % Descarga fluvial em m3/s;
A=1000.0; % Área uniforme da secção transversal
xw=3000.0 % Posição longitudinal do ponto de lançamento;
k=2.5e-003; % Coeficiente de atrito;
W=20.0; % Transporte inicial do efluente na posição xw=3000.0 m
fo=0.65; % Fração de água doce no estuário, na posição 
         % de lançamento do efluente, determinada pela distribuição
         % estacionária longitudinal da salinidade. 

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)

% Cálculo da velocidade gerada pela descarga fluvial uf

uf=Qf/A 

% Cálculo da concentração inicial "co"

co=(W/Qf)*fo

% Cálculo da função "psi" 

psi=(4.0*Kxc*k)/(uf)^2

% Cálculo de variação longitudinal da concentração C=C(x)
% Estuário abaixo, no sentido da desembocadura do estuário
% Localizada em x=0
% x maior ou igual a -xw

xw=-3000.0
x=0.0:-100.0:-3000.0

coef1=uf/(2*Kxc);
arg=1.-(sqrt(1+psi))*(x-xw)
Cab=co*exp(coef1*arg)

figure

plot(x,Cab,'*')
xlabel('Distância longitudinal (m)','fontsize',14)
ylabel('Concentração (kg/m^3)','fontsize',14)

 gtext('Estuário abaixo','fontsize',14)	
 
   
% Cálculo de variação longitudinal da concentração C=C(x)
% Estuário acima, no sentido da ZR do estuário
% x menor ou igual a -xw

%xx=-3000.:-100.0:-6000.0
xx=-3000.:-100.0:-6000.0
arg1=1.+(sqrt(1+psi))*(xx-xw)
Cab=co*exp(coef1*arg1)

figure

plot(xx,Cab,'*')
xlabel('Distância longitudinal (m)','fontsize',14)
ylabel('Concentração (kg/m^3)','fontsize',14)

gtext('Estuário acima','fontsize',14)
      
print -dbitmap f:\Academica42\cap07\concentracao



 


