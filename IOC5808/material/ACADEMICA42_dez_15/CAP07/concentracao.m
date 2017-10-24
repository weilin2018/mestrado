%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA concentracao.m                                      %
% LUIZ bRUNER DE MIRANDA                                       %
% Simula��o da concentra��o de subst�ncias n�o conservativas   %
% Teoria de Stommel (1953) e Fischer et al. (1979)             %
% Maiz detaalhes em Miranda et al. (2002, 2012 2nd edi��o      %
% (eqs. 7.146, 7.147 e 7.148 p. 282-283)                       %
% ESTU�RIO UNIDIMENSIONAL E TIPO BEM MISTURADO                 % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%COEFICIENTES, PAR�METROS (geometria, descarga fluvial)

nr=1; % input('Numero da simulacao:');
co=0.13; % Concentra��o inicial em kg/m3 
Kxc=60.0e+0002; % Coef. de difus�o turbulenta horizotal m.m/s - Kxc');
Qf=100.0; % Descarga fluvial em m3/s;
A=1000.0; % �rea uniforme da sec��o transversal
xw=3000.0 % Posi��o longitudinal do ponto de lan�amento;
k=2.5e-003; % Coeficiente de atrito;
W=20.0; % Transporte inicial do efluente na posi��o xw=3000.0 m
fo=0.65; % Fra��o de �gua doce no estu�rio, na posi��o 
         % de lan�amento do efluente, determinada pela distribui��o
         % estacion�ria longitudinal da salinidade. 

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)

% C�lculo da velocidade gerada pela descarga fluvial uf

uf=Qf/A 

% C�lculo da concentra��o inicial "co"

co=(W/Qf)*fo

% C�lculo da fun��o "psi" 

psi=(4.0*Kxc*k)/(uf)^2

% C�lculo de varia��o longitudinal da concentra��o C=C(x)
% Estu�rio abaixo, no sentido da desembocadura do estu�rio
% Localizada em x=0
% x maior ou igual a -xw

xw=-3000.0
x=0.0:-100.0:-3000.0

coef1=uf/(2*Kxc);
arg=1.-(sqrt(1+psi))*(x-xw)
Cab=co*exp(coef1*arg)

figure

plot(x,Cab,'*')
xlabel('Dist�ncia longitudinal (m)','fontsize',14)
ylabel('Concentra��o (kg/m^3)','fontsize',14)

 gtext('Estu�rio abaixo','fontsize',14)	
 
   
% C�lculo de varia��o longitudinal da concentra��o C=C(x)
% Estu�rio acima, no sentido da ZR do estu�rio
% x menor ou igual a -xw

%xx=-3000.:-100.0:-6000.0
xx=-3000.:-100.0:-6000.0
arg1=1.+(sqrt(1+psi))*(xx-xw)
Cab=co*exp(coef1*arg1)

figure

plot(xx,Cab,'*')
xlabel('Dist�ncia longitudinal (m)','fontsize',14)
ylabel('Concentra��o (kg/m^3)','fontsize',14)

gtext('Estu�rio acima','fontsize',14)
      
print -dbitmap f:\Academica42\cap07\concentracao



 


