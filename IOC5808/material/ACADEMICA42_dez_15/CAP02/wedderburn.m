%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA wedderburn.m                               %
% C�LCULO DO N�MERO DE WEDDERBURN                     %
% TEORIA em Miranda et al. (2002, 2012 2nd edition)   %
% Equa��es 2.40 e 2.42 p. 87 e 88                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C�lculo do N�mero de Wedderburn em fun��o da  %
% velocidade do vento                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COEFICIENTES, PAR�METROS (geometria, descarga fluvial)
% S�o dados 
% A velocidade do vento: V
% A diferen�a longitudinal de densidade 
% entre a �gua do mar e da cabeceira do estu�rio: deltaro

%Constantes no SI
%densidade do ar: roar
%Comprimento do canal estuarino: L
%Profundidade: h - L/h "raz�o de aspecto"

function wedderburn(deltaro,V)

 diary e:\ACADEMICA42_maio_14\CAP02\diariowedderburn.txt
  
 V=6.0; % velocidade do vento em m/s
 roar=1.2; % densidade do ar em kg/m3
 g=9.8; % acelera��o da gravidade
 L=10000.0; % dimens�o longitudinal em m
 h=3.0; % profundidade m�dia do estu�rio 
 cd=1.5e-003; %coeficiente de arrasto adimensional
 rob= 1030.0; % densidade na boca do estu�rio
 roc= 1000.0; % densidade na cabeceira do estu�rio
 
% C � o eixo das abcissas (intensidade do vento)

% C�lculo da diferen�a longitudinal de densidade deltaro 

deltaro=rob-roc

% C�lculo do n�mero de Wedderburn para diferentes valores
% da intensidade do vento

C=[0:30];

figure

we=(g*deltaro*h^2)/(L*roar*cd*V^2)

semilogy(C,we,'b-')
            
hold on

V=2.0;

we=(g*deltaro*h^2)/(L*roar*cd*V^2);

semilogy(V,we,'*');

xx=[0,30];
yy=[1,1];
line (xx,yy);

hold on

V=4.0;

we=(g*deltaro*h^2)/(L*roar*cd*V^2);

semilogy(V,we,'*');

V=6.0;

we=(g*deltaro*h^2)/(L*roar*cd*V^2);

semilogy(V,we,'*');

V=8.0;

we=(g*deltaro*h^2)/(L*roar*cd*V^2);

semilogy(V,we,'*');

V=10.0;

we=(g*deltaro*h^2)/(L*roar*cd*V^2);

semilogy(V,we,'*');

V=12.0;

we=(g*deltaro*h^2)/(L*roar*cd*V^2);

semilogy(V,we,'*');

V=14.0;

we=(g*deltaro*h^2)/(L*roar*cd*V^2);

semilogy(V,we,'*');

V=16.0;

we=(g*deltaro*h^2)/(L*roar*cd*V^2);


semilogy(V,we,'*');

xlabel('Intensidade do vento (ms^{-1})','fontsize',14)
ylabel('N�mero de Wedderburn','fontsize',14)
legend('Baroclinicidade x Tens�o do Vento',1)

print -dbitmap e:\ACADEMICA42\CAP02\wedderburn.bmp



