%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA wedderburn.m                               %
% CÁLCULO DO NÚMERO DE WEDDERBURN                     %
% TEORIA em Miranda et al. (2002, 2012 2nd edition)   %
% Equações 2.40 e 2.42 p. 87 e 88                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cálculo do Número de Wedderburn em função da  %
% velocidade do vento                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COEFICIENTES, PARÂMETROS (geometria, descarga fluvial)
% São dados 
% A velocidade do vento: V
% A diferença longitudinal de densidade 
% entre a água do mar e da cabeceira do estuário: deltaro

%Constantes no SI
%densidade do ar: roar
%Comprimento do canal estuarino: L
%Profundidade: h - L/h "razão de aspecto"

function wedderburn(deltaro,V)

 diary e:\ACADEMICA42_maio_14\CAP02\diariowedderburn.txt
  
 V=6.0; % velocidade do vento em m/s
 roar=1.2; % densidade do ar em kg/m3
 g=9.8; % aceleração da gravidade
 L=10000.0; % dimensão longitudinal em m
 h=3.0; % profundidade média do estuário 
 cd=1.5e-003; %coeficiente de arrasto adimensional
 rob= 1030.0; % densidade na boca do estuário
 roc= 1000.0; % densidade na cabeceira do estuário
 
% C é o eixo das abcissas (intensidade do vento)

% Cálculo da diferença longitudinal de densidade deltaro 

deltaro=rob-roc

% Cálculo do número de Wedderburn para diferentes valores
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
ylabel('Número de Wedderburn','fontsize',14)
legend('Baroclinicidade x Tensão do Vento',1)

print -dbitmap e:\ACADEMICA42\CAP02\wedderburn.bmp



