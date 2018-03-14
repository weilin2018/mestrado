%  Simulacao de perfis de velocidade transversal (v=v(y,Z) estacion�rios 
%  Foi utilizada a equa��o de Miranda et al. (2002) (eq. 10.82 p. 363)ou 
%  Miranda et al. (2012, 2nd ed.) (eq. 10.82 p. 371)
%  Programa��o de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Estu�rios parcialmente misturados lateralmente estratificado 
%  Equil�brio dos componentes barotr�pico e barocl�nico com o atrito interno
%  com atrito m�ximo no fundo
%  For�antes e parcelas que dependem da profundidade (Z) tiveram o sinal trocado
%  em rela��o � orienta��o utilizada no livro (Miranda et al., 2012), logo,
%  a profundidade adimensional Z est� orientada positivamente para cima com origem na superf�cie -1<Z<0
%  Tens�o de cisalhamento do vento tau>0 e tau<0 estu�rio abaixo e acima

clear
clc
close all

%COEFICIENTES, PAR�METROS (geometria, descarga fluvial)

nr=1; %input('Numero da simulacao:');
nz=0.02; %coeficiente cinem�tico de viscosidade; 
h=-10.0; %input('Prof. local (m) - H:');
roa=1021.5; %input a densidade na boca do estu�rio;
rob=1020.0; %input a densidade na cabeceira do estu�rio;
deltay=200; %dist�ncia transversal;
g=9.8; % aceleracao da gravidade (m/s2)
ro=1005.0; %densidade da �gua do mar em kg/m3

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    


% Gradiente de densidade transversal 

delro=(roa-rob);

roy=(delro./deltay)/ro

% Formula��o matem�tica - Inclina��o da superf�cie livre (etax)
% Equa��o (10.80) p. 362

etay=-0.375*h*roy

% Formula��o matem�tica - perfil de velocidade - (v=v(y,Z));
% Equa��o 10.82 p. 363 

Z=0.0:-0.1:-1.0

coef1=-[(g*h^3)/(nz)]*roy

vzv=(coef1)*(-0.167*Z.^3-0.188*Z.^2+0.021)

vmedia=mean(vzv)

% Resultados

x=[0,0]
X=[0,-1]

figure

f1 = plot(vzv,Z,'k');
set(f1,'linewidth',3);

line(x,X)

legend('Barocl�nico','Atrito m�ximo no fundo',4)
xlabel('Componente transversal, v(m s^{-1})','fontsize',14);
ylabel('Profundidade, Z','fontsize',14);

%Salvando os perfis verticais com e sem vento

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   



diary off


