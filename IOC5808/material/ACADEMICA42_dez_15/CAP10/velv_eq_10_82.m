%  Simulacao de perfis de velocidade transversal (v=v(y,Z) estacionários 
%  Foi utilizada a equação de Miranda et al. (2002) (eq. 10.82 p. 363)ou 
%  Miranda et al. (2012, 2nd ed.) (eq. 10.82 p. 371)
%  Programação de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Estuários parcialmente misturados lateralmente estratificado 
%  Equilíbrio dos componentes barotrópico e baroclínico com o atrito interno
%  com atrito máximo no fundo
%  Forçantes e parcelas que dependem da profundidade (Z) tiveram o sinal trocado
%  em relação à orientação utilizada no livro (Miranda et al., 2012), logo,
%  a profundidade adimensional Z está orientada positivamente para cima com origem na superfície -1<Z<0
%  Tensão de cisalhamento do vento tau>0 e tau<0 estuário abaixo e acima

clear
clc
close all

%COEFICIENTES, PARÂMETROS (geometria, descarga fluvial)

nr=1; %input('Numero da simulacao:');
nz=0.02; %coeficiente cinemático de viscosidade; 
h=-10.0; %input('Prof. local (m) - H:');
roa=1021.5; %input a densidade na boca do estuário;
rob=1020.0; %input a densidade na cabeceira do estuário;
deltay=200; %distância transversal;
g=9.8; % aceleracao da gravidade (m/s2)
ro=1005.0; %densidade da água do mar em kg/m3

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    


% Gradiente de densidade transversal 

delro=(roa-rob);

roy=(delro./deltay)/ro

% Formulação matemática - Inclinação da superfície livre (etax)
% Equação (10.80) p. 362

etay=-0.375*h*roy

% Formulação matemática - perfil de velocidade - (v=v(y,Z));
% Equação 10.82 p. 363 

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

legend('Baroclínico','Atrito máximo no fundo',4)
xlabel('Componente transversal, v(m s^{-1})','fontsize',14);
ylabel('Profundidade, Z','fontsize',14);

%Salvando os perfis verticais com e sem vento

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   



diary off


