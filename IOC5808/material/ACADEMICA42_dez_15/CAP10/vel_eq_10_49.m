%  Simulacao de perfis de velocidade estacionários 
%  Foi utilizada a equação de Miranda et al. (2002) (eq. 10.49 - p. 351) ou
%  Miranda et al. 2012 2nd. ed. - p. 359
%  Programação de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Estuários bem misturados 
%  Equilíbrio dos componentes barotrópico e baroclínico com o atrito interno
%  com as forçantes vento e descarga fluvial e atrito no fundo moderado
%  Forçantes e parcelas que dependem da profundidade (Z) tiveram o sinal trocado
%  em relação à orientação utilizada no livro (Miranda et al., 2002, 2012), logo,
%  a profundidade adimensional Z está orientada positivamente para cima com origem na superfície -1<Z<0
%  Tensão de cisalhamento do vento tau>0 e tau<0 estuário abaixo e acima

clear
clc
close all

%COEFICIENTES, PARÂMETROS (geometria, descarga fluvial)
%Para o coeficiente cinemático de viscosidade nz foi utilizada
%a aproximação nz=k.Uo.h, com Uo denotando a amplitude de acordo com
%(Bowden, 1953 e Prandle, 1982)
%da velocidade gerada pela maré

nr=1; %input('Numero da simulacao:');
tau=-0.05; % -0.5 tensão de cisalhamento do vento em Pa, as parcelas dessa
uf=0.1; %velocidade gerada pela descarga fluvial;
uo=0.5; %Amplitude da velocidade da maré; 
k=2.5e-003; % Coeficiente de Bowden e Prandle;
h=-10.0; % -10 input('Prof. local (m) - H:');
rob=1025.0; %input a densidade na boca do estuário;
roc=1000.0; %input a densidade na cabeceira do estuário;
deltax=1.0e+003; %1.0e+003 distância longitudinal;
%deltax=1.0e+002; %1.0e+002 distância longitudinal;
beta=7.0e-004;% Coeficiente de contração salina;
g=9.8; % aceleracao da gravidade (m/s2)
ro=1020.0 %densidade da água do mar em kg/m3

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    


% Gradiente de densidade 

delro=(rob-roc)

rox=(delro./deltax)/ro

% Formulação matemática - Inclinação da superfície livre (etax)
% Equação (10.45) p. 350

cb1=2.212*((k*uo*uf)/(h*g));
cb2=-(0.592)*h*rox;
cb3=(0.631)*tau/(ro*g*h);

etax=-(cb1+cb2+cb3)

% Formulação matemática - perfil de velocidade - (u=u(x,Z));
% Equação 10.49 p. 351 

Z=0.0:-0.1:-1.0

coef1=[(g*h^2/k*uo)*(rox)]

coef2=[tau/(ro*k*uo)]


%Quando o vento é desprezado "uzv"

uzv=(coef1)*(-0.167*Z.^3-0.296*Z.^2+0.058)+1.5*uf*(1.106*Z.^2+0.630)

vmedia=mean(uzv)

%Quando vento não é despresado "uzvv"

uzvv=(coef1)*(-0.167*Z.^3-0.296*Z.^2+0.058)+1.5*uf*(1.106*Z.^2+0.630)+(coef2)*(0.316*Z.^2+Z+0.395)

vmedia=mean(uzvv)


% Resultados

x=[0,0]
X=[0,-1]

figure

f1 = plot(uzv,Z,'k');
set(f1,'linewidth',3);
hold
f2 = plot(uzvv,Z,'k');
set(f2,'linewidth',2);

line(x,X)

legend('Barotrópico+Baroclínico+atrito moderado','Com vento estuário acima',4)
xlabel('Componente, u(m s^{-1})','fontsize',14);
ylabel('Profundidade, Z','fontsize',14);

%Salvando os perfis verticais com e sem vento

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

 
diary off

