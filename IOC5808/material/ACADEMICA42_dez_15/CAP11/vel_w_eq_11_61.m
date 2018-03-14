%  Simulacao de perfis do componente vertical (w) de velocidade estacionários
%  Estuários estratificados. 
%  Equilíbrio dos componentes barotrópico e baroclínico com o atrito interno
%  com as forçantes vento e descarga fluvial e atrito máximo no fundo"
%  (Andutta & Miranda, 2006)
%  Foi utilizada a equação de Miranda et al. (2002) (eq. 11.61 p. 381)
%  Foi trocada a orientação do eixo Oz


clear
clc
close all

nr=1; %input('Numero da simulacao:');
nz=1.5e-001; %input('Coef. de viscosidade turbulenta vertical m.m/s - Nz:');
tau=-0.8;%tensão de cisalhamento do vento em Pa, as parcelas dessa

% forçante que dependem da profundidade devem ter o sinal trocado
% quando a tensão de cisalhamento é estuário abaixo;

uf=0.1;% Velocidade gerada pela descarga fluvial;
uo=0.5;% Velocidade gerada pela maré;
k=2.5e-003; %Coeficiente de Prandle
h=-8.0; %input('Prof. local (m) - H:');
rob=1025.0;%input a densidade na boca do estuário;
roc=1000.0;%input a densidade na cabeceira do estuário;
deltax=1.5e+003;%distância longitudinal;

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

%Z orientado positivamente para cima com origem na superfície -1<Z<0

Z=0.0:-0.1:-1.0

%Z=0:-.1:-1.0

% Parâmetros e cálculo de quantidades físicas

g=9.8;% aceleracao da gravidade (m/s2)

ro=1020;%densidade da água do mar em kg/m3

delro=(rob-roc);

rox=(delro./deltax)/ro

d2rox=rox/deltax

dhdx=h./deltax


% Formulação matemática - perfil da velocidade vertical - (w=w(x,Z));
% Equação 11.61 p. 381 
% Condições de contorno superior e inferior:w(x,0)=w(x,-1)=0

coef1=-[(2*g*h*dhdx)/(k*uo)*(rox)];
coef2=-[(g*h.^3)/(k*uo)]*(1/ro)*d2rox;

wzx=(coef1+coef2)*(-0.0417*Z.^4-0.083*Z.^3+0.041*Z)

wmedia=mean(wzx)


% Resultados

x=[0,0]
X=[0,-1]

figure

f1 = plot(wzx,Z,'k');
set(f1,'linewidth',3);
%hold
%f2 = plot(wzvv,Z,'k');
%set(f2,'linewidth',2);

line(x,X)

legend('Baroclínico+atrito',4)
xlabel('Componente vertical - w(m s^{-1})','fontsize',14);
ylabel('Profundidade, Z','fontsize',14);

%Salvando os perfis verticais com e sem vento

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

 
diary off


 


