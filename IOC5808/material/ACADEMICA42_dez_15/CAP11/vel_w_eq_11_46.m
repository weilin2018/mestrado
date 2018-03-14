%  Simulacao de perfis de velocidade estacion�rios w=w(x,Z)
%  Foi utilizada a equa��o de Miranda et al. (2011) (eq. 11.46 p. 377)
%  Programa��o de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Estu�rios parcialmente misturados  
%  Equil�brio dos componentes barotr�pico e barocl�nico com o atrito interno
%  com as for�antes vento e descarga fluvial e atrito m�ximo no fundo"
%  For�antes e parcelas que dependem da profundidade (Z) tiveram o sinal trocado
%  em rela��o � orienta��o utilizada no livro (Miranda et al., 2011), logo,
%  a profundidade adimensional Z est� orientada positivamente para cima com origem na superf�cie -1<Z<0
%  Tens�o de cisalhamento do vento tau>0 e tau<0 estu�rio abaixo e acima


clear
clc
close all

%COEFICIENTES, PAR�METROS (geometria, descarga fluvial)
nr=1; %input('Numero da simulacao:');
nz=1.5e-001; %input('Coef. de viscosidade turbulenta vertical m.m/s - Nz:');
tau=-0.8; %tens�o de cisalhamento do vento em Pa, as parcelas dessa
uf=0.1; %velocidade gerada pela descarga fluvial;
h=-8.0; %input('Prof. local (m) - H:');
rob=1025.0; %input a densidade na boca do estu�rio;
roc=1000.0; %input a densidade na cabeceira do estu�rio;
deltax=1.5e+002; %dist�ncia longitudinal;
B=500; %Foi feita a aproxima��o em que a largura (B=cte.)foi considerada constante
g=9.8; % aceleracao da gravidade (m/s2)
ro=1020; %densidade da �gua do mar em kg/m3

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

% Gradiente longitudinal de densidade 

delro=(rob-roc);

rox=(delro./deltax)/ro


%Z=0:-.1:-1.0

d2rox=rox/deltax

dhdx=h./deltax


% Formula��o matem�tica - perfil da velocidade vertical - (w=w(x,Z));
% Equa��o 11.46 p. 377 

%Z orientado positivamente para cima com origem na superf�cie -1<Z<0

Z=0.0:-0.1:-1.0

coef1=-[(g./nz)*(d2rox)]*h.^4;
coef2=[tau/h./(ro*nz)]*2*dhdx

%Quando o vento � desprezado "uzv"

wzx=coef1*(-0.0417*Z.^4-0.0625*Z.^3+0.0208*Z)

wmedia=mean(wzx)

%Quando vento n�o � despresado "uzvv"

wzvv=wzx+coef2*(0.25*Z.^3+0.5*Z.^2+0.25*Z)

wmedia=mean(wzvv)


% Resultados

x=[0,0]
X=[0,-1]

figure

f1 = plot(wzx,Z,'k');
set(f1,'linewidth',3);
hold
f2 = plot(wzvv,Z,'k');
set(f2,'linewidth',2);

line(x,X)

legend('Barocl�nico','Com vento estu�rio acima',4)
xlabel('Componente vertical, w(m s^{-1})','fontsize',14);
ylabel('Profundidade, Z','fontsize',14);

%Salvando os perfis verticais com e sem vento

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

 
diary off


 


