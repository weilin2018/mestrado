
%  Simula��o de perfis de velocidade transversais estacion�rios 
%  Foram utilizadas as equa��es de Miranda et al. (2002) (eq. 10.88 e 10.89 p. 365)
%  ou Miranda et al. (2012, 2nd ed.) (eq. 10.88 e 10.89 p. 373)
%  Programa��o de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Estu�rios parcialmente misturados e lateralmente estratificados 
%  Equil�brio do componente barotr�pico e barocl�nico com o atrito interno
%  com atrito de fundo moderado 
%  For�antes e parcelas que dependem da profundidade (Z) tiveram o sinal trocado
%  em rela��o � orienta��o utilizada no livro (Miranda et al., 2012), logo,
%  a profundidade adimensional Z est� orientada positivamente para cima com origem na superf�cie -1<Z<0

clear
clc
close all

%COEFICIENTES, PAR�METROS (geometria, descarga fluvial)

nr=1; %input('Numero da simulacao:');
nz=4.5e-002; %Amplitude da velocidade da mar�; 
k=2.5e-003; %Coeficiente de Prandle;
h=-10.0; %input('Prof. local (m) - H:');
roa=1021.5; %input a densidade na boca do estu�rio;
rob=1020.0; %input a densidade na cabeceira do estu�rio;
deltay=200; %dist�ncia transversal;
uo=0.5; %Amplitude da velocidade da mar�; 
g=9.8; %aceleracao da gravidade (m/s2)
ro=1020.0; %densidade da �gua do mar em kg/m3

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    


% Gradiente transversal de densidade 

delro=(roa-rob);

roy=(delro./deltay)/ro

% Formula��o matem�tica - Inclina��o da superf�cie livre (etay)
% Equa��o (10.87) p. 365


etay=-(9./16.)*h*roy


% Formula��o matem�tica - perfil de velocidade - (v=v(y,Z));

Z=0.0:-0.1:-1.0

% Com o coeficiente de viscosidade nz

coef1=-[(g*h^3/nz)*(roy)]

vzv=(coef1)*(-0.167*Z.^3-0.281*Z.^2+0.052)

vmedia=mean(vzv)

% Com o coeficiente de viscosidade nz
% com a aproxima��o de Bowden (1953) e Prandle (1982)
 
coef2=[(g*h^2/k*uo)*(roy)]

vzvv=(coef2)*(-0.167*Z.^3-0.281*Z.^2+0.052)

vmedia=mean(vzvv)


% Resultados

x=[0,0]
X=[0,-1]

figure

f1 = plot(vzv,Z,'k');
set(f1,'linewidth',3);
hold
f2 = plot(vzvv,Z,'k');
set(f2,'linewidth',2);

line(x,X)

legend('Barocl�nico + atrito (nz)','baroclinico+atrio (k)',4)
xlabel('Componente transversal, v(m s^{-1})','fontsize',14);
ylabel('Profundidade, Z','fontsize',14);

%Salvando os perfis verticais com e sem vento

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

 
diary off

%

