%  Simulacao de perfis de velocidade estacion�rios u=u(x,Z)
%  Foi utilizada a equa��o de Miranda et al. (2002) (eq. 11.45 p. 377)
%  ou Miranda et al. (2012) (eq. 11.45 p. 385) e Officer (1977, p. 15)
%  Programa��o de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Estu�rios parcialmente misturados ou bem misturados 
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
kz=4.0e-002; %input('Coef. de difus�o turbulenta vertical m.m/s - Kz:');
SM=35.0; %input salinidade na boca do estu�rio SM; 
SR=1.0; %input salinidade na zona de mar� do rio SR; 
tau=-0.08; %tens�o de cisalhamento do vento em Pa, as parcelas dessa
uf=0.1;%velocidade gerada pela descarga fluvial;
h=8.0; %input(8.'Prof. local (m) - H:');
rob=1025.0;%input a densidade na boca do estu�rio;
roc=1000.0;%input a densidade na cabeceira do estu�rio;
deltax=1.5e+004; %dist�ncia longitudinal;
g=9.8; % aceleracao da gravidade (m/s2)
ro=1020.0; %densidade da �gua do mar em kg/m3
tide=2.0; %altura da mar�;
salm=35.; %salinidade m�dia na coluna;


Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

% Gradiente longitudinal de densidade e salinidade

delro=(rob-roc);

rox=(delro./deltax)/ro;

delsal=(SM-SR)*1.e-003;

salx=(delsal./deltax);

% C�lculo da velocidade na superf�cie us e no fundo ub(eq. 1.5 e 1.6 Officer)
% e da velocidade pela equa��o 1.4 de Officer "uof"

%Z orientado positivamente para baixo (+g) com origem na superf�cie 0<Z<1

% C�lculo do perfil vertical de velocidade equacoes 1.4 e 1.6
% e salinidade (equa��o 1.7) de acordo com Officer (1977)

% Perfil de VELOCIDADE "velof"

Z=0.0:0.1:1.0

eta=(tide/deltax);

us=((g*h.^2)/(6*nz))*(eta-0.25*rox*h)+uf

ub=-(g*h.^2)/(3*nz)*(eta-0.75*rox*h)+uf

velof=us*(1-9*Z.^2+8*Z.^3)+ub*(-3*Z.^2+4*Z.^3)+uf*(12*Z.^2-12*Z.^3)

vmedia=mean(velof)

% Perfil de SALINIDADE "salof"

 coefsal=((h.^2)/(kz))*salx
 cs1=(us*(0.4*Z.^5-0.75*Z.^4+.5*Z.^2-1./12.))
 cs2=(ub*(.20*Z.^5-0.25*Z.^4+1./60.))
 cs3=(uf*(-0.6*Z.^5+1.*Z.^4-0.5*Z.^2+1./15.))
 
 salof=salm+coefsal*(cs1+cs2+cs3)
 
 salmedia=mean(salof)
 

% Resultados

x=[0,0]
X=[0,1]

figure

Z=1.0:-0.1:0.0;

plot(velof,Z,'--')

line(x,X)

print -dbitmap e:\ACADEMICA42_maio_14\CAP11\vel

xlabel('U-Component u(m s^{-1})','fontsize',14);
ylabel('Depth, Z','fontsize',14);



%hold

%f1 = plot(velof,Z,'k');
%set(f1,'linewidth',3);
%hold
%f2 = plot(salof,Z,'k');
%set(f2,'linewidth',2);

%plot(Z,salof,'k')

%print -dbitmap e:\ACADEMICA42_maio_14\CAP11\sal
%xlabel('Salinity - psu','fontsize',14);
%ylabel('Depth, Z','fontsize',14);

%line(x,X)

%legend('Barotr�pico+Barocl�nico',4)
%xlabel('U-Component u(m s^{-1})','fontsize',14);
%ylabel('Depth, Z','fontsize',14);

%Salvando os perfis verticais com e sem vento

%salva=['print -dbitmap simulacao',Y,'.bmp'];
%eval(salva)  

%figure

%plot(xx,Cab,'*')
%xlabel('Dist�ncia longitudinal (m)','fontsize',14)
%ylabel('Concentra��o (kg/m^3)','fontsize',14)

%gtext('Estu�rio acima','fontsize',14)
      
%print -dbitmap f:\Academica42\cap07\concentracao


 
diary off


 


