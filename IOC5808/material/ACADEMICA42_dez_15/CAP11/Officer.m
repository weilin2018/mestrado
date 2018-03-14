%  Simulacao de perfis de velocidade estacionários u=u(x,Z)
%  Foi utilizada a equação de Miranda et al. (2002) (eq. 11.45 p. 377)
%  ou Miranda et al. (2012) (eq. 11.45 p. 385) e Officer (1977, p. 15)
%  Programação de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Estuários parcialmente misturados ou bem misturados 
%  Equilíbrio dos componentes barotrópico e baroclínico com o atrito interno
%  com as forçantes vento e descarga fluvial e atrito máximo no fundo"
%  Forçantes e parcelas que dependem da profundidade (Z) tiveram o sinal trocado
%  em relação à orientação utilizada no livro (Miranda et al., 2011), logo,
%  a profundidade adimensional Z está orientada positivamente para cima com origem na superfície -1<Z<0
%  Tensão de cisalhamento do vento tau>0 e tau<0 estuário abaixo e acima


clear
clc
close all

%COEFICIENTES, PARÂMETROS (geometria, descarga fluvial)

nr=1; %input('Numero da simulacao:');
nz=1.5e-001; %input('Coef. de viscosidade turbulenta vertical m.m/s - Nz:');
kz=4.0e-002; %input('Coef. de difusão turbulenta vertical m.m/s - Kz:');
SM=35.0; %input salinidade na boca do estuário SM; 
SR=1.0; %input salinidade na zona de maré do rio SR; 
tau=-0.8; %tensão de cisalhamento do vento em Pa, as parcelas dessa
uf=0.1;%velocidade gerada pela descarga fluvial;
h=-8.0; %input('Prof. local (m) - H:');
rob=1025.0;%input a densidade na boca do estuário;
roc=1000.0;%input a densidade na cabeceira do estuário;
deltax=1.5e+004; %distância longitudinal;
g=9.8; % aceleracao da gravidade (m/s2)
ro=1020.0; %densidade da água do mar em kg/m3
tide=2.0; %altura da maré);


Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

% Gradiente longitudinal de densidade e salinidade

delro=(rob-roc);

rox=(delro./deltax)/ro;

delsal=(SM-SR)*1.e-003;

salx=(delsal./deltax);

% Formulação matemática - Inclinação da superfície livre (etax)
% Equação (11.52) p. 379

cb1=-0.375*(h*rox);
cb2=3.*(nz/g*h.^2)*uf;
cb3=1.5*tau/(ro*g*h);

etax=cb1+cb2+cb3

% Cálculo da velocidade na superfície us e no fundo ub(eq. 1.5 e 1.6 Officer)
% e da velocidade pela equação 1.4 de Officer "uof"

%Z orientado positivamente para cima com origem na superfície -1<Z<0

Z=0.0:-0.1:-1.0;

eta=(tide/deltax);

us=((g*h.^2)/(6*nz))*(eta-0.25*rox*h)+uf

ub=-(g*h.^2)/(3*nz)*(eta-0.75*rox*h)+uf

velof=us*(1-9*Z.^2-8*Z.^3)+ub*(3*Z.^2+4*Z.^3)+uf*(12*Z.^2+12*Z.^3)

% Cálculo do perfil vertical de salinidade Officer (1977) eq. 1.7

 coefsal=((h.^2)/(kz))*salx
 cs1=(us*(0.4*Z.^5-0.75*Z.^4-.5*Z.^2-1./12.))
 cs2=(ub*(.20*Z.^5+0.25*Z.^4+1./60.))
 cs3=(uf*(-0.6*Z.^5-0.5*Z.^2+1./15.))
 
 salof=coefsal*(cs1+cs2+cs3)
 

% Formulação matemática - perfil longitudinal de velocidade - (u=u(x,Z));
% Equação 11.45 p. 377 

%Z orientado positivamente para cima com origem na superfície -1<Z<0

Z=0.0:-0.1:-1.0

coef1=-[(g*h.^3/nz)*(rox)]

coef2=[tau*h./(ro*nz)]


%Quando o vento é desprezado "uzv"


uzv=(coef1)*(-0.167*Z.^3-0.188*Z.^2+0.0208)+uf*(-1.5*Z.^2+1.5)

vmedia=mean(uzv)

%Quando vento não é despresado "uzvv"

uzvv=(coef1)*(-0.167*Z.^3-0.188*Z.^2+0.0208)+uf*(-1.5*Z.^2+1.5)-(coef2)*(0.75*Z.^2+Z+0.25)

%uzvv=uzv+coef2*(-0.75*Z.^2-Z-0.25)

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

legend('Barotrópico+Baroclínico','Com vento estuário acima',4)
xlabel('Componente, u(m s^{-1})','fontsize',14);
ylabel('Profundidade, Z','fontsize',14);

%Salvando os perfis verticais com e sem vento

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

 
diary off


 


