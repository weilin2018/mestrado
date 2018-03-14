%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA vel_eq_10_23.m                             %
% AJUSTE PERFIL DE VELOCIDADE                         %
% ALESSANDRO LUVIZON BÉRGAMO                          %
% LUIZ bRUNER DE MIRANDA                              %
% ALESSANDRO LUVIZON BÉRGAMO                          %
% MARÇO DE 2001/2012                                  % 
% TEORIA Miranda et al. (2002) - (eq. 10.23 p. 345)   %
% ou Miranda et al. (2012, 2nd ed.), p. 353           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  Simulacao de perfis de velocidade estacionários 
%  Foi utilizada a equação de Miranda et al. (2011) (eq. 10.23 p. 345)
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

nr=1; % input('Numero da simulacao:');
nz=1.0e-002; % 5.0e-002 input('Coef. de viscosidade turbulenta vertical m.m/s - Nz:');
tau=0.01 % tensão de cisalhamento do vento em Pa
uf=0.001; % 0.1 velocidade gerada pela descarga fluvial;
h=-10.0; % -10 m input('Prof. local (m) - H:');
rob=1020.0; % input a densidade na boca do estuário;
roc=1000.0; % input a densidade na cabeceira do estuário;
%deltax=1.0e+004; % deltax=1.0e+004 distância longitudinal # "bem misturado";
deltax=1.0e+004; % distância longitudinal # "parcialmente misturado";
g=9.8; % aceleracao da gravidade (m/s2)
ro=1020.0; %densidade da água do mar em kg/m3

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

% Gradiente longitudinal de densidade 

delro=(rob-roc);

rox=(delro./deltax)/ro

% Formulação matemática - Inclinação da superfície livre (etax)
% Equação (10.19) p. 344


cb1=3*(nz*uf)/(g*h^2)

cb2=-(0.375)*h*rox

cb3=(1.5)*(tau)/(ro*g*h)

etax=-(cb1+cb2+cb3)

% Formulação matemática - perfil de velocidade - (u=u(x,Z)) Equação 10.23 # p.345

Z=0.0:-0.1:-1.0

coef1=-[(g*h.^3/nz)]*(rox)

coef2=-(tau*h./(ro*nz))


%Quando o vento é desprezado "uzv"

uzv=(coef1)*(-0.167*Z.^3-0.188*Z.^2+0.0208)-uf*(1.5*Z.^2-1.5)

vmedia=mean(uzv)

%Quando vento não é despresado "uzvv"

uzvv=uzv+coef2*(0.75*Z.^2+1.0*Z+0.25)

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

legend('Barotropic+Baroclinic','Wind',4)
xlabel('u-Component, u(m s^{-1})','fontsize',14);
ylabel('Depth (Z)','fontsize',14);

%Salvando os perfis verticais com e sem vento

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

 
diary off



 


