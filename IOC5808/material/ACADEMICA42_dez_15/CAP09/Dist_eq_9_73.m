%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA vel_eq_9_73.m                           %
% DIST�NCIA DE PENETRA��O DA CUNHA SALINA          %
% LUIZ BRUNER DE MIRANDA                           %
% TEORIA Miranda et al. (2002, 2012 2nd)           %
% (eq. 9.73 p. 333 ou p. 341)                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  Simula��o da dist�ncia de penetra��o da cunha salina 
%  em um estu�rio tipo cunha salina 
%  Foi utilizada a equa��o de Miranda et al. (2011, 2nd) (eq. 9.73 p. 333)
%  Programa��o de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Equil�brio dos componentes barotr�pico, descarga fluvial, tens�o
%  interfacial de atrito e o atrito interno na cunha salina e atrito m�ximo no fundo"
%  Tens�o de cisalhamento do vento foi considerada desprez�vel

clear
clc
close all

%COEFICIENTES, PAR�METROS (geometria, descarga fluvial)

nr=1; % input('Numero da simulacao:');
Qf=1040.; % Descarga fluvial em m3/s;
Ho=11.5; % espessura da  camada na cabeceira do estu�rio (m) - Ho:');
h1=4.0; % espessura da  camada sobrejacente � cunha salina (m) - h1:');
h2=6.0; % espessura da  camada da cunha salina (m) - h2' na posi��o x);
B=1000.0; % Largura do estu�rio na posi��o x=7,8 km da boca do estu�rio;
hm=7.3; % Profundidade da cunha salina na boca do estu�rio;
k=2.5e-003; % Coeficiente de atrito;
g=9.8; % Acelera��o da gravidade;
ro1=1000.0; % Densidade da camada sobrejacente � interface da cunha salina
ro2=1020.0; % Densidade da cunha salina

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

% C�lculo da profundidade adimensional Hm

Hm=hm/Ho

% C�lculo da velocidade gerada pela descarga fluvial "u1=uf" na camada
% sobrejacente � cunha salina

uf=Qf/(h1*B)

% C�lculo da gravidade espec�fica "delta"

delta=g*(ro2-ro1)/ro2

% C�lculo do coeficiente "gama" (eq. 9.70 - p. 332)com a nota��o "gm"

gm=(uf^2)/(g*delta*Ho)

% C�lculo do argumento do logar�tmo neperiano "arg"

arg=(3-Hm)/3

% C�lculo da dist�ncia de penetra��o da cunha salina Xc (eqs. 9.72 ou 9.73)

coef1=(2*Ho)/(k*gm)

parc1=((3./2.)*(Hm)^2+(1./4.)*(Hm)^4)

% parc2=8*(3*log(arg)+Hm)
parc2=8*(3*log(arg)+Hm)

Xc=(coef1*(parc1+parc2))/1000.0

% Xc = Dist�ncia de penetra��o da cunha salina em "km"

figure

gtext('Dist�ncia de penetra��o da cunha salina (km)=','fontsize',14)
 
I=num2str(Xc);
texto=['gtext(','I',')'];
eval(texto)

diary off




