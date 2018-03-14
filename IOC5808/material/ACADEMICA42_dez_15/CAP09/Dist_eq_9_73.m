%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA vel_eq_9_73.m                           %
% DISTÂNCIA DE PENETRAÇÃO DA CUNHA SALINA          %
% LUIZ BRUNER DE MIRANDA                           %
% TEORIA Miranda et al. (2002, 2012 2nd)           %
% (eq. 9.73 p. 333 ou p. 341)                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  Simulação da distância de penetração da cunha salina 
%  em um estuário tipo cunha salina 
%  Foi utilizada a equação de Miranda et al. (2011, 2nd) (eq. 9.73 p. 333)
%  Programação de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Equilíbrio dos componentes barotrópico, descarga fluvial, tensão
%  interfacial de atrito e o atrito interno na cunha salina e atrito máximo no fundo"
%  Tensão de cisalhamento do vento foi considerada desprezível

clear
clc
close all

%COEFICIENTES, PARÂMETROS (geometria, descarga fluvial)

nr=1; % input('Numero da simulacao:');
Qf=1040.; % Descarga fluvial em m3/s;
Ho=11.5; % espessura da  camada na cabeceira do estuário (m) - Ho:');
h1=4.0; % espessura da  camada sobrejacente à cunha salina (m) - h1:');
h2=6.0; % espessura da  camada da cunha salina (m) - h2' na posição x);
B=1000.0; % Largura do estuário na posição x=7,8 km da boca do estuário;
hm=7.3; % Profundidade da cunha salina na boca do estuário;
k=2.5e-003; % Coeficiente de atrito;
g=9.8; % Aceleração da gravidade;
ro1=1000.0; % Densidade da camada sobrejacente à interface da cunha salina
ro2=1020.0; % Densidade da cunha salina

Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

% Cálculo da profundidade adimensional Hm

Hm=hm/Ho

% Cálculo da velocidade gerada pela descarga fluvial "u1=uf" na camada
% sobrejacente à cunha salina

uf=Qf/(h1*B)

% Cálculo da gravidade específica "delta"

delta=g*(ro2-ro1)/ro2

% Cálculo do coeficiente "gama" (eq. 9.70 - p. 332)com a notação "gm"

gm=(uf^2)/(g*delta*Ho)

% Cálculo do argumento do logarítmo neperiano "arg"

arg=(3-Hm)/3

% Cálculo da distância de penetração da cunha salina Xc (eqs. 9.72 ou 9.73)

coef1=(2*Ho)/(k*gm)

parc1=((3./2.)*(Hm)^2+(1./4.)*(Hm)^4)

% parc2=8*(3*log(arg)+Hm)
parc2=8*(3*log(arg)+Hm)

Xc=(coef1*(parc1+parc2))/1000.0

% Xc = Distância de penetração da cunha salina em "km"

figure

gtext('Distância de penetração da cunha salina (km)=','fontsize',14)
 
I=num2str(Xc);
texto=['gtext(','I',')'];
eval(texto)

diary off




