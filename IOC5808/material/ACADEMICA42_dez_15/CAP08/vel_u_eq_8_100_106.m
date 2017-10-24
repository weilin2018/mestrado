%  Simulação da inclinação da superfície livre 
%  e do perfil de velocidade estacionário 
%  gerados pela tensão de cisalhamento do vento
%  Foi utilizada a equação de Miranda et al. (2002) (eq. 8.100 e 8.106, p. 312-314)
%  ou Miranda et al. (2012, 2nd ed.) (eq. 8.100 e 8.106,  p. 316-318)
%  Equilíbrio do componente barotrópico e o atrito interno
%  com atritos máximo e moderado no fundo. tendo como forçantes 
%  o vento longitudinal (TauV) 
%  O eixo Oz está orientado em sentido oposto ao da teoria (-1<Z<0)

clear
clc
close all

%COEFICIENTES, PARÂMETROS (geometria, descarga fluvial)

nr=1; %input('Numero da simulacao:');
nz=2.0e-003; %coeficiente cinemático de viscosidade;
k=2.5e-003 %input coeficiente de atrito
Ho=10.0; %input('Prof. local (m) - Ho:');
B=500.0; % input largura do estuário (m)
g=9.8; % aceleracao da gravidade (m/s2)
TauW=0.05 % input tensão de cisalhamento do vento em Pa
ro=1.020; % input a densidade da água


Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario) 

% Formulação matemática - Inclinação da superfície livre (EtaxW)
% Equação (8.97) p. 316

 EtaxW=(1.5*TauW)/(ro*g*Ho)
 
%Profundidade adimensional Z

Z=0.0:.1:1.0

 %Z=0.0:-0.1:-1.0
 %Z=[.0;-.1;-.2;-.3;-.4;-.5;-.6;-.7;-.8,;-.9;-1.0] 
 
% Z1=-Z
% ZZ=1+Z
 
 
 % Cálculo da velocidade gerada 
 % pela tensão de cisalhamento do vento %
 % com a eq. 8.99 - p. 317 (2nd ed.),%
 
 % Iniciando pelo cálculo da inclinação da superfície livre Etax
 
 Etax=1.5*(TauW)/(ro*g*nz)
  
 % Com atrito de fundo máximo (eq. 8.100, p. 317)
 
 coef1=(TauW*Ho)/(ro*nz)
 
 uTau1=(coef1)*(0.75*Z.^2-0.5*Z)
 
 uRel1=uTau1/coef1
 
 % Salvando o perfil vertical da velocidade 
 % relativa gerada pela tensão do vento
 % Com atrito de fundo máximo (eq. 8.100, p. 317)

 x=[0,0]
 X=[0,1]

 figure
 
 
 f1 = plot(uRel1,Z,'--');
 set(f1,'linewidth',3);
 legend('Vento+Atrito máximo no fundo',4)
 xlabel('Componente longitudinal, uRelativa)','fontsize',14);
 ylabel('Profundidade, Z','fontsize',14);

 line(x,X)
 
 salva=['print -dbitmap simulacao',Y,'.bmp'];
 eval(salva)   
 
 figure
 
 % Formulação matemática - perfil de velocidade - (uTau)
 % gerado pelo vento Equação 8.105 e 8.106 p. 318 
 % com atrito moderado no fundo
 
 % Formulação matemática - Inclinação da superfície livre (EtaxW)
 % Equação (8.104) p. 314
 
 EtaxWM=1.149*TauW/(ro*g*Ho)
 
 coef2=(TauW*Ho)/(ro*nz)
 
 uTau2=coef2*(0.574*Z.^2-0.149*Z-0.117)
 
 uRel2=uTau2/coef2

 % Salvando o perfil vertical da velocidade 
 % relativa gerada pela tensão do vento
 % Com atrito de fundo moderado (eq. 8.106, p. 318)

 
 f2 = plot(uRel2,Z,'--');
 set(f2,'linewidth',3);
 legend('Vento+Atrito moderado no fundo',4)
 xlabel('Componente longitudinal, uRelativa','fontsize',14);
 ylabel('Profundidade, Z','fontsize',14);

 line(x,X)
 
 salva=['print -dbitmap simulacao',Y,'.bmp'];
 eval(salva)   
 