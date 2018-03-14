%  Simulação da inclinação da superfície livre Etaxqf
%  e do perfil de velocidade estacionário (uqf) gerados pela
%  pela descarga fluvial (Qf).
%  Foi utilizada a equação de Miranda et al. (2002) (eq. 8.90 e 8.92, p. 311)
%  ou Miranda et al. (2012, 2nd ed.) (eq. 8.90 e 8.92,  p. 315)
%  Equilíbrio do componente barotrópico e o atrito interno
%  com atrito máximo moderado no fundo. tendo como forçantes a descarga fluvial
%  O eixo Oz está orientado em sentido oposto ao da teoria (-1<Z<0)


clear
clc
close all

%COEFICIENTES, PARÂMETROS (geometria, descarga fluvial)

nr=1; %input('Numero da simulacao:');
nz=2.0e-003; %coeficiente cinemático de viscosidade; 
Ho=10.0; %input('Prof. local (m) - Ho:');
Uo=0.5; %input a amplitude da velocidade da maré
k=2.5e-003 %input coeficiente de atrito
B=500.0; % input largura do estuário (m) (m/s)
Qf=100.0; % input descarga fluvial em m3/s 
g=9.8; % aceleracao da gravidade (m/s2)


Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    


% Formulação matemática - Inclinação da superfície livre (Etaxqf)

 Etaxqf=-0.89*k*Qf*Uo/(g*B*Ho^2)

 
% Formulação matemática - perfil de velocidade - (uqf)
% Profundidade adimensional Z
% O sentido e o sinal do eixo Z 
% estão invertidos em relação a teoria 


 Z=[.0;-.1;-.2;-.3;-.4;-.5;-.6;-.7;-.8,;-.9;-1.0] 
 
 % Troca de sinal do eixo da profundidade adimensional Z
 
 Z1=-Z
 ZZ=1+Z

 % uf= velocidade média na seção transversal e estacionária
 % gerada pela descarga fluvial Qf
 
 uf=Qf/(B*Ho)
 
 
 % Cálculo do perfil uQf=f(Z)
 % Fórmula Alternativa
 
 %coef1=g*Etaxqf/(2*nz)
 
 %coef2=(pi*nz)/(2*k*Ho*Uo)
 
 %uQf=coef1*Ho^2*(1.0*Z.^2-2.0*Z-coef2)
 
 % Fórmula direta com a eq. 8.92 p. 315 (2nd Ed.)
 
 uQf=0.89*uf*(-0.5*Z.^2+1.0*Z+pi/4.)
 
 % Cálculo da velocidade uQf e da velocidade relativa uQf/uf
   
 uRel=uQf/uf
 
 % Salvando os perfis verticais 
 % Incluindo uma linha no zero de velocidade

 x=[0,0]
 X=[0,1]

 figure
 
 f1 = plot(uQf,ZZ,'--');
 set(f1,'linewidth',3);
 legend('Atrito moderado no fundo',2)
 xlabel('Componente longitudinal, uQf (m s^{-1})','fontsize',14);
 ylabel('Profundidade, Z','fontsize',14);

 line(x,X)
 
 salva=['print -dbitmap simulacao',Y,'.bmp'];
 eval(salva)   
 
 hold
 
 figure
 
 f2 = plot(uRel,ZZ,'--');
 set(f2,'linewidth',3);
 legend('Atrito moderado no fundo',2)
 xlabel('Componente longitudinal, uQf/uf ','fontsize',14);
 ylabel('Profundidade, Z','fontsize',14);
 
 
 line(x,X)
 
 salva=['print -dbitmap simulacao',Y,'.bmp'];
 eval(salva)   


 diary off



