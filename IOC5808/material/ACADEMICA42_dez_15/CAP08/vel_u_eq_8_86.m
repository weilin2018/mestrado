%  Simulação da inclinação da superfície livre Etaxqf
%  e do perfil de velocidade estacionário (uqf) gerados pela
%  pela descarga fluvial (Qf).
%  Foi utilizada a equação de Miranda et al. (2002) (eq. 8.85 e 8.86, p. 308-309)
%  ou Miranda et al. (2012, 2nd ed.) (eq. 8.85 e 8.86,  p. 312-313)
%  Equilíbrio do componente barotrópico e o atrito interno
%  com atrito máximo no fundo. tendo como forçantes a descarga fluvial
%  O eixo Oz está orienado em sentido oposto ao da teoria (-1<Z<0)

clear
clc
close all

%COEFICIENTES, PARÂMETROS (geometria, descarga fluvial)

nr=1; %input('Numero da simulacao:');
nz=2.0e-003; %coeficiente cinemático de viscosidade; 
Ho=10.0; %input('Prof. local (m) - Ho:');
B=500,0; % input largura do estuário (m)
Qf=100.0; % input descarga fluvial em m3/s 
g=9.8; % aceleracao da gravidade (m/s2)


Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    


% Formulação matemática - Inclinação da superfície livre (Etaxqf)
% Equação (8.85) p. 312

 Etaxqf=-3.0*nz*Qf/(g*B*Ho^3)

 
% Formulação matemática - perfil de velocidade - (uqf)
% gerado pela descarga fluvial Equação 8.85 p. 312 

% Profundidade adimensional Z
% O sentido e o sinal do eixo Z 
% estão invertidos em relação ao livro 

 %Z=0.0:-0.1:-1.0
 Z=[.0;-.1;-.2;-.3;-.4;-.5;-.6;-.7;-.8,;-.9;-1.0] 
 
 Z1=-Z
 ZZ=1+Z
 
 uf=Qf/(B*Ho)
 
 %uqf=uf*[(-3/2)*Z.^2+3*Z.]
 
 uQf=-uf*(1.5*Z.^2+3.*Z)
 
 uRel=uQf/uf

% %Salvando os perfis verticais sem tensão do vento

 %x=[0,0]
 %X=[0,-1]
 
 %x=[0,0]
 %X=[uf,-1]

 figure
 z=flipud(Z)
 f1 = plot(uQf,z,'--');
 set(f1,'linewidth',3);
 legend('Barotrópico','Atrito máximo no fundo',4)
 xlabel('Componente longitudinal, uQf (m s^{-1})','fontsize',14);
 ylabel('Profundidade, Z','fontsize',14);

 %line(x,X)
 
 salva=['print -dbitmap simulacao',Y,'.bmp'];
 eval(salva)   
 
 figure
 
 f2 = plot(uRel,z,'--');
 set(f2,'linewidth',3);
  legend('Barotrópico-Atrito máximo no fundo',4)
 xlabel('Componente longitudinal, uQf/uf )','fontsize',14);
 ylabel('Profundidade, Z','fontsize',14);

 %line(x,X)
 
 salva=['print -dbitmap simulacao',Y,'.bmp'];
 eval(salva)   


 diary off

