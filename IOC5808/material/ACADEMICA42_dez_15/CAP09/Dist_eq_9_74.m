%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA vel_eq_9_74.m                                    %
% DISTÂNCIA RELATIVA DA PENETRAÇÃO DA CUNHA SALINA (x/Xc)   %
% LUIZ BRUNER DE MIRANDA                                    %
% TEORIA Miranda et al. (2002, 2012 2nd)                    %
% (eq. 9.74 p. 333 ou p. 341)                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Simulação da distância relativa de penetração da cunha salina 
%  em um estuário tipo cunha salina (x/Xc) 
%  Foi utilizada a equação de Miranda et al. (2011, 2nd) (eq. 9.74 p. 333)
%  Programação de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Equilíbrio dos componentes barotrópico, descarga fluvial, tensão
%  interfacial de atrito e o atrito interno na cunha salina e atrito máximo no fundo"
%  Tensão de cisalhamento do vento foi considerada desprezível

clear
clc
close all

nr=1; % input('Numero da simulacao:');
Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

%COEFICIENTES, PARÂMETROS (geometria, descarga fluvial)

Ho=12.0; % espessura da  camada na cabeceira do estuário (m) - Ho:');
hm=12.0; % Profundida da cunha salina na boca do estuário

% x/Xc=Xrel; Xrel denota a distância relativa 
% H(x)/Hm=Hrel; Hrel denota a espessura relativa da cunha salina

% Cálculo da profundidade adimensional Hm

Hm=hm/Ho

% Cálculo da espessura relativa da cunha salina

H=0:0.1:1.0

% Cálculo do argumento do logarítmo neperiano "arg"

arg1=(3.0-H)/3.

arg2=(3.0-Hm)/3.

% Cálculo da distância relativa "Xrel" de penetração da cunha salina

coef1=(H./Hm).^2; % Coeficiente da p. 333 - eq.9.74. Observe que no livro
                  % esse coeficiente está com o símbolo de "ao quadrado"
                  % na posição errada. H e Hm ambos devem estar ao
                  % quadrado.

parc1=(3/2)+(H.^2)/4.
parc2=8*(3*log(arg1)+H)./H.^2
parc3=(3/2)+(Hm^2)/4.
parc4=8*(3*log(arg2)+Hm)/(Hm^2)

Xrel=coef1.*[(parc1+parc2)/(parc3+parc4)]

figure

plot(Xrel,H,'*')

xlabel('Distância relativa (x/Xc)','fontsize',14);
ylabel('Espessura relativa (H/Hm)','fontsize',14);

hold
plot (0,0,'*')

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

 
diary off



