%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA vel_eq_9_74.m                                    %
% DIST�NCIA RELATIVA DA PENETRA��O DA CUNHA SALINA (x/Xc)   %
% LUIZ BRUNER DE MIRANDA                                    %
% TEORIA Miranda et al. (2002, 2012 2nd)                    %
% (eq. 9.74 p. 333 ou p. 341)                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Simula��o da dist�ncia relativa de penetra��o da cunha salina 
%  em um estu�rio tipo cunha salina (x/Xc) 
%  Foi utilizada a equa��o de Miranda et al. (2011, 2nd) (eq. 9.74 p. 333)
%  Programa��o de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Equil�brio dos componentes barotr�pico, descarga fluvial, tens�o
%  interfacial de atrito e o atrito interno na cunha salina e atrito m�ximo no fundo"
%  Tens�o de cisalhamento do vento foi considerada desprez�vel

clear
clc
close all

nr=1; % input('Numero da simulacao:');
Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

%COEFICIENTES, PAR�METROS (geometria, descarga fluvial)

Ho=12.0; % espessura da  camada na cabeceira do estu�rio (m) - Ho:');
hm=12.0; % Profundida da cunha salina na boca do estu�rio

% x/Xc=Xrel; Xrel denota a dist�ncia relativa 
% H(x)/Hm=Hrel; Hrel denota a espessura relativa da cunha salina

% C�lculo da profundidade adimensional Hm

Hm=hm/Ho

% C�lculo da espessura relativa da cunha salina

H=0:0.1:1.0

% C�lculo do argumento do logar�tmo neperiano "arg"

arg1=(3.0-H)/3.

arg2=(3.0-Hm)/3.

% C�lculo da dist�ncia relativa "Xrel" de penetra��o da cunha salina

coef1=(H./Hm).^2; % Coeficiente da p. 333 - eq.9.74. Observe que no livro
                  % esse coeficiente est� com o s�mbolo de "ao quadrado"
                  % na posi��o errada. H e Hm ambos devem estar ao
                  % quadrado.

parc1=(3/2)+(H.^2)/4.
parc2=8*(3*log(arg1)+H)./H.^2
parc3=(3/2)+(Hm^2)/4.
parc4=8*(3*log(arg2)+Hm)/(Hm^2)

Xrel=coef1.*[(parc1+parc2)/(parc3+parc4)]

figure

plot(Xrel,H,'*')

xlabel('Dist�ncia relativa (x/Xc)','fontsize',14);
ylabel('Espessura relativa (H/Hm)','fontsize',14);

hold
plot (0,0,'*')

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

 
diary off



