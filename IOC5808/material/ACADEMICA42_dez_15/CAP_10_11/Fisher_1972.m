%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programa para aplica��o de modelos bidimensionais            %
% estacion�rios com c�lculo dos componente longitudinal (u) e  %
% vertical (w) de velocidade e compara��o com dados            %
% experimentais.                                               %
% A simula��o considera atrito m�ximo e moderado               %
% (com escorregamento de fundo - "u" de fundo � dif. de 0)     %
% O modelo utilizado foi proposto por Fischer et al. (1972).   %
% Mais detalhes da teoria em Miranda et al. 2002 - Cap. 10     %
% p. 340 ou Miranda et al. 2012 - Cap. 10 - p. 348             %  
% Programa��o: Marcos Eduardo Cordeiro Bernardes - 2001        % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nota��o utilizada:                                             %
% g: Acelera��o da gravidade (m/s2)                                   %
% h: Profundidade m�dia no tempo de amostragem (m)               %
%    (j� corrigido o efeito da mar� e da topog. de fundo)        %
% B: Largura do canal (m)                                        %
% Rho: Densidade utilizada no calculo de Rhox (5 kg/m3)          %
%  ESTE VALOR DE Rho REPRESENTA A VARIACAO DE DENSIDADE          %
%  OBSERVADA POR MIRANDA (1996) PARA BERTIOGA                    %
% x: Comprimento do canal (m)                                    %
% Rho0: Densidade de refer�ncia (kg/m3)                          %
% Beta: Coeficiente de contra��o salina (1/partes por mil)       %
% Rhox: Gradiente longitudinal de densidade (kg/m4)              %
% Rhox2: Derivada segunda do grad. long. de densidade (kg/m4)    %
% TCV: Tens�o de cisalhamento do vento (N/m2)                    %
% TCF: Tensao de cisalhamento do fundo (N/m2)                    %
% uf: Velocidade da descarga de �gua fluvial (m/s)               %
% u0: Velocidade de superf�cie (m/s)                             %
% Ro: Densidade da coluna de �gua (kg/m3)                        %
%     (calculada a partir da equa��o de estado simplificada)     %
% Az: Coeficiente de viscosidade turbulenta (m2/s)               %
% UD: Velocidade te�rica gerada pelo gradiente de densidade (m/s)%
% UR: Velocidade te�rica gerada pela descarga de �gua doce (m/s) %
% UV: Velocidade te�rica gerada pelo vento (m/s)                 %
% UF: Velocidade te�rica gerada pelo atrito com o fundo (m/s)    %
% U: Velocidade resultante em u - comp. longitudinal (m/s)       %
% a: Amplitude m�dia da mar� (m)                                 %
% k1: Constante adimensional para c�lculo do coef. viscosidade    %
%    turbulenta                                                  %
% C0: Celeridade da onda (m/s)                                   %
% WD: Velocidade te�rica gerada pelo gradiente de densidade (m/s)%
% WV: Velocidade te�rica gerada pelo vento (m/s)                 %
% W: Velocidade resultante em w - comp. vertical (m/s)           %
% Bh2: Prod. da larg. do canal pela prof. quadrada  do mesmo (m2)%
% Bh4: Prod. da larg. do canal pela prof. a quarta do mesmo (m4) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diary cananeia.txt

g=9.87
h=10.6
Rho0=1000;
Beta=0.0007;
B=1700
x=6000
Sref=3.7
Sx=Sref*0.001/x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sx e multiplicado por 0.001 para tornar a salinidade%
% numa grandeza adimensional (g/g)                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rho=2.6
Rhox=Rho/x
TCV=0
k1=0.0003
a=1.022
Kz=0.0000003
t5=['Kz= ',num2str(Kz),' m^{2} s^{-1}'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Carregando arquivo-fonte com os dados experimenais %
% u=velocidade - s=salinidade e z= profundidade      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Cananeia.dat
z=Cananeia(:,1);
z1=-z;
z4=1+z;
u=Cananeia(:,2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Devido ao fato das esta��es realizadas no Mar de Cananeia       %
% apresentarem uma orienta��o de eixos diferente do convencional  %
% (no Mar de Cananeia, valores positivos de veloc. representariam %
% correntes estu�rio acima e nao abaixo, como comumente conven-   %
% cionado), foi feita uma padronizacao. Assim, os valores de      %
% veloc. para tais situacoes foram invertidos (u=-u3). Portanto,  %
% os gradientes de salinidade e de densidade serao positivos.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uf=mean(u);
u0=u(1);
s=Cananeia(:,3);
usz4=[u s z4]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equa��o simplificada do estado linear %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ro=1000*(1+(Beta.*s));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C�lculo da amplitude da velocidade gerada           %
% pela mar� U0, com base na celeridade da onda        %
% assumindo que a propaga��o da mar� � progressiva-   %
% e sem atrito Miranda et al. 2012  - eq. 2.20 -      %
% p. 76. A amplitude da velocidade gerada pela mar�   %
% � calculada em fun��o da profundidade h             %
% eq. 2.23 - p. 77                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C0=sqrt(g*h);
U0=(a*C0)/h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C�lculo do coeficiente de viscosidade turbulenta Az  %
% e seu valor m�dio                                    %
% Miranda et al. 2002 - equa��o 10.29 p. 347           %
% ou Miranda et al. 2012 - equa��o 10.29 p. 355        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Az=(1/(2*u0-3*uf))*(((0.0416*g*h^3*Rhox)/Rho0)+((0.5*TCV*h)./Ro));
Az1=mean(Az)
t1=['Az= ',num2str(Az1), ' m^{2} s^{-1}'];
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C�lculo da tens�o de cisalhamento do fundo (TCF)      %
% em fun��o do gradiente longitudinal de densidade      %
% da amplitude da velocidade gerada pela mar� Uo, da    %
% descarga fluvial e da tens�o de cisalhamento do       %
% vento- Miranda et al. 2012 - Cap. 11 p. 389, eq. 11.62% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TCF=Ro.*((-0.092*g*h^2*Rhox)/Rho0)+Ro.*(2.212*k1*U0*uf)-0.369*TCV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estratificacao verticalmente homog�nea         %
% Solu��o com escorregamento                     %
% Miranda et al. 2002 - equa��o 10.29 p. 351     %
% ou Miranda et al. 2012 - equa��o 10.49 p. 359  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

umistd=((g*(h^2)*Rhox)/(k1*U0*Rho0))*(0.167*z1.^3-0.296*z1.^2+0.058);
umistr=uf*(1.106*z1.^2+0.63);
umistv=(TCV/(k1*U0))*(0.316*z1.^2-z1+0.395).*(1./Ro);
umistot=umistd+umistr+umistv
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C�lculo do perfil te�rico da velocidade u        % 
% Solu��o sem escorregamento                       %
% Miranda et al. 2002 - equa��o 10.33 p. 348       %
% ou Miranda et al. 2012 - equa��o 10.33 9 p. 356  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Roaz=Ro.*Az;

umistd1=((g*(h^3)*Rhox)/Rho0)*(0.167*z1.^3-0.188*z1.^2+0.0208).*(1./Az);
umistr1=-uf*(1.5*z1.^2-1.5);
umistv1=TCV*h*(0.75*z1.^2-z1+0.25).*(1./Roaz);
umistot1=umistd1+umistr1+umistv1

plot(u,z4,'k');
hold on
plot(umistot,z4,'k*-');
plot(umistot1,z4,'ko-');
plot([0 0],[0 1],'k');
hold off
title('Baia do Trapande  27/08/98 - Quadratura ')
xlabel('Velocidade Longitudinal (m s^{-1})')
ylabel('Profundidade Adimensional, Z')
%print -dpcx256 /user/angra/home2/bernardes/simcana/fiscan2708umist.pcx
legend('u','umistot','umistot1',4)
figure


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estratificacao parcialmente misturada %
% Solu��o com escorregamento            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C�lculo do perfil te�rico da velocidade u  % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


UD=((g*h^2*Rhox)/(k1*U0*Rho0))*(0.167*z1.^3-0.25*z1.^2+0.0417);
UR=uf;
UV=(TCV/(k1*U0))*(0.5*z1.^2-z1+0.333).*(1./Ro);
UF=(TCF./(k1*U0)).*(0.5*z1.^2-0.0417).*(1./Ro);
UTOT=UD+UR+UV+UF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estratificacao parcialmente misturada %
% Solucao sem escorregamento            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C�lculo do perfil te�rico da velocidade u  % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UD1=((g*h^3*Rhox)/Rho0)*(0.167*z1.^3-0.188*z1.^2+0.0208).*(1./Az);
UR1=uf*(-1.5*z1.^2+1.5);
UV1=(TCV*h)*(0.75*z1.^2-z1+0.25).*(1./Roaz);
UTOT1=UD1+UR1+UV1


plot(u,z4,'k');
hold on
plot(UTOT,z4,'k*-');
plot(UTOT1,z4,'ko-');
plot([0 0],[0 1],'k');
hold off
title('Baia do Trapande  27/08/98 - Quadratura ')
xlabel('Velocidade Longitudinal (m s^{-1})');
ylabel('Profundidade Adimensional, Z');
legend('u','UTOT','UTOT1',4);
%print -dpcx256 /user/angra/home2/bernardes/simcana/fiscan2708u.pcx
figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estratifica��o parcialmente misturada %
% Solu��o com escorregamento            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C�lculo do perfil te�rico do componente vertical %
% da velocidade w                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Rhox2=(24.5-2*12)/(x/2)^2;
Bh4=B*h^4;
Bh4x=Bh4/x;
Bh2=B*h^2;
Bh2x=Bh2/x;
gh2=g*h^2;
gh2x=gh2/x;

WD=((h*Rhox*gh2x)/(k1*U0*Rho0)+((g*h^3*Rhox2)/(k1*U0*Rho0)))*(-0.0417*z1.^4+0.083*z1.^3-0.0417*z1);
W=WD*100000

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estratificacao parcialmente misturada %
% Solu��o sem escorregamento            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C�lculo do perfil te�rico da velocidade w  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WD1=((g/(B*Rho0))*(Bh4x*Rhox+Bh4*Rhox2)*(-0.0417*z1.^4+0.0625*z1.^3-0.0208*z1)).*(1./Az);
WV1=(((TCV*Bh2x)/B)*(-0.25*z1.^3+0.5*z1.^2-0.25*z1)).*(1./Roaz);
W1=(WD1+WV1)*100000

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% As velocidades verticais est�o com sinal negativo p/%
% que os referenciais utilizados sejam padronizados   %
% (modelo original o referencial cresce para baixo,   %
% enquanto os resultados est�o sendo apresentados     %
% com referencial crescendo para cima.  		      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


plot(-W,z4,'k*-');
hold on
plot(-W1,z4,'ko-');
plot([0 0],[0 1],'k');
hold off
title('Baia do Trapande  27/08/98 - Quadratura ')
xlabel('Velocidade Vertical (10^{-5} m s^{-1})');
ylabel('Profundidade Adimensional, Z');
legend('W','W1',3);

figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C�lculo da salinidade te�rica %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solu��o com escorregamento    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

smed=mean(s);
smed2=[smed smed smed smed smed smed smed smed smed smed]';
SD=smed+((u0*h^2*Sx)/Kz)*(0.144*z1.^5-0.425*z1.^4+0.5*z1.^2-0.105);
SR=uf*(-0.091*z1.^5+0.36*z1.^4-0.057);
ST=SD+SR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solu��o com atrito m�ximo de fundo %
% ou seja, velocidade nula fundo     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SD1=smed+((u0*h^2*Sx)/Kz)*(0.4*z1.^5-0.75*z1.^4+0.5*z1.^2-0.083);
SR1=uf*(-0.6*z1.^5+z1.^4-0.1);
ST1=SD1+SR1


plot(s,z4,'k');
hold on
plot(ST1,z4,'ko-');
plot(ST,z4,'k*-');
plot([smed2 smed2],[0 1],'k');
hold off
title('Baia do Trapande  27/08/98 - Quadratura ')
xlabel('Salinidade');
ylabel('Profundidade Adimensional, Z');
legend('s','ST1','ST',3);


