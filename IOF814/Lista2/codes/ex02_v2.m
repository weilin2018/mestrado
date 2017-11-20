% SOLUCAO DO SISTEMA DE EQUACOES HIDRODINAMICAS 1D LINEAR
% COM ESQUEMA DE SEGUNDA ORDEM NO TEMPO.
% CONDICAO DE CONTORNO: ELEVACAO DE -0.2M NO PONTO INICIAL(SUL) E DE 0.2 NO FINAL(NORTE) DA GRADE
% EFEITOS DE DIFUSAO E DECAIMENTO

% CONFIGURACAO DA GRADE
% U E U E U E U E U E U E U E U E U E U E
% - + - + - + - + - + - + - + - + - + - +
% (1) (2) ... j-2 j-1 (j) j+1 j+2 ... (jmax)
clear all; close all; clc

% Parametros do Modelo e da grade
jmax          = 50;               % numero de pontos da grade
nmax          = 400;              % numero de passos de tempo
dx            = 100;              % espacamento de grade
dt            = 5;                % passo de tempo
H             = 3;                % profundidade media
amp           = 0.2;              % amplitude da elevacao
per           = 120;              % periodo da elevacao
freqplot      = 40;               % frequencia de plotagem
D             = 0.002;             % coeficiente de difusao
fric          = 0.004;            % coeficiente de decaimento
g             = 9.8;
omega         = 2*pi/per;
xgrid         = ((1:jmax)-1)*dx;

% Criacao das matrizes dos dados
eant          = zeros(jmax,1);
eatu          = zeros(jmax,1);
eren          = zeros(jmax,1);
vant          = zeros(jmax,1);
vatu          = zeros(jmax,1);
vren          = zeros(jmax,1);

% Calculo dos coeficientes das equacoes discretizadas
qu            = g*(dt/dx);
qe            = dt*H/dx;
qd            = 2*D*dt/(dx*dx);
qr            = fric*2*dt;
amp2          = 2*amp;
cor2          = 2;

% CONDICOES INICIAIS
eant(1)       = amp*sin(omega*dt);
eatu(1)       = amp*sin(omega*dt*2);

% LOOP NO TEMPO
kplot=1;

for n=2:nmax
     tempo=n*dt;
     kplot=kplot+1;

     % calculo de v e eta renovados
     vren(2:jmax-1) = vant(2:jmax-1) - qu*(eatu(2:jmax-1)-eatu(1:jmax-2))...
      + qd*(vant(3:jmax) - 2*vant(2:jmax-1) + vant(1:jmax-2)) - qr*(vant(2:jmax-1));

     eren(2:jmax-1) = eant(2:jmax-1) - qe*(vatu(3:jmax) - vatu(2:jmax-1));

     % definicao das condicoes de contorno para cada variavel
     eren(1) = -amp;
     eren(jmax-1) = amp;
     eren(jmax) = 2*eren(jmax-1) - eren(jmax-2);

     vren(1) = 2*vren(2)-vren(3);
     vren(jmax) = 2*vren(jmax-1)-vren(jmax-2);

     % plotagem dos resultados

     if(kplot==freqplot)
        kplot = 0;

        vplot=vren;
        vplot(2:jmax-1) = (vplot(3:jmax) + vplot(2:jmax-1))/2;

        subplot(2,1,1)
        plot(xgrid,eren,'LineWidth',2)
        axis([xgrid(1) xgrid(jmax) -amp2 amp2]);
        grid on
        title(['Cana de Sao Sebastiao com ventos de sudoeste - tempo = ',...
            num2str(tempo/60),' minutos'],'fontsize',12)
        ylabel('Elevacao (m)','fontsize',12)
        subplot(2,1,2)
        plot(xgrid,vplot,'LineWidth',2)
        axis([xgrid(1) xgrid(jmax) -cor2 cor2]);
        grid on
        xlabel('DISTANCIA NA GRADE (m)','fontsize',12)
        ylabel('Corrente (m/s)','fontsize',12)
        out = ['../outputs/ex02/ssb_circ_',num2str(tempo)];
        grafico=['print -djpeg ', out];
        eval(grafico);

     end

     % evolucao no tempo das variaveis
     eant   = eatu;
     eatu   = eren;
     vant   = vatu;
     vatu   = vren;

end
