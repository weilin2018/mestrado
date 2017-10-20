%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programa desenvolvivo para a Questão 4 da Lista 1 de exercícios da disciplina
% IOF814 - Modelos Numéricos Aplicados a Processos Costeiros e Estuarinos
% ministrada pelo Prof Joseph Harari, do Instituto Oceanográfico da USP.
%
% Para mais detalhes do desenvolvimento da discretização, pode ser conferido
% na solução da Lista, em IOF814/Lista1/outputs/Lista1_IOF814.pdf
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

% definicao dos parametros do modelo, como tamanho da grade, passo de tempo e espaco
% e outras informacoes
%ler parametros do modelo
jmax=250;      %numero de pontos da grade em x
nmax=720;      %numero de passos de tempo
dt=10;         %passo de tempo, em segundos
dx=100;        %espacamento da grade, em m
c=2;           %velocidade da corrente em x, em m/s
ci=1;          %valor da concentracao inicial
jin=25;        %indice do ponto inicial da grade com ci
jfi=50;        %indice do ponto final da grade com ci
freqplot=60;   %frequencia de plotagem
conc=2;        %limite maximo do eixo da concentracao no grafico
freqfilt=30;   %frequencia de filtragem

%calculo das constantes do modelo
q=c*dt/dx;
fatu=zeros(jmax,1);
fren=zeros(jmax,1);
xgrid=((1:jmax)-1)*dx;

%condicao inicial
fatu(jin:jfi)=ci;
fcin=fatu;
kplot=1;

for n=2:nmax
  tempo=n*dt;
  kplot=kplot+1;

  % funcoes de recorrencia
  fren(2:jmax-1)=fatu(2:jmax-1)- q*(fatu(2:jmax-1)-fatu(1:jmax-2));

  if(kplot==freqplot)
    kplot=0;
    %kfig=kfig+1;
    figure(1)
    plot(xgrid,fcin,'r','LineWidth',2)
    hold
    plot(xgrid,fren,'LineWidth',2)
    axis([xgrid(1) xgrid(jmax) -conc conc]);
    title(['Advec sinal retangular (esq. 4a. ordem) - tempo ',...
        num2str(tempo/3600),'h'],'fontsize',12)
    xlabel('DISTANCIA NA GRADE (m)','fontsize',12)
    ylabel('CONCENTRACAO','fontsize',12)
    grid on
    XX=num2str(tempo);
    % saving graph: please change the directory for some like result/.
    % this structure Im using is for my github acc
    grafico=['print -djpeg ../outputs/Q04_ord01/q04_ord01_',XX];
    eval(grafico);
    pause(1)
    hold
end
fatu=fren;


end
