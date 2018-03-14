%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programa desenvolvivo para a Questão 4 da Lista 1 de exercícios da disciplina
% IOF814 - Modelos Numéricos Aplicados a Processos Costeiros e Estuarinos
% ministrada pelo Prof Joseph Harari, do Instituto Oceanográfico da USP.
%
% Para mais detalhes do desenvolvimento da discretização, pode ser conferido
% na solução da Lista, em IOF814/Lista1/outputs/Lista1_IOF814.pdf
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

clear all;
close all;

% definicao dos parametros do modelo, como tamanho da grade, passo de tempo e espaco
% e outras informacoes
%ler parametros do modelo
jmax=250;      %numero de pontos da grade em x
nmax=7200/2;      %numero de passos de tempo
dt=1;         %passo de tempo, em segundos
dx=100;        %espacamento da grade, em m
c=0.2;           %velocidade da corrente em x, em m/s
ci=1;          %valor da concentracao inicial
jin=25;        %indice do ponto inicial da grade com ci
jfi=50;        %indice do ponto final da grade com ci
freqplot=100;   %frequencia de plotagem
conc=2;        %limite maximo do eixo da concentracao no grafico
freqfilt=30;   %frequencia de filtragem

%calculo das constantes do modelo
dx12=12*dx;
q=c*dt/dx;
fatu=zeros(jmax,1);
fren=zeros(jmax,1);
ffilt=zeros(jmax,1);
xgrid=((1:jmax)-1)*dx;

%condicao inicial
fatu(jin:jfi)=ci;
fcin=fatu;
kplot=0;
contfilt=2;

for n=2:nmax
  tempo=n*dt;
  kplot=kplot+1;

  % FORMULA DE RECORRENCIA
  fren(4:jmax-3)=fatu(4:jmax-3) - q*(fatu(2:jmax-5) - 8*fatu(3:jmax-4) + 8*fatu(5:jmax-2) - fatu(6:jmax-1));
  %% aplicando formula de recorrencia em cada ponto extremo da grade
  % formula para 1o e 2o ponto de grade
  fren(2)=fatu(1)-q*(-25*fatu(1)+48*fatu(2)-36*fatu(3)+16*fatu(4)-3*fatu(5));
  fren(3)=fatu(2)-q*(-3*fatu(2)-10*fatu(3)+18*fatu(4)-6*fatu(5)+fatu(6));
  % formula para penultimo e ultimo pontos de grade
  fren(jmax-2)=fatu(jmax-2)-q*(-fatu(jmax-5)+6*fatu(jmax-4)-18*fatu(jmax-3)+10*fatu(jmax-2)+3*fatu(jmax-1));
  fren(jmax-1)=fatu(jmax-1)-q*(3*fatu(jmax-5)-16*fatu(jmax-4)+36*fatu(jmax-3)-48*fatu(jmax-2)+25*fatu(jmax-1));

  if(kplot==freqplot)
    kplot=0;
    %kfig=kfig+1;
    %figure(kfig)
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
    grafico=['print -djpeg ../outputs/Q04_ord04/q04_ord04_',XX];
    eval(grafico);
    pause(0.1)
    hold
end
fatu=fren;


end

toc
