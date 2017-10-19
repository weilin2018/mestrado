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
q=c*dt/(12*dx);
fatu=zeros(jmax,1);
fant=zeros(jmax,1)
fren=zeros(jmax,1);
ffilt=zeros(jmax,1);
xgrid=((1:jmax)-1)*dx;

%condicao inicial
fatu(jin:jfi)=ci;
fcin=fatu;

% FORMULAS DE RENOVACAO
% formula geral
fren(3:jmax-2)=fatu(3:jmax-2) - q*(fatu(1:jmax-4) - 8*fatu(2:jmax-3) + 8*fatu(4:jmax-1) - fatu(5>jmax));
% formula para 1o e 2o ponto de grade
fren(1)=fatu(1)-q*(-25*fatu(1)+48*fatu(2)-36*fatu(3)+16*fatu(4)-3*fatu(5));
fren(2)=fatu(2)-q*(-3*fatu(2)-10*fatu(2)+18*fatu(2)+18*fatu(3)-g*fatu(4)+fatu(5));
% formula para penultimo e ultimo
