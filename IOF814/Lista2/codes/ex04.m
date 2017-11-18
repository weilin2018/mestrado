%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IOF814 - Prof Joseph Harari                                        %
%  Aluno: Danilo Augusto Silva         nusp: 7279456                  %
%                                                                     %
%                         QUESTÃO 04 - LISTA 02                       %                                                                  %
%   Os dados de saidas sao salvos no diretorio ../outputs/Q05         %
%                                                                     %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;clc

% aplicacao do metodo dos minimo quadrados

% criar grade conforme dada na lista

kmax = 10;
jmax = 10;

xgrid = linspace(1,10,1);
ygrid = linspace(1,10,1);

temp = ones(kmax,jmax)*NaN;

% preencher grade com os valores
temp(3,2) = 22.0;
temp(6,2) = 22.5;
temp(9,4) = 24.5;
temp(2,6) = 22.4;
temp(4,8) = 22.5;
temp(9,9) = 25.5;

%% metodo dos minimos quadrados pela inpaints_nans
T = inpaint_nans(temp,1);

figure(1)
contourf(T);
colorbar;
title('Condicao Inicial: Grade com Temperatura Superficial Interpolada','fontsize',15);
xlabel('Distancia (EW)','fontsize',12);
ylabel('Distancia (NS)','fontsize',12);

figure(3)
pcolor(temp);
title('Grade com Temperatura Superficial - Dados Originais','fontsize',15);
xlabel('Distancia (EW)','fontsize',12);
ylabel('Distancia (NS)','fontsize',12);