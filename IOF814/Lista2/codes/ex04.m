%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IOF814 - Prof Joseph Harari                                        %
%  Aluno: Danilo Augusto Silva         nusp: 7279456                  %
%                                                                     %
%                         QUESTÃO 04 - LISTA 02                       %                                                                  %
%   Os dados de saidas sao salvos no diretorio ../outputs/Q04         %
%                                                                     %
%  Metodo dos Minimos Quadrados atraves de uma solucao de             %
%           sistema linear para determinacao dos coeficientes         %
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

%% metodo dos minimos quadrados para estabelecimento de condicao inicial de modelos
T     = temp;             % duplicando variavel para poder trabalhar
[n,m] = size(T);          % salvar as dimensoes da matriz de entrada
T     = T(:);             % todas as linhas
nm    = n*m;              % total de valores que serao necessarios
k     = isnan(T(:));      % indices de elementos NaN

nan_list  = find(k);      % lista com indices dos pontos NaN
know_list = find(~k);     % lista de valores conhecidos

% converter os indices (linear) de valores NaN para o formato [linha,coluna]
[nr,nc]   = ind2sub([n,m],nan_list);

nan_list  = [nan_list,nr,nc]; % coluna 1: indices em linha, col 2: indices de linha e col 3: indices de colunas

% utilizando diferencas finitas, computa-se a segunda derivada das linhas e colunas
% primeiro para as linhas
[i,j]     = ndgrid(2:(n-1),1:m);
ind       = i(:)+(j(:)-1)*n;
np        = (n-2)*m;
fda       = sparse(repmat(ind,1,3), [ind-1,ind,ind+1], repmat([1 -2 1], np, 1),n*m,n*m);
% agora para as colunas
[i,j]     = ndgrid(1:n,2:(m-1));
ind       = i(:)+(j(:)-1)*n;
np        = n*(m-2);
fda       = fda+sparse(repmat(ind,1,3), [ind-n,ind,ind+n], repmat([1 -2 1],np,1),nm,nm);

rhs = -fda(:,know_list)*T(know_list);   % eliminar os valores conhecidos
k         = find(any(fda(:,nan_list),2));

% resolver o sistema linear usando o mldivide do matlab
B         = T;
B(nan_list(:,1)) = fda(k,nan_list(:,1))\rhs(k);

% retornando a matriz interpolada para o formato da matriz original
T         = reshape(B,n,m);

%% plotar os dados para conferir o resultado
figure(1)
contourf(T);
colorbar;
caxis([21 26]);
title({'Condicao Inicial: Grade com Temperatura','Superficial Interpolada'},'fontsize',14);
xlabel('Distancia (EW)','fontsize',12);
ylabel('Distancia (NS)','fontsize',12);
out = ['../outputs/ex04/gradeInterpolada'];
grafico=['print -djpeg ', out];
eval(grafico);

figure(2)
pcolor(temp);
colorbar;
caxis([21 26]);
title({'Grade com Temperatura Superficial','Dados Originais'},'fontsize',14);
xlabel('Distancia (EW)','fontsize',12);
ylabel('Distancia (NS)','fontsize',12);
out = ['../outputs/ex04/gradeOriginal'];
grafico=['print -djpeg ', out];
eval(grafico);
