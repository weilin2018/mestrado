%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programa desenvolvivo para a Questão 2 da Lista 1 de exercícios da disciplina
% IOF814 - Modelos Numéricos Aplicados a Processos Costeiros e Estuarinos
% ministrada pelo Prof Joseph Harari, do Instituto Oceanográfico da USP.
%
% Para mais detalhes do desenvolvimento da discretização, pode ser conferido
% na solução da Lista, em IOF814/Lista1/outputs/Lista1_IOF814.pdf
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

% definindo os parametros iniciais do modelo
jmax=20;         % numero de celulas em x
kmax=20;         % numero de celulas em y
nmax=100;         % periodo de simulacao em segundos
pmax=kmax*jmax;  % quantidade de pontos total
u=0.3;           % componente u da velocidade
v=0.8;           % componente v da velocidade

dx=5;            % passo de espaco em x
dy=5;            % passo de espaco em y
dt=1;            % passo de tempo em s

freqplo=5;       % frequencia de plot

xgrid=((1:jmax)-1)/dx;
ygrid=((1:kmax)-1)/dy;

kdesp=2:2;
jdesp=2:2;      % posicao de despejo
conc=1;         % concentracao de despejo total

concorte=0.1;  % concentracao limite de corte

% CALCULOS INICIAIS
qu = u/(4*dx);
qv = v/(4*dy);

% Definicao das condicoes iniciais e inicializacao das matrizes de dado
fren=zeros(pmax,1);         % funcao de nivel n+1 (renovado)
fatu=zeros(pmax,1);         % funcao de nível n (atual)
fdes=zeros(kmax,jmax);      % funcao de despejo
faux=zeros(pmax,1);         % funcao auxiliar
grid=zeros(kmax,jmax);      % grade
d=zeros(pmax,1);            % definicao do vetor d para os termos f(n+1)
A=zeros(pmax);              % definicao da matriz A para os coeficientes
fdes(kdesp,jdesp)=conc;     % condicao inicial de despejo

% como fdesp e 2D, precisamos passar para 1D
for ii = 1:kmax
   faux(jmax*(ii-1)+1:jmax*ii) = fdes(ii,:);
end
fatu=fatu+faux;
fant=fatu;

% para popular a matriz A com as diagonais principais preenchidas pelos
% coeficientes determinados durante a discretizacao, fazemos:

% determinando os indices do contorno
inferiorBoundary=2:jmax-1;
superiorBoundary=pmax-jmax+2:pmax-1;
leftBoundary=find(rem(1:pmax,jmax)==0);
rightBoundary=find(rem(0:pmax-1,jmax)==0);

% concatenando todos os indices referentes ao contorno
contorno=cat(2,inferiorBoundary,superiorBoundary,leftBoundary,rightBoundary);

% pegando todos os indices que nao sao referentes aos contornos
interior=setdiff(1:pmax,contorno);

for ii=1:pmax
  if ii~=contorno
      A(ii,ii)=1/2*dt;      % coef de f^{n+1}_{j,k}
      A(ii,ii+1)=qu;        % coef de f^{n+1}_{j+1,k}
      A(ii,ii-1)=-qu;       % coef de f^{n+1}_{j-1,k}
      A(ii,ii+jmax)=qv;     % coef de f^{n+1}_{j,k+1}
      A(ii,ii-jmax)=-qv;    % coef de f^{n+1}_{j,k-1}
  else
      A(ii,ii)=1;           % condicao de contorno rigida
  end
end

contplo=0;                  % contador de plotagem
contFig=0;                  % contador de figuras para salvar o nome
% LOOP NO TEMPO
% CONDICOES DE CONTORNO
% FORMULA DE RECORRENCIA
% PLOTAGEM (PRESSIONE ENTER PARA EVOLUIR NO TEMPO)
% EVOLUCAO NO TEMPO DAS VARIAVEIS
for n=2:nmax
  tempo=n*dt;
  contplo=contplo+1;        % iteracao no contador

  % formula de recorrencia
  d(interior) = 0.5*fant(interior)/dt - qu*fatu(interior+1) + ...
                qu*fatu(interior-1) - qv*fatu(interior+jmax) + ...
                qv*fatu(interior-jmax);

  % resolver o sistema linear
  fren=A\d; % note que a barra invertida serve para resolver um sistema linear. Mais em https://www.mathworks.com/help/matlab/ref/mldivide.html

  if(contplo==freqplo)
    contplo=0;
    contFig=contFig+1;

    maxconc=max(fren);

    % convertendo o vetor para matriz
    plt=zeros(kmax,jmax);
    for ii=1:kmax
        plt(ii,:)=fren(jmax*(ii-1)+1:jmax*ii);
    end
    contourf(xgrid,ygrid,plt,[concorte:0.01:maxconc]);
    grid;
    axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)]);
    colorbar
    title(['Adveccao Bidimensional (semi implic, 2a ordem) - tempo ',...
        num2str(tempo),' segundos'],'fontsize',12);
    xlabel('DISTANCIA NA GRADE(m)','fontsize',12);
    ylabel('DISTANCIA NA GRADE(m)','fontsize',12);
    numberFig=sprintf('%03d',contFig);
    grafico=['print -djpeg ../outputs/Q02/q02_', numberFig];
    eval(grafico);

    hold off;
    pause(0.01);
  end

  fant=fatu;
  fatu=fren;

end
