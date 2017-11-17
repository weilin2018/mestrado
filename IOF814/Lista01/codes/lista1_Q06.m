%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IOF814 - Prof Joseph Harari                                        %
%  Aluno: Danilo Augusto Silva         nusp: 7279456                  %
%                                                                     %
%                         QUESTÃO 06 - LISTA 01                       %
%                                                                     %
%   Advecção: explicito, centrado no tempo                            %
%   Difusão: explicito, avancado no tempo                             %
%   decaimento: implicito                                             %
%                                                                     %
%   Informacoes quanto a simulacao realizada:                         %
%     Despejo continuo (emissario) em 46.35oW 24.01oS                 %
%     Despejo instantaneo (derrame) em 46.35oW 24.10oS                %
%                                                                     %
%                                                                     %
%   Os dados de saidas sao salvos no diretorio ../outputs/Q06         %
%                                                                     %
%   Em caso de esquecer como determinar as componentes de velocidade  %
%   a partir de uma velocidade, consultar:                            %
%   https://www.khanacademy.org/science/physics/two-dimensional-motion/two-dimensional-projectile-mot/a/what-are-velocity-components
%                                                                     %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;clc

%% Parte I - definindo parâmetros do modelo (coeficientes, grades)
nmax=1800;           % tempo de simulacao
jmax=150;            % tamanho da grade em x
kmax=150;            % tamanho da grade em y
dt=10;               % passo de tempo
kx=10;               % coeficiente de difusao
ky=10;               % coeficiente de difusao
freqplot=60;         % frequencia de plotagem
concorte=0.0000001;     % limite de corte da concentracao
cderrame=20;         % concentracao do derramente de oleo
cdespejo=5;          % concentracao de despejo do emissario
r=0.5e-4;            % fator de decaimento
lon_desp=-46.35;     % ponto X (lon) de despejo do emissario
lat_desp=-24.01;     % ponto Y (lat) de despejo do emissario
lon_derr=-46.35;     % ponto X (lon) de derrame
lat_derr=24.10;      % ponto Y (lat) de derrame

% definicao da grade batimetrica baseada no arquivo .xyz em anexo
batFile = load('../data/etopo1.xyz');
lon=batFile(:,1);
lat=batFile(:,2);
bat=batFile(:,3);

% extrair os valores maximos e minimos dos dados carregados
lonmin=min(lon);
lonmax=max(lon);
latmin=min(lat);
latmax=max(lat);
batmin=min(bat);
batmax=max(bat);

% gerar grade
nlat = linspace(latmin, latmax, kmax);
nlon = linspace(lonmin, lonmax, jmax);
[llat, llon] = meshgrid(nlat, nlon); % gridar os pontos de grade

% baseado na grade criada, vamos determinar os pontos da grade
% despejo do emissario.
% Importante: encontrar utilizando a menor distancia entre o ponto dado e grade
% nao estava dando certo, entao foi definido manualmente os pontos de despejo,
% determnando-se visualmente as localizacoes.
xdesp=75;
ydesp=118;
% despejo do derrame
xderr=75;
yderr=38;

% interpolar os dados de batimetria para a grade gerada acima
nbat=griddata(lon,lat,bat,llon,llat);

%% definindo a resolucao do modelo baseados em xy do arquivo carregado

% para determinar a resolucao da grade baseado nos pontos obtidos
% no arquivo .xyz e nas definicoes no inicio do modelo para kmax e jmax,
% devemos converter os dados de llat e llon para metros, utilizando a
% funcao do matlab deg2rad

theta = (latmax+latmin)/2;       %  latitude central da regicao modelada
Rt=6371e3;                       %raio da Terra
Ct=2*pi*Rt;                      %comprimento da Terra
dlat=Ct/360;                     %latitude em metro
dlon=dlat*cos(deg2rad(theta));   %longitude em metro

% comprimento da grade em graus convertendo para metros
latL = (latmax - latmin)*dlat;
lonL = (lonmax - lonmin)*dlon;

% definindo a resolucao da grade
dx=lonL/jmax;
dy=latL/kmax;

fprintf('Gerando figura 1: batimetria utilizada\n');
prange=linspace(batmin,60,30);
batPlot = nbat;
batPlot(batPlot>0)=60;

figure(1)
contourf(llon,llat,batPlot,prange)
title('Batimetria (m)')
xlabel('Longitude')
ylabel('Latitude')
colormap('winter')
colorbar
% saving figure
grafico=['print -djpeg ../outputs/Q06/batimetria'];
eval(grafico);
%close(1)

%% Definicao das chaves para terra (0) e mar (1)
kmar=nbat;
kmar(kmar>0)=0;
kmar(kmar~=0)=1;

fprintf('Gerando figura 2: chave da regiao modelada\n')
figure(2)
contourf(llon,llat,kmar)
title('Chave da regiao modelada (1 mar, 0 terra)')
xlabel('Longitude')
ylabel('Latitude')
map = [210,180,140; 135,206,250]/255;
colormap(map)
colorbar
% saving figure
grafico=['print -djpeg ../outputs/Q06/chaves'];
eval(grafico);
%close(2)

%% Loop no tempo para verificar a evolucao das plumas

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          CENARIO 1 - VENTO DE SUL                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc

fprintf('Cenario 1: ventos de 1m/s para Norte\n');
% definicao das matrizes de dados
fant=zeros(kmax,jmax);    % matriz de valores anteriores (n-1)
fatu=zeros(kmax,jmax);    % matriz de valores atuais (n)
fren=zeros(kmax,jmax);    % matriz de valores renovados (n+1)
fdes=zeros(kmax,jmax);    % matriz de despejo

% definicao do campo de velocidade (vento de sul)
theta=90;
angle=(theta*3.14)/180;
speed=1.;
u=speed*cos(angle);
v=speed*sin(angle);

% criando o campo de velocidade
v=ones(kmax,jmax)*u;
u=ones(kmax,jmax)*v;

% condicoes iniciais de concentracao

% derramente, por ser instantaneo, fica fora do loop do tempo
fdes(xderr,yderr)=cderrame;
fant=fdes;
fatu=fdes;

quadv=dt/dx;
qvadv=dt/dy;
qudif=2*dt*kx/dx/dx;
qvdif=2*dt*ky/dy/dy;
rdec=1+2*dt*r;

kplot=2;
contFig=0;    % contador de figuras para nomear as imagens salvas
for n=3:nmax
  % determinando o despejo continuo
  fant(xdesp,ydesp)=cdespejo;
  fatu(xdesp,ydesp)=cdespejo;

  % formula de recorrencia

  fren(2:kmax-1,2:jmax-1)=(fant(2:kmax-1,2:jmax-1)...
      -kmar(2:kmax-1,3:jmax).*kmar(2:kmax-1,1:jmax-2).*...
      u(2:kmax-1,2:jmax-1)*quadv.*(fatu(2:kmax-1,3:jmax)-fatu(2:kmax-1,1:jmax-2))...
      -kmar(3:kmax,2:jmax-1).*kmar(1:kmax-2,2:jmax-1).*...
       v(2:kmax-1,2:jmax-1)*qvadv.*(fatu(3:kmax,2:jmax-1)-fatu(1:kmax-2,2:jmax-1))...
      +qudif*kmar(2:kmax-1,3:jmax).*kmar(2:kmax-1,1:jmax-2).*...
      (fant(2:kmax-1,3:jmax)-2*fant(2:kmax-1,2:jmax-1)+fant(2:kmax-1,1:jmax-2))...
      +qvdif*kmar(3:kmax,2:jmax-1).*kmar(1:kmax-2,2:jmax-1).*...
      (fant(3:kmax,2:jmax-1)-2*fant(2:kmax-1,2:jmax-1)+fant(1:kmax-2,2:jmax-1)))/rdec;

  % atualizando a matriz de contaminantes
  fdes = fdes+fren;

  fren=fren.*kmar;
  fren(fren<concorte)=0;

  kplot=kplot+1;
  if(kplot==freqplot)
      contFig=contFig+1;
      kplot=0;
      maximo=max(max(fren));
      figure(5)
      contour(llon,llat,nbat,[0.1 0.2 0.3],'LineWidth',2);
      hold;
      %plot(nlon(xderr),nlat(yderr),'xm','LineWidth',2)
      contourf(llon,llat,fren,[concorte:maximo])
      colorbar
      title(['Vento para Norte: Adv,dif e dec, conc em t=',num2str(n*dt),'seg'],'fontsize',12)
      xlabel('Longitude')
      ylabel('Latitude')
      grid on
      % saving figure
      numberFig=sprintf('%03d',contFig);
      grafico=['print -djpeg ../outputs/Q06/windFromSouth/windS_', numberFig];
      eval(grafico);

      hold off;
      pause(0.2)
  end

  fant=fatu;
  fatu=fren;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     CENARIO 2 - VENTO DE SUDOESTE                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc

fprintf('Cenario 1: ventos de 1m/s para Nordeste\n');
% definicao das matrizes de dados
fant=zeros(kmax,jmax);    % matriz de valores anteriores (n-1)
fatu=zeros(kmax,jmax);    % matriz de valores atuais (n)
fren=zeros(kmax,jmax);    % matriz de valores renovados (n+1)
fdes=zeros(kmax,jmax);    % matriz de despejo

% definicao do campo de velocidade (vento de sudoeste)
theta=45;
angle=(theta*3.14)/180;
speed=1.;
u=speed*cos(angle);
v=speed*sin(angle);

% criando o campo de velocidade
u=ones(kmax,jmax)*u;
v=ones(kmax,jmax)*v;

% condicoes iniciais de concentracao

% derramente, por ser instantaneo, fica fora do loop do tempo
fdes(xderr,yderr)=cderrame;
fant=fdes;
fatu=fdes;

quadv=dt/dx;
qvadv=dt/dy;
qudif=2*dt*kx/dx/dx;
qvdif=2*dt*ky/dy/dy;
rdec=1+2*dt*r;

kplot=2;
contFig=0;    % contador de figuras para nomear as imagens salvas
for n=3:nmax
  % determinando o despejo continuo
  fant(xdesp,ydesp)=cdespejo;
  fatu(xdesp,ydesp)=cdespejo;

  % formula de recorrencia

  fren(2:kmax-1,2:jmax-1)=(fant(2:kmax-1,2:jmax-1)...
      -kmar(2:kmax-1,3:jmax).*kmar(2:kmax-1,1:jmax-2).*...
      u(2:kmax-1,2:jmax-1)*quadv.*(fatu(2:kmax-1,3:jmax)-fatu(2:kmax-1,1:jmax-2))...
      -kmar(3:kmax,2:jmax-1).*kmar(1:kmax-2,2:jmax-1).*...
       v(2:kmax-1,2:jmax-1)*qvadv.*(fatu(3:kmax,2:jmax-1)-fatu(1:kmax-2,2:jmax-1))...
      +qudif*kmar(2:kmax-1,3:jmax).*kmar(2:kmax-1,1:jmax-2).*...
      (fant(2:kmax-1,3:jmax)-2*fant(2:kmax-1,2:jmax-1)+fant(2:kmax-1,1:jmax-2))...
      +qvdif*kmar(3:kmax,2:jmax-1).*kmar(1:kmax-2,2:jmax-1).*...
      (fant(3:kmax,2:jmax-1)-2*fant(2:kmax-1,2:jmax-1)+fant(1:kmax-2,2:jmax-1)))/rdec;

  % atualizando a matriz de contaminantes
  fdes = fdes+fren;

  fren=fren.*kmar;
  fren(fren<concorte)=0;

  kplot=kplot+1;
  if(kplot==freqplot)
      contFig=contFig+1;
      kplot=0;
      maximo=max(max(fren));
      figure(6)
      contour(llon,llat,nbat,[0.1 0.2 0.3],'LineWidth',2);
      hold;
      %plot(nlon(xderr),nlat(yderr),'xm','LineWidth',2)
      contourf(llon,llat,fren,[concorte:maximo])
      colorbar
      title(['Vento para Nordeste: Adv,dif e dec, conc em t=',num2str(n*dt),'seg'],'fontsize',12)
      xlabel('Longitude')
      ylabel('Latitude')
      grid on
      % saving figure
      numberFig=sprintf('%03d',contFig);
      grafico=['print -djpeg ../outputs/Q06/windFromSW/windSW_', numberFig];
      eval(grafico);

      hold off;
      pause(0.2)
  end

  fant=fatu;
  fatu=fren;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     CENARIO 3 - VENTO DE SUDESTE                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc

fprintf('Cenario 1: ventos de 1m/s para Noroeste\n');
% definicao das matrizes de dados
fant=zeros(kmax,jmax);    % matriz de valores anteriores (n-1)
fatu=zeros(kmax,jmax);    % matriz de valores atuais (n)
fren=zeros(kmax,jmax);    % matriz de valores renovados (n+1)
fdes=zeros(kmax,jmax);    % matriz de despejo

% definicao do campo de velocidade (vento de sudoeste)
theta=135;
angle=(theta*3.14)/180;
speed=-1.;
u=speed*cos(angle);
v=speed*sin(angle);

% criando o campo de velocidade
u=ones(kmax,jmax)*u;
v=ones(kmax,jmax)*v;

% condicoes iniciais de concentracao

% derramente, por ser instantaneo, fica fora do loop do tempo
fdes(xderr,yderr)=cderrame;
fant=fdes;
fatu=fdes;

quadv=dt/dx;
qvadv=dt/dy;
qudif=2*dt*kx/dx/dx;
qvdif=2*dt*ky/dy/dy;
rdec=1+2*dt*r;

kplot=2;
contFig=0;    % contador de figuras para nomear as imagens salvas
for n=3:nmax
  % determinando o despejo continuo
  fant(xdesp,ydesp)=cdespejo;
  fatu(xdesp,ydesp)=cdespejo;

  % formula de recorrencia

  fren(2:kmax-1,2:jmax-1)=(fant(2:kmax-1,2:jmax-1)...
      -kmar(2:kmax-1,3:jmax).*kmar(2:kmax-1,1:jmax-2).*...
      u(2:kmax-1,2:jmax-1)*quadv.*(fatu(2:kmax-1,3:jmax)-fatu(2:kmax-1,1:jmax-2))...
      -kmar(3:kmax,2:jmax-1).*kmar(1:kmax-2,2:jmax-1).*...
       v(2:kmax-1,2:jmax-1)*qvadv.*(fatu(3:kmax,2:jmax-1)-fatu(1:kmax-2,2:jmax-1))...
      +qudif*kmar(2:kmax-1,3:jmax).*kmar(2:kmax-1,1:jmax-2).*...
      (fant(2:kmax-1,3:jmax)-2*fant(2:kmax-1,2:jmax-1)+fant(2:kmax-1,1:jmax-2))...
      +qvdif*kmar(3:kmax,2:jmax-1).*kmar(1:kmax-2,2:jmax-1).*...
      (fant(3:kmax,2:jmax-1)-2*fant(2:kmax-1,2:jmax-1)+fant(1:kmax-2,2:jmax-1)))/rdec;

  % atualizando a matriz de contaminantes
  fdes = fdes+fren;

  fren=fren.*kmar;
  fren(fren<concorte)=0;

  kplot=kplot+1;
  if(kplot==freqplot)
      contFig=contFig+1;
      kplot=0;
      maximo=max(max(fren));
      figure(7)
      hold;
      %plot(nlon(xderr),nlat(yderr),'xm','LineWidth',2)
      contourf(llon,llat,fren,[concorte:maximo])
      colorbar
      title(['Vento para Noroeste: Adv,dif e dec, conc em t=',num2str(n*dt),'seg'],'fontsize',12)
      xlabel('Longitude')
      ylabel('Latitude')
      grid on
      % saving figure
      numberFig=sprintf('%03d',contFig);
      grafico=['print -djpeg ../outputs/Q06/windFromSE/windSE_', numberFig];
      eval(grafico);

      hold off;
      pause(0.2)
  end

  fant=fatu;
  fatu=fren;

end
