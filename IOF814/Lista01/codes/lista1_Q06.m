%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%                         QUESTÃO 06 - LISTA 01                       %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

% carregar arquivo com batimetria e linha de costa
file = load('../data/etopo1.xyz');

%extraindo latitude, longitude e batimetria
%coordenadas em utm
lon = file(:,1)*60*1852;
lat = file(:,2)*60*1852;
bat = file(:,3);

lonmin=min(lon);
lonmax=max(lon);
latmin=min(lat);
latmax=max(lat);
profmin=min(bat);
profmax=max(bat);

%Quantidade de pontos de grade e parametros iniciais
dx=151;%espacamento em x em metros
dy=151;%espacamento em y em metros
dt=30;%intervalo de tempo em segundos
kx=10; %difusão em x em m/s
ky=10;%difusão em y em m/s
kmax = round((abs(latmin - latmax)/dx));%definição de kmax (em y)
jmax = round((abs(lonmin - lonmax)/dx)); %definição de jmax (em x)
nmax=3600;%passos de tempo
freqplot=100;% frequencia de plotagem
concorte=0.0001;
xgrid=((1:jmax)-1)*dx;
ygrid=((1:kmax)-1)*dy;

%Localizacao dos pontos de descarte
jcont=-46.35*60*1852; % longitude do despejo continuo
kcont=-24.01*60*1852; % latitude do despejo continuo
ccont=0.5;            % concentracao do despejo continuo

jdesp=-46.35*60*1852; % longitude do despejo pontual
kdesp=-24.10*60*1852; % latitude do despejo pontual
cdesp=10.0;           % concentracao do despejo pontual

% Plot da batimetria
latplot=latmin:dy:latmax;
lonplot=lonmin:dx:lonmax;
[LAT,LON]=meshgrid(latplot,lonplot);
%BAT=griddata(lon,lat,bat,LON,LAT);

load ../data/danilo.mat;

%%%%% plotar batimetria
figure (1)
contourf(LON,LAT,BAT)
title('BATIMETRIA (m)')
xlabel('LONGITUDE')
ylabel('LATITUDE')
colorbar

% gridar os pontos de grade criado
[XBAT,YBAT]=meshgrid(xgrid,ygrid);

%%%%% Definindo uma chave booleana para pontos de grade: terra (1) e oceano (0)
kmar=BAT;
kmar(kmar>0)=0;
kmar(kmar~=0)=1;

[XKMAR,YKMAR]=meshgrid(1:jmax,1:kmax);
figure(2)
contourf(XKMAR,YKMAR,kmar','LineWidth',2);
colorbar;
title(['Chaves da regiao: 1-mar, 0-terra'], 'fontsize',12);
axis([1 jmax 1 kmax]);
xlabel('Indices na grade - EW','fontsize',12);
ylabel('Indices na grade - NS','fontsize',12);

%% Secao 1: corrente de 1m/s para Norte

%%%%% Definindo as condições iniciais
fant=zeros(kmax,jmax);
fatu=zeros(kmax,jmax);
fren=zeros(kmax,jmax);
fdesp=zeros(kmax,jmax);
% pontos aleatorios para teste
kdesp=120; jdesp=120;
cdesp=0.5;

% componente zonal e meridional da velocidade: 1m/s para norte
uNorth=1;
vNorth=0;

quadv=(dt/dx);
qvadv=(dt/dy);
qudif=2*dt*kx/dx/dx;
qvdif=2*dt*ky/dy/dy;
rdec=1+2*dt*r;

kplot=2;
for n=3:nmax
  %fren


end
