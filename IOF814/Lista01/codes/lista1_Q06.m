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
kmax = round((abs(latmin - latmax) /dx));%definição de kmax (em y)
jmax = round((abs(lonmin - lonmax)/dx)); %definição de jmax (em x)
nmax=3600;%passos de tempo
freqplot=100;% frequencia de plotagem
concorte=0.0001;
xgrid=((1:jmax)-1)*dx;
ygrid=((1:kmax)-1)*dy;

% Plot da batimetria
latplot=latmin:dy:latmax;
lonplot=lonmin:dx:lonmax;
[LAT,LON]=meshgrid(latplot,lonplot);
BAT=griddata(lon,lat,bat,LON,LAT);

figure (1)
contourf(LON,LAT,BAT)
title('BATIMETRIA (m)')
xlabel('LONGITUDE')
ylabel('LATITUDE')
colorbar
