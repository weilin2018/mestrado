k%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODELO NUMERICO DE ADVECCAO E DIFUSAO 2D PARA A REGIÃO COSTEIRA DE      %
% VITÓRIA - ES(40.5W-40W; 20S-20.5S), COM DISPERSAO DE CONTAMINANTES NA   %
% SUPERFICIE DESPEJADOS EM UM PONTO SUJEITO A CORRENTE DE 0.5M/S.         %
% A BATIMETRIA ORIGINAL, POSSUI GRANDEZA GRAU DECIMAL, QUE PRONTAMENTE FOI%
% TRANFORMADA EM UTM                                                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

%vendo a cara do arquivo nc
%ncdisp('GEBCO_2014_2D_-55.0_-40.0_-30.0_-20.0.nc');
load ../data/danilo.mat;

%Carregamento do arquivo com dados de batimetria do litoral de Vitória - ES
file = load('../data/etopo1.xyz');
%
% %extraindo latitude, longitude e batimetria
% %coordenadas em utm
% lon = file(:,1)*60*1852;
% lat = file(:,2)*60*1852;
% bat = file(:,3);
%
% lonmin=min(lon);
% lonmax=max(lon);
% latmin=min(lat);
% latmax=max(lat);
% profmin=min(bat);
% profmax=max(bat);
%
% %Quantidade de pontos de grade e parametros iniciais
% dx=149;%espacamento em x em metros
% dy=150;%espacamento em y em metros
% dt=30;%intervalo de tempo em segundos
% kx=10; %difusão em x em m/s
% ky=10;%difusão em y em m/s
% kmax = round((abs(latmin - latmax) /dx));%definição de kmax (em y)
% jmax = round((abs(lonmin - lonmax)/dx)); %definição de jmax (em x)
% nmax=3600;%passos de tempo
% freqplot=100;% frequencia de plotagem
% concorte=0.0001;
% xgrid=((1:jmax)-1)*dx;
% ygrid=((1:kmax)-1)*dy;
%
% % Plot da batimetria
% latplot=latmin:dy:latmax;
% lonplot=lonmin:dx:lonmax;
% [LAT,LON]=meshgrid(latplot,lonplot);
% BAT=griddata(lon,lat,bat,LON,LAT);

figure (1)
contourf(LON,LAT,BAT)
title('BATIMETRIA (m)')
xlabel('LONGITUDE')
ylabel('LATITUDE')
colorbar

[XBAT,YBAT]=meshgrid(xgrid,ygrid);

% Definicao de chave para os pontos de grade
% 1 maritimo, 0 terrestre
kmar=BAT;
kmar(kmar>0)=0;
kmar(kmar~=0)=1;

[XKMAR,YKMAR]=meshgrid(1:jmax,1:kmax);
%plotando chaves
figure(2)
contourf(XKMAR,YKMAR,kmar','LineWidth',2);
colorbar;
title(['Chaves da regiao modelada (1 mar, 0 terra)'],'fontsize',12)
axis([1 jmax 1 kmax])
xlabel('Indices na grade - EW','fontsize',12)
ylabel('Indices na grade - NS','fontsize',12)
%print -djpeg fig_kmar

% Detectando se o contaminante atinge areas costeiras
% 1 terrestre e mar prof., 0 area costeira
mask=BAT;
mask(mask>0)=0;
mask(mask<-10)=0;
mask(mask~=0)=1;
figure(3)
contourf(XKMAR,YKMAR,mask')
colorbar
title('Regiao costeira')

% Separando a área costeira em regioes
costa=zeros(kmax,jmax);
ne=costa; c=costa; so=costa;

% 1-120 121-240 241-368
ne(241:368,1:end)=1;
c(121:240,1:end)=1;
so(1:120,1:end)=1;

% grafico das regioes
figure(4)
subplot(3,1,1); contourf(ne); hold
contour(XKMAR,YKMAR,BAT',[0.1 0.2 0.3]','LineWidth',2,'LineColor','k'); hold off
subplot(3,1,2); contourf(c); hold
contour(XKMAR,YKMAR,BAT',[0.1 0.2 0.3]','LineWidth',2,'LineColor','k'); hold off
subplot(3,1,3); contourf(so); hold
contour(XKMAR,YKMAR,BAT',[0.1 0.2 0.3]','LineWidth',2,'LineColor','k'); hold off



%-------------------------------CORRENTE SUL------------------------------%

%ponto de despejo
kdesp=300;jdesp=300;  % local do despejo continuo
despejo=0.5;        % lancamento continuo

%inicializacao das matrizes:
fant=zeros(kmax,jmax);
fatu=zeros(kmax,jmax);
fren=zeros(kmax,jmax);
fdespejo=zeros(kmax,jmax);

%Condicao inicial
fdespejo(kdesp,jdesp)=despejo;

u=-0.50;      %componente zonal da velocidade (m/s)
v=0;          %componente meridional da velocidade (m/s)

quadv=(dt/dx);
qvadv=(dt/dy);
qudif=2*dt*kx/dx/dx;
qvdif=2*dt*ky/dy/dy;


%%%%%Equacao da adveccao - difusao 2D
kplot=2;
for n=3:nmax
fren(2:kmax-1,2:jmax-1)=(fant(2:kmax-1,2:jmax-1)...
    -kmar(2:kmax-1,3:jmax).*kmar(2:kmax-1,1:jmax-2).*...
    u*quadv.*(fatu(2:kmax-1,3:jmax)-fatu(2:kmax-1,1:jmax-2))...
    -kmar(3:kmax,2:jmax-1).*kmar(1:kmax-2,2:jmax-1).*...
     v*qvadv.*(fatu(3:kmax,2:jmax-1)-fatu(1:kmax-2,2:jmax-1))...
    +qudif*kmar(2:kmax-1,3:jmax).*kmar(2:kmax-1,1:jmax-2).*...
    (fant(2:kmax-1,3:jmax)-2*fant(2:kmax-1,2:jmax-1)+fant(2:kmax-1,1:jmax-2))...
    +qvdif*kmar(3:kmax,2:jmax-1).*kmar(1:kmax-2,2:jmax-1).*...
    (fant(3:kmax,2:jmax-1)-2*fant(2:kmax-1,2:jmax-1)+fant(1:kmax-2,2:jmax-1)));

fren=fren.*kmar;
ind=find(fren<concorte);
fren(ind)=0;


% plot de resultados
kplot=kplot+1;
if(kplot==freqplot)
kplot=0;
maximo=max(max(fren));
figure(5)
contour(XBAT,YBAT,BAT',[0.1 0.2 0.3],'LineWidth',2);
hold
plot(xgrid(kdesp),ygrid(jdesp),'xm','LineWidth',2)
contourf(xgrid,ygrid,fren',[concorte:maximo])
colorbar
%axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
tempo=n*dt;
title(['Adv,dif conc em t=',num2str(tempo),'seg'],'fontsize',12)
xlabel('DISTANCIA NA GRADE (m)','fontsize',12)
ylabel('DISTANCIA NA GRADE (m)','fontsize',12)
grid on
hold off
%pause
end

fren=fren+fdespejo;
fant=fatu;
fatu=fren;
end


%-------------------------------CORRENTE SUDOESTE-------------------------%

%ponto de despejo
kdesp=300;jdesp=300;  % local do despejo continuo
despejo=0.5;        % lancamento continuo


%inicializacao das matrizes:
fant=zeros(kmax,jmax);
fatu=zeros(kmax,jmax);
fren=zeros(kmax,jmax);
fdespejo=zeros(kmax,jmax);

%Condicao inicial
fdespejo(kdesp,jdesp)=despejo;

u=-0.50;      %componente zonal da velocidade (m/s)
v=-0.50;          %componente meridional da velocidade (m/s)

quadv=(dt/dx);
qvadv=(dt/dy);
qudif=2*dt*kx/dx/dx;
qvdif=2*dt*ky/dy/dy;

x=0;

%%%%%Equacao da adveccao - difusao 2D
kplot=2;
for n=3:nmax
fren(2:kmax-1,2:jmax-1)=(fant(2:kmax-1,2:jmax-1)...
    -kmar(2:kmax-1,3:jmax).*kmar(2:kmax-1,1:jmax-2).*...
    u*quadv.*(fatu(2:kmax-1,3:jmax)-fatu(2:kmax-1,1:jmax-2))...
    -kmar(3:kmax,2:jmax-1).*kmar(1:kmax-2,2:jmax-1).*...
     v*qvadv.*(fatu(3:kmax,2:jmax-1)-fatu(1:kmax-2,2:jmax-1))...
    +qudif*kmar(2:kmax-1,3:jmax).*kmar(2:kmax-1,1:jmax-2).*...
    (fant(2:kmax-1,3:jmax)-2*fant(2:kmax-1,2:jmax-1)+fant(2:kmax-1,1:jmax-2))...
    +qvdif*kmar(3:kmax,2:jmax-1).*kmar(1:kmax-2,2:jmax-1).*...
    (fant(3:kmax,2:jmax-1)-2*fant(2:kmax-1,2:jmax-1)+fant(1:kmax-2,2:jmax-1)));


%determinando se o contaminate chega a costa e em que perido de tempo
tempo=n*dt;
fren=fren.*kmar;
tdesp=fren.*mask;
tdesp=max(max(tdesp));
if x==0
	if tdesp>0; t=tempo;
	x=1;
	end
end

ind=find(fren<concorte);
fren(ind)=0;


% plot de resultados
kplot=kplot+1;
if(kplot==freqplot)
kplot=0;
maximo=max(max(fren));
figure(6)
contour(XBAT,YBAT,BAT',[0.1 0.2 0.3],'LineWidth',2);
hold
plot(xgrid(kdesp),ygrid(jdesp),'xm','LineWidth',2)
contourf(xgrid,ygrid,fren',[concorte:maximo])
colorbar
%axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
tempo=n*dt;
title(['Adv,dif conc em t=',num2str(tempo),'seg'],'fontsize',12)
xlabel('DISTANCIA NA GRADE (m)','fontsize',12)
ylabel('DISTANCIA NA GRADE (m)','fontsize',12)
grid on
hold off
%pause
end
fren=fren+fdespejo;
fant=fatu;
fatu=fren;
end

% Removendo concentracoes fora da linha de costa
ndesp=fren.*mask;
maximo=max(max(ndesp));
figure(7)
contourf(XKMAR,YKMAR,ndesp',linspace(concorte,maximo,10)); hold
contour(XKMAR,YKMAR,BAT',[0.1 0.2 0.3]','LineWidth',2,'LineColor','k'); hold off
colorbar
title('concentracao na regiao costeira')

% Contaminate na área costeira em regioes

ne=max(max(ndesp.*ne));
c=max(max(ndesp.*c));
so=max(max(ndesp.*so));

% define a região
t = num2str(t);
str = sprintf('O contaminante atingiu em t= %s ',t);
if ne > 0; disp([str, 'seg a regiao nordeste']); end;
if c > 0; disp([str, 'seg a regiao centro']); end;
if so > 0; disp([str, 'seg a regiao sudeste']); end;


%-------------------------------CORRENTE OESTE------------------------------%

%ponto de despejo
kdesp=300;jdesp=200;  % local do despejo continuo
despejo=0.5;        % lancamento continuo

%inicializacao das matrizes:
fant=zeros(kmax,jmax);
fatu=zeros(kmax,jmax);
fren=zeros(kmax,jmax);
fdespejo=zeros(kmax,jmax);
%Condicao inicial
fdespejo(kdesp,jdesp)=despejo;

u=-0.50;      %componente zonal da velocidade (m/s)
v=-0.50;          %componente meridional da velocidade (m/s)

quadv=(dt/dx);
qvadv=(dt/dy);
qudif=2*dt*kx/dx/dx;
qvdif=2*dt*ky/dy/dy;

x=0;

%%%%%Equacao da adveccao - difusao 2D
kplot=2;
for n=3:nmax
fren(2:kmax-1,2:jmax-1)=(fant(2:kmax-1,2:jmax-1)...
    -kmar(2:kmax-1,3:jmax).*kmar(2:kmax-1,1:jmax-2).*...
    u*quadv.*(fatu(2:kmax-1,3:jmax)-fatu(2:kmax-1,1:jmax-2))...
    -kmar(3:kmax,2:jmax-1).*kmar(1:kmax-2,2:jmax-1).*...
     v*qvadv.*(fatu(3:kmax,2:jmax-1)-fatu(1:kmax-2,2:jmax-1))...
    +qudif*kmar(2:kmax-1,3:jmax).*kmar(2:kmax-1,1:jmax-2).*...
    (fant(2:kmax-1,3:jmax)-2*fant(2:kmax-1,2:jmax-1)+fant(2:kmax-1,1:jmax-2))...
    +qvdif*kmar(3:kmax,2:jmax-1).*kmar(1:kmax-2,2:jmax-1).*...
    (fant(3:kmax,2:jmax-1)-2*fant(2:kmax-1,2:jmax-1)+fant(1:kmax-2,2:jmax-1)));

%determinando se o contaminate chega a costa e em que perido de tempo
tempo=n*dt;
fren=fren.*kmar;
tdesp=fren.*mask;
tdesp=max(max(tdesp));
if x==0
	if tdesp>0; t=tempo;
	x=1;
	end
end

ind=find(fren<concorte);
fren(ind)=0;


% plot de resultados
kplot=kplot+1;
if(kplot==freqplot)
kplot=0;
maximo=max(max(fren));
figure(8)
contour(XBAT,YBAT,BAT',[0.1 0.2 0.3],'LineWidth',2);
hold
plot(xgrid(kdesp),ygrid(jdesp),'xm','LineWidth',2)
contourf(xgrid,ygrid,fren',[concorte:maximo])
colorbar
%axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
tempo=n*dt;
title(['Adv,dif conc em t=',num2str(tempo),'seg'],'fontsize',12)
xlabel('DISTANCIA NA GRADE (m)','fontsize',12)
ylabel('DISTANCIA NA GRADE (m)','fontsize',12)
grid on
hold off
%pause
end
fren=fren+fdespejo;
fant=fatu;
fatu=fren;
end

% Removendo concentracoes fora da linha de costa
ndesp=fren.*mask;
maximo=max(max(ndesp));
figure(9)
contourf(XKMAR,YKMAR,ndesp',linspace(concorte,maximo,10)); hold
contour(XKMAR,YKMAR,BAT',[0.1 0.2 0.3]','LineWidth',2,'LineColor','k'); hold off
colorbar
title('concentracao na regiao costeira')

% Contaminate na área costeira em regioes

ne=max(max(ndesp.*ne));
c=max(max(ndesp.*c));
so=max(max(ndesp.*so));

% define a região
t = num2str(t);
str = sprintf('O contaminante atingiu em t= %s ',t);
if ne > 0; disp([str, 'seg a regiao nordeste']); end;
if c > 0; disp([str, 'seg a regiao centro']); end;
if so > 0; disp([str, 'seg a regiao sudeste']); end;
