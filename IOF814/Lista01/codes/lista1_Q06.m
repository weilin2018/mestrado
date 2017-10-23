%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%                         QUESTÃO 06 - LISTA 01                       %
%                                                                     %
%   Advecção: explicito, centrado no tempo                            %
%   Difusão: explicito, avancado no tempo                             %
%   decaimento: implicito                                             %
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
dx=100;%espacamento em x em metros
dy=100;%espacamento em y em metros
dt=30;%intervalo de tempo em segundos
kx=10; %difusão em x em m/s
ky=10;%difusão em y em m/s
r=0.5e-3; % coeficiente de decaimento
kmax = round((abs(latmin - latmax)/dx));%definição de kmax (em y)
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
%'
colorbar;
title(['Chaves da regiao: 1-mar, 0-terra'], 'fontsize',12);
axis([1 jmax 1 kmax]);
xlabel('Indices na grade - EW','fontsize',12);
ylabel('Indices na grade - NS','fontsize',12);

%%%%%Condicoes iniciais
fant=zeros(kmax,jmax);
fatu=zeros(kmax,jmax);
fren=zeros(kmax,jmax);

%Localizacao dos pontos de descarte
jcont=-46.35*60*1852; % longitude do despejo continuo
kcont=-24.01*60*1852; % latitude do despejo continuo
ccont=0.5;            % concentracao do despejo continuo

xderr=-46.35*60*1852; % longitude do derrame
yderr=-24.10*60*1852; % latitude do derrame
cderr=10.0;           % concentracao do derrame

%[xder,yder]=meshgrid(xderr,yderr);
xder=120;
yder=120;
fant(yder,xder)=cderr;
fatu(yder,xder)=cderr;

% campo de velocidade
u=ones(kmax,jmax)*0.00;
v=ones(kmax,jmax)*1.00;

quadv=dt/dx;
qvadv=dt/dy;
qudif=2*dt*kx/dx/dx;
qvdif=2*dt*ky/dy/dy;
rdec=1+2*dt*r;

%%%%%Equacao da adveccao - difusao - decaimento 2D
kplot=2;
for n=3:nmax
fren(2:kmax-1,2:jmax-1)=(fant(2:kmax-1,2:jmax-1)...
    -u(2:kmax-1,2:jmax-1)*quadv.*(fatu(2:kmax-1,3:jmax)-fatu(2:kmax-1,1:jmax-2))...
    -v(2:kmax-1,2:jmax-1)*qvadv.*(fatu(3:kmax,2:jmax-1)-fatu(1:kmax-2,2:jmax-1))...
    +qudif*(fant(2:kmax-1,3:jmax)-2*fant(2:kmax-1,2:jmax-1)+fant(2:kmax-1,1:jmax-2))...
    +qvdif*(fant(3:kmax,2:jmax-1)-2*fant(2:kmax-1,2:jmax-1)+fant(1:kmax-2,2:jmax-1)))/rdec;
ind=find(fren<concorte);
fren(ind)=0;
kplot=kplot+1;
if(kplot==freqplot)
kplot=0;
%%%%%calculo da soma das concentracoes e seu maximo
soma=sum(sum(fren));
maximo=max(max(fren));
figure(1)
plot(xgrid(xder),ygrid(yder),'xm','LineWidth',2)
hold
contourf(xgrid,ygrid,fren,[concorte:maximo])
colorbar
axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
tempo=n*dt;
%%%%%calculo da soma das concentracoes e seu maximo
soma=sum(sum(fren));
maximo=max(max(fren));
title(['Adv,dif,dec conc em t=',num2str(tempo),'seg soma ',num2str(soma),' max ',num2str(maximo)],'fontsize',12)
xlabel('DISTANCIA NA GRADE (m)','fontsize',12)
ylabel('DISTANCIA NA GRADE (m)','fontsize',12)
grid on
hold off
%pause
end
fant=fatu;
fatu=fren;
end
