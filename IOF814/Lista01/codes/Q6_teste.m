
clear all; close all; clc

fprintf(['Modelos Numericos Aplicados a Processos Costeiros e Estuarinos - IOF814\r',...
    '1a. Lista de exercicios\r',...
    'Nome: Octavio Ambrosio Pereira Jr\r\r',...
    'Opcoes:\r',...
    '1 - velocidade para sul\r',...
    '2 - velocidade para sudoeste\r',...
    '3 - velocidade para oeste\r\r',...
    '4 - despejo superior direito\r',...
    '5 - despejo central direito\r',...
    '6 - despejo inferior direito\r\r',...
    ]);

vel=input('Digite o numero da velocidade: ');
des=input('Digite o numero do despejo: ');

% Definicao de variaveis
nmax=1800;
jmax=150;
kmax=150;
dt=10;
kx=10;
ky=10;
freqplot=60;
concorte=0.0001;
cderr=5;
r=5e-3;

% Grade batimetrica
h_file=load('../data/etopo1.xyz');
lon=h_file(:,1);
lat=h_file(:,2);
bat=h_file(:,3);

lonmin=min(lon);
lonmax=max(lon);
latmin=min(lat);
latmax=max(lat);
batmin=min(bat);
batmax=max(bat);

% Grade numerica
nlat=linspace(latmin,latmax,kmax);
nlon=linspace(lonmin,lonmax,jmax);
[dLAT,dLON]=meshgrid(nlat,nlon);

% Interpolacao
nbat=griddata(lon,lat,bat,dLON,dLAT);
nbat=nbat'; dLON=dLON'; dLAT=dLAT';

% Fator de conversao: graus - metros
latc=(latmax+latmin)/2;       %latitude central
Rt=6371e3;                    %raio da Terra
Ct=2*pi*Rt;                   %comprimento da Terra
dlat=Ct/360;                  %latitude em metro
dlon=dlat*cos(deg2rad(latc)); %longitude em metro

% comprimento da grade em graus
Llat=latmax-latmin;
Llon=lonmax-lonmin;
% comprimento da grade em metros
Ly=Llat*dlat;
Lx=Llon*dlon;
% resolucao da grade em metros
dy=Ly/kmax;
dx=Lx/jmax;

% plotando batimetria
pmax=60; prange=linspace(batmin,pmax,30);
pbat=nbat; pbat(pbat>0)=pmax;

figure (1)
contourf(dLON,dLAT,pbat,prange)
title('Batimetria (m)')
xlabel('Longitude')
ylabel('Latitude')
colormap('jet')
colorbar

% Definicao de chave para os pontos de grade
% 1 maritimo, 0 terrestre
kmar=nbat;
kmar(kmar>0)=0;
kmar(kmar~=0)=1;

% Plotando chave
figure (2)
contourf(dLON,dLAT,kmar)
title('Chave da regiao modelada (1 mar, 0 terra)')
xlabel('Longitude')
ylabel('Latitude')
map = [210,180,140; 135,206,250]/255;
colormap(map)
colorbar

% Definicao de chave para area costeira
% 1 area costeira
kcost=nbat;
kcost(kcost>0)=0;
kcost(kcost<-10)=0;
kcost(kcost~=0)=1;

% Plotando chave
figure(3)
contourf(dLON,dLAT,kcost)
title('Chave da regiao costeira')
xlabel('Longitude')
ylabel('Latitude')
colorbar



% matriz do somatorio do contaminante
scon=zeros(kmax,jmax);

%Condicoes iniciais
t0=-1;
fant=zeros(kmax,jmax);
fatu=zeros(kmax,jmax);
fren=zeros(kmax,jmax);

if vel==1
    u=0; v=-0.5;
elseif vel==2
    u=-0.35; v=-0.35;
elseif vel==3
    u=-0.5;  v=0;
else
    fprintf(['Opcao nao identificada\r',...
        'velocidade 3 selecionada\r'])
    u=-0.5;  v=0;
end;

%Campo de velocidades
u=ones(kmax,jmax)*u;
v=ones(kmax,jmax)*v;

if des==4
    xderr=ceil(jmax*0.75);
    yderr=ceil(kmax*0.75);
elseif des==5
    xderr=ceil(jmax*0.75);
    yderr=ceil(kmax*0.50);
elseif des==6
    xderr=ceil(jmax*0.75);
    yderr=ceil(kmax*0.25);
else
    fprintf(['Opcao nao identificada\r',...
        'velocidade 3 selecionada\r'])
    xderr=ceil(jmax*0.75);
    yderr=ceil(kmax*0.25);
end;

[xder,yder]=meshgrid(xderr,yderr);

quadv=dt/dx;
qvadv=dt/dy;
qudif=2*dt*kx/dx/dx;
qvdif=2*dt*ky/dy/dy;
rdec=1+2*dt*r;
%'
kplot=2;
for n=3:nmax
    fant(yder,xder)=cderr;
    fatu(yder,xder)=cderr;

    fren(2:kmax-1,2:jmax-1)=fant(2:kmax-1,2:jmax-1)...
        -kmar(2:kmax-1,3:jmax).*kmar(2:kmax-1,1:jmax-2).*...
        u(2:kmax-1,2:jmax-1)*quadv.*(fatu(2:kmax-1,3:jmax)-fatu(2:kmax-1,1:jmax-2))...
        -kmar(3:kmax,2:jmax-1).*kmar(1:kmax-2,2:jmax-1).*...
         v(2:kmax-1,2:jmax-1)*qvadv.*(fatu(3:kmax,2:jmax-1)-fatu(1:kmax-2,2:jmax-1))...
        +qudif*kmar(2:kmax-1,3:jmax).*kmar(2:kmax-1,1:jmax-2).*...
        (fant(2:kmax-1,3:jmax)-2*fant(2:kmax-1,2:jmax-1)+fant(2:kmax-1,1:jmax-2))...
        +qvdif*kmar(3:kmax,2:jmax-1).*kmar(1:kmax-2,2:jmax-1).*...
        (fant(3:kmax,2:jmax-1)-2*fant(2:kmax-1,2:jmax-1)+fant(1:kmax-2,2:jmax-1));

    % registrando regioes contaminadas
    scon=scon+fren;

    % registrando tempo que o contaminante atinge areas costeiras
    if t0<0
        concost=fren.*kcost;
        concost=max(max(concost));
        if concost>concorte
            t0=n*dt;
        end
    end

    fren=fren.*kmar;
    fren(fren<concorte)=0;

    kplot=kplot+1;
    if(kplot==freqplot)
        kplot=0;
        maximo=max(max(fren));
        figure(5)
        contour(dLON,dLAT,nbat,[0.1 0.2 0.3],'LineWidth',2);
        hold;
        plot(nlon(xder),nlat(yder),'xm','LineWidth',2)
        contourf(dLON,dLAT,fren,[concorte:maximo])
        colorbar
        title(['Adv,dif, conc em t=',num2str(n*dt),'seg'],'fontsize',12)
        xlabel('Longitude')
        ylabel('Latitude')
        grid on
        hold off;
        pause(0.2)
    end

    fant=fatu;
    fatu=fren;
end

% 2 - removendo concentracoes fora da linha de costa
mdec=scon.*kcost;

% Determinando regioes
c=zeros(kmax,jmax);
n=c; s=c; l=c; o=c; ne=c; no=c; se=c; so=c;

rlat=ceil(linspace(1,jmax,4));
rlat=sort([rlat rlat(2:3)+1]);
rlon=ceil(linspace(1,kmax,4));
rlon=sort([rlon rlon(2:3)+1]);

% linha (y) - coluna (x)
c(rlon(3):rlon(4),rlat(3):rlat(4))=1;
n(rlon(5):rlon(6),rlat(3):rlat(4))=1;
s(rlon(1):rlon(2),rlat(3):rlat(4))=1;
l(rlon(3):rlon(4),rlat(5):rlat(6))=1;
o(rlon(3):rlon(4),rlat(1):rlat(2))=1;
ne(rlon(5):rlon(6),rlat(5):rlat(6))=1;
no(rlon(5):rlon(6),rlat(1):rlat(2))=1;
se(rlon(1):rlon(2),rlat(5):rlat(6))=1;
so(rlon(1):rlon(2),rlat(1):rlat(2))=1;

% grafico das regioes
figure(4)
subplot(3,3,1); contourf(no);
subplot(3,3,3); contourf(ne);
subplot(3,3,4); contourf(o);
subplot(3,3,5); contourf(c);
subplot(3,3,6); contourf(l);
subplot(3,3,7); contourf(so);
subplot(3,3,8); contourf(s);
subplot(3,3,9); contourf(se);
subplot(3,3,2); contourf(n);
title('regioes do modelo')

mno=max(max(mdec.*no));
mn=max(max(mdec.*n));
mne=max(max(mdec.*ne));
mo=max(max(mdec.*o));
mc=max(max(mdec.*c));
ml=max(max(mdec.*l));
mso=max(max(mdec.*so));
ms=max(max(mdec.*s));
mse=max(max(mdec.*se));

% resultado da regi�o
s_string = 'O contaminante atingiu a regiao ';
if mno > 0; disp([s_string, 'noroeste']); end;
if mn > 0; disp([s_string, 'norte']); end;
if mne > 0; disp([s_string, 'nordeste']); end;
if mo > 0; disp([s_string, 'oeste']); end;
if mc > 0; disp([s_string, 'central']); end;
if ml > 0; disp([s_string, 'leste']); end;
if mso > 0; disp([s_string, 'sudoeste']); end;
if ms > 0; disp([s_string, 'sul']); end;
if mse > 0; disp([s_string, 'sudeste']); end;

% tempo em que o contaminante atingiu a costa
if t0>0
    disp(['O contaminante atingiu a costa no tempo ',num2str(t0),'s'])
else
    disp(['O contaminante n�o atingiu a costa'])
end
