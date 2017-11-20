clear all; close all; clc

% definindo espacamento (que faz mais sentido)

%% controle de plotagem, para facilitar a vida durante desenvolvimento.
% 0 = True: plotar
% 1 = False: nao plotar
DEPLOY=0;

%% carregar batimetria e dados de latitude/longitude
% definicao da grade batimetrica baseada no arquivo .xyz em anexo
file = load('../data/etopo1.xyz');
lon       = file(:,1)*60*1852; % conversao para UTM
lat       = file(:,2)*60*1852;
bat       = file(:,3);

% extrair os valores maximos e minimos dos dados carregados
lonmin    = min(lon);
lonmax    = max(lon);
latmin    = min(lat);
latmax    = max(lat);
batmin    = min(bat);
batmax    = max(bat);

%% Parametros numericos do modelo
nmax      = 1200;             % numero de passos de tempo maximo
dt        = 120;              % passo de tempo
dx        = 2000;             % espacamento de grade em x
dy        = 2000;             % espacamento de grade em y
kx        = 10;               % coeficiente de difusao
ky        = 10;               % coeficiente de difusao
freqplot  = 30;               % frequencia de plotagem

% parametros baseado na grade batimetrica utilizada
jmax      = round((abs(min(lon) - max(lon)))/dx);    % nro de pontos da grade em x
kmax      = round((abs(min(lat) - max(lat)))/dy);    % nro de pontos da grade em y
jmax2     = jmax*2;           % pontos de grade em x duplicado (grade alternada)
kmax2     = kmax*2;           % pontos de grade em y duplicado (grade alternada)
xgrid     = ((1:jmax)-1)*dx;
ygrid     = ((1:kmax)-1)*dy;
[X,Y]     = meshgrid(xgrid,ygrid);  % grid com pontos de plotagem tipo eta
xgrid2    = ((1:jmax2)-1)*dx/2;
ygrid2    = ((1:kmax2)-1)*dy/2;
[XBAT,YBAT] = meshgrid(xgrid2,ygrid2); % grid com pontos de plotagem tipo u,v

%% Parametros hidrodinamicos do modelo
dens      = 1024;                  % densidade media da agua do mar
latid     = 25*pi/180;             % latitude (em graus, para rad)
fco       = 2*7.292E-5*sin(latid); % parametro de Coriolis
skipWind  = 5;                     % reduzir densidade de vetores de vento
g         = 9.8;                   % aceleracao da gravidade
rfric     = 0.02;                  % coeficiente de friccao no fundo
etaContorno   = 0.00025;           % elevacao nos contornos abertos

%% Definindo uma nova batimetria para a grade determinada
% gerar grade
nlat         = latmin:dy/2:latmax;
nlon         = lonmin:dx/2:lonmax;
[llat, llon] = meshgrid(nlat, nlon); % gridar os pontos de grade
llon2 = llon(1:2:end,1:2:end-1);
llat2 = llat(1:2:end,1:2:end-1);

% interpolar os dados de batimetria para a grade gerada acima
nbat=griddata(lon,lat,bat,llon,llat);

fprintf('Generating Figure 1: bathymetry\n');
prange              = linspace(batmin,5,40);
batPlot             = nbat;
batPlot(batPlot>0)  = 60;

if DEPLOY==0
  figure(1)
  contourf(llon,llat,batPlot,prange);
  title('Bathymetry (m)')
  xlabel('DISTANCE (m) EW', 'fontsize', 12)
  ylabel('DISTANCE (m) NS', 'fontsize', 12)
  colormap('winter')
  colorbar
  % saving figure
  out = ['../outputs/ex05/bathymetry'];
  grafico=['print -djpeg ', out];
  eval(grafico);
end

%% Define keys for land (0) and sea (1)
kmar          = nbat;
kmar(kmar>0)  = 0;
kmar(kmar~=0) = 1;

%determinacao da batimetria e chaves maritimas em u, v e eta:
% nota: usando somente metade da grade
bate          = kmar(1:2:end-1,2:2:end)';     % batimetria em pontos do tipo eta
batu          = kmar(1:2:end-1,1:2:end-1)';   % batimetria em pontos do tipo u
batv          = kmar(2:2:end,2:2:end)';       % batimetria em pontos do tipo v %'
kmare         = bate*0;
kmare(bate>0) = 1;                  % valores positivos indicam terra (1)
kmaru         = batu*0;
kmaru(batu>0) = 1;
kmarv         = batv*0;
kmarv(batv>0) = 1;

[XKMAR,YKMAR]=meshgrid(1:jmax,1:kmax);
if DEPLOY==0
  % plotando chaves
  figure(2)
  contourf(XKMAR,YKMAR,kmarv,'LineWidth',2);
  colorbar; grid;
  title(['Chaves da regiao modelada (1 mar, 0 terra)'],'fontsize',15)
  axis([1 XKMAR(end) 1 YKMAR(end)])
  xlabel('Indices na grade - EW','fontsize',12)
  ylabel('Indices na grade - NS','fontsize',12)
  % saving figure
  out = ['../outputs/ex05/keys'];
  grafico=['print -djpeg ', out];
  eval(grafico);
end

%% calculo dos coeficientes e demais termos das equacoes discretizadas
D2y       = kx./dy/dy;
D2x       = ky./dx/dx;
dt2       = dt*2;
denbatu   = dens.*batu;
denbatv   = dens.*batv;
coef      = 1+rfric*dt2;          % termo de decaimento

% Condicoes iniciais de repouso (0-valores anteriores, 1-atuais, 2-renovados)
eta0      = zeros(kmax,jmax);
U0        = zeros(kmax,jmax);
V0        = zeros(kmax,jmax);
eta1      = zeros(kmax,jmax);
U1        = zeros(kmax,jmax);
V1        = zeros(kmax,jmax);
eta2      = zeros(kmax,jmax);
U2        = zeros(kmax,jmax);
V2        = zeros(kmax,jmax);

% forcante do sistema: vento
rhoAr=1.25;         % air density
fric=2.6*1E-3;      % drag coefficient

uwind(1:kmax,1:jmax)=7.0711;          % zonal component
vwind(1:kmax,1:jmax)=7.0711;          % meridional component
wwind=sqrt(uwind.^2 + vwind.^2);      % speed
taux=fric*rhoAr.*uwind.*wwind;        % wind stress x
tauy=fric*rhoAr.*vwind.*wwind;        % wind stress y

ventomax=max(max(wwind));

if DEPLOY==0
  fprintf('Generating figure 2: wind field\n');
  figure(3)
    contour(llon,llat,batPlot,'LineWidth',1.5,'color','k');
  title('Bathymetry (m)')
  xlabel('DISTANCE (m) EW', 'fontsize', 12)
  ylabel('DISTANCE (m) NS', 'fontsize', 12)
  hold on
  quiver(llon2',llat2',uwind,vwind,'LineWidth',2);%''
  title(['Wind - Maximum Intensity ',...
    num2str(ventomax),' m/s'],'fontsize',12)
  axis([llon(1) llon(end) llat(1) llat(end)])
  %axis equal
  xlabel('DISTANCE (m) EW','fontsize',12)
  ylabel('DISTANCE (m) NS','fontsize',12)
  % saving figure
  out = ['../outputs/ex05/windField'];
  grafico=['print -djpeg ', out];
  eval(grafico);
end

%% loop no tempo

kplot   = 1;
for n=2:nmax
  tempo=n*dt;
  kplot=kplot+1;
  % print to control all process of the program, printing mean u,v and eta
  fprintf('Calculando timestep %i/%i - umean=%2.5f, vmean=%2.5f, elevmean=%2.5f\n',...
                             n,nmax,mean(mean(U2)),mean(mean(V2)),mean(mean(eta2)));

  % condicao de contorno aberto de elevacao
  etamet=ones(kmax,jmax).*etaContorno;

  for j=1:jmax
      eta2(1,j)=eta1(1,j) + etamet(1,j)*kmare(1,j);
      eta2(2,j)=eta1(2,j) + etamet(2,j)*kmare(2,j);
      eta2(kmax-1,j)=eta1(kmax-1,j) + etamet(kmax-1,j)*kmare(kmax-1,j);
      eta2(kmax,j)=eta1(kmax,j) + etamet(kmax,j)*kmare(kmax,j);
  end
  for k=1:kmax
      eta2(k,1)=eta1(k,1) + etamet(k,1)*kmare(k,1);
      eta2(k,2)=eta1(k,2) + etamet(k,2)*kmare(k,2);
      eta2(k,jmax-1)=eta1(k,jmax-1) + etamet(k,jmax-1)*kmare(k,jmax-1);
      eta2(k,jmax)=eta1(k,jmax) + etamet(k,jmax)*kmare(k,jmax);
  end



  % eq da continuidade
  for j=2:jmax-1
      for k=2:kmax-1
         if kmare(k,j)>0
         forcx=(batu(k,j+1).*U1(k,j+1)-batu(k,j).*U1(k,j))/dx;
         forcy=(batv(k,j).*V1(k,j)-batv(k-1,j).*V1(k-1,j))/dy;
         eta2(k,j)=eta0(k,j)-dt2*(forcx+forcy);
         end
      end
  end

% eq do movimento em x
  for j=2:jmax-1
      for k=2:kmax-1
          if (kmaru(k,j)*kmare(k,j)*kmare(k,j-1))>0
            vmedu=(V1(k,j)+V1(k,j-1)+V1(k-1,j-1)+V1(k-1,j))/4;
            umed=(U1(k,j-1)+2*U1(k,j)+U1(k,j+1))/4;
            advx=(umed./dx*2)*(U1(k,j+1)-U1(k,j-1))+(vmedu./dy*2)*(U1(k+1,j)-U1(k-1,j));
            difusx=(D2x)*(U0(k,j+1)-2*U0(k,j)+U0(k,j-1))+(D2y)*(U0(k+1,j)-2*U0(k,j)+U0(k-1,j));
            forc=fco*vmedu-g.*(eta1(k,j)-eta1(k,j-1))./dx+taux(k,j)./denbatu(k,j);
            U2(k,j)=(U0(k,j)+(forc-advx+difusx)*dt2)/coef;
          end
       end
  end

  % eq do movimento em y
  for j=2:jmax-1
      for k=2:kmax-1
          if (kmarv(k,j)*kmare(k,j)*kmare(k+1,j))>0
            umedv=(U1(k,j)+U1(k+1,j)+U1(k+1,j+1)+U1(k,j+1))/4;
            vmed=(V1(k-1,j)+2*V1(k,j)+V1(k+1,j))/4;
            advy=(umedv./dx*2)*(V1(k,j+1)-V1(k,j-1))+(vmed./dy*2)*(V1(k+1,j)-V1(k-1,j));
            difusy=(D2x)*(V0(k,j+1)-2*V0(k,j)+V0(k,j-1))+(D2y)*(V0(k+1,j)-2*V0(k,j)+V0(k-1,j));
            forc=-fco*umedv-g.*(eta1(k+1,j)-eta1(k,j))./dy+tauy(k,j)./denbatv(k,j);
            V2(k,j)=(V0(k,j)+(forc-advy+difusy)*dt2)/coef;
          end
      end
  end

  % condicao de contorno (u e v) com extrapolacao linear
  for j=1:jmax
      %eta2(1,j)=(2*eta2(2,j)-eta2(3,j))*kmare(1,j);
      %eta2(kmax,j)=(2*eta2(kmax-1,j)-eta2(kmax-2,j))*kmare(kmax,j);
      U2(1,j)=(2*U2(2,j)-U2(3,j))*kmaru(1,j);
      V2(1,j)=(2*V2(2,j)-V2(3,j))*kmarv(1,j);
      U2(kmax,j)=(2*U2(kmax-1,j)-U2(kmax-2,j))*kmaru(kmax,j);
      V2(kmax,j)=(2*V2(kmax-1,j)-V2(kmax-2,j))*kmarv(kmax,j);
  end
  for k=1:kmax
      %eta2(k,1)=(2*eta2(k,2)-eta2(k,3))*kmare(k,1);
      %eta2(k,jmax)=(2*eta2(k,jmax-1)-eta2(k,jmax-2))*kmare(k,jmax);
      U2(k,1)=(2*U2(k,2)-U2(k,3))*kmaru(k,1);
      V2(k,1)=(2*V2(k,2)-V2(k,3))*kmarv(k,1);
      U2(k,jmax)=(2*U2(k,jmax-1)-U2(k,jmax-2))*kmaru(k,jmax);
      V2(k,jmax)=(2*V2(k,jmax-1)-V2(k,jmax-2))*kmarv(k,jmax);
  end

  % renovando as variaveis no tempo
  eta0      = eta1;
  eta1      = eta2;
  U0        = U1;
  V0        = V1;
  U1        = U2;
  V1        = V2;

  % plotar resultados
  if kplot==freqplot
      kplot = 0 ;           % resetamos o contatos de plotagem

      % definicao de u e v para pontos de grade do tipo eta:
      uplot = U2;
      vplot = V2;
      for j=1:jmax-1
         for k=1:kmax
            if kmare(k,j)>0
              uplot(k,j)=(U2(k,j)+U2(k,j+1))/2;
            end
         end
      end
      for j=1:jmax
         for k=2:kmax
            if kmare(k,j)>0
              vplot(k,j)=(V2(k,j)+V2(k-1,j))/2;
            end
         end
      end

      % plotar elevacao e corrente
      velo=sqrt(uplot.^2+vplot.^2);
      velomax=max(max(velo));
      etamax = max(max(eta2(:,:)));
      etamin = min(min(eta2(:,:)));

      figure(5)
      contourf(X,Y,eta2,'LineWidth',0.1);
      colorbar;
      title(['Elev (m) - time ',num2str(tempo/60),...
            ' min. Limits ',num2str(etamin),' a ',num2str(etamax),' m'],'fontsize',12)
      axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
      xlabel('DISTANCE (m) EW','fontsize',12)
      ylabel('DISTANCE (m) NS','fontsize',12)
      % saving figure
      out = ['../outputs/ex05/elev_',num2str(tempo)];
      grafico=['print -djpeg ', out];
      eval(grafico);

      figure(6)
      quiver(X,Y,uplot,vplot,'LineWidth',0.8);
      title(['Currents (m/s) - time ',num2str(tempo/60),...
          ' min - intens max ',num2str(velomax),' m/s'],'fontsize',12)
      %axis equal
      axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
      xlabel('DISTANCE (m) EW','fontsize',12)
      ylabel('DISTANCE (m) NS','fontsize',12)
      % saving figure
      out = ['../outputs/ex05/curr_',num2str(tempo)];
      grafico=['print -djpeg ', out];
      eval(grafico);
  end

  pause(.2)
end
