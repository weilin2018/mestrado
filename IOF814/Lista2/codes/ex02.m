%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ainda falta conferir se as equacoes estao corretas                    %
%   ainda falta colocar o vento forcando a circulacao tambem            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Problema: em k+1 ha problema de continuacao da grade quado eu     %
% testo um condicao pra saber se os vizinhos sao pontos de agua .       %
% Ver como resolver isso.                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all;clc

%% carregar batimetria e tratar informacoes
% import bathymetry file
load ('../data/batimetria_litoralnorte.dat');
lon=batimetria_litoralnorte(:,1);
lat=batimetria_litoralnorte(:,2);
bat=batimetria_litoralnorte(:,3);

% pelos limites especificados no exercicio:
latmin = -23.9;
latmax = -23.7;
lonmin = -45.5;
lonmax = -45.3;

% recortar a batimetria carregada para a regiao de interesse

% comecando pela latitude:
ilat = find(lat<=latmax & lat>=latmin);
cutlat = lat(ilat);
cutlon = lon(ilat);
cutbat = bat(ilat);

% agora pela longitude
ilon = find(cutlon<=lonmax & cutlon>=lonmin);
lat = cutlat(ilon);
lon = cutlon(ilon);
bat = cutbat(ilon);

clear ilon ilat cutlat cutlon cutbat batimetria_litoralnorte

%% parametro constantes do modelo
nmax = 900;                   % tempo de simulacao em timestep
kmax = 150;                   % tamanho da grade em x
jmax = 150;                   % tamanho da grade em y
jmax2 = jmax/2;
kmax2 = kmax/2;
dt = 30;                      % passo de tempo (em segundos)
Dx = 10;                      % coeficiente de difusao em x
Dy = 10;                      % coeficiente de difusao em y
freqplot = 10;                % frequencia de plotagem
etaContorno = 0.2;            % condicao de contorno de elevacao

% parametros de modelagem hidrodinamica
dens = 1024;                  % densidade media da agua do mar
latRef = -25*pi/180;          % latitude de referencia em radianos
fco  = 2*7.292E-5*sin(latRef);% frequencia de coriolis
skipWind = 10;                % para melhorar visualizacao
g = 9.8;                      % aceleracao da gravidade
rfric = 0.02;                 % coeficiente de friccao no fundo

%% geracao da grade

% gerar grade
nlon = linspace(lonmin, lonmax, jmax);
nlat = linspace(latmin, latmax, kmax);
[llat, llon] = meshgrid(nlat,nlon);

% interpolar batimetria para a grade gerada
nbat = griddata(lon,lat,bat,llon,llat);

% definindo a resolucao do modelo baseados em xy do arquivo carregado

% para determinar a resolucao da grade baseado nos pontos obtidos
% no arquivo .xyz e nas definicoes no inicio do modelo para kmax e jmax,
% devemos converter os dados de llat e llon para metros, utilizando a
% funcao do matlab deg2rad

theta = (latmax+latmin)/2;       % central latitude of the region
Rt=6371e3;                       % earths radius
Ct=2*pi*Rt;                      % earths length
dlat=Ct/360;                     % latitude in meters
dlon=dlat*cos(deg2rad(theta));   % longitude in meters

% Converting grids Length to meters
latL = (latmax - latmin)*dlat;
lonL = (lonmax - lonmin)*dlon;

% Define grid resolution
dx=lonL/jmax;
dy=latL/kmax;

xgrid=((1:jmax2)-1)*dx;
ygrid=((1:kmax2)-1)*dy;

figure(1)
contourf(llon,llat,nbat);
title('Batimetria [m] da regiao');
xlabel('Longitude','LineWidth',14);
ylabel('Latitude', 'LineWidth', 14);
axis equal
colorbar;
% saving figure
out = ['../outputs/ex02/batimetria'];
grafico=['print -djpeg ', out];
eval(grafico);

% criando chaves para mar (1) e terra (0)
A=nbat(1:end-1,1:end-1);
inderr=find(A>=-5);
A(inderr)=NaN;
kmar1=(isnan(A));

B=kmar1;
inderrA=find(kmar1>0);
inderrB=find(kmar1<=0);
B(inderrA)=0;
B(inderrB)=1;
kmar2=B;
A=B';%'
kmar=A;

kmar(1,1:30)=1;
kmar(1:40,1)=1;

% limpar ambiente
clear inderr inderrA inderrB A B kmar2 kmar1

figure(2)
contourf(llon(1:end-1,1:end-1),llat(1:end-1,1:end-1),kmar)
title('Chave da regiao modelada (1 mar, 0 terra)')
xlabel('Longitude')
ylabel('Latitude')
map = [210,180,140; 135,206,250]/255;
colormap(map)
colorbar
% saving figure
out = ['../outputs/ex02/chaves'];
grafico=['print -djpeg ', out];
eval(grafico);


% geracao das matrizes de batimetria para pontos do tipo u,v e eta

%determinacao da batimetria e chaves maritimas em u, v e eta:
% important note: working only with half of the grid (skipping each two points)
bate=kmar(1:2:end-1,2:2:end);     % bathymetry on elevation points
batu=kmar(1:2:end-1,1:2:end-1);   % bathymetry on u points
batv=kmar(2:2:end,2:2:end);       % bathymetry on v points
kmare=bate*0;                     % creating new matrix fill with zeros for eta
kmare(bate>0)=1;                  % positive values for bathymetry => sea
kmaru=batu*0;                     % creating new matrix fill with zeros for u
kmaru(batu>0)=1;
kmarv=batv*0;                     % creating new matrix fill with zeros for v
kmarv(batv>0)=1;
[X,Y]=meshgrid(xgrid,ygrid);      % gridding lat and lon

indterra = find(kmare==0);          % land indexes

% calculo dos coeficientes (constantes) das equacoes discretizadas

denbatu = dens.*batu;
denbatv = dens.*batv;

coef    = 1+rfric*dt*2;
difus2x = Dx/dx/dx;
difus2y = Dy/dy/dy;
dx2     = 2*dx;
dy2     = 2*dy;
dt2     = 2*dt;

% Condicoes iniciais de repouso (0-valores anteriores, 1-atuais, 2-renovados)
eta0=zeros(kmax,jmax);
U0=zeros(kmax,jmax);
V0=zeros(kmax,jmax);
eta1=zeros(kmax,jmax);
U1=zeros(kmax,jmax);
V1=zeros(kmax,jmax);
eta2=zeros(kmax,jmax);
U2=zeros(kmax,jmax);
V2=zeros(kmax,jmax);

%Condicoes de vento e calculo das tensoes de cisalhamento na superficie
dens_ar=1.25;
fric=2.6*1E-3;

uwind(1:kmax,1:jmax)=7.0711;
vwind(1:kmax,1:jmax)=7.0711;
wwind=sqrt(uwind.^2+vwind.^2);
taux=fric*dens_ar.*uwind.*wwind;
tauy=fric*dens_ar.*vwind.*wwind;

vento=sqrt(uwind.^2+vwind.^2);
ventomax=max(max(vento));


figure(3)
contourf(llon(1:end-1,1:end-1),llat(1:end-1,1:end-1),kmar)
map = [210,180,140; 135,206,250]/255;
colormap(map)
colorbar
title('Bathymetry (m)')
xlabel('DISTANCE (m) EW', 'fontsize', 12)
ylabel('DISTANCE (m) NS', 'fontsize', 12)
hold on
quiver(llon(1:skipWind:end,1:skipWind:end),llat(1:skipWind:end,1:skipWind:end),uwind(1:skipWind:end,1:skipWind:end),vwind(1:skipWind:end,1:skipWind:end),'LineWidth',1.5,'color','k')
title(['Wind - Maximum Intensity ',...
    num2str(ventomax),' m/s'],'fontsize',12)
xlabel('Longitude')
ylabel('Latitude')
% saving figure
out = ['../outputs/ex02/windField'];
grafico=['print -djpeg ', out];
eval(grafico);

clc;
kmax=(kmax/2)-1;
jmax=(jmax/2)-1;

kplot=1;

for n=2:nmax
  tempo=n*dt;
  kplot=kplot+1;
  fprintf('Calculating timestep %i/%i - umean=%2.5f, vmean=%2.5f, elevmean=%2.5f\n',...
                             n,nmax,mean(mean(U2)),mean(mean(V2)),mean(mean(eta2)));

  %Desnivel devido ao vento e oscilações no contorno: norte e sul
  for j=1:jmax
      eta2(1,j)=(-etaContorno)*kmare(1,j);
      eta2(2,j)=(-etaContorno)*kmare(2,j);
      eta2(kmax-1,j)=(etaContorno)*kmare(kmax-1,j);
      eta2(kmax,j)=(etaContorno)*kmare(kmax,j);
  end
  for k=1:kmax
      eta2(k,1)=(-etaContorno)*kmare(k,1);
      eta2(k,2)=(-etaContorno)*kmare(k,2);
      eta2(k,jmax)=(etaContorno)*kmare(k,jmax);
      eta2(k,jmax-1)=(etaContorno)*kmare(k,jmax-1);
  end


      %Eq da Continuidade (10.64)
      for j=2:jmax-1
          for k=2:kmax-1
             if kmare(k,j)>0 % se for ponto marítimo
             forcx=(batu(k,j+1).*U1(k,j+1)-batu(k,j).*U1(k,j))/dx;
             forcy=(batv(k,j).*V1(k,j)-batv(k-1,j).*V1(k-1,j))/dy;
             eta2(k,j)=eta0(k,j)-dt2*(forcx+forcy);
             end
          end
      end

    %Eq. do movimento em x
      for j=2:jmax-1
          for k=2:kmax-1
              % ponto marítimo e vizinho de pontos marítimos
              if (kmaru(k,j)*kmare(k,j)*kmare(k,j-1))>0
              vmedu=(V1(k,j)+V1(k,j-1)+V1(k-1,j-1)+V1(k-1,j))/4;
              umed=(U1(k,j-1)+U1(k,j)*2+U1(k,j+1))/4;
              forc=fco*vmedu-g.*(eta1(k,j)-eta1(k,j-1))./dx...
                 +taux(k,j)./denbatu(k,j) ...
                 -umed*(U1(k,j+1)-U1(k,j-1))/dx2-vmedu*(U1(k+1,j)-U1(k-1,j))/dy2 ...
                 +difus2x*(U0(k,j+1)-2*U0(k,j)+U0(k,j-1))+difus2y*(U0(k+1,j)-2*U0(k,j)+U0(k-1,j));
              U2(k,j)=(U0(k,j)+forc*dt2)/coef;
              end
           end
      end

      %Eq. do movimento em y
      for j=2:jmax-1
          for k=2:kmax-1
              % da mesma forma, o ponto deve ser marítimo e seus vizinhos
              % cima/baixo também devem
              if (kmarv(k,j)*kmare(k,j)*kmare(k+1,j))>0
              umedv=(U1(k,j)+U1(k+1,j)+U1(k+1,j+1)+U1(k,j+1))/4;
              vmed=(V1(k-1,j)+2*V1(k,j)+V1(k+1,j))/4;
              forc=-fco*umedv-g.*(eta1(k+1,j)-eta1(k,j))./dy...
                 +tauy(k,j)./denbatv(k,j)...
                 -umedv*(V1(k,j+1)-V1(k,j-1))/dx2-vmed*(V1(k+1,j)-V1(k-1,j))...
                 +difus2x*(V0(k,j+1)-2*V0(k,j)+V0(k,j-1))+difus2y*(V0(k+1,j)-2*V0(k,j)+V0(k-1,j));
              V2(k,j)=(V0(k,j)+forc*dt2)/coef;
              end
          end
      end

  % Condicoes de contorno com extrapolacao linear
  for j=1:jmax
      eta2(1,j)=(2*eta2(2,j)-eta2(3,j))*kmare(1,j);
      eta2(kmax,j)=(2*eta2(kmax-1,j)-eta2(kmax-2,j))*kmare(kmax,j);
      U2(1,j)=(2*U2(2,j)-U2(3,j))*kmaru(1,j);
      V2(1,j)=(2*V2(2,j)-V2(3,j))*kmarv(1,j);
      U2(kmax,j)=(2*U2(kmax-1,j)-U2(kmax-2,j))*kmaru(kmax,j);
      V2(kmax,j)=(2*V2(kmax-1,j)-V2(kmax-2,j))*kmarv(kmax,j);
  end
  for k=1:kmax
      eta2(k,1)=(2*eta2(k,2)-eta2(k,3))*kmare(k,1);
      eta2(k,jmax)=(2*eta2(k,jmax-1)-eta2(k,jmax-2))*kmare(k,jmax);
      U2(k,1)=(2*U2(k,2)-U2(k,3))*kmaru(k,1);
      V2(k,1)=(2*V2(k,2)-V2(k,3))*kmarv(k,1);
      U2(k,jmax)=(2*U2(k,jmax-1)-U2(k,jmax-2))*kmaru(k,jmax);
      V2(k,jmax)=(2*V2(k,jmax-1)-V2(k,jmax-2))*kmarv(k,jmax);
  end

  % Renovando as variaveis no tempo
  eta0=eta1;
  U0=U1;
  V0=V1;
  eta1=eta2;
  U1=U2;
  V1=V2;


  % Plotagem de resultados
  if(kplot==freqplot)
     kplot=0;

        %definir u e v nos pontos tipo eta para plotagem
        uplot=U2;
        vplot=V2;
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

        % Plotando elevacoes e correntes
        velo=sqrt(uplot.^2+vplot.^2);
        veloma=max(velo);
        velomax=max(veloma);
        etama=max(eta2(:,:));
        etami=min(eta2(:,:));
        etamax=max(etama);
        etamin=min(etami);

        figure(6)
        contourf(llon,llat,eta2,'LineWidth',2);
        colorbar;
        title(['Elev (m) - tempo ',num2str(tempo/60),...
              ' min. Limites ',num2str(etamin),' a ',num2str(etamax),' m'],'fontsize',12)
        %axis equal
        %axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
        xlabel('DISTANCIA (m) EW','fontsize',12)
        ylabel('DISTANCIA (m) NS','fontsize',12)
        % print -djpeg fig_elev

        figure(7)
        quiver(llon(1:5:end,1:5:end),llat(1:5:end,1:5:end),uplot(1:5:end,1:5:end),vplot(1:5:end,1:5:end),'LineWidth',2);
        title(['Correntes (m/s) - tempo ',num2str(tempo/60),...
            ' min - intens max ',num2str(velomax),' m/s'],'fontsize',12)
        %axis equal
        %axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
        xlabel('DISTANCIA (m) EW','fontsize',12)
        ylabel('DISTANCIA (m) NS','fontsize',12)
        % print -djpeg fig_corr
     end

  pause(1)

end
