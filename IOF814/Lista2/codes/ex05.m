%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IOF814 - Prof Joseph Harari                                        %
%  Aluno: Danilo Augusto Silva         nusp: 7279456                  %
%                                                                     %
%                         QUESTÃO 05 - LISTA 02                       %
%   Os dados de saidas sao salvos no diretorio ../outputs/Q06         %
%                                                                     %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;clc

% author note:
% uncomment the next line only if you want to save images for report. Otherwise,
% keep this commented to avoid conflict in the system
%genpath('/home/tparente/Dropbox/TOOLBOX_MATLAB/new_codes/m_map/');

%% Parte I - definicao dos parametros do modelo
nmax=90;                % tempo de simulacao
jmax=150;               % tamanho da grade em x
kmax=150;               % tamanho da grade em y
jmax2=jmax/2;           % reduzindo a grade pela metade
kmax2=kmax/2;           % reduzindo a grade pela metade
dt=30;                  % passo de tempo
kx=10;                  % coeficiente de difusao
ky=10;                  % coeficiente de difusao
freqplot=10;            % frequencia de plotagem
boundaryElevation=0.00025;  % open boundary condition of elevation

% parametros especificos para modelagem hidrodinamica
dens=1024;                  % densidade media da agua do mar
latid=25*pi/180;            % latitude (em graus, para rad)
fco=2*7.292E-5*sin(latid);  % parametro de Coriolis
skipWind=10;                % parametro para plotar vento com menor densidade
                            % de vetores
g=9.8;                      % aceleracao da gravidade
rfric=0.02;     % coeficiente de friccao no fundo

% definicao da grade batimetrica baseada no arquivo .xyz em anexo
file = load('../data/etopo1.xyz');
lon=file(:,1);
lat=file(:,2);
bat=file(:,3);

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

% defining a location to control free surface elevation
indLon = 38;
indLat = 55;
control = []; % matrix with nmax dimension to store all data


% interpolar os dados de batimetria para a grade gerada acima
nbat=griddata(lon,lat,bat,llon,llat);

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

fprintf('Generating Figure 1: bathymetry\n');
prange=linspace(batmin,5,40);
batPlot = nbat;
batPlot(batPlot>0)=60;

figure(1)
contourf(llon,llat,batPlot,prange)
title('Bathymetry (m)')
xlabel('DISTANCE (m) EW', 'fontsize', 12)
ylabel('DISTANCE (m) NS', 'fontsize', 12)
colormap('winter')
colorbar
% saving figure
out = ['../outputs/ex05/bathymetry'];
grafico=['print -djpeg ', out];
eval(grafico);


%% Define keys for land (0) and sea (1)
kmar=nbat;
kmar(kmar>0)=0;
kmar(kmar~=0)=1;

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

indterra=find(kmare==0);          % land indexes

denbatu=dens.*batu;
denbatv=dens.*batv;

%% creating matrix for hydrodynamics data

% update kmax,jmax for half of the grid and lon,lat
kmax=kmax2; jmax=jmax2;
llon2=llon(1:2:end,1:2:end);
llat2=llat(1:2:end,1:2:end);

% Initial conditions - rest (1 for currently values and 2 for new values)
eta1=zeros(kmax,jmax);
U1=zeros(kmax,jmax);
V1=zeros(kmax,jmax);
eta2=zeros(kmax,jmax);
U2=zeros(kmax,jmax);
V2=zeros(kmax,jmax);

% Wind conditions and constants
rhoAr=1.25;         % air density
fric=2.6*1E-3;      % drag coefficient

uwind(1:kmax,1:jmax)=7.0711;          % zonal component
vwind(1:kmax,1:jmax)=7.0711;          % meridional component
wwind=sqrt(uwind.^2 + vwind.^2);      % speed
taux=fric*rhoAr.*uwind.*wwind;        % wind stress x
tauy=fric*rhoAr.*vwind.*wwind;        % wind stress y

ventomax=max(max(wwind));

fprintf('Generating figure 2: wind field\n');
figure(2)
contour(llon,llat,batPlot);%,[0.5],'LineWidth',2,'color','k')
title('Bathymetry (m)')
xlabel('DISTANCE (m) EW', 'fontsize', 12)
ylabel('DISTANCE (m) NS', 'fontsize', 12)
hold on
quiver(llon2(1:skipWind:end,1:skipWind:end),llat2(1:skipWind:end,1:skipWind:end),uwind(1:skipWind:end,1:skipWind:end),vwind(1:skipWind:end,1:skipWind:end),'LineWidth',2)
title(['Wind - Maximum Intensity ',...
    num2str(ventomax),' m/s'],'fontsize',12)
%axis equal
xlabel('DISTANCE (m) EW','fontsize',12)
ylabel('DISTANCE (m) NS','fontsize',12)
% saving figure
out = ['../outputs/ex05/windField'];
grafico=['print -djpeg ', out];
eval(grafico);

%kmax=kmax/2;
%jmax=jmax/2;

clc;

%Loop no tempo
kplot=1;
for n=2:nmax
   tempo=n*dt;
   kplot=kplot+1;
   % print to control all process of the program, printing mean u,v and eta
   fprintf('Calculating timestep %i/%i - umean=%2.5f, vmean=%2.5f, elevmean=%2.5f\n',...
                              n,nmax,mean(mean(U2)),mean(mean(V2)),mean(mean(eta2)));

   % open boundary condition: elevation given by user
    etamet=ones(kmax,jmax)*boundaryElevation;

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

    %continuity eq
    for j=2:jmax-1
        for k=2:kmax-1
           if kmare(k,j)>0
           forcx=(batu(k,j+1).*U1(k,j+1)-batu(k,j).*U1(k,j))/dx;
           forcy=(batv(k,j).*V1(k,j)-batv(k-1,j).*V1(k-1,j))/dy;
           eta2(k,j)=eta1(k,j)-dt*(forcx+forcy);
           end
        end
    end

  %momentum eq (x)
    for j=2:jmax-1
        for k=2:kmax-1
            if (kmaru(k,j)*kmare(k,j)*kmare(k,j-1))>0
            vmedu=(V1(k,j)+V1(k,j-1)+V1(k-1,j-1)+V1(k-1,j))/4;
            forc=fco*vmedu-g.*(eta2(k,j)-eta2(k,j-1))./dx...
               +taux(k,j)./denbatu(k,j)-rfric*U1(k,j);
            U2(k,j)=U1(k,j)+forc*dt;
            end
         end
    end

    %momentum eq (y)
    for j=2:jmax-1
        for k=2:kmax-1
            if (kmarv(k,j)*kmare(k,j)*kmare(k+1,j))>0
            umedv=(U2(k,j)+U2(k+1,j)+U2(k+1,j+1)+U2(k,j+1))/4;
            forc=-fco*umedv-g.*(eta2(k+1,j)-eta2(k,j))./dy...
               +tauy(k,j)./denbatv(k,j)-rfric*V1(k,j);
            V2(k,j)=V1(k,j)+forc*dt;
            end
        end
    end

    % boundary condition - linear extrapolation
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

    % updating variables in time
    eta1=eta2;
    U1=U2;
    V1=V2;

    % save eta to control
    control = [control, eta2(indLon,indLat)];

    % plotting results
     if(kplot==freqplot)
        kplot=0;

        %define u and v in eta points
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

        % plotting elevation and current
        velo=sqrt(uplot.^2+vplot.^2);
        veloma=max(velo);
        velomax=max(veloma);
        etama=max(eta2(:,:));
        etami=min(eta2(:,:));
        etamax=max(etama);
        etamin=min(etami);

        figure(5)
        contourf(llon2,llat2,eta2,'LineWidth',0.1);
        colorbar;
        title(['Elev (m) - time ',num2str(tempo/60),...
              ' min. Limits ',num2str(etamin),' a ',num2str(etamax),' m'],'fontsize',12)
        %axis equal
        %axis([llon(1) llon(jmax) llat(1) llat(kmax)])
        xlabel('DISTANCE (m) EW','fontsize',12)
        ylabel('DISTANCE (m) NS','fontsize',12)
        % saving figure
        out = ['../outputs/ex05/elev_',num2str(tempo)];
        grafico=['print -djpeg ', out];
        eval(grafico);

        figure(6)
        quiver(X(1:5:end,1:5:end)',Y(1:5:end,1:5:end)',uplot(1:5:end,1:5:end),vplot(1:5:end,1:5:end),'LineWidth',2);
        title(['Currents (m/s) - time ',num2str(tempo/60),...
            ' min - intens max ',num2str(velomax),' m/s'],'fontsize',12)
        axis equal
        axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
        xlabel('DISTANCE (m) EW','fontsize',12)
        ylabel('DISTANCE (m) NS','fontsize',12)
        % saving figure
        out = ['../outputs/ex05/curr_',num2str(tempo)];
        grafico=['print -djpeg ', out];
        eval(grafico);
     end

     pause(1)

end

% plot eta control in all simulation
figure(7)
plot(control,'LineWidth',2);
xlabel('Timestep (dt) in seconds');
ylabel('Sea Surface Elevation in meters ');
title(['Time Series of Sea Surface Elevation at ', num2str(llon2(indLon,indLat)), ', ', num2str(llat2(indLon,indLat))]);
%title('Time Series of Sea Surface Elevation');
grid
% saving figure
out = ['../outputs/ex05/elevationTimeseries'];
grafico=['print -djpeg ', out];
eval(grafico);
