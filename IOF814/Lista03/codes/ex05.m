clear all; close all; clc

% definindo espacamento (que faz mais sentido)

%% controle de plotagem, para facilitar a vida durante desenvolvimento.
% 0 = True: plotar
% 1 = False: nao plotar
DEPLOY=1;

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
nmax      = 1440;             % numero de passos de tempo maximo
dt        = 30;               % passo de tempo (segundos)
dx        = 1000;             % espacamento de grade em x
dy        = 1000;             % espacamento de grade em y
dz        = 1;
kx        = 10;               % coeficiente de difusao
ky        = 10;               % coeficiente de difusao
freqplot  = 10;               % frequencia de plotagem
freqForc  = (3*60*60)/dt;          % frequencia de atualizada da forcante (a cada 3h)
difus     = 0.1;              % coeficiente de difusao
difz      = 0.001;            % coeficiente de difusao na vertical

% parametros baseado na grade batimetrica utilizada
jmax      = round((abs(min(lon) - max(lon)))/dx);    % nro de pontos da grade em x
kmax      = round((abs(min(lat) - max(lat)))/dy);    % nro de pontos da grade em y
lmax      = abs(min(bat));                           % nro de pontos da grade em z
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
ampl      = 1.0;                   % elevacao nos contornos abertos
period    = 3600;

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
kmar = -kmar;

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

dt2=dt*2;
dx2=dx*2;
dy2=dy*2;
difus2x=difus/dx/dx;
difus2y=difus/dy/dy;
difzdz=difz/dz;
denbatu=dens.*batu;
denbatv=dens.*batv;

amplit=ones(kmax,jmax)*ampl;
omega=2*pi/period;

% Condicoes iniciais de repouso (0 - valores anteriores, 1 - valores atuais, 2 - renovados)
eta0=zeros(kmax,jmax);
u0=zeros(kmax,jmax,lmax);
v0=zeros(kmax,jmax,lmax);
eta1=zeros(kmax,jmax);
u1=zeros(kmax,jmax,lmax);
v1=zeros(kmax,jmax,lmax);
eta2=zeros(kmax,jmax);
u2=zeros(kmax,jmax,lmax);
v2=zeros(kmax,jmax,lmax);

%Condicoes de vento e calculo das tensoes de cisalhamento na superficie
dens_ar=1.25;
fric=2.6*1E-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   comeca com mare + vento de sul
% após 10800s: nova mare com vento de sudeste
% após +10800s: nova mare com vento de leste
% após +10800s: nova mare com vento de sudoeste
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% variaveis para controlar as forcantes do modelo durante o processamento
nforc = 0;            % controlador para especificar o tipo de forcante (vento)
contf = freqForc-1;   % contador para inserir novo vento e mare

load wind_data.dat;
etamet=ones(kmax,jmax)*ampl; % elevacao de 1m
uwind(1:kmax,1:jmax)=wind_data(1,1);
vwind(1:kmax,1:jmax)=wind_data(1,2);
wwind=sqrt(uwind.^2+vwind.^2);
taux=fric*dens_ar.*uwind.*wwind;
tauy=fric*dens_ar.*vwind.*wwind;

%Loop no tempo
kplot=1;
for n=1:nmax
   tempo=n*dt;
   kplot=kplot+1;
   contf=contf+1;
   % print to control all process of the program, printing mean u,v and eta
   fprintf('Calculando timestep %i/%i  \n\n\n',n,nmax);

   % teste para usar as primeiras forcantes
   if contf==freqForc+1
      contf = 0;
      nforc = nforc +1;

      etamet=ones(kmax,jmax)*ampl; % elevacao de 1m
      uwind(1:kmax,1:jmax)=wind_data(nforc,1);
      vwind(1:kmax,1:jmax)=wind_data(nforc,2);
      wwind=sqrt(uwind.^2+vwind.^2);
      taux=fric*dens_ar.*uwind.*wwind;
      tauy=fric*dens_ar.*vwind.*wwind;

     % plotando os ventos
      ventomax=max(max(wwind));

      figure(2)
      quiver(X,Y,uwind,vwind,'LineWidth',2)
      title(['Vento - intensidade maxima ',...
          num2str(ventomax),' m/s'],'fontsize',12)
      axis equal
      axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
      xlabel('DISTANCIA (m) EW','fontsize',12)
      ylabel('DISTANCIA (m) NS','fontsize',12)

      % desnivel na borda devido ao vento e oscilacoes periodicas no contorno
      for j=1:jmax
         eta2(1,j)=(etamet(1,j)+amplit(1,j)*cos(omega*tempo))*kmare(1,j);
         eta2(2,j)=(etamet(2,j)+amplit(2,j)*cos(omega*tempo))*kmare(2,j);
         eta2(kmax-1,j)=(etamet(kmax-1,j)+amplit(kmax-1,j)*cos(omega*tempo))*kmare(kmax-1,j);
         eta2(kmax,j)=(etamet(kmax,j)+amplit(kmax,j)*cos(omega*tempo))*kmare(kmax,j);
      end
      for k=1:kmax
         eta2(k,1)=(etamet(k,1)+amplit(k,1)*cos(omega*tempo))*kmare(k,1);
         eta2(k,2)=(etamet(k,2)+amplit(k,2)*cos(omega*tempo))*kmare(k,2);
         eta2(k,jmax-1)=(etamet(k,jmax-1)+amplit(k,jmax-1)*cos(omega*tempo))*kmare(k,jmax-1);
         eta2(k,jmax)=(etamet(k,jmax)+amplit(k,jmax)*cos(omega*tempo))*kmare(k,jmax);
      end
    end

      %Eq da Continuidade: 13.49
      for j=2:jmax-1
          for k=2:kmax-1
             if kmare(k,j)>0
                 forcx=0;
                 forcy=0;
                 for l=2:lmax-1
                      forcx=(u1(k,j+1,l)-u1(k,j,l))*dz/dx+forcx;
                      forcy=(v1(k,j,l)-v1(k-1,j,l))*dz/dy+forcy;
                 end
                 eta2(k,j)=eta0(k,j)-dt2*(forcx+forcy);
             end
          end
      end

      % filtragem
      etafil=eta1;
      etafil=(eta0+2*eta1+eta2)/4.*kmare;
      eta1=etafil;

    %Eq. do movimento em x
      for j=2:jmax-1
          for k=2:kmax-1
              if (kmaru(k,j)*kmare(k,j)*kmare(k,j-1))>0
                bat3=[batu(k,j) bate(k,j) bate(k,j-1)];
                batmin=min(bat3);
                lfim=batmin/dz-1;
                for l=2:lfim
                    dudzinf=(u0(k,j,l+1)-u0(k,j,l))*difzdz;
                    if (l==2) dudzsup=-taux(k,j)/dens; end
                    if (l~=2) dudzsup=(u0(k,j,l)-u0(k,j,l-1))*difzdz; end
                    umed=(u1(k,j-1,l)+u1(k,j,l)*2+u1(k,j+1,l))/4;
                    vmedu=(v1(k,j,l)+v1(k,j-1,l)+v1(k-1,j-1,l)+v1(k-1,j,l))/4;
                    forc=fco*vmedu-g.*(eta1(k,j)-eta1(k,j-1))./dx...
                     +(dudzinf-dudzsup)/dz-rfric*u0(k,j,l)+...
                     -umed*(u1(k,j+1,l)-u1(k,j-1,l))/dx2-vmedu*(u1(k+1,j,l)-u1(k-1,j,l))/dy2+...
                     difus2x*(u0(k,j+1,l)-u0(k,j,l)*2+u0(k,j-1,l))+...
                     difus2y*(u0(k+1,j,l)-u0(k,j,l)*2+u0(k-1,j,l));
                    u2(k,j,l)=u0(k,j,l)+forc*dt2;
                end
              end
           end
      end

      % filtragem u
      ufil=u1;
      for j=2:jmax-1
          for k=2:kmax-1
              ufil(k,j,:)=(u0(k,j,:)+2*u1(k,j,:)+u2(k,j,:))/4.*kmaru(k,j);
          end
      end

      u1=ufil;

      %Eq. do movimento em y
      for j=2:jmax-1
          for k=2:kmax-1
              if (kmarv(k,j)*kmare(k,j)*kmare(k+1,j))>0
                bat3=[batv(k,j) bate(k,j) bate(k+1,j)];
                batmin=min(bat3);
                lfim=batmin/dz-1;
                for l=2:lfim
                  dvdzinf=(v0(k,j,l+1)-v0(k,j,l))*difzdz;
                  if (l==2) dvdzsup=-tauy(k,j)/dens; end
                  if (l~=2) dvdzsup=(v1(k,j,l)-v0(k,j,l-1))*difzdz; end
                  vmed=(v1(k-1,j,l)+v1(k,j,l)*2+v1(k+1,j,l))/4;
                  umedv=(u2(k,j,l)+u2(k+1,j,l)+u2(k+1,j+1,l)+u2(k,j+1,l))/4;
                  forc=-fco*umedv-g.*(eta1(k+1,j)-eta1(k,j))./dy...
                     +(dvdzinf-dvdzsup)/dz-rfric*v0(k,j,l)+...
                   -umedv*(v1(k,j+1,l)-v1(k,j-1,l))/dx2-vmed*(v1(k+1,j,l)-v1(k-1,j,l))/dy2+...
                     difus2x*(v0(k,j+1,l)-v0(k,j,l)*2+v0(k,j-1,l))+...
                     difus2y*(v0(k+1,j,l)-v0(k,j,l)*2+v0(k-1,j,l));
                  v2(k,j,l)=v0(k,j,l)+forc*dt2;
                end
              end
          end
      end

      % filtragem de v
      vfil=v1;
      for j=2:jmax-1
          for k=2:kmax-1
              vfil(k,j,:)=(v0(k,j,:)+2*v1(k,j,:)+v2(k,j,:))/4.*kmarv(k,j);
          end
      end

      v1=vfil;

  % Condicoes de contorno com extrapolacao linear
  for j=1:jmax
      eta2(1,j)=(2*eta2(2,j)-eta2(3,j))*kmare(1,j);
      eta2(kmax,j)=(2*eta2(kmax-1,j)-eta2(kmax-2,j))*kmare(kmax,j);
      u2(1,j,2:lmax-1)=(2*u2(2,j,2:lmax-1)-u2(3,j,2:lmax-1))*kmaru(1,j);
      v2(1,j,2:lmax-1)=(2*v2(2,j,2:lmax-1)-v2(3,j,2:lmax-1))*kmarv(1,j);
      u2(kmax,j,2:lmax-1)=(2*u2(kmax-1,j,2:lmax-1)-u2(kmax-2,j,2:lmax-1))*kmaru(kmax,j);
      v2(kmax,j,2:lmax-1)=(2*v2(kmax-1,j,2:lmax-1)-v2(kmax-2,j,2:lmax-1))*kmarv(kmax,j);
  end
  for k=1:kmax
      eta2(k,1)=(2*eta2(k,2)-eta2(k,3))*kmare(k,1);
      eta2(k,jmax)=(2*eta2(k,jmax-1)-eta2(k,jmax-2))*kmare(k,jmax);
      u2(k,1,2:lmax-1)=(2*u2(k,2,2:lmax-1)-u2(k,3,2:lmax-1))*kmaru(k,1);
      v2(k,1,2:lmax-1)=(2*v2(k,2,2:lmax-1)-v2(k,3,2:lmax-1))*kmarv(k,1);
      u2(k,jmax,2:lmax-1)=(2*u2(k,jmax-1,2:lmax-1)-u2(k,jmax-2,2:lmax-1))*kmaru(k,jmax);
      v2(k,jmax,2:lmax-1)=(2*v2(k,jmax-1,2:lmax-1)-v2(k,jmax-2,2:lmax-1))*kmarv(k,jmax);
  end

  % Renovando as variaveis no tempo
  eta0=eta1;
  eta1=eta2;
  u0=u1;
  u1=u2;
  v0=v1;
  v1=v2;

  if kplot==freqplot
    kplot = 0;

   %definir u e v nos pontos tipo eta para plotagem
   uplot=u2;
   vplot=v2;
   for j=1:jmax-1
      for k=1:kmax
         if kmare(k,j)>0
           uplot(k,j,2:lmax-1)=(u2(k,j,2:lmax-1)+u2(k,j+1,2:lmax-1))/2;
         end
      end
   end
   for j=1:jmax
      for k=2:kmax
         if kmare(k,j)>0
           vplot(k,j,2:lmax-1)=(v2(k,j,2:lmax-1)+v2(k-1,j,2:lmax-1))/2;
         end
      end
   end

   % Plotando elevacoes e correntes
   etama=max(eta2(:,:));
   etami=min(eta2(:,:));
   etamax=max(etama);
   etamin=min(etami);
   velo=sqrt(uplot.^2+vplot.^2);
   velomax2=max(max(velo(:,:,2)));
   velomax3=max(max(velo(:,:,3)));
   velomax7=max(max(velo(:,:,7)));
   velomax8=max(max(velo(:,:,8)));

   figure(3)
   contourf(X,Y,eta2,'LineWidth',2);
   colorbar;
   title(['Elev (m) - tempo ',num2str(tempo/60),...
         ' min. Limites ',num2str(etamin),' a ',num2str(etamax),' m'],'fontsize',12)
   axis equal
   axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
   xlabel('DISTANCIA (m) EW','fontsize',12)
   ylabel('DISTANCIA (m) NS','fontsize',12)
   % print -djpeg fig_elev

     figure(4)
     subplot(1,2,1)
     quiver(X,Y,uplot(:,:,6),vplot(:,:,6));
     title(['Prof 05m - max ',num2str(velomax2),' m/s'])
     axis equal
     axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
     xlabel('DISTANCIA (m) EW','fontsize',12)
     ylabel('DISTANCIA (m) NS','fontsize',12)
     subplot(1,2,2)
     quiver(X,Y,uplot(:,:,11),vplot(:,:,11));
     title(['Prof 10m - max ',num2str(velomax3),' m/s'])
     axis equal
     axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
     xlabel('DISTANCIA (m) EW','fontsize',12)
     ylabel('DISTANCIA (m) NS','fontsize',12)

  end

  pause(0.1)
end
