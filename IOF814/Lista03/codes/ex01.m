%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IOF814 - Prof Joseph Harari                                        %
%  Aluno: Danilo Augusto Silva         nusp: 7279456                  %
%                                                                     %
%                         QUESTÃO 01 - LISTA 03                       %                                                                  %
%   Os dados de saidas sao salvos no diretorio ../outputs/Q01         %
%                                                                     %
%                             Modelo 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;clc

%% Parametros numericos

dx=1000;                                    %espacamento em x em metros
dy=1000;                                    %espacamento em y em metros
dy2=dy*2;
dx2=dx*2;
dt=30;                                      %intervalo de tempo em segundos
dz=1;                                       % espacamento de grade em z
difus=0.1;                                  % coeficiente de difusao
difz=0.001;                                 % coeficiente de difusao na vertical
kmax = 11;                                  % definicao de kmax (em y)
jmax = 11;                                  % definicao de jmax (em x)
jmax2=jmax*2;
kmax2=kmax*2;
lmax=100;
nmax=360;                                   %passos de tempo
rfric=0.02;                                 %coeficiente de friccao no fundo
freqplot=20;                                % frequencia de plotagem
concorte=0.0001;
xgrid=((1:jmax)-1)*dx;
ygrid=((1:kmax)-1)*dy;
dens=1024;                                  %densidade media da agua do mar
latid=45*pi/180;                            %latitude(de graus para rad)
fco=2*7.29E-5*sin(latid);                   %parametro de coriolis
g=9.8;
H=100;                                      % profundidade media

%% batimetria e grade
bat(1:kmax2,1:jmax2) = H;

% matrizes de batimetria para os 3 tipos de grades
bate=bat(1:2:end-1,2:2:end);
batu=bat(1:2:end-1,1:2:end-1);
batv=bat(2:2:end,2:2:end);

[X,Y]=meshgrid(xgrid,ygrid);

%% constantes de discretizacao
difzdz=difz/dz;
dens2=dens;

%% definicao das condicoes iniciais: 2 niveis no tempo
eta1=zeros(kmax,jmax);
u1=zeros(kmax,jmax,lmax);
v1=zeros(kmax,jmax,lmax);
eta2=zeros(kmax,jmax);
u2=zeros(kmax,jmax,lmax);
v2=zeros(kmax,jmax,lmax);

%Condicoes de vento e calculo das tensoes de cisalhamento na superficie
dens_ar=1.025;
fric=2.6*1E-3;
uwind(1:kmax,1:jmax)=-5.65682;
vwind(1:kmax,1:jmax)=-5.65682;
wwind=sqrt(uwind.^2+vwind.^2);
taux=fric*dens_ar.*uwind.*wwind;
tauy=fric*dens_ar.*vwind.*wwind;

vento=sqrt(uwind.^2+vwind.^2);
ventoma=max(vento);
ventomax=max(ventoma);

%plotando os ventos
figure(1)
quiver(X,Y,uwind,vwind,'LineWidth',2)
title(['Vento - intensidade maxima ',...
    num2str(ventomax),' m/s'],'fontsize',12)
axis equal
axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
xlabel('DISTANCIA (m) EW','fontsize',12)
ylabel('DISTANCIA (m) NS','fontsize',12)
%kfig=kfig+1; eval(['print -dpng figuras1a/fig_vento_' num2str(kfig)])

%Loop no tempo
kplot=2;
kfig=0;


for n=3:nmax
   tempo=n*dt;
   kplot=kplot+1;

   %Eq. do movimento em x
    for j=2:jmax-1
        for k=2:kmax-1
            batmin=batu;
            for l=2:lmax-1
            dudzinf=(u1(k,j,l+1)-u1(k,j,l))*difzdz;
            if (l==2) dudzsup=-taux(k,j)/dens; end
            if (l~=2) dudzsup=(u1(k,j,l)-u1(k,j,l-1))*difzdz; end
            vmedu=(v1(k,j,l)+v1(k,j-1,l)+v1(k-1,j-1,l)+v1(k-1,j,l))/4;
            forc=fco*vmedu+(dudzinf-dudzsup)/dz;
            u2(k,j,l)=u1(k,j,l)+forc*dt;
            end
            end
         end

    %Eq. do movimento em y
    for j=2:jmax-1
        for k=2:kmax-1

            for l=2:lmax-1
            dvdzinf=(v1(k,j,l+1)-v1(k,j,l))*difzdz;
            if (l==2) dvdzsup=-tauy(k,j)/dens; end
            if (l~=2) dvdzsup=(v1(k,j,l)-v1(k,j,l-1))*difzdz; end
            umedv=(u2(k,j,l)+u2(k+1,j,l)+u2(k+1,j+1,l)+u2(k,j+1,l))/4;

            forc=-fco*umedv+(dvdzinf-dvdzsup)/dz;
            v2(k,j,l)=v1(k,j,l)+forc*dt;
            end
            end
        end

% Condicoes de contorno com extrapolacao linear
for j=1:jmax
    u2(1,j,2:lmax-1)=(2*u2(2,j,2:lmax-1)-u2(3,j,2:lmax-1));
    v2(1,j,2:lmax-1)=(2*v2(2,j,2:lmax-1)-v2(3,j,2:lmax-1));
    u2(kmax,j,2:lmax-1)=(2*u2(kmax-1,j,2:lmax-1)-u2(kmax-2,j,2:lmax-1));
    v2(kmax,j,2:lmax-1)=(2*v2(kmax-1,j,2:lmax-1)-v2(kmax-2,j,2:lmax-1));
end
for k=1:kmax
    u2(k,1,2:lmax-1)=(2*u2(k,2,2:lmax-1)-u2(k,3,2:lmax-1));
    v2(k,1,2:lmax-1)=(2*v2(k,2,2:lmax-1)-v2(k,3,2:lmax-1));
    u2(k,jmax,2:lmax-1)=(2*u2(k,jmax-1,2:lmax-1)-u2(k,jmax-2,2:lmax-1));
    v2(k,jmax,2:lmax-1)=(2*v2(k,jmax-1,2:lmax-1)-v2(k,jmax-2,2:lmax-1));
end

% Renovando as variaveis no tempo
u1=u2;
v1=v2;

% Plotagem de resultados
   if(kplot==freqplot)
   kplot=0;

%definir u e v  para plotagem
uplot=u2;
vplot=v2;


% Plotando elevacoes e correntes na superficie, a 5 metros, a 30 metros e a
% 50 metros
velo=sqrt(uplot.^2+vplot.^2);
velomax2=max(max(velo(:,:,2)));
velomax6=max(max(velo(:,:,6)));
velomax31=max(max(velo(:,:,31)));
velomax51=max(max(velo(:,:,51)));


if(n>freqplot) close (1); end
figure(1)
subplot(2,2,1)
quiver(X,Y,uplot(:,:,2),vplot(:,:,2));
title(['Prof 1m - max ',num2str(velomax2),' m/s, tempo ',num2str(tempo/60),...
      ' min.'])
axis equal
axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
xlabel('DISTANCIA (m) EW','fontsize',12)
ylabel('DISTANCIA (m) NS','fontsize',12)
subplot(2,2,2)
quiver(X,Y,uplot(:,:,6),vplot(:,:,6));
title(['Prof 05m - max ',num2str(velomax6),' m/s'])
axis equal
axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
xlabel('DISTANCIA (m) EW','fontsize',12)
ylabel('DISTANCIA (m) NS','fontsize',12)
subplot(2,2,3)
quiver(X,Y,uplot(:,:,31),vplot(:,:,31));
title(['Prof 30m - max ',num2str(velomax31),' m/s'])
axis equal
axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
xlabel('DISTANCIA (m) EW','fontsize',12)
ylabel('DISTANCIA (m) NS','fontsize',12)
subplot(2,2,4)
quiver(X,Y,uplot(:,:,51),vplot(:,:,51));
title(['Prof 50m - max ',num2str(velomax51),' m/s'])
axis equal
axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
xlabel('DISTANCIA (m) EW','fontsize',12)
ylabel('DISTANCIA (m) NS','fontsize',12)

   end

pause(1)

end
