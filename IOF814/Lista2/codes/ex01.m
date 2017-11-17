%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IOF814 - Prof Joseph Harari                                        %
%  Aluno: Danilo Augusto Silva         nusp: 7279456                  %
%                                                                     %
%                         QUESTÃO 01 - LISTA 02                       %                                                                  %
%   Os dados de saidas sao salvos no diretorio ../outputs/Q01         %
%                                                                     %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;clc

g=9.8;          % aceleracao da gravidade
nmax=360;       % niveis de tempo
dt=60;          % passo de tempo
freqplot=10;    % frequencia de plotagem

% definir variaveis
Ly = 40000000;                  % comprimento meridional [m]
Lx = Ly;                        % comprimento zonal da bacia [m]
H  = 200;                       % profundidade media [m]
B  = 10e-11;                    % beta [m^-1s]
r  = 0.2;                       % taxa de decaimento [kg/m²/s]
To = 0.1;                       % tau0 [N/m^2]
rho= 1024;                      % densidade do oceano

%Tox = To*cos((pi*y)/Ly);          % vento zonal com variacao meridional

% criacao da grade
dx = 100000;                    % espacamento da grade em x
dy = 100000;                    % espacamento da grade em y

kmax = Ly/dy;                   % comprimento da grade em y
jmax = Lx/dx;                   % comprimento da grade em x
jmax2=jmax*2;                   % devido a grade alternada
kmax2=kmax*2;

xgrid=((1:jmax)-1)*dx;          % quantidade de pontos em cada direcao
ygrid=((1:kmax)-1)*dy;
[X,Y]=meshgrid(xgrid,ygrid);    % gridar

% geracao das matrizes de armazenamento de dados
eta1 = zeros(kmax,jmax);
U1   = zeros(kmax,jmax);
V1   = zeros(kmax,jmax);
eta2 = zeros(kmax,jmax);
U2   = zeros(kmax,jmax);
V2   = zeros(kmax,jmax);
psi  = zeros(kmax,jmax);

% calcular a tensao de cisalhamento variando conforme a posicao em y
for k=1:1:kmax
    for j=1:1:jmax
        taux(k,j) = To*cos(pi*k/Ly);
    end
end

% calcular psi usando a equacao discretizada no exercicio
for j=2:jmax-1
  for k=2:kmax-1
    termo1 = 1/(2*rho*H*dy) * (taux(k+1,j+1) + taux(k+1,j-1) + taux(k-1,j+1) + taux(k-1,j-1));
    termo2 = B/(4*dx) * (psi(k+1,j+1) + psi(k,j+1) - psi(k-1,j-1) - psi(k,j-1));
    termo3 = r/(rho*H) * ( (psi(k,j+1) + psi(k,j-1))/dx/dx + (psi(k+1,j) + psi(k-1,j))/dy/dy );
    termo4 = (2*r)/(rho*H) * (1/dx/dx + 1/dy/dy);
    psi(k,j) = (termo1 + termo2 + termo3)./termo4;
  end
end

psi0=psi;

%Loop no tempo
kplot=1;
for n=2:nmax
   tempo=n*dt;
   kplot=kplot+1;

    %Eq da Continuidade
    for j=2:jmax-1
        for k=2:kmax-1
           forcx=(psi0(k,j+1)-2*psi0(k,j)+psi0(k,j-1))/dx/dx;
           forcy=(psi0(k+1,j)-2*psi0(k,j)+psi0(k-1,j))/dy/dy;
           eta2(k,j)=forcx+forcy;
        end
    end

  %Eq. do movimento em x
    for j=2:jmax-1
        for k=2:kmax-1
            U2(k,j)=-((psi0(k+1,j)+psi0(k+1,j-1))./2-(psi0(k-1,j)+psi0(k-1,j-1))./2)./2*dy;
         end
    end

    %Eq. do movimento em y
    for j=2:jmax-1
        for k=2:kmax-1
            V2(k,j)=((psi0(k+1,j+1)+psi0(k,j+1))./2-(psi0(k+1,j-1)+psi0(k,j-1))./2)./2*dx;
        end
    end

% Condicoes de contorno com extrapolacao linear
for j=1:jmax
    eta2(1,j)=(2*eta2(2,j)-eta2(3,j));
    eta2(kmax,j)=(2*eta2(kmax-1,j)-eta2(kmax-2,j));
    U2(1,j)=(2*U2(2,j)-U2(3,j));
    V2(1,j)=(2*V2(2,j)-V2(3,j));
    U2(kmax,j)=(2*U2(kmax-1,j)-U2(kmax-2,j));
    V2(kmax,j)=(2*V2(kmax-1,j)-V2(kmax-2,j));
end
for k=1:kmax
    eta2(k,1)=(2*eta2(k,2)-eta2(k,3));
    eta2(k,jmax)=(2*eta2(k,jmax-1)-eta2(k,jmax-2));
    U2(k,1)=(2*U2(k,2)-U2(k,3));
    V2(k,1)=(2*V2(k,2)-V2(k,3));
    U2(k,jmax)=(2*U2(k,jmax-1)-U2(k,jmax-2));
    V2(k,jmax)=(2*V2(k,jmax-1)-V2(k,jmax-2));
end

% Renovando as variaveis no tempo
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
        uplot(k,j)=(U2(k,j)+U2(k,j+1))/2;
   end
end
for j=1:jmax
   for k=2:kmax
        vplot(k,j)=(V2(k,j)+V2(k-1,j))/2;
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
quiver(X,Y,uplot,vplot,'LineWidth',2);
title(['Correntes (m/s) - tempo ',num2str(tempo/60),...
    ' min - intens max ',num2str(velomax),' m/s'],'fontsize',12)
axis equal
axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
xlabel('DISTANCIA (m) EW','fontsize',12)
ylabel('DISTANCIA (m) NS','fontsize',12)
% print -djpeg fig_corr
   end

pause(1)

end
