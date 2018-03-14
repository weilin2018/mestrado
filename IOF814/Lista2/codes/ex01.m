clear all; close all; clc

% constantes iniciais do modelo
g       = 9.8 ;             % aceleracao da gravidade
nmax    = 50;               % nivei de tempo de processamento
dx      = 1;                % espacamento zonal da grade
dy      = 1;                % espacamento meridional da grade
dt      = 60;               % passo de tempo
freqplot= 10;               % frequencia de outputs
dens    = 1024;             % densidade media da agua do mar

% informacoes dadas pelo exercicio
Ly      = 4e6;              % comprimento meridional da bacia [m]
tau0    = 0.1;              % vento zonal [N/m^2]
tau     = tau0*cos((pi*11)/Ly);
H       = 200;              % profundidade media [m]
r       = 0.2;              % coeficiente de friccao no fundo [kg/m^2/s]
beta    = 10e-11;           % beta [s/m]

% informacoes da grade: tamanho, batimetria e outros
kmax    = 20;            % comprimento da grade em y
jmax    = 20;            % comprimento da grade em x
kmax2   = kmax*2;           % devido a grade alternada
jmax2   = jmax*2;

bat(1:kmax2,1:jmax2)=H;      % batimetria uniforme

% criacao da batimetria para pontos do tipo eta,u,v
bate=bat(1:2:end-1,1:2:end);
batu=bat(1:2:end-1,1:2:end-1);
batv=bat(1:2:end-1,1:2:end);

% criacao da grade apenas para generalizar o codigo para qualquer batimetria
kmare           = bate*0;
kmare(bate>0)   = 1;
kmaru           = batu*0;
kmaru(batu>0)   = 1;
kmarv           = batv*0;
kmarv(batv>0)   = 1;

xgrid2          = ((1:jmax2)-1)*dx/2;
ygrid2          = ((1:kmax2)-1)*dy/2;
[XBAT,YBAT]     = meshgrid(xgrid2,ygrid2);
xgrid           = ((1:jmax)-1)*dx;
ygrid           = ((1:kmax)-1)*dy;
[X,Y]           = meshgrid(xgrid,ygrid);

% verificar se a batimetria e uniforme ou nao. Se nao, plotar
if max(max(bat))~=min(min(bat))
  figure(1)
  contourf(XBAT,YBAT,bat,'LineWidth',2);
  colorbar;
  title(['Batimetria da regiao modelada (m)'],'fontsize',12)
  axis equal
  axis([xgrid2(1) xgrid2(jmax2) ygrid2(1) ygrid2(kmax2)])
  xlabel('DISTANCIA (m) EW','fontsize',12)
  ylabel('DISTANCIA (m) NS','fontsize',12)
end

%% determinacao das condicoes iniciais [em repouso]
eta1    = zeros(kmax,jmax);
U1      = zeros(kmax,jmax);
V1      = zeros(kmax,jmax);
eta2    = zeros(kmax,jmax);
U2      = zeros(kmax,jmax);
V2      = zeros(kmax,jmax);

psi     = ones(kmax,jmax);
taux    = ones(kmax,jmax)*tau;

% utilizar a equacao da vorticidade discretizada para determinar
% a funcao de corrente e, entao, modelar u,v,eta

for j=2:jmax-1
    for k=2:kmax-1
        % aqui realizamos uma condicao para calcular psi somente em pontos de agua
        if kmare(k,j)>0
            % calculo da unica forcante considerada: vento
            termo1 = 1/(2*dens*H*dy) * (taux(k+1,j+1) + taux(k+1,j-1) + taux(k-1,j+1) + taux(k-1,j-1));
            termo2 = beta/(4*dx) * (psi(k+1,j+1) + psi(k,j+1) - psi(k-1,j-1) - psi(k,j-1));
            termo3 = r/(dens*H) * ( (psi(k,j+1) + psi(k,j-1))/dx/dx + (psi(k+1,j) + psi(k-1,j))/dy/dy );
            termo4 = (2*r)/(dens*H) * (1/dx/dx + 1/dy/dy);
            psi(k,j) = (termo1 + termo2 + termo3)./termo4;
        end
    end
end

%Loop no tempo
kplot=1;
for n=2:nmax
    tempo=n*dt;
    kplot=kplot+1;

    %Eq da Continuidade
    for j=2:jmax-1
        for k=2:kmax-1
           if kmare(k,j)>0
           forcx=(batu(k,j+1).*psi(k,j+1)-batu(k,j).*2*psi(k,j)+batu(k,j-1).*psi(k,j-1))/dx/dx;
           forcy=(batu(k+1,j).*psi(k+1,j)-batu(k,j).*2*psi(k,j)+batu(k-1,j).*psi(k-1,j))/dy/dy;
           eta2(k,j)=(forcx+forcy)/-100;
           end
        end
    end

    %Eq. do movimento em x
    for j=2:jmax-1
        for k=2:kmax-1
            if (kmaru(k,j)*kmare(k,j)*kmare(k,j-1))>0 %o ponto e os vizinhos tem que ser mar (kmares)
            U2(k,j)=-((psi(k+1,j)+psi(k+1,j-1))./2-(psi(k-1,j)+psi(k-1,j-1))./2)./2*dy;
            end
         end
    end

    %Eq. do movimento em y
    for j=2:jmax-1
        for k=2:kmax-1
            if (kmarv(k,j)*kmare(k,j)*kmare(k+1,j))>0 %o ponto e os vizinhos tem que ser mar (kmares)
            V2(k,j)=((psi(k+1,j+1)+psi(k,j+1))./2-(psi(k+1,j-1)+psi(k,j-1))./2)./2*dx;
            end
        end
    end

    eta2(1,:) = 0;
    eta2(kmax,:) = 0;
    U2(1,:) = 0;
    V2(1,:) = 0;
    U2(kmax,:) = 0;
    V2(kmax,:) = 0;

    eta2(:,1) = 0;
    U2(:,1) = 0;
    V2(:,1) = 0;
    eta2(:,jmax) = 0;
    U2(:,jmax) = 0;
    V2(:,jmax) = 0;


    % Renovando as variaveis no tempo
    eta1=eta2;
    U1=U2;
    V1=V2;

end

% como nao ha variacao no tempo da corrente/elevacao, plotamos somente o ultimo instante modelado

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

figure(3)
contourf(X,Y,eta2,'LineWidth',.5);
colorbar;
title(['Elev (m) - tempo ',num2str(tempo/60),...
       ' min. Limites ',num2str(etamin),' a ',num2str(etamax),' m'],'fontsize',12)
%axis equal
axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
xlabel('DISTANCIA (m) EW','fontsize',12)
ylabel('DISTANCIA (m) NS','fontsize',12)
% saving figure
out = ['../outputs/ex01/elevacao_',num2str(tempo)];
grafico=['print -djpeg ', out];
eval(grafico);

figure(4)
quiver(X,Y,uplot,vplot,'LineWidth',1.);
title(['Correntes (m/s) - tempo ',num2str(tempo/60),...
    ' min - intens max ',num2str(velomax),' m/s'],'fontsize',12)
%axis equal
axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)])
xlabel('DISTANCIA (m) EW','fontsize',12)
ylabel('DISTANCIA (m) NS','fontsize',12)
% saving figure
out = ['../outputs/ex01/circulacao_',num2str(tempo)];
grafico=['print -djpeg ', out];
eval(grafico);
