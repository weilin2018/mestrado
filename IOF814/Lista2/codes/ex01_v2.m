clear all; close all; clc

% constantes iniciais do modelo
g       = 9.8 ;             % aceleracao da gravidade
nmax    = 360;              % nivei de tempo de processamento
dx      = 100000;           % espacamento zonal da grade
dy      = 100000;           % espacamento meridional da grade
dt      = 60;               % passo de tempo
freqplot= 10;               % frequencia de outputs
dens    = 1024;             % densidade media da agua do mar

% informacoes dadas pelo exercicio
tau0    = 0.1;              % vento zonal [N/m^2]
Ly      = 4e6;              % comprimento meridional da bacia [m]
H       = 200;              % profundidade media [m]
r       = 0.2;              % coeficiente de friccao no fundo [kg/m^2/s]
beta    = 10e-11;           % beta [s/m]

% informacoes da grade: tamanho, batimetria e outros
kmax    = Ly/dy;            % comprimento da grade em y
jmax    = Ly/dx;            % comprimento da grade em x
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

%% determinacao das condicoes iniciais [em repouso]
eta1    = zeros(kmax,jmax);
U1      = zeros(kmax,jmax);
V1      = zeros(kmax,jmax);
eta2    = zeros(kmax,jmax);
U2      = zeros(kmax,jmax);
V2      = zeros(kmax,jmax);

psi     = ones(kmax,jmax);
taux    = ones(kmax,jmax)*tau0;

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


for n=2:10
    %Eq da Continuidade
    for j=2:jmax-1
        for k=2:kmax-1
           if kmare(k,j)>0
           forcx=(batu(k,j+1).*psi(k,j+1)-batu(k,j).*2*psi(k,j)+batu(k,j-1).*psi(k,j-1))/dx/dx;
           forcy=(batu(k+1,j).*psi(k+1,j)-batu(k,j).*2*psi(k,j)+batu(k-1,j).*psi(k-1,j))/dy/dy;
           eta2(k,j)=forcx+forcy;
           end
        end
    end
    % calcular u como teste
    for j=2:jmax-1
        for k=2:kmax-1
            if (kmaru(k,j)*kmare(k,j)*kmare(k,j-1))>0 % testa o ponto e seus vizinhos
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

    % 0 nos contornos
    eta2(1,:)       =   0;
    eta2(kmax,:)    =   0;
    U2(1,:)         =   0;
    V2(1,:)         =   0;
    U2(kmax,:)      =   0;
    V2(kmax,:)      =   0;

    eta2(:,1)       =   0;
    eta2(:,jmax)    =   0;
    U2(:,1)         =   0;
    V2(:,1)         =   0;
    U2(:,jmax)      =   0;
    V2(:,jmax)      =   0;
%    % Condicoes de contorno com extrapolacao linear
%    for j=1:jmax
%        eta2(1,j)=(2*eta2(2,j)-eta2(3,j))*kmare(1,j);
%        eta2(kmax,j)=(2*eta2(kmax-1,j)-eta2(kmax-2,j))*kmare(kmax,j);
%        U2(1,j)=(2*U2(2,j)-U2(3,j))*kmaru(1,j);
%        V2(1,j)=(2*V2(2,j)-V2(3,j))*kmarv(1,j);
%        U2(kmax,j)=(2*U2(kmax-1,j)-U2(kmax-2,j))*kmaru(kmax,j);
%        V2(kmax,j)=(2*V2(kmax-1,j)-V2(kmax-2,j))*kmarv(kmax,j);
%    end
%    for k=1:kmax
%        eta2(k,1)=(2*eta2(k,2)-eta2(k,3))*kmare(k,1);
%        eta2(k,jmax)=(2*eta2(k,jmax-1)-eta2(k,jmax-2))*kmare(k,jmax);
%        U2(k,1)=(2*U2(k,2)-U2(k,3))*kmaru(k,1);
%        V2(k,1)=(2*V2(k,2)-V2(k,3))*kmarv(k,1);
%        U2(k,jmax)=(2*U2(k,jmax-1)-U2(k,jmax-2))*kmaru(k,jmax);
%        V2(k,jmax)=(2*V2(k,jmax-1)-V2(k,jmax-2))*kmarv(k,jmax);
%    end

    % renovacao das variaveis no tempo
    eta1    = eta2;
    U1      = U2;
    V1      = V2;

end

% definicao de u e v nos mesmos pontos do tipo eta para visualizacao
% basicamente e feita a media entre dois pontos vizinhos da grade
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

quiver(X,Y,uplot*10,vplot*10);
axis([xgrid(1) xgrid(jmax) ygrid(1) ygrid(kmax)]);
