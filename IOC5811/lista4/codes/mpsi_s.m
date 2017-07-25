
clear;hold off;
close all

%*******************************************************
%                    PROGRAM MPSI_S
%
% calcula numericamente a solução de Stommel para 
% campos de rotacional do vento e de topografia 
% de fundo teoricos especificados
%
% resolução por esquema iterativo em forma nao-dimensional
% bacia oceanica quadrada de dimensao L=1e6
%
%
%             por Ilson da Silveira 
%
%******************************************************


tol=1e-3;  % tolerancia para erro medio durante a iteracao
 
% parametros do modelo

L=1e6; 
H=5e3;     
f0=7.3e-5;
beta=2e-11;
hE=sqrt(2e-2/f0) % espessura da camada de Ekman teórica [Av= 1e-2 m2/s]
r= f0*hE/H;
tau0=1e-4        % tau é a tensao de cisalhamento/densidade [m2/s2] 

% a nao-dimensionalizacao é feita assumindo que a escala da
% funcao de corrente geostrofica psi é dada pela rel. de Sverdrup:
%                psi= tau0/(rho*H*beta)



r=r/(beta*L);     % parametro de friccao normalizado
ff=f0/(beta*L);   % parametro de Coriolis normalizado


% constroi grade

nmax=101;
dx=0.01; 

nsel1=[1:nmax-2];
nsel2=[3:nmax];
nsel=[2:nmax-1];


x=[0:nmax-1]*dx;
y=x';
[xg,yg]=meshgrid(x,y);
xg=flipud(xg);
yg=flipud(yg);

yL=max(y);

% campo de rotacional do vento nao-dimensional

%curl=zeros(size(xg));                    % sem vento
curl=-pi/yL*sin(pi*yg./yL);              % p/giro unico
%curl=-pi/yL*sin(2*pi*yg./yL);            % p/giros duplos 
%curl=-pi/yL*exp(-(xg-yL/2).*(xg-yL/2));  % campo gaussiano
%curl(:,1:51)=-pi/yL*ones(101,51);        % campo em funcao de Heaviside 
%curl=-pi/yL*ones(size(xg));              % campo uniforme

% campo de topografia de fundo nao-dimensional-- assumido H=5000 m. 
% ou seja hb=b/H


hb=zeros(size(xg));                                                    % fundo plano
%hb=-0.275*yg;                                                         % fundo inclinado linearmente em y;
%hb=-3*0.275*xg;                                                       % fundo inclinado linearmente em x;
%hb=0.1*exp(-30*(xg-yL/2).*(xg-yL/2));                                 % cordilheira meso-atlantica
%hb= 0.1*exp(-30*(xg-yL/2).*(xg-yL/2)).*exp(-30*(yg-yL/2).*(yg-yL/2)); % banco submarino


% constroi campo de gradiente de VP ambiente 

fac=1.0                     % **use fac=0.0 para solucao no plano f **
                            % caso contraio, use fac=1.0
 
y0=0.5;                     % latitude central
B=fac*(yg-y0)+ff*hb;        % funcao de VP ambiente

[Bx,By]=gradient(B,dx,dx);  % calculo do gradiente de VP ambiente 
By=-By;                     % ambient PV


gray2=jet(256);
%gray2=gray2(50:256-50,:);
colormap(gray2);

% plotando VP ambiente

	    b1=0.1*floor(min(min(B))*10);
            b2=0.1*ceil(max(max(B))*10);

 lb=b1:.05:b2; 
 figure(1) 
 cs=contourf(x,y,flipud(B),lb,'w'); colorbar
% clabel(cs,lp)
 axis('square')
 xlabel('x [em 10^6 m]')
 ylabel('y [em 10^6 m]')
 title('A Vorticidade Potencial Ambiente')

%pause

% plotando topografia

% em 2D 

 lt=-1:.05:1; 
 figure(2) 

 subplot(121)
 cs=contour(x,y,flipud(hb),lt,'k'); %colorbar
% clabel(cs,lp)
 axis('square')
 xlabel('x [em 10^6 m]')
 ylabel('y [em 10^6 m]')
 title('A topografia de fundo em 2D')

% em 3D

 subplot(122)
 cs=mesh(x,y,4e2*flipud(hb)); % colorbar
% clabel(cs,lp)
 axis('square')
 xlabel('x [em 10^6 m]')
 ylabel('y [em 10^6 m]')
 zlabel('z [em m]')
 title('A topografia de fundo em 3D')

%pause

% condicoes de contorno e iniciais

psi=zeros(size(xg));


% Processo de iteracao (loop) 

  for niter=1:2000

 niter
  
  psi_old=psi;

  av=psi(nsel1,nsel)+psi(nsel2,nsel)+psi(nsel,nsel1)+psi(nsel,nsel2);
  F=1/r*(curl(nsel,nsel)*dx*dx - psi(nsel,nsel2).*(By(nsel,nsel))*dx ...
  +psi(nsel2,nsel).*(Bx(nsel,nsel))*dx );
  psi(nsel,nsel)=(av-F)./(4+(dx/r)*(By(nsel,nsel)-Bx(nsel,nsel) ) );
  crit= max(max(abs(psi-psi_old)))
  if crit <= tol,break,end

  end

% calcula vorticidade relativa

  zeta=(psi(nsel1,nsel)+psi(nsel2,nsel)+psi(nsel,nsel1)+psi(nsel,nsel2)- ...
  4*psi(nsel,nsel))./dx./dx;;


% plotando funcao de corrente

 orient portrait

	    p1=0.1*floor(min(min(psi))*10);
            p2=0.1*ceil(max(max(psi))*10);

 figure(3) 
 lp=p1:.1:p2; 
 cs=contourf(x,y,flipud(psi),lp,'w'); colorbar
% clabel(cs,lp)
 axis('square')
 xlabel('x [em 10^6 m]')
 ylabel('y [em 10^6 m]')
 title('A Solucao de Stommel')









