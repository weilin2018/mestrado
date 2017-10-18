clear all; close all;

jmax=200;
kmax=200;
nmax=90;
u=0.1;           % componente u da velocidade
v=0.1;           % componente v da velocidade

dx=10;
dy=10;
dt=4;
pol=100;
posini=95;
posfim=105;
freqplo=5;
xgrid=((1:jmax)-1)/dx;

% CALCULOS INICIAIS: calcular q!!!
qu = (dt*u)/(2*dx);
qv = (dt*v)/(2*dy);

% CONDICOES INICIAIS
fatu=zeros(jmax,kmax);
fren=zeros(jmax,kmax);
fant=zeros(jmax,kmax);

fatu(posini:posfim)=pol;
fcin=fatu;

contplo=2;
pol050=0.2*pol;
pol150=1.2*pol;

% criando os indices
sj=zeros(kmax,jmax);
pj=zeros(kmax,jmax);
dj=zeros(kmax,jmax);

aj=(dt*u)/(2*dx);
bj=1;
cj=(-dt*u)/(2*dx);

sk=zeros(kmax,kmax);
pk=zeros(kmax,kmax);
dk=zeros(kmax,kmax);

ak=-(2*v)/(2*dy);
bk=0;
ck=(2*v)/(2*dy);


% LOOP NO TEMPO
% CONDICOES DE CONTORNO
% FORMULA DE RECORRENCIA
% PLOTAGEM (PRESSIONE ENTER PARA EVOLUIR NO TEMPO)
% EVOLUCAO NO TEMPO DAS VARIAVEIS
for n=3:nmax
   tempo=n*dt;
   % calcular dj e dk
   dj(2:jmax-1,2:kmax-1)=fant(2:jmax-1,2:kmax-1) - qu*(fatu(3:jmax,2:kmax-1) - fatu(1:jmax-2,2:kmax-1));
   dk(2:jmax-1,2:kmax-1)=-qv*(fatu(2:jmax-1,3:jmax) - fatu(2:jmax-1, 1:kmax-2));

   %varredura ascendente
   for j=2:jmax-1
     sj(j)=-cj/(bj+aj*sj(j-1));
     pj(j)=(dj(j)-aj*pj(j-1))/(bj+aj*sj(j-1));
   end
   for k=2:kmax-1
     sj(k)=-ck/(bk+ak*sk(k-1));
     pj(k)=(dk(k)-ak*pk(k-1))/(bk+ak*sk(k-1));
   end

  % varredura descendente
  % considerando que jmax=kmax, entao fazemos somente um loop para calcular fren
  for j=jmax-1:-1:2
    k=j; % facilitar a leitura da formula seguinte:
    fren(j,k) = sj(j)*fren(j+1,k) + pj(j) + sk(k)*fren(j,k+1) + pk(k);
  end

  contplo=contplo+1;
  if(contplo==freqplo)
    contplo=0;
    figure (1)
    plot(xgrid,fcin,'r','LineWidth',2)
    hold
    plot(xgrid,fren,'LineWidth',2)
    axis([xgrid(1) xgrid(jmax) -pol050 pol150]);
    title(['Adveccao de sinal retangular (semi implic, 2a ordem) - tempo ',...
        num2str(tempo),' segundos'],'fontsize',12)
    xlabel('DISTANCIA NA GRADE(m)','fontsize',12)
    ylabel('conc','fontsize',12)
    grid on
    pause
    hold off
  end

  fant=fatu;
  fatu=fren;

end
