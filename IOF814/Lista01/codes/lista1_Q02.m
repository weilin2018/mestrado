%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programa desenvolvivo para a Questão 2 da Lista 1 de exercícios da disciplina
% IOF814 - Modelos Numéricos Aplicados a Processos Costeiros e Estuarinos
% ministrada pelo Prof Joseph Harari, do Instituto Oceanográfico da USP.
%
% Para mais detalhes do desenvolvimento da discretização, pode ser conferido
% na solução da Lista, em IOF814/Lista1/outputs/Lista1_IOF814.pdf
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
fant=fatu;
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

sk=zeros(kmax,jmax);
pk=zeros(kmax,jmax);
dk=zeros(kmax,jmax);

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
   % lembrando que: dj usa 2 niveis de tempo (fant e fatu) e dk usa somente 1 (fatu)

   for j=2:1:jmax-1
        for k=2:1:jmax-1
            dj(k,j)=fant(k,j) - qu*(fatu(k,j+1)-fatu(k,j-1));
            dk(k,j)=-qv*(fatu(k+1,j)-fatu(k-1,j));
        end
   end

   %varredura ascendente
   for j=2:jmax-1
     sj(j)=-cj/(bj+aj*sj(j-1));
     pj(j)=(dj(j)-aj*pj(j-1))/(bj+aj*sj(j-1));
   end
   for k=2:kmax-1
     sk(k)=-ck/(bk+ak*sk(k-1));
     pk(k)=(dk(k)-ak*pk(k-1))/(bk+ak*sk(k-1));
   end

  % varredura descendente
  % considerando que jmax=kmax, entao fazemos somente um loop para calcular fren
  for j=jmax-1:-1:2
    k=j; % facilitar a  leitura da formula seguinte:
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
    pause(0.01)
    hold off
  end

  fant=fatu;
  fatu=fren;

end
