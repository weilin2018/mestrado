% INICIALMENTE SOMENTE LEITURA DOS DADOS OBSERVADOS E TRATAMENTO ESTATÍSTICO DO MESMO

close all;
clear all;
clc;

addpath(genpath('/media/danilo/Danilo/mestrado/ventopcse/data/ECOSAN/morais2016/exp_sECOM/ocean_toolbox')); % adicionando pastas e subpastas ao toolbox

fontes;


cd /media/danilo/Danilo/mestrado/ventopcse/data/ECOSAN/morais2016/exp_sECOM;
load /media/danilo/Danilo/mestrado/ventopcse/data/ECOSAN/Dados_mexidos/ES0802.Dat;

%% Estabelecendo vari�veis

data=ES0802;
v=data(:,1);%velocidade meridional
u=data(:,2);%velocidade zonal
velcor=data(:,3);%resultante da velocidades
dircor=data(:,4);%direcao, em graus, das correntes
hora=data(:,6); %horas em que foram feitas as medi��es
minu=data(:,7); %min em que foram feitas as medi��es
seg=data(:,8); %segundos em que foram feitas as medi��es
mes=data(:,9); %mes em que foram feitas as medi��es
dia=data(:,10); %dias em que foram feitas as medi��es
ano=data(:,11); %anos em que foram feitas as medi��es
ang=55; % angulo de rota��o do vetor
dec=24; % correcao da declinacao magnetica
[ur, vr]=rot(u, v, ang+dec); % rotacao dos vetores para as componentes paralelas e perpendiculares (u, v, angulo de rota��o)
tempo=julian([ano mes dia hora minu seg]); %cria dia juliano

%% Calculando media de 3 minutos, a cada 30 minutos
% pois o tempo estava espa�ado em 3 medi��es, a cada 30 minutos

% ur
cont=0;
for j=1:length(tempo);
 if data(:,7)==2;
  if ((tempo(j-1)-tempo(j))>-7e-4)&((tempo(j+1)-tempo(j))<7e-4);
   cont=cont+1;
   ur_m(cont)=mean(ur(j-1:j+1));
   tempo_m(cont)=tempo(j);
  end;
  end;
 if data(j,7)==32;
  if ((tempo(j-1)-tempo(j))>-7e-4)&((tempo(j+1)-tempo(j))<7e-4);
   cont=cont+1;
   ur_m(cont)=mean(ur(j-1:j+1));
   tempo_m(cont)=tempo(j);
  end;
  end;
  end; clear j;

  % vr
cont=0;
for j=1:length(tempo);
 if data(:,7)==2;
  if ((tempo(j-1)-tempo(j))>-7e-4)&((tempo(j+1)-tempo(j))<7e-4);
   cont=cont+1;
   vr_m(cont)=mean(vr(j-1:j+1));
   tempo_m(cont)=tempo(j);
  end;
  end;
 if data(j,7)==32;
  if ((tempo(j-1)-tempo(j))>-7e-4)&((tempo(j+1)-tempo(j))<7e-4);
   cont=cont+1;
   vr_m(cont)=mean(vr(j-1:j+1));
   tempo_m(cont)=tempo(j);
  end;
  end;
  end; clear j;

%% Interpolacao das medias para cada hora

tempo_m=tempo_m-min(tempo_m);% para comecar do tempo 'zero'
tempo_m=tempo_m*24;% transformando de dias para horas
tempo_i=0:6:max(tempo_m);% criando novo vetor de tempo em horas

ind=find(isnan(ur_m)==0);
for j=1:length(tempo_i);
    ur_m_6h(j)=interp1(tempo_m(ind),ur_m(ind),tempo_i(j),'spline');
    vr_m_6h(j)=interp1(tempo_m(ind),vr_m(ind),tempo_i(j),'spline');
end;

tempo=julian(ano(1),mes(1),dia(1),hora(1))+tempo_i./24;%transformando para julian (data inicial + horas./24)_________________MUDAR DATA PARA CADA UM!
gregorian(tempo(1))
gregorian(tempo(end))

% close all
plot(tempo,ur_m_6h)
hold on
gregaxd(tempo,7)

ur_m_6h=ur_m_6h(251:370);
vr_m_6h=vr_m_6h(251:370);
tempo=tempo(251:370);
gregorian(tempo(1))
gregorian(tempo(2))
gregorian(tempo(end-1))
gregorian(tempo(end))

%% Controle de qualidade

u1=ur_m_6h;
v1=vr_m_6h;
limu=mean(u1)+3*std(u1);%limite superior da media + 3 desvios padrao para u
limv=mean(v1)+3*std(v1);%limite superior da media + 3 desvios padrao para v
i=1;
limu_=mean(u1)-3*std(u1);%limite inferior da media - 3 desvios padrao para u
limv_=mean(v1)-3*std(v1);%limite inferior da media - 3 desvios padrao para v
j=1;

ind=find(isnan(u1)==0);

    while i<=length(u1);
    if u1(i)>limu;
        u1(i)=nan;
    end
    if v1(i)>limv;
        v1(i)=nan;
    end
    i=i+1;
    end;

    while j<=length(v1);
    if u1(j)<limu_;
        u1(j)=nan;
    end

    if v1(j)<limv_;
        v1(j)=nan;
    end
    j=j+1;
    end;

%%  interpolação dos NaN

u1=interp1(tempo,u1,tempo,'pchip');
v1=interp1(tempo,v1,tempo,'pchip');

%% VELOCIDADES DAS COMPONENTES NO TEMPO

fig1=figure;
make_nice(fig1);
maximize(fig1);

subplot(3,1,1);
h=plot(tempo,u1,'k*-');% dados sem filtro hanning no tempo
set(h,'linewidth',2,'color',[.7 .7 .7]);% plotar com maior espessura e em escala de cinza
hold on;
ylabel('Velocidades cm.s^{-1}'); %colocando um rotulo no eixo y
gregaxd(tempo,7);%transformando para calendário gregoriano, com intervalos de 10 dias
xlim([min(tempo)-1 max(tempo)+1]);% definindo limites do eixo 'x'          PARA f_40304 O LIMITE É % 2.447560958333334e+06 e para o restante é max(tempo)+2
grid minor;
legend('Componente Perpendicular');%,'Location','SouthEast');
title(['ES0802 pp - Perpendicular (u)','       mean(u)=',num2str(nanmean(u1)),'       std(u)=',num2str(nanstd(u1))]);
set(gca,'FontSize',9);

subplot(3,1,2);
hh=plot(tempo,v1,'b*-');% dados sem filtro hanning no tempo
set(hh,'linewidth',2,'color',[.7 .7 .7]);% plotar com maior espessura e em escala de cinza
hold on;
ylabel('Velocidades cm.s^{-1}'); %colocando um rotulo no eixo y
gregaxd(tempo,7);%transformando para calendário gregoriano, com intervalos de 10 dias
xlim([min(tempo)-1 max(tempo)+1]);% definindo limites do eixo 'x'          PARA f_40304 O LIMITE É % 2.447560958333334e+06 e para o restante é max(tempo)+2
grid minor;
title(['ES0802 pp - Paralela (v)','       mean(v)=',num2str(nanmean(v1)),'       std(v)=',num2str(nanstd(v1))]);
legend('Componente Paralela');%,'Location','SouthEast');
set(gca,'FontSize',9);

subplot(3,1,3);
hh=plot(tempo,u1,'k',tempo,v1,'r');
set(hh,'linewidth',2.5);
gregaxd(tempo,7);%transformando para calendário gregoriano, com intervalos de 10 dias
set(gca,'fontsize',9);
ylabel('Velocidades cm.s^{-1}'); %colocando um rotulo no eixo y
xlim([min(tempo)-1 max(tempo)+1]);% definindo limites do eixo 'x'          PARA f_40304 O LIMITE É % 2.447560958333334e+06 e para o restante é max(tempo)+2
grid minor;
title(['ES0802 pp','       min(u)=',num2str(min(u1)),'       max(u)=',num2str(max(u1)),'       min(v)=',num2str(min(v1)),'       max(v)=',num2str(max(v1))]);
legend('Perpendicular','Paralela');%,'Location','SouthEast')
set(gca,'FontSize',9);

eval(['print -dpng -r300 /home/ocpm/Dropbox/MESTRADO_USP/DADOS_PCI/ecosan_pci/figures_ecosan/ES0802_pp_fil_subplot_serietemp_u_v_validacao.png']);

%% SALVANDO VARIAVEIS PARA CORRELACIONAR

u11=u1;
v11=v1;
save('/home/ocpm/Dropbox/MESTRADO_USP/exp_sECOM/u_v_ecosan_100.mat','u11','v11');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load u_v_ecosan_100.mat;
load valida_ecosan_100.mat;

%% AJUSTANDO FASE DAS VARIAVEIS

% urh_ecosan_ES0802=urh_ecosan_ES0802(5:end);
% vrh_ecosan_ES0802=vrh_ecosan_ES0802(5:end);
% wv_ecosan_ES0802=wv_ecosan_ES0802(5:end);
% wu_ecosan_ES0802=wu_ecosan_ES0802(5:end);
% u11=u11(1:end-4);
% v11=v11(1:end-4);
% tempo0=tempo0(1:length(u11));
zero=zeros(size(tempo0));

%% SKILL

willmott_perp = skill_willmott(u11,urh_ecosan_ES0802);
willmott_para = skill_willmott(v11,vrh_ecosan_ES0802);

gregorian(tempo0(1))
gregorian(tempo0(end))
%%

fig5=figure;
make_nice(fig5);maximize(fig5);

subplot(3,1,1);
h=plot(tempo0,wv_ecosan_ES0802,'r')
set(h,'linewidth',2);
hold on;
plot(tempo0,zero,'k');
ylabel('m.s^{-1}'); %colocando um rotulo no eixo y
xlim([min(tempo0)-1 max(tempo0)+1]);% definindo limites do eixo 'x'          PARA f_40304 O LIMITE É % 2.447560958333334e+06 e para o restante é max(tempo)+2
t=title(['V3 - Vento - Componente Paralela - 10 m de altitude']);
set(t,'fontsize',14);
set(gca, 'xtick', []);
l=legend('Paralela','location','southeast');
set(l,'fontsize',10);
grid minor;
box off;

subplot(3,1,2);
h=plot(tempo0,urh_ecosan_ES0802,'k',tempo0,u11,'k:')
set(h,'linewidth',2);
hold on;
plot(tempo0,zero,'k');
ylabel('cm.s^{-1}'); %colocando um rotulo no eixo y
xlim([min(tempo0)-1 max(tempo0)+1]);% definindo limites do eixo 'x'          PARA f_40304 O LIMITE É % 2.447560958333334e+06 e para o restante é max(tempo)+2
t=title(['V3 - Corrente - Componente Perpendicular - Isóbata de 100 m - Prof. real de 56 m','       skill=',num2str(willmott_perp,'%.2f')]);
set(t,'fontsize',14);
set(gca, 'xtick', []);
l=legend('Modelo','Real','location','southeast');
set(l,'fontsize',10);
grid minor;
box off;

subplot(3,1,3);
h=plot(tempo0,vrh_ecosan_ES0802,'r',tempo0,v11,'r:')
set(h,'linewidth',2);
hold on;
plot(tempo0,zero,'k');
ylabel('cm.s^{-1}'); %colocando um rotulo no eixo y
xlim([min(tempo0)-1 max(tempo0)+1]);% definindo limites do eixo 'x'          PARA f_40304 O LIMITE É % 2.447560958333334e+06 e para o restante é max(tempo)+2
set(gca, 'xtick', []);
l=legend('Modelo','Real','location','southeast');
set(l,'fontsize',10);
t=title(['V3 - Corrente - Componente Paralela - Isóbata de 100 m - Prof. real de 56 m','       skill=',num2str(willmott_para,'%.2f')]);
set(t,'fontsize',14);
grid minor;
box off;
gregaxd(tempo0,3);
xlim([min(tempo0)-1 max(tempo0)+1]);% definindo limites do eixo 'x'          PARA f_40304 O LIMITE É % 2.447560958333334e+06 e para o restante é max(tempo)+2
