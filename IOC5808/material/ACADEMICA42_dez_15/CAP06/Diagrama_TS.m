%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA Diagrama_TS.m                           %
% PLOTA DADOS DE TEMPERATURA E DE SALINIDADE       % 
% NO CLÁSSICO DIAGRAMA TS DE HELLAND HANSEN        %      
% GERANDO O DIAGRAMA DE ESTADO NA FORMA DE IMAGENS % 
% DE PARES DE PONTOS (T,S) ESPALHADOS              %
% PREPARADO POR FERNANDO PINHEIRO ANDUTTA EM 2009  %
% MAIS DETALHES SOBRE APLICAÇÃO PARA MASSAS        % 
% DE ÁGUA ESTUARINAS EM Miranda et al. (2002, 2012)% 
% p. 233 e p. 237-238                              %
% ATENÇÃO: O PROGRAMA USA A SUBROTINA "SW_DENS0"   %
% DO MORGAN (1994)- CSIRO MARINE LABORATORIES,     %
% AUSTRÁLIA, 222, 28 P. PARA CALCULAR A DENSIDADE  %
% DA ÁGUA ESTUARINA À PRESSÃO HIDROSTÁTICA ZERO    % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  close all
  clear all


  esp_title = 13;
  esp_xlabel = 13;
  esp_ylabel = 13;
  line = 2;

  % ESCOLHA DO IDIOMA: Português =1; Inglês=0
  % PARA O PRODUTO FINAL: FIGURA DO DIAGRAMA T-S
  % A BARRA DE CORES DO DIAGRAMA INDICA O PARÂMETRO SIGMA-T
  
  language = 1;
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % DEVIDO A RÁPIDA VARIAÇÃO NO TEMPO DAS MASSAS DE ÁGUA ESTUARINAS  %
 % O PROGRAMA FOI PREPARADO PARA PLOTAR O DIAGRAMA COM              %
 % DADOS AMOSTRADOS DURANTE MARÉS DE QUADRATURA E DE SIZÍGIA        %
 % OS ARQUIVOS DESSES DADOS SÃO OS SEGUINTES:                       %
 % "quadratura.dat" e "sizígia.dat"                                 % 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  load quadratura.dat
  
  data1 = quadratura;        
  prof1 = data1(:,2);
  temp1 = data1(:,5);
  sal1 = data1(:,4);
  dens1 = data1(:,6);

  load sizigia.dat
  
  data2=sizigia;        
  prof2 = data2(:,2);
  temp2 = data2(:,5);
  sal2 = data2(:,4);
  dens2 = data2(:,6);
  
  dens2=dens2+1.e+003
  
  min_sal1 = min(sal1);
  min_temp1 = min(temp1);
  min_sal2 = min(sal2);
  min_temp2 = min(temp2);

  max_sal1 = max(sal1);
  max_temp1 = max(temp1);
  max_sal2 = max(sal2);
  max_temp2 = max(temp2);

%%%%Criando os vetores que formarão a grade para o cálculo das isopcnais

if min_sal1<min_sal2
    minS = min_sal1;
else
    minS = min_sal2;
end
if max_sal1<max_sal2
    maxS = max_sal2;
else
    maxS = max_sal1;
end

if min_temp1<min_temp2
    minT = min_temp1;
else
    minT = min_temp2;
end
if max_sal1<max_temp2
    maxT = max_temp2;
else
    maxT = max_temp1;
end

minT = floor(minT);
maxT = ceil(maxT);
minS = floor(minS);
maxS = ceil(maxS);

dT = (maxT-minT)/4;
dS = (maxS-minS)/4;
T = (minT:dT:maxT);  %% atention delta T must fill well from minT to maxT
S = (minS:dS:maxS);
[GrS,GrT] = meshgrid(S,T); %% atention delta S must fill well from minS to maxS

%%%%Builting isopcinais lines
%%%%Criando o grafico das isopcnais
dens = (sw_dens0(GrS,GrT)) - 1000;

figure;

% A BARRA DE CORES DO DIAGRAMA INDICA O PARÂMETRO SIGMA-T

hold on;
min_dens = floor(min(min(dens)));
max_dens = ceil(max(max(dens)));

dD = (max_dens-min_dens)/5;

c = contourf(GrS,GrT,dens,min_dens-3:dD:max_dens+3);
colormap('winter');
shading flat
axis([minS maxS minT maxT])
%clabel(c,min_dens:1:max_dens);
%text(36.8,20,'AT'); 

if language==1
xlabel('Salinidade','fontsize',esp_xlabel); 
ylabel('Temperatura (ºC)','fontsize',esp_ylabel);
colorbar
%title('Diagrama TS espalhado','fontsize',esp_title)
else 
xlabel('Salinity','fontsize',esp_xlabel); 
ylabel('Temperature (ºC)','fontsize',esp_ylabel);
colorbar
%title(' TS Diagram scattered','fontsize',esp_title)    
end

   hold on
  plot(sal1,temp1,'W.');
  %p1 = plot(sal1,temp1,'W.');
  %set(p1,'linewidth',line);
  hold on
  plot(sal2,temp2,'k.');
  %p2 = plot(sal2,temp2,'k.');
  %set(p1,'linewidth',line);
  
 print -dbitmap f:\academica42\cap06\Diagrama_color

 clabel(c,min_dens:0.5:max_dens);
 grid



