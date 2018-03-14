
% Editada por B.Sc. Danilo Augusto Silva
% Lab. Hidrodinâmica Costeira - IOUSP
% Fevereiro/2018
% rotina original de M.Sc. Carine de Godoi Rezende Costa
% atualizada por B.Sc. Rafaela Farias do Nascimento(agosto/2015)

%%%%%%%%%%%%%%%%%%%%%%%%% O QUE ESTE PROGRAMA FAZ %%%%%%%%%%%%%%%%%%%%%%%%%

% Importante: model_grid deve estar completo, com os pontos de terra todos
% 
% Em necessidade de baixar os dados, utlizar a rotina em python 'download_ncep.py'
% 
% 1) Lê os dados dos arquivos do diretorio selecionado
%     . como os dados são horarios, lê-se os dados a cada 6 horas
%     . interpola para a grade do modelo
%     . cria arquivo ASCII de entrada para o wind.f: wind_ncep
% OBS: a convenção das componentes do vento é vetorial, ou seja, o vetor
% apontar na direção que o vento vai

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variáveis
% load longitude and latitude from model_grid

clear all; close all; clc
addpath /home/tparente/danilo/matlab_programs/nctoolbox/
setup_nctoolbox

load /home/tparente/danilo/matlab_programs/linha_costaJR.mat

fid=fopen('/home/tparente/danilo/mestrado/github/artigoTG/grade/model_grid_comterra');

grade=textscan(fid,'%5d%5d%10.2f%10.2f%10.2f%10.2f%10.5f%10.5f%5.1f%5s',...
    'headerlines',20);
fclose(fid);

% clocar II e JJ da grade utilizada
II=94;
JJ=297;
i=II-2;
j=JJ-2;

lat=grade{7};
lon=grade{8};
 
lon=reshape(lon,i,j);
lonsecom=lon';
lat=reshape(lat,i,j);
latsecom=lat';

%% Cria arquivo de saida da interpolação
wf='wind_castelhanos';
fid=fopen(wf,'wt');

%% Leitura dos arquivos grib2 que serão utilizados na interpolação

arq = dir('/home/tparente/danilo/mestrado/artigo_tg/dados_vento/CSFR/castelhanos/wnd10m*');
% pr = dir('/home/tparente/danilo/mestrado/artigo_tg/dados_vento/CSFR/prmsl*');

% aqui posso colocar um código pra baixar direto do ftp: ftp://nomads.ncdc.noaa.gov/CFSR/HP_time_series/

%% loop de interpolação de cada arquivo listado

% vetor de tempo
tempo = [0:6:2184]; % contabilizando apenas 111 dias (reduzindo o aquecimento para 8 dias)
cont_tempo = 1;

disp('Inicio da leitura dos arquivos')

%% fazer para o primeiro arquivo 
k=1;

   % leitura do arquivo grib2
    disp(arq(k).name)

    nc = ncgeodataset(['/home/tparente/danilo/mestrado/artigo_tg/dados_vento/CSFR/castelhanos/' arq(k).name]);
    % extrair variáveis
    lat = nc.data('lat');
    lon = nc.data('lon');
    lon=lon-360;
    U = nc.variable('u-component_of_wind_height_above_ground');
    V = nc.variable('v-component_of_wind_height_above_ground');
    time = nc.data('time');
    pre = 1024;

    % criar loop aqui para ler os dados de 6 em 6 horas
    lat     =   double(lat);
    lon     =   double(lon);
    pre     =   double(pre);
    % lendo a partir do 10o dia do mês (2400 horas)
    for step=1:6:time(end)
        % ler o dado
        us = squeeze(U.data(step,1,1:numel(lat),1:numel(lon)));
        vs = squeeze(V.data(step,1,1:numel(lat),1:numel(lon)));
        
        us      =   double(us);
        vs      =   double(vs);

        % recortar 
        pla     =   find(lat>=min(latsecom(:))-2 & lat<=max(latsecom(:))+2);
        plo     =   find(lon>=min(lonsecom(:))-2 & lon<=max(lonsecom(:))+2);
        la      =   lat(pla);
        lo      =   lon(plo);
        u       =   us(pla,plo);
        v       =   vs(pla,plo);
        %p=pre(pla,plo);

        % vetoriza preservando a correspondência
        lonc    =   repmat(lo,numel(la),1);  % roda todas as lons primeiro
        latc    =   repmat(la,1,numel(lo));  % repete a mesma lat primeiro
        latc    =   latc';                   % transpõe para
        latc    =   latc(:);                 % vetorizar rodando todas as lons para uma mesma lat
        uc      =   u';
        uc      =   uc(:);
        vc      =   v';
        vc      =   vc(:);
        %pc=p';
        %pc=pc(:);

        % interpolar
        ui      =   griddata(lonc,latc,uc,lonsecom,latsecom,'linear');
        vi      =   griddata(lonc,latc,vc,lonsecom,latsecom,'linear');
        % multiplicando o vento por 2
        ui      = ui * 2;
        vi      = vi * 2;
        %pi=griddata(lonc,latc,pc,lonsecom,latsecom,'linear');

        % gravar 
        % criar as linhas e colunas de saída
        usecom  =   nan(JJ,II);
        vsecom  =   nan(JJ,II);
        psecom  =   nan(JJ,II);
        
        % copia os contornos adjacentes
        usecom([1 JJ],2:II-1)=ui([1 JJ-2],:);
        usecom(2:JJ-1,[1 II])=ui(:,[1 II-2]);
        % copia o cantinho com o vizinho
        usecom([1 JJ],1)=usecom([1 JJ],2);
        usecom([1 JJ],II)=usecom([1 JJ],II-1);
        % preenche o interior
        usecom(2:JJ-1,2:II-1)=ui;
        
        vsecom([1 JJ],2:II-1)=vi([1 JJ-2],:);
        vsecom(2:JJ-1,[1 II])=vi(:,[1 II-2]);
        vsecom([1 JJ],1)=vsecom([1 JJ],2);
        vsecom([1 JJ],II)=vsecom([1 JJ],II-1);
        vsecom(2:JJ-1,2:II-1)=vi;
        
        psecom([1 JJ],2:II-1)=    pre              %pi([1 JJ-2],:);
        psecom(2:JJ-1,[1 II])=    pre              %pi(:,[1 II-2]);
        psecom([1 JJ],1)=     pre              %psecom([1 JJ],2);
        psecom([1 JJ],II)=    pre              %psecom([1 JJ],II-1);
        psecom(2:JJ-1,2:II-1)=    pre              %pi;
        
        % cria arquivos de entrada do modelo 

        % time=arq(k).name(18:19);
        time=tempo(cont_tempo);
        cont_tempo = cont_tempo + 1;
        fprintf(fid,'%10.5f\n',time);
        for j=1:size(usecom,1)
            for i=1:size(usecom,2)
                fprintf(fid,'%5.0f%5.0f%10.3f%10.3f%10.3f\n',...
                    [i j usecom(j,i) vsecom(j,i) psecom(j,i)]);
            end 
        end

    end

%% fazer para Janeiro, Fevereiro e Março

for k=2:numel(arq)
    % leitura do arquivo grib2
    disp(arq(k).name)

    nc = ncgeodataset(['/home/tparente/danilo/mestrado/artigo_tg/dados_vento/CSFR/castelhanos/' arq(k).name]);
    % extrair variáveis
    lat = nc.data('lat');
    lon = nc.data('lon');
    lon=lon-360;
    U = nc.variable('u-component_of_wind_height_above_ground');
    V = nc.variable('v-component_of_wind_height_above_ground');
    time = nc.data('time');
    pre = 1024;

    % criar loop aqui para ler os dados de 6 em 6 horas
    lat     =   double(lat);
    lon     =   double(lon);
    pre     =   double(pre);

    for step=1:6:time(end)
        % ler o dado
        us = squeeze(U.data(step,1,1:numel(lat),1:numel(lon)));
        vs = squeeze(V.data(step,1,1:numel(lat),1:numel(lon)));
        
        us      =   double(us);
        vs      =   double(vs);

        % recortar 
        pla     =   find(lat>=min(latsecom(:))-2 & lat<=max(latsecom(:))+2);
        plo     =   find(lon>=min(lonsecom(:))-2 & lon<=max(lonsecom(:))+2);
        la      =   lat(pla);
        lo      =   lon(plo);
        u       =   us(pla,plo);
        v       =   vs(pla,plo);
        %p=pre(pla,plo);

        % vetoriza preservando a correspondência
        lonc    =   repmat(lo,numel(la),1);  % roda todas as lons primeiro
        latc    =   repmat(la,1,numel(lo));  % repete a mesma lat primeiro
        latc    =   latc';                   % transpõe para
        latc    =   latc(:);                 % vetorizar rodando todas as lons para uma mesma lat
        uc      =   u';
        uc      =   uc(:);
        vc      =   v';
        vc      =   vc(:);
        %pc=p';
        %pc=pc(:);

        % interpolar
        ui      =   griddata(lonc,latc,uc,lonsecom,latsecom,'linear');
        vi      =   griddata(lonc,latc,vc,lonsecom,latsecom,'linear');
        % multiplicando o vento por 2
        ui      = ui * 2;
        vi      = vi * 2;
        %pi=griddata(lonc,latc,pc,lonsecom,latsecom,'linear');

        % gravar 
        % criar as linhas e colunas de saída
        usecom  =   nan(JJ,II);
        vsecom  =   nan(JJ,II);
        psecom  =   nan(JJ,II);
        
        % copia os contornos adjacentes
        usecom([1 JJ],2:II-1)=ui([1 JJ-2],:);
        usecom(2:JJ-1,[1 II])=ui(:,[1 II-2]);
        % copia o cantinho com o vizinho
        usecom([1 JJ],1)=usecom([1 JJ],2);
        usecom([1 JJ],II)=usecom([1 JJ],II-1);
        % preenche o interior
        usecom(2:JJ-1,2:II-1)=ui;
        
        vsecom([1 JJ],2:II-1)=vi([1 JJ-2],:);
        vsecom(2:JJ-1,[1 II])=vi(:,[1 II-2]);
        vsecom([1 JJ],1)=vsecom([1 JJ],2);
        vsecom([1 JJ],II)=vsecom([1 JJ],II-1);
        vsecom(2:JJ-1,2:II-1)=vi;
        
        psecom([1 JJ],2:II-1)=    pre              %pi([1 JJ-2],:);
        psecom(2:JJ-1,[1 II])=    pre              %pi(:,[1 II-2]);
        psecom([1 JJ],1)=     pre              %psecom([1 JJ],2);
        psecom([1 JJ],II)=    pre              %psecom([1 JJ],II-1);
        psecom(2:JJ-1,2:II-1)=    pre              %pi;
        
        % cria arquivos de entrada do modelo 

        % time=arq(k).name(18:19);
        time=tempo(cont_tempo);
        cont_tempo = cont_tempo + 1;
        fprintf(fid,'%10.5f\n',time);
        for j=1:size(usecom,1)
            for i=1:size(usecom,2)
                fprintf(fid,'%5.0f%5.0f%10.3f%10.3f%10.3f\n',...
                    [i j usecom(j,i) vsecom(j,i) psecom(j,i)]);
            end 
        end

    end

end

%% complementar o arquivo repetindo o ultimo instante de tempo   

% vetor de tempo
tempo = [4368:6:4380];
%tempo(182) = tempo(182) + 1;
cont_tempo = 1;

disp('Inicio da leitura dos arquivos')
k = numel(arq);

% repetir o ultimo instante de tempo 2x

for repeat=1:1:2
    % leitura do arquivo grib2
    disp(arq(k).name)

    nc = ncgeodataset(['/home/tparente/danilo/mestrado/artigo_tg/dados_vento/CSFR/castelhanos/' arq(k).name]);
    
    % extrair variáveis
    lat = nc.data('lat');
    lon = nc.data('lon');
    lon=lon-360;
    U = nc.variable('u-component_of_wind_height_above_ground');
    V = nc.variable('v-component_of_wind_height_above_ground');
    time = nc.data('time');
    pre = 1024;

    % criar loop aqui para ler os dados de 6 em 6 horas
    lat     =   double(lat);
    lon     =   double(lon);
    pre     =   double(pre);


	step = time(end);

        % ler o dado
        us = squeeze(U.data(step,1,1:numel(lat),1:numel(lon)));
        vs = squeeze(V.data(step,1,1:numel(lat),1:numel(lon)));
        
        us      =   double(us);
        vs      =   double(vs);

        % recortar 
        pla     =   find(lat>=min(latsecom(:))-2 & lat<=max(latsecom(:))+2);
        plo     =   find(lon>=min(lonsecom(:))-2 & lon<=max(lonsecom(:))+2);
        la      =   lat(pla);
        lo      =   lon(plo);
        u       =   us(pla,plo);
        v       =   vs(pla,plo);
        %p=pre(pla,plo);

        % vetoriza preservando a correspondência
        lonc    =   repmat(lo,numel(la),1);  % roda todas as lons primeiro
        latc    =   repmat(la,1,numel(lo));  % repete a mesma lat primeiro
        latc    =   latc';                   % transpõe para
        latc    =   latc(:);                 % vetorizar rodando todas as lons para uma mesma lat
        uc      =   u';
        uc      =   uc(:);
        vc      =   v';
        vc      =   vc(:);
        %pc=p';
        %pc=pc(:);

        % interpolar
        ui      =   griddata(lonc,latc,uc,lonsecom,latsecom,'linear');
        vi      =   griddata(lonc,latc,vc,lonsecom,latsecom,'linear');
        %pi=griddata(lonc,latc,pc,lonsecom,latsecom,'linear');

        % gravar 
        % criar as linhas e colunas de saída
        usecom  =   nan(JJ,II);
        vsecom  =   nan(JJ,II);
        psecom  =   nan(JJ,II);
        
        % copia os contornos adjacentes
        usecom([1 JJ],2:II-1)=ui([1 JJ-2],:);
        usecom(2:JJ-1,[1 II])=ui(:,[1 II-2]);
        % copia o cantinho com o vizinho
        usecom([1 JJ],1)=usecom([1 JJ],2);
        usecom([1 JJ],II)=usecom([1 JJ],II-1);
        % preenche o interior
        usecom(2:JJ-1,2:II-1)=ui;
        
        vsecom([1 JJ],2:II-1)=vi([1 JJ-2],:);
        vsecom(2:JJ-1,[1 II])=vi(:,[1 II-2]);
        vsecom([1 JJ],1)=vsecom([1 JJ],2);
        vsecom([1 JJ],II)=vsecom([1 JJ],II-1);
        vsecom(2:JJ-1,2:II-1)=vi;
        
        psecom([1 JJ],2:II-1)=    pre              %pi([1 JJ-2],:);
        psecom(2:JJ-1,[1 II])=    pre              %pi(:,[1 II-2]);
        psecom([1 JJ],1)=     pre              %psecom([1 JJ],2);
        psecom([1 JJ],II)=    pre              %psecom([1 JJ],II-1);
        psecom(2:JJ-1,2:II-1)=    pre              %pi;
        
        % cria arquivos de entrada do modelo 

        % time=arq(k).name(18:19);
        time=tempo(cont_tempo);
        cont_tempo = cont_tempo + 1;
        fprintf(fid,'%10.5f\n',time);
        for j=1:size(usecom,1)
            for i=1:size(usecom,2)
                fprintf(fid,'%5.0f%5.0f%10.3f%10.3f%10.3f\n',...
                    [i j usecom(j,i) vsecom(j,i) psecom(j,i)]);
            end 
        end

end

%% parametros numericos


clc
disp('Parametros numericos e de output para run_data');

% imprimir quantos dias e os parametros necessarios no run_data
dias = 0;

for k=1:numel(arq)
    nc = ncgeodataset(['/home/tparente/danilo/mestrado/artigo_tg/dados_vento/CSFR/castelhanos/' arq(k).name]);
    
    time = nc.data('time');
    s = size(time);
    s = s - 1;
    dias = dias + s(1)/24;
    
end

% parametros numéricos
dti         = 15; %seg
nstep       = (dias/dti)*3600*24;

% parametros de output/netcdf
saving_time = 6*60; % minutos
avge        = (saving_time*60)/dti;
jhm         = nstep/avge;

avge_print = sprintf('#       AVGE: %i', avge);
nstep_print = sprintf('#      NSTEP: %i', nstep);
jhm_print = sprintf('#      JHM: %i', jhm);

disp('#-----------------------#');
disp(nstep_print);
disp(avge_print);
disp(jhm_print);
disp('#        DTI: 15s');
disp('#-----------------------#');
