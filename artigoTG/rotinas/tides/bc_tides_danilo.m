addpath('/home/tparente/Dropbox/TOOLBOX_MATLAB/new_codes/tmd_toolbox/');

tpxo = '/home/tparente/Dropbox/TOOLBOX_MATLAB/new_codes/tmd_toolbox/DATA/Model_tpxo7.2';

grd  = '/home/tparente/Downloads/shelfwaves/model_grid_comterra';

II=94;
JJ=297;
IKB=15;

fi=fopen(grd);
grade=textscan(fi,'%5d%5d%10.2f%10.2f%10.2f%10.2f%10.5f%10.5f%5.1f%5s',...
    'headerlines',IKB+5);
fclose(fi);

I=grade{1}; % pontos de células na perpendicular
J=grade{2}; % pontos de células longitudinais
H=grade{5}; % Profundidade
lat=grade{7}; % Latitude 
lon=grade{8}; % Longitude 
i=II-2; % retirando as últimas células
j=JJ-2; % retirando as últimas células
grdx=reshape(lon,i,j);
grdy=reshape(lat,i,j);
depth=reshape(H,i,j);
Ii=reshape(I,i,j);
Jj=reshape(J,i,j);

%% ===== ACHANDO OS CONTORNOS ============
 
id=find(H>0);
lat_water=lat(id);lon_water=lon(id);
I_water=Ii(id);J_water=Jj(id);
% ----- LON / LAT / DO CONTORNO LESTE ------
id=find(I_water==I_water(end)); % leste
lon_w1=lon_water(id);lat_w1=lat_water(id);
 
% ----- LON / LAT / DO CONTORNO NORTE
id=find(J_water==J_water(end)); % contorno norte
lat_w2=lat_water(id);lon_w2=lon_water(id);
id=find(lon==lon_w2(1));
I_norte=I(id);
 
% ----- LON / LAT / DO CONTORNO SUL
id=find(J_water==J_water(1)); % contorno sul
lat_w3=lat_water(id);lon_w3=lon_water(id);
id=find(lon==lon_w3(1));
I_sul=I(id);

%% ===== ACHANDO AMPLITUDE E FASE PARA PONTOS DO CONTORNO ===
 
% ------ CONTORNO LESTE -------
 
LESTEamp=[];LESTEfase=[];
for kk=1:length(lat_w1)
    [a,f,p,c]=tmd_extract_HC(tpxo,lat_w1(kk),lon_w1(kk),'z');
    LESTEamp(:,kk)=a; LESTEfase(:,kk)=f;
end
 
% ------ CONTORNO NORTE -------
NORTEamp=[];NORTEfase=[];
for kk=1:length(lat_w2)
    [a,f,p,c]=tmd_extract_HC(tpxo,lat_w2(kk),lon_w2(kk),'z');
    NORTEamp(:,kk)=a; NORTEfase(:,kk)=f;
end
 
% ------ CONTORNO SUL -------
SULamp=[];SULfase=[];
for kk=1:length(lat_w3)
    [a,f,p,c]=tmd_extract_HC(tpxo,lat_w3(kk),lon_w3(kk),'z');
    SULamp(:,kk)=a; SULfase(:,kk)=f;
end
 
%% ====== INDICES DOS HARMONICOS ===============
% M2 = 1 / S2=2 / N2=3 / K1=5 / O1=6 / P1=7
 
%% ====== COLOCANDO NO FORMATO DO RUN_DATA ==================
 
ngrid1=0; % CONTADOR PARA O LOOP
 
% Offshore boundary
veti=[II-1 II-2];
vetj=3:JJ-2;
for k=1:length(vetj)
    east(k,:)=[veti(1) vetj(k) veti(2) vetj(k)];
    ngrid1=ngrid1+1;
end
emeaneast(1:length(vetj))=0;
 
%% North boundary
veti=I_norte:II-1;
vetj=[JJ-1 JJ-2];
for k=1:length(veti)
    north(k,:)=[veti(k) vetj(1) veti(k) vetj(2)];
    ngrid1=ngrid1+1;
end
emeannorth(1:length(veti))=0;
 
%% South boundary
veti=I_sul:II-1;
vetj=[2 3];
for k=1:length(veti)
    south(k,:)=[veti(k) vetj(1) veti(k) vetj(2)];
    ngrid1=ngrid1+1;
end
emeansouth(1:length(veti))=0;
 
fprintf(1,'ngrid1 = %10.0f',ngrid1);
 
%% ===== ESCREVENDO =====
ngrid2=0
out='/home/tparente/danilo/mestrado/github/artigoTG/rotinas/tides/bc_tide'
fid=fopen(out,'w')
 
% --- ORDEM DO RUN_DATA -----
% S2 // M2 // N2 // K1 // P1 // O1
harm=[2,1,3,5,7,6];
 
for ik=1:length(south)
    fprintf(fid,'%5.0f%5.0f%5.0f%5.0f%10.5f\n',south(ik,:),emeansouth(ik))
    fprintf(fid,'%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n',SULamp(harm,ik))
    fprintf(fid,'%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n',SULfase(harm,ik))
    ngrid2=ngrid2+1;
end
 
for ik=1:length(east)
    fprintf(fid,'%5.0f%5.0f%5.0f%5.0f%10.5f\n',east(ik,:),emeaneast(ik))
    fprintf(fid,'%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n',LESTEamp(harm,ik))
    fprintf(fid,'%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n',LESTEfase(harm,ik))
    ngrid2=ngrid2+1;
end
 
for ik=1:length(north)
    fprintf(fid,'%5.0f%5.0f%5.0f%5.0f%10.5f\n',north(ik,:),emeannorth(ik))
    fprintf(fid,'%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n',NORTEamp(harm,ik))
    fprintf(fid,'%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n',NORTEfase(harm,ik))
    ngrid2=ngrid2+1;
end

fclose(fid)
fprintf(1,'ngrid2 = %10.0f',ngrid2);