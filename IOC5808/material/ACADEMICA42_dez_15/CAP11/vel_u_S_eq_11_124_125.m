  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% PROGRAMA vel_u_S_eq_11_124_125.m                        %
% ALESSANDRO LUVIZON BÉRGAMO, LUIZ BRUNER DE MIRANDA      %
% E FERNANDO PINHEIRO ANDUTTA                             % 
% SIMULAÇÃO DE PERFIS DE VELOCIDADE e de SALINIDADE       %  
% TEORIA x DADOS EXPERIMENTAIS DE MARÉS DE QUADRATURA     %
% E DE SIZÍGIA                                            %
% OS DADOS ESTÃO NOS ARQUIVOS Qres_pia.dat e Siz_pia.dat  %
% ESSES ARQUIVOS SÃO LIDOS NO INÍCIO DO PROCESSAMENTO     %  
% TEORIA Hansen & Rattray (1965) E DETALHES               %
% NO TRABALHO ORIGINAL OU EM Miranda et al. (2012)        %      
% (eqs. 11.124 e % 11.125, p. 408)                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dinâmica gerada pelas forçantes: descarga fluvial e efeito baroclínico;
% Dissipação de energia por atrito interno e atrito máximo no fundo;
% Orientação do eixo NS # v>0 e v<0 para sul e norte, respectivamente; 
% vento (>0 para o sul)e gradiente longitudinal de salinidade (densidade):
% Unidades do SI
% Validação da teoria com o método Skill (Wilmott, 1981)
% Adaptado por Fernando Pinheiro Andutta para perfis verticais

clear
clc
close all


quadra_sizigia= menu('Wich station to study?','quadratura','sizigia');


%%% loading files
if quadra_sizigia==1;
load Qres_pia.txt
M = Qres_pia;
Z=M(:,1);
S=M(:,2);
V=M(:,3);
else
load Sres_pia.txt
M = Sres_pia;
Z=M(:,1);
S=M(:,2);
V=M(:,3);
end
    
    
 %%%% Setting free parameters for simulation   
 %%%% Indicando do parâmetros livres para simulação
 %%%% E escolhendo os dados de maré de quadratura ou de
 %%%% sizígia (Qres_pia.dat e Siz_pia.dat, respectivamente)
 
if quadra_sizigia==1;
    
ar= 0.009;%input('Vel. da descarga fluvial (m/s)- Ua:');
br= 0.004;%input('Coef. de viscosidade turbulenta vertical m.m/s - Nz:');
kh=1000;%input('Coef. de difusão turbulenta horizontal m.m/s - Kh:');
kz= 0.00015;%input('Coef. de difusão turbulenta vertical m.m/s - Kz:');
cr=0.02;%input('Tensao de cisalhamento do vento (N/m2) - Tw:');
dr=33;%input('Maximo de salinidade (boca) - S_bc:');
er=1;%input('Minimo de salinidade (cabeceira) - S_cb:');
fr=20000;%input('Delta x (m) - D:');
gr=11;%input('Prof. media local (m) - H:');
so=26.5;%input('Salinidade cc:');
xx=17000;%input('distância longitudinal(m) - dist:');
nu=0.85;%input('coeficiente constante, ni, adimencional:');
else
ar=0.003;%input('Vel. da descarga fluvial (m/s)- Ua:');
br= 0.009;%input('Coef. de viscosidade turbulenta vertical m.m/s - Nz:');
kh=1000;%input('Coef. de difusão turbulenta horizontal m.m/s - Kh:');
kz= 0.000025;%input('Coef. de difusão turbulenta vertical m.m/s - Kz:');
cr=0.02;%input('Tensao de cisalhamento do vento (N/m2) - Tw:');
dr=33;%input('Maximo de salinidade (boca) - S_bc:');
er=1;%input('Minimo de salinidade (cabeceira) - S_cb:');
fr=20000;%input('Delta x (m) - D:');
gr=11;%input('Prof. media local (m) - H:');
so=27.8;%input('Salinidade cc:');
xx=17000;%input('distância longitudinal(m) - dist:');
nu=0.85;%input('coeficiente constante, ni, adimencional:');
end

%Z orientado positivamente para cima com origem na superfície -1<Z<0

z=zeros(11,1);

g=9.8;                      % aceleracao da gravidade (m/s2)

% calculo do perfil e do grad long de densidade usando: ro=ro0*(1+B*S)

   B=0.77*10^(-3);          % coef. de contracao salina
   ro0=1000;                % densidade
   
ro_bc=ro0*(1+(B*dr));       % densidade na boca
ro_cb=ro0*(1+(B*er));       % densidade na cabeceira

Dro=(ro_bc-ro_cb);          % gradiente long. de densidade

ro=ro0*(1+(B*S));           % perfil medio de densidade



% Formulação matemática - perfil de velocidade  
% longitudinal (v) - pois neste exemplor o canal estuarino
% está orientado no sentido Norte-Sul
% e orientado positivamente para o Sul 

Drodx=Dro./(fr); 


Ug=nu./48.*(g.*(gr).^(3)./(br)).*((1./ro).*(Dro./(fr))).*(1-9.*(Z.^2)-8.*(Z.^3))

Udf=(3./2.*(ar).*(1-1.*(Z.^2)))

Uv=(1./4.)*((cr.*gr)./(ro.*br).*(1+(4.*Z)+3.*(Z.^2)))

%U=Ug+Udf+Uv;
%Soma dos modos baroclínico, descarga fluvial e da tensão do vento

U=(Ug+Udf+Uv)
ures=mean(U)

% Formulação matemática - perfil da salinidade - (Sal)
% parâmetros 


coef1=(nu*ar*xx)/kh

coef2=nu*(ar*gr)^2/(kh*kz)

rox=(g.*(gr).^(3)./(br.*ar)).*((1./ro).*(Dro./(fr)))

S1=(-Z-0.5)-1/2.*((Z.^2)-1./3.)-1/2.*(-2.*Z-3./2.*(Z.^2)+1./4.*(Z.^4))
S2=(nu./48.)*(rox).*(1/2.*(Z.^2)-3./4.*(Z.^4)-2./5.*(Z.^5))
%S2=nu./48.*(g.*(gr).^(3)./(br.*ar)).*((1./ro).*(Dro./(fr))).*(1/2.*(Z.^2)-3./4.*(Z.^4)...
%   -2./5.*(Z.^5))
S3=(1./4.)*((cr)*(gr)./ro.*(br.*ar)).*(1/2.*(Z.^2)+2./3.*(Z.^3)+1./4.*(Z.^4))
S4=((11./40.)-(1./80.)*((cr).*(gr)./ro.*(br.*ar))-nu./576.*(g.*(gr).^(3)./(br*ar)).*...
   ((1./ro).*(Dro./(fr))))

%Calculo do perfil teórico da salinidade 
%Salinidade=Sal

Sal=(so)*(1+coef1+(coef2)*(S1+S2+S3+S4))
Sres=mean (Sal)

% Validação: Diferenças teoria - experimento

deltaS=Sal-S
deltau=U-V

% parametros

Vel_desc_Ua=ar;
Coef_visc_turb_Nz=br;
Coef_visc_turb_Kh=kh
Tens_cis_vento_Tw=cr;
Max_sal_Sbc=dr;
Min_sal_Scb=er;
Dist_D=fr;
Prof_media_H=gr;


% resultados
Ug=Ug;
Udf=Udf;
Uv=Uv;
U=U;


%%% Calculating the skill
%%% Cálculo do parâmetro de validação
%%% pelo método Skill

U_cal = U;
U_obs = V;
S_cal = Sal;
S_obs = S;

mean_U_obs = nanmean(U_obs);
mean_U_obs = mean_U_obs*ones(11,1);
mean_S_obs = nanmean(S_obs);
mean_S_obs = mean_S_obs*ones(11,1);

sum_vel = nansum(abs(U_cal-U_obs).^2);
diff_vel1 = abs(U_cal-mean_U_obs);
diff_vel2 = abs(U_obs-mean_U_obs);
diff_vel = nansum((diff_vel1+diff_vel2).^2);

sum_sal = nansum(abs(S_cal-S_obs).^2);
diff_sal1 = abs(S_cal-mean_S_obs);
diff_sal2 = abs(S_obs-mean_S_obs);
diff_sal = nansum((diff_sal1+diff_sal2).^2);

matrix_one = ones(11,1);

skill_vel = nanmean(matrix_one-(sum_vel./diff_vel));
skill_sal = nanmean(matrix_one-(sum_sal./diff_sal));

skill_vel = 100*skill_vel;
skill_vel = round(skill_vel);
skill_vel = skill_vel/100;

skill_sal = 100*skill_sal;
skill_sal = round(skill_sal);
skill_sal = skill_sal/100;

% Salvando e arquivando comparativamente os resultados
% teóricos e experimentais


data = [Z U_obs U_cal S_obs S_cal];
if quadra_sizigia==1;
fnome=['save -ASCII out_quadra_Uobs_Ucal_Sobs_Scal.dat data'];
eval(fnome);
else
fnome=['save -ASCII out__sizigia_Uobs_Ucal_Sobs_Scal.dat data'];
eval(fnome);
end


if quadra_sizigia==1;
fnome=['save -ASCII Skill_vel_quadra.dat skill_vel'];
eval(fnome);
fnome=['save -ASCII Skill_sal_quadra.dat skill_sal'];
eval(fnome);
else
fnome=['save -ASCII Skill_vel_sizigia.dat skill_vel'];
eval(fnome);
fnome=['save -ASCII Skill_sal_sizigia.dat skill_sal'];
eval(fnome);

end

%%%% velocidades Ug (modo baroclínico)
%%%% Udf (modo da descarga fluvial) 
%%%% Uv (modo gerado pela tensão do vento)
%%%% U denota a velocidade resultante
%%%% U=Ug+Udf+Uv

figure
plot(Ug,Z,'c','linewidth',2)
hold on
plot(Udf,Z,'r','linewidth',2)
hold on
plot(Uv,Z,'g','linewidth',2)
hold on
plot(U,Z,'k','linewidth',2)
line(z,Z)
legend('Vg','Vdf','Vv','v_{c}',4)
xlabel('Componente, v (m s^{-1})','fontsize',16);
ylabel('Profundidade, Z','fontsize',16);

if quadra_sizigia==1
print -djpeg Fig_01_quadra.bmp
else
print -djpeg Fig_01_sizigia.bmp  
end

%%%% velocidades U(teórica) e V(experimental)

figure
plot(U_cal,Z,'k--','linewidth',2)
hold on
plot(U_obs,Z,'k','linewidth',2)
line(z,Z)
legend('Teoria','Observação',4)
xlabel('Componente longitudinal(m s^{-1})','fontsize',16);
ylabel('Profundidade, Z','fontsize',16);
gtext('Skill =')
gtext(num2str(skill_vel));


if quadra_sizigia==1
print -djpeg Fig_02_quadra.bmp
else
print -djpeg Fig_02_sizigia.bmp  
end

%%%%%%%%%%%% Sal (theoretical) and S (experimental)
figure
plot(S_cal,Z,'k--','linewidth',2)
hold on
plot(S_obs,Z,'k','linewidth',2)
%line(z,Z)
legend('Teoria','Observação',3)
xlabel('Salinidade (^o /_o_o)','fontsize',16);
ylabel('Profundidade, Z','fontsize',16);
gtext('Skill =')
gtext(num2str(skill_sal));


if quadra_sizigia==1
print -djpeg Fig_03_quadra.bmp
else
print -djpeg Fig_03_sizigia.bmp  
end



