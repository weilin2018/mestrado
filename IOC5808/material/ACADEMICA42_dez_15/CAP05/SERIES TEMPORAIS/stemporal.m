function Y = stemporal
load serie1.dat
%serie1 = serieC-2-jan; clear serieC-2-jan;

% ESTE PROGRAMA UTILIZA AS SUBROTINAS AUXILIARES QUE ESTÃO 
% GRAVADAS NESTA PASTA: "ROTMAJAX", "DVPVEC  e "SMOOTH"

clc
str= ...                                                             
    ['                                                     '
     '                                                     '
     ' O arq. serie1.dat está divido em 02 trechos:        '  
     '                                                     '
     ' 1) início: 22/08 às 08,0 h, término: 22/08 às 20,5 h'
     ' 2) início: 28/08 às 04,5 h, término: 28/08 às 17,0 h'
     '                                                     '
    ];
disp(str);                                               
str= ...                                                             
    ['                                           '
     '                                           '
     ' Deseja analisar algum trecho específico?  '
     '                                           '
    ];
disp(str);                                               
R = input(' sim(s)  ou  não(n)? ','s');

while R=='s'
str=['   '];
disp(str);                                               
Resp = input(' Qual trecho (de 1 a 2)?  ');
if Resp== 1
iniser=1;finser=286;
diaini=22; mesini=08; horaini=08.0; diafin=22; mesfin=08; horafin=20.5;
elseif Resp == 2
iniser=287;finser=572;
diaini=28; mesini=08; horaini=04.5; diafin=28; mesfin=08; horafin=17.0;
end

serie(:,1)=serie1(iniser:finser,1);%time
serie(:,2)=-serie1(iniser:finser,2);%prof e atenção pela mudança do sinal (+ ou -) na profundidade na série 2
serie(:,3)=serie1(iniser:finser,3);%temp
serie(:,4)=serie1(iniser:finser,4);%sali
serie(:,5)=serie1(iniser:finser,5);%sigt
serie(:,6)=serie1(iniser:finser,6);%cpNS é a mais impoprtante
serie(:,7)=serie1(iniser:finser,7);%cpEW 

analise_da_serie(serie,diaini,mesini,horaini,diafin,mesfin,horafin)
clear serie;
clc
str= ...                                                             
    ['                                           '
     '                                           '
     ' Deseja analisar outro trecho?             '
     '                                           '
    ];
disp(str);                                               
R = input(' sim(s)  ou  não(n)?','s');
end

function analise_da_serie(serie,diaini,mesini,horaini,diafin,mesfin,horafin)
time=serie(:,1);
prof=serie(:,2);
temp=serie(:,3);
sali=serie(:,4);
sigt=serie(:,5);
cpNS=serie(:,6);
cpEW=serie(:,7);
Npto = max(size(serie));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Escolha do Nível (profundidade) adimensional            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str=['   ']; disp(str); 
nivel = input(' Qual o nível (profundidade) desejada? ');

%Profundidade adimensiona

k=0;
for i = 1 : Npto
if prof(i) == nivel
k = k+1;   
sernivel(k,1)= k;
sernivel(k,2)= time(i);
sernivel(k,3)= prof(i);
sernivel(k,4)= temp(i);
sernivel(k,5)= sali(i);
sernivel(k,6)= sigt(i);
sernivel(k,7)= cpNS(i);
sernivel(k,8)= cpEW(i);
end
end

Nkpto = max(size(sernivel));

indice=sernivel(:,1);
time  =sernivel(:,2);
prof  =sernivel(:,3);
temp  =sernivel(:,4);
sali  =sernivel(:,5);
sigt  =sernivel(:,6);
cpNS  =sernivel(:,7);
cpEW  =sernivel(:,8);

sigt=sigt-1000.;

nuobs=size(indice,1);
nuobs2=nuobs/2;
dt_amos = 1;
famos=60/dt_amos;
%t_hor=indice/famos;
%t_dias=t_hor/24;

time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Séries Temporais %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%    Velocidade    %%%%%%%%%%%
yzero(1:nuobs)=0;

%thetad=-23.0*(pi./180.)

[cpNS,cpEW,thetad]=ROTMAJAX(cpNS,cpEW) ;

figure

cpNS = smooth(cpNS,6);plot(time,cpNS,'b');

%subplot(2,1,1)
plot(time,cpNS)
gtext('u','fontsize',14')

hold

plot(time,yzero,'k')
titulo=[' u and v-Components , Depth: Z=', num2str(nivel)];
title(titulo)
ylabel('Velocity (m s^{-1})','fontsize',14)
axis([min(time) max(time) min(cpNS) max(cpNS)])

cpEW = smooth(cpEW,6);plot(time,cpEW,'b');

%subplot(2,1,2)
plot(time,cpEW,'r')
gtext('v','fontsize',14')
%hold
%plot(time,yzero,'k')
%hold
%titulo=[' Transversal Component, Depth: ', num2str(nivel),'m'];
%title(titulo)
xlabel('Time (hour)','fontsize',14)
%ylabel('Velocity (m s^{-1})')
%axis([min(time) max(time) min(cpNS) max(cpNS)]) % NS > EW

%ext='stmp';
%arquifig = [nomedoarquivo num2str(nivel) ext]; 
%eval(['print -dbitmap ' arquifig]); 

pause

%print -dbitmap e:\estuario-24-FEV\programas\cap05\series temporais\Serie_u_v.bmp


%%%%%%%%%%   S,T e SigmaT    %%%%%%%%%%%
figure
subplot(3,1,1)
temp= smooth(temp,6);plot(time,temp,'b');
plot(time,temp)
titulo=['Depth: Z=',num2str(nivel)];
title(titulo)
ylabel('Temperature (^{o}C)','fontsize',14)
axis([min(time) max(time) min(temp) max(temp)])

subplot(3,1,2)
sali= smooth(sali,6);plot(time,sali,'b');
plot(time,sali)
ylabel('Salinity','fontsize',14)
axis([min(time) max(time) min(sali) max(sali)])

subplot(3,1,3)
sigt= smooth(sigt,6);plot(time,sigt,'b');
plot(time,sigt)
ylabel('Sigma-t (kg m^{-3})','fontsize',14)
xlabel('Time (hour)','fontsize',14)
axis([min(time) max(time) min(sigt) max(sigt)])

%ext='stmp';
%arquifig = [nomedoarquivo num2str(nivel) ext]; 
%eval(['print -dbitmap ' arquifig]); 
pause

%print -dbitmap e:\estuario-24-FEV\programas\cap05\series temporais\Seriehidro.bmp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  INTERPOLAÇAO %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
str=['']; disp(str); 
Resp2 = input(' Deseja interpolar, sim(s) ou não(n)? ','s');
if Resp2=='s'
timedt=min(time):.1:max(time);
tempdt=interp1(time,temp,timedt);
salidt=interp1(time,sali,timedt);
sigtdt=interp1(time,sigt,timedt);
tempdt=interp1(time,temp,timedt);
cpNSdt=interp1(time,cpNS,timedt);
cpEWdt=interp1(time,cpEW,timedt);

figure
subplot(2,1,1)
plot(timedt,cpNSdt,'g')
hold
plot(time,cpNS,'+r')
plot(time,yzero,'k')
hold
titulo=[' Longitudinal Component, Depth: Z=', num2str(nivel)];
title(titulo)
xlabel('Time (hour)')
ylabel('Velocity (m s^{-1})')
axis([min(time) max(time) min(cpNS) max(cpNS)]) % NS>EW

subplot(2,1,2)
plot(timedt,cpEWdt,'g')
hold
plot(time,cpEW,'+r')
plot(time,yzero,'k')
hold
titulo=['Longitudinal Component, Depth: Z=', num2str(nivel)];
title(titulo)
ylabel('Velocity (m s^{-1})')
axis([min(time) max(time) min(cpNS) max(cpNS)])

%ext='stmp';
%arquifig = [nomedoarquivo num2str(nivel) ext]; 
%eval(['print -dbitmap ' arquifig]); 
pause

%%%%%%%%%%   S,T e SigmaT    %%%%%%%%%%%
clc
figure
subplot(3,1,1)
plot(timedt,tempdt,'g')
hold
plot(time,temp,'+r')
hold
titulo=['Depth: ', num2str(nivel)];
title(titulo)
ylabel('Temperature (^{o}C)')
axis([min(time) max(time) min(temp) max(temp)])

subplot(3,1,2)
plot(timedt,salidt,'g')
hold
plot(time,sali,'+r')
hold
ylabel('Salinity')
axis([min(time) max(time) min(sali) max(sali)])

subplot(3,1,3)
plot(timedt,sigtdt,'g')
hold
plot(time,sigt,'+r')
hold
ylabel('Sigma-t (kg m{-3})')
xlabel('Time (hour)')
axis([min(time) max(time) min(sigt) max(sigt)])

%ext='stmp';
%arquifig = [nomedoarquivo num2str(nivel) ext]; 
%eval(['print -dbitmap ' arquifig]); 
pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Diagrama stikplot  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
figure
feather(cpEW,cpNS)
%axis equal

titulo=[' "stick-plot" Diagram, Depth: Z=',num2str(nivel)];
title(titulo)

xlabel('Time (hour)','fontsize',14)
ylabel('(m s^{-1})','fontsize',14)

%ext='stkp';
%arquifig = [nomedoarquivo num2str(nivel) ext]; 
%eval(['print -dbitmap ' arquifig]); 
pause

print -dbitmap e:\Academica42\cap05\series temporais\stick.bmp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Diagrama rosa das correntes %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
figure
%subplot(1,2,1)
compass(cpEW,cpNS)

titulo=['Current Rose, velocity (m s^{-1}), Depth: Z=',num2str(nivel)];
title(titulo)

%subplot(1,2,2)
%dir=(atan2(cpNS,cpEW))*180/pi;
%rose(dir)
%title('histograma angular das correntes')

%ext='rosa';
%arquifig = [nomedoarquivo num2str(nivel) ext]; 
%eval(['print -dbitmap ' arquifig]); 
pause

print -dbitmap e:\Academica42\cap05\series temporais\rose.bmp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Diagrama vetorial progressivo %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
figure
dvpvec(cpEW,cpNS,60)

%titulo=['vetor progressivo da corrente, (+) início e (o) a cada 24h, prof.: ',num2str(nivel),'m'];
titulo=['Progressive Vectorial Diagram, Depth: Z=',num2str(nivel)];
title(titulo)

%ext='dvp';
%arquifig = [nomedoarquivo num2str(nivel) ext]; 
%eval(['print -dbitmap ' arquifig]); 
pause

print -dbitmap e:\Academica42\cap05\series temporais\vetorialProgr.bmp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Histogramas %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
figure
hist(cpNS,20)

titulo=['Longitudinal Component, Depth: Z=',num2str(nivel)];
title(titulo)

xlabel('Velocity (m s^{-1})','fontsize',14)
ylabel('Occurrence Numbers','fontsize',14)

%ext='hsty';
%arquifig = [nomedoarquivo num2str(nivel) ext]; 
%eval(['print -dbitmap ' arquifig]); 
pause

print -dbitmap e:\Academica42\cap05\series temporais\histoU.bmp

figure
hist(cpEW,20)

titulo=['Transversal Component, Depth: Z=',num2str(nivel)];
title(titulo)

xlabel('Velocity (m s^{-1})','fontsize',14)
ylabel('Occurrence Numbers','fontsize',14)

%ext='hstx';
%arquifig = [nomedoarquivo num2str(nivel) ext]; 
%eval(['print -dbitmap ' arquifig]); 
pause

print -dbitmap e:\Academica42\cap05\series temporais\histoV.bmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Estatísicas das séries           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stat(1,1)=min(temp);
stat(2,1)=max(temp);
stat(3,1)=mean(temp);
stat(4,1)=median(temp);
stat(5,1)=std(temp);

stat(1,2)=min(sali);
stat(2,2)=max(sali);
stat(3,2)=mean(sali);
stat(4,2)=median(sali);
stat(5,2)=std(sali);

stat(1,3)=min(sigt);
stat(2,3)=max(sigt);
stat(3,3)=mean(sigt);
stat(4,3)=median(sigt);
stat(5,3)=std(sigt);


stat(1,5)=min(cpNS);
stat(2,5)=max(cpNS);
stat(3,5)=mean(cpNS);
stat(4,5)=median(cpNS);
stat(5,5)=std(cpNS);

stat(1,4)=min(cpEW);
stat(2,4)=max(cpEW);
stat(3,4)=mean(cpEW);
stat(4,4)=median(cpEW);
stat(5,4)=std(cpEW);

%save  estatist.dat  stat -ascii

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


























