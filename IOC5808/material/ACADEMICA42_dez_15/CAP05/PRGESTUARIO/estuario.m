%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA estuario.m                                         %
% TEM COMO PROGRAMAS AUXILIARES OS PROGRAMAS "hidro.m         % 
% E "decomp.m" J� APRESENTADOS E QUE GERAM OS ARQUIVOS        %
% DE DADOS HIDROGR�FICOS "ts_.dat "vel_.dat" e zl_.dat        % 
% J� GRAVADOS ANTERIORMENTE E DEVEM SER CARREGADOS COM        % 
% COM O COMANDO "LOAD" NO IN�CIO DO PROCESSAMENTO             % 
% CALCULA VALORES M�DIOS NO TEMPO E NO ESPA�O                 %
% COM INTEGRA��ES NUM�RICAS COM A TEORIA                      %
% APRESENTADA EM Miranda et al. (2002, 2012 2nd ed.)          %
% Equa��es 5.19 - p. 169 e 5.23 p. - 171                      % 
% Calcula os componentes do transporte advectivo de sal       %  
% Mais detalhes da teoria em Miranda et al. (2012) p.184 a 190%
% PARA MAIS DETALHES, CONSULTAR                               %
% B�rgamo, Miranda & Corr�a et al. (2002)-                    % 
% Relat.t�c. inst. oceanogr. S�o Paulo,                       % 
% no. 49 pp.1-16.                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	% Adapta��o do programa dinest.m, para o c�lculo dos % 
	%    componentes do transporte e par�metros do     %
	%             sistema estuarino                    %
    %      Alunos da IOF-827 do ano de 1996            %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Exemplo das vari�veis do programa estuario%%%%%%%%%%%%%% 
	%estuario(vl_22_08.dat,ts_22_08.dat,zl_22_08.dat,26,11) %
    %26=n�mero de esta��es, 11=n�mero de n�veis             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	function  Y=estuario(vel,ts,h,nest,npro)
	global va sa vs ss sss vss pa vdoce
	
   diary G:\Academica42_maio_14\cap05\prgestuario\estuario.txt

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%   vel   = matriz de velocidade longitudinal          %
	%   ts   = matriz de temperatura e salinidade          %
	%   h   = matriz-linha de prof. local instant�nea      %	
	%   nest = n�mero de esta��es                          %
	%   npro = n�mero de profundidades em cada esta��o     %
    %   Esses resulados s�o obtidos com os programas       %
    %   hidro.m e decomp.m                                 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%   Dados de Temperatura e Salinidade     %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	prof=-ts(:,2); % Profundidade adimensional sinal (-) devido a orienta��o de OZ
	temp=ts(:,3); % Temperatura em graus Celsius
	salt=ts(:,4); % Salinidade em EPS
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	% ATEN��O: O PROGRAMA USA A SUBROTINA "SW_DENS"    %
    % DO MORGAN (1994)- CSIRO MARINE LABORATORIES,     %
    % AUSTR�LIA, 222, 28 P. PARA CALCULAR A DENSIDADE  %
    % DA �GUA ESTUARINA � PRESS�O HIDROST�TICA ZERO    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	dens=sw_dens(salt,temp,0);         % Densidade em (kg/m^3)
	%dens=sw_dens(salt,temp,prof);     % Densidade em (kg/m^3)

	sigma_t=(dens-1000);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%C�lculo da Prof. Adm. Z=z/h %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	np=max(size(prof));
	k=1;

	for i=1:npro:np-(npro-1);
		maxp=max(abs(prof(i:i+(npro-1))));
 		paa=prof(i:i+(npro-1))/(maxp);
  	for j=1:npro
    		pats(k)=paa(j);
    		k=k+1;
  	end
    end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%     Velocidade                            %
    %Usa os dados de velocidade                 %
    %processados pelo decomp.m                  %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   u=vel(:,2); % velocidade longitudinal em (m/s)
   prov=-vel(:,1);
   %prov=vel(:,1);

   
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% C�lculo da velocidade nas prof. Adm. Z        %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	nv=max(size(u));
	k=1;
	for i=1:(npro):np-(npro-1);
 		maxpv=max(abs(prov(i:i+(npro-1))));
 		paav=prov(i:i+(npro-1))/(maxpv);
   	for j=1:(npro)
    		pauv(k)=paav(j);
    	k=k+1;
   	end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Interpola��o Pelo M�todo "Cubic Spline" todas  % 
    % vari�veis S, T, densidade e velocidade         % 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	pa=0:.1:1;
	k=0;
 	for i=1:(npro):np-(npro-1);
 		k=k+1;
  		ta=spline(pats(i:i+((npro-1))),temp(i:i+((npro-1))),pa);
  		sa=spline(pats(i:i+((npro-1))),salt(i:i+((npro-1))),pa);
  		da=spline(pats(i:i+((npro-1))),dens(i:i+((npro-1))),pa);
  		siga=spline(pats(i:i+((npro-1))),sigma_t(i:i+((npro-1))),pa);
  		va=spline(pauv(i:i+((npro-1))),u(i:i+((npro-1))),pa);
  
   	for j=1:11
     		estac(j,k)=k;
     		tta(j,k)=ta(j);		%temperatura
     		ssa(j,k)=sa(j);		%salinidade
		    daa(j,k)=da(j);		%densidade
     		sigaa(j,k)=siga(j);	%sigma-t
     		vva(j,k)=va(j);		%velocidade longitudinal
   	end
 	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%  hidro: Esta��o  Prof.Adm.  Temp.  Salin. Dens.  Sigma_T   Veloc. %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	k=0;
	for i=1:nest
 	for j=1:11
  		k=k+1;
  		psva(k,1)=estac(j,i);
  		psva(k,2)=-pa(j);
  		psva(k,3)=tta(j,i);
  		psva(k,4)=ssa(j,i);
		psva(k,5)=daa(j,i);
  		psva(k,6)=sigaa(j,i);
  		psva(k,7)=vva(j,i);
 	end
	end

	linha=zeros(1,7);

	hidro=[psva;linha];


	% C�lculo da diferen�a entre densidade de superficie e fundo (Delta-ro)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   %p=hidro(:,2);
    p=-hidro(:,2);
	d=hidro(:,5);

		x=find(p==0);  
	
   tam=length(x);      
   
   for i=1:tam        
     if i~=tam       
	I=num2str(i);
		

	droi=['dro',I];
	delta=[droi,'=d(x(',I,'+1)-1)-d(x(',I,'));'];
	eval(delta);

      end
    end	

	% Organiza uma matriz para delta_ro %
	


    for i=1        
 
   	DRO=['Dro=[dro',num2str(i),' dro',num2str(i+1),' dro',num2str(i+2),' dro',num2str(i+3)...
,' dro',num2str(i+4),' dro',num2str(i+5),' dro',num2str(i+6),' dro',num2str(i+7),' dro',num2str(i+8)...
,' dro',num2str(i+9),' dro',num2str(i+10),' dro',num2str(i+11),' dro',num2str(i+12),' dro',num2str(i+13)...
,' dro',num2str(i+14),' dro',num2str(i+15),' dro',num2str(i+16),' dro',num2str(i+17),' dro',num2str(i+18)...
,' dro',num2str(i+19),' dro',num2str(i+20),' dro',num2str(i+21),' dro',num2str(i+22),' dro',num2str(i+23)...
,' dro',num2str(i+24),' dro',num2str(i+25),' ];'];
	eval(DRO);
  									     
     end


	Delta_ro=DRO



	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% C�lculo das m�dias por integra��o num�rica        %
    % ao longo da coluna de �gua                        %
    % Miranda et al. (2002, 2011 2nd) eq. 5.19 - p. 169 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
	%    M�dia no espa�o     % 
	%%%%%%%%%%%%%%%%%%%%%%%%%%

	for j=1:(nest)
 		somas=0;
 		somav=0;
 		somavs=0;
		somast=0;
  	for i=2:10
   		somas=somas+ssa(i,j);
   		somav=somav+vva(i,j);
   		somavs=somavs+vva(i,j)*ssa(i,j);
		somast=somast+sigaa(i,j);
  	end
 		sme(j)=(somas+ssa(1,j)/2+ssa(11,j)/2)*(0.1);
 		vme(j)=(somav+vva(1,j)/2+vva(11,j)/2)*(0.1);
 		vsme(j)=(somavs + vva(1,j)*ssa(1,j)/2 + vva(11,j)*ssa(11,j)/2)*(0.1);
		rome(j)=(somast+daa(1,j)/2+daa(11,j)/2)*(0.1);
	end

	media_esp_S=sme
	media_esp_U=vme
	media_esp_US=vsme
	media_esp_ro=rome
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% C�lculo das m�dias por integra��o num�rica           %
    % no intervalo de tempo de um ou mais ciclos de mar�   %
    % Miranda et al. (2002, 2011 2nd) eq. 5.23 - p. 171    %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	somas=0;
	somav=0;
	somah=0;
	somavs=0;
  	for i=2:(nest-1)
   		somas=somas+sme(i);
   		somav=somav+vme(i);
   		somah=somah+h(i);
   		somavs=somavs+vsme(i)*h(i);
  	end
		sa=(somas+sme(1)/2+sme(nest)/2)/(nest-1);
		va=(somav+vme(1)/2+vme(nest)/2)/(nest-1);
		ha=(somah+h(1)/2+h(nest)/2)/(nest-1);

	media_temp_Z=ha % m�dia temporal da profundidade = ha  
	media_esptemp_U=va % m�dia temporal e espacial da velocidade = sa
	media_esptemp_S=sa % m�dia temporal e espacial da salinidade = sa

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% C�lculo do transporte m�dio de volume resultante %
    % por unidade de largura do canal estuarino        %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	transp_medio_resultante=(somavs+vsme(1)*h(1)/2+vsme(nest)*h(nest)/2)/(nest-1);


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% C�lculo dos componentes do transporte advectivo de sal %
    % Mais detalhes em Miranda et al. (2002, 2011 2nd) -     % 
    % eqs. 5.37 e seguintes - p. 186                         %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% m�dia temporal e espacial da velocidade (residual) aprox.= Qf/Area: va %                         
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% componente barotr�pica associada � mar� e, portanto, peri�dica: vt %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	vt=vme-va;
	st=sme-sa;
	ht=h-ha;

	U_barotrop=vt
	S_barotrop=st
	oscil_mare=ht

	
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% componente barocl�nica - circula��o gravitacional, estacion�ria: vs %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for i=1:11
 		somas=0;
 		somav=0;
  	for j=2:(nest-1)
  		 somas=somas+ssa(i,j);
  		 somav=somav+vva(i,j);
  	end
 		vss(i)=(somav+vva(i,1)/2+vva(i,nest)/2)/(nest-1);
 		sss(i)=(somas+ssa(i,1)/2+ssa(i,nest)/2)/(nest-1);
	end
		vs=vss-va;
		ss=sss-sa;

	U_baroclin=vs
	S_baroclin=ss

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%    flutua��es turbulentas: s�o calculadas na parte do transporte   %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Parcelas componentes do transporte total %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



	% Descarga de �gua doce %
	%%%%%%%%%%%%%%%%%%%%%%%%%

	descarga_doce=va*ha*sa;


	% Deriva de Stokes %
	%%%%%%%%%%%%%%%%%%%%

	deriva_stokes=sa*(sum(ht(2:(nest-1)).*vt(2:(nest-1)))+ht(1)*vt(1)/2+ht(nest)*vt(nest)/2)/(nest-1);



	% Correla��o de mar� %
	%%%%%%%%%%%%%%%%%%%%%%

	correl_mare=ha*(sum(st(2:(nest-1)).*vt(2:(nest-1)))+st(1)*vt(1)/2+st(nest)*vt(nest)/2)/(nest-1);



	% Circula��o gravitacional %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	circ_gravit=ha*(sum(vs(2:10).*ss(2:10))+vs(1)*ss(1)/2+vs(11)*ss(11)/2)/(10);



	% Difus�o turbulenta %
	%%%%%%%%%%%%%%%%%%%%%%

	vxvx=vva-va;
	sxsx=ssa-sa;
	for i=1:11
	for j=1:nest
		vxv(i,j)=vxvx(i,j)-vt(j);
		sxs(i,j)=sxsx(i,j)-st(j);
	end
	end

	for j=1:nest
	for i=1:11
		vlinha(i,j)=vxv(i,j)-vs(i);
		slinha(i,j)=sxs(i,j)-ss(i);
	end
	end


	%  M�dia no espa�o para a difus�o turbulenta
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for j=1:(nest)
 		soma=0;
  	for i=2:10
   		soma=soma+vlinha(i,j)*slinha(i,j);
  	end
		mdte(j)=(soma+vlinha(1,j)*slinha(1,j)/2+vlinha(11,j)*slinha(11,j)/2)*(0.1);
	end


	%  M�dia total para a difus�o turbulenta
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	soma=0;
  	for i=2:(nest-1)
   		soma=soma+mdte(i);
  	end
		mdtet=(soma+mdte(1)/2+mdte(nest)/2)/(nest-1);

	
	difusao_turb=mdtet*ha;



	% Dispers�o da mar� %
	%%%%%%%%%%%%%%%%%%%%%

	sddd=sum(vt(2:(nest-1)).*st(2:(nest-1)).*ht(2:(nest-1)));
	ssdd=vt(1)*st(1)*ht(1)/2 + vt(nest)*st(nest)*ht(nest)/2;

	dispersao_mare=(sddd+ssdd)/(nest-1);



	% Circula��o residual %
	%%%%%%%%%%%%%%%%%%%%%%%

	cddd=sum(st(2:(nest-1)).*ht(2:(nest-1)));
	ccdd=st(1)*ht(1)/2 + st(nest)*ht(nest)/2;

	circ_resid=va*(cddd+ccdd)/(nest-1);


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Transporte total de sal por unidade da largura do canal estuarino %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  transp_total=descarga_doce+deriva_stokes+correl_mare+circ_gravit+difusao_turb+dispersao_mare+circ_resid;


	descarga_doce
	deriva_stokes
	correl_mare
	circ_gravit
	difusao_turb
	dispersao_mare
	circ_resid

	transp_total
	transp_medio_resultante

	%erro=(transp_medio_resultante-transp_total)/transp_medio_resultante


	
	% C�lculo do n�mero de Richardson por camada  RiL   %
    % Teoria em Miranda et al. (2002, 2012, 2nd ed)     %
    % p. 85 
	

	RiL=(9.8.*(abs(h))./vme.^2).*(Dro./rome)	%valores adimensionais
	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Resultados                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



	% Tempo ajustado de acordo com os hor�rios de coleta dos experimentos
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	  tempo=5.5:.5:18.; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Velocidade                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




   % Perfis verticais m�dios            
   %%%%%%%%%%%%%%%%%%%%%%%%%%
   
   figure
   
	plot(vss,-pa,'k--')	% vss = vel. media temporal (<u>)	
	hold			% vs  = comp. barocl. da vel. (us)	     	
	%plot(vs,-pa,'k:')	% va  = vel. de descarga de agua doce (ua)				
					
	for i=1:11
		vdoce(i)=va;
	end			% vs=vss-va : us=<u>-ua

   xx=[0,0];
   yy=[0,-1];
   plot(xx,yy,'k-');


	plot(vdoce,-pa,'k-.')

	gtext('u_{a}','fontsize',14)	% linha tra�o e ponto

	%gtext('us')	% linha pontilhada

	gtext('<u(Z)>','fontsize',14)	% linha tracejada
   
    gtext('u_{a}(m s^{-1})=','fontsize',14)
    I=num2str(va);
    texto=['gtext(','I',')'];
    eval(texto)
   
   %gtext('us=<u>-ua')

   %title('Gravitational circulation # spring tide')

   %title('Gravitational circulation # neap tide')

   ylabel('Depth, Z','fontsize',14)
   xlabel('u-Component (m s^{-1})','fontsize',14)
	%hold
	
   %pause

	%print -dbitmap c:\testecuri\perfilv_si.bmp

%    print -dbitmap e:\Academica42_maio_14\cap05\prgestuario\perfilu_qu.bmp


    % Perfis Eulerianos de velocidade
	% Varia��o no tempo da velocidade longitudinal na coluna d'�gua
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   figure 
   plot(-1,0)
   
      hold

	for i=1:nest
		plot(vva(:,i),-pa,'k')
        plot(xx,yy,'k-');
      
   end
   
     
   xlabel('u-Component(m s^{-1})','fontsize',14)
   ylabel('Depth, Z','fontsize',14)
   
	%title('Hourly variability of velocity # spring tide')
   
   %title('Hourly variability of velocity # neap tide')

   %hold

   %pause

	   
%    print -dbitmap e:\Academica42_maio_14\cap05\prgestuario\varu_qu.bmp

   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Salinidade                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	

	% Perfis verticais m�dios 
	%%%%%%%%%%%%%%%%%%%%%%%%%%
   figure
	plot(sss,-pa,'k')
	hold
   plot(sss,-pa,'k.')
   gtext('<S(Z)>','fontsize',14')
   gtext('S_{a}=','fontsize',14)
   
   Sa=sa*1000;
   Sar=round(Sa);
   Sa=Sar/1000;
      
   I=num2str(Sa);
   texto=['gtext(','I',')'];
   eval(texto)

	%title('Perfil m�dio temporal da salinidade')
	ylabel('Depth, Z','fontsize',14)
	xlabel('Salinity','fontsize',14)
   hold
   
   

   title('Mean salinity profile # spring tide')
   
   %title('Mean salinity profile # neap tide')

   
	ylabel('Depth, Z','fontsize',14)
	xlabel('Salinity','fontsize',14)
	hold

   %pause

   %print -dbitmap c:\testecuri\perfilS_si.bmp

%    print -dbitmap e:\Academica42_maio_14\cap05\prgestuario\perfilS_qu.bmp

    % Perfis Eulerianos de salinidade                  %                
	% varia��o no tempo da salinidade na coluna d'�gua %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   figure
	plot(min(min(ssa)),0)
	hold

	for i=1:nest
		plot(ssa(:,i),-pa,'k.')
		plot(ssa(:,i),-pa,'k')
	end

   %title('Hourly variability of salinity # spring tide')

   %title('Hourly variability of salinity # neap tide')

	xlabel('Salinity','fontsize',14)
	ylabel('Depth, Z','fontsize',14)
	
	%print -dbitmap c:\testecuri\varsal_si.bmp

%   print -dbitmap e:\Academica42_maio_14\cap05\prgestuario\varsal_qu.bmp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estratifica��o da coluna d'�gua - N�mero de Richardson por Camada,    %
% Mais detalhes em Miranda et al. (2002, 2012 2nd ed.) eq. 2.35 - p. 85 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure


dois=ones(1,26)*2;
vinte=ones(1,26)*20;
   
   semilogy(tempo,RiL)	
   
   hold on
   
   semilogy(tempo,RiL,'.')
   semilogy(tempo,dois,':')
   semilogy(tempo,vinte,':')
   
   %title('Richardson number, Ri_{L}# spring tide')
   
   title('Richardson number, Ri_{L}# neap tide')

      
   ylabel('Ri_{L}','fontsize',14)
   xlabel('Time (h)','fontsize',14)
   
	hold off
   

	%print -dbitmap c:\testecuri\riL_si.bmp
   
%    print -dbitmap e:\Academica42_maio_14\cap05\prgestuario\riL_qu.bmp


% A seguir o programa apresenta a varia��o temporal de v�rias propriedades:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Mar� x Velocidade m�dia na coluna d'�gua  = DERIVA DE STOKES           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   figure
   
	plotyy(tempo,ht,tempo,vme)
	
	% h=h(t)
	% vel. media esp.
	
   %title('Hourly tide variation and mean velocity in water column')
   xlabel('Time(h)','fontsize',14)  
   ylabel('Tide (m) - Mean velocity (m s^{-1})','fontsize',14)
   
   	gtext('h(t)','fontsize',14)	% azul
   	gtext('u(t)','fontsize',14)	% verde

	%print -dbitmap c:\testecuri\hv_si.bmp
   
%    print -dbitmap e:\Academica42_maio_14\cap05\prgestuario\hu_qu.bmp

   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Mar� x salinidade m�dia na coluna d'�gua             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   figure
   
	plotyy(tempo,ht,tempo,sme)
	
	% h=h(t)
	% sal. media esp.
	
   %title('Hourly tide variation and mean salinity in water column')
	xlabel('Time (h)','fontsize',14)
    ylabel('Tide(m) - Mean salinity','fontsize',14)
   
   	gtext('h(t)','fontsize',14) % azul
   	gtext('S(t)','fontsize',14)	  % verde

%     print -dbitmap e:\Academica42_maio_14\cap05\prgestuario\hS_si.bmp

   %print -dbitmap g:\Caravelas-A-EST\hS_qu.bmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  velocidade m�dia na coluna d'�gua x salinidade m�dia na coluna d'�gua = DIFUS�O DA MAR�  %         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   figure
   
	plotyy(tempo,vme,tempo,sme)
	
	% h=h(t)
	% sal. media esp.
	
   %title('Varia��o da velocidade m�dia na coluna de �gua e da salinidade m�dia na coluna de �gua')
	xlabel('Time(h)','fontsize',14)
   ylabel('esq:u(t)(m/s) - dir:sal.m�dia esp.','fontsize',14)
   
  	gtext('u(t)','fontsize',14) % azul
  	gtext('S(t)','fontsize',14) % verde

% 	print -dbitmap G:\Academica42_maio_14\cap05\prgestuario\uS_qu.bmp

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Componentes do transporte total de sal %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	tt=[ descarga_doce deriva_stokes correl_mare circ_gravit difusao_turb dispersao_mare circ_resid 0 transp_total transp_medio_resultante];
   
   figure	
   
   x=1:1:10;
	
    	bar(x,tt)
    
    	hold on
    
    	y=zeros(1,10);
    
    	plot(x,y)    	
        
   %title('Salt transport # spring tide')
       
   %title('Salt transport # neap tide')
    
	ylabel('Transport (Kg m^{-1} s^{-1})','fontsize',14)
    xlabel('Salt transport components','fontsize',14)
   
	hold on
   
%    print -dbitmap G:\Academica42_maio_14\cap05\prgestuario\trans_si.bmp

	%print -dbitmap e:\proces_A_JAN-08\quadratura\trans_qu.bmp


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%  C�lculo dos par�metros circula��o pc  %
    %  e estratifica��o pe                %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	deltas=sss(11)-sss(1);
	szero=sa;

	pe=deltas/szero;

	pc=abs(vss(1)/va);

	
	% Par�metro de estratifica��o %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	pe

	
	% Par�metro de circula��o %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	pc

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diagrama Circula��o x Estratifica��o    %
	% comando para chamar a rotina diagrama.m %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	diagrama_e_c(pc,pe)	
   
    diary off
	
