%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA Richcamada.m   "Número de Richardson por camada"       %
% O PROGRMA CALCULA A VARIAÇAO NO TEMPO DO NÚMERO DE RICHARDSON   %
% TEM COMO PROGRAMAS AUXILIARES OS PROGRAMAS "hidro.m             % 
% E "decomp.m", QUE SÃO APRESENTADOS NOS PROGRAMAS DO CAP. 05     %
% QUE GERAM  OS ARQUIVOS  DE DADOS HIDROGRÁFICOS "ts.dat"         %
% E "vel.dat" e zl.dat E DEVEM SER CARREGADOS COM                 % 
% COM O COMANDO "load" NO INÍCIO DO PROCESSAMENTO. O PROGARAMA    %
% PARA OS DADOS DESTE PROCESSAMENTO "nest=26" e "npro=11".        % 
% CALCULA VALORES MÉDIOS NO TEMPO E NO ESPAÇO                     %
% COM INTEGRAÇÕES NUMÉRICAS COM A TEORIA                          %
% APRESENTADA EM Miranda et al. (2002, 2012(2nd ed.))             %
% Equação 2.35 - p. 85                                            %
% PARA MAIS DETALHES, CONSULTAR:                                  %
% Bérgamo, Miranda e Corrêa (2002), Relat.téc.                    %
% inst. oceanogr. São Paulo,  no. 49 pp.1-16.                     %
% Quando  RiL >20, RiL<20 ou RiL<2 a coluna de água é estável,   %
% a turbulência no fundo pode ocasionar mistura vertical e na     %
% última condição a mistura turbulenta torna a coluna de água     %
% instável de acordo com Dyer & New (1986).                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

		
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	% Adaptação do programa dinest.m e estuario.m          % 
	% Programado pelos Alunos da IOF-827 do ano de 1996    %
    % Os argumentos da "function" são os arquivos de dados:%
    % vel.dat, ts.dat, zl.dat,nest=26 e npro=11            %
    % Obs. h=zl.dat                                        %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	function  Y=richcamada(vel,ts,h,nest,npro)
	%global va sa vs ss sss vss pa vdoce
	
   %diary E:\ACADEMICA42_maio_14\CAP02\Richcamada\Richardson.txt

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%   vel   = matriz de velocidade longitudinal          %
	%   ts   = matriz de temperatura e salinidade          %
	%   h   = matriz-linha de prof. local instantânea      %	
	%   nest = número de estações                          %
	%   npro = número de profundidades em cada estação     %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%   Dados de Temperatura e Salinidade     %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	prof=-ts(:,2); % Profundidade adimensional sinal (-) devido a orientação de OZ
	temp=ts(:,3); % Temperatura em graus Celsius
	salt=ts(:,4); % Salinidade em EPS
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
	% ATENÇÃO: O PROGRAMA USA A SUBROTINA "SW_DENS"    %
    % DO MORGAN (1994)- CSIRO MARINE LABORATORIES,     %
    % AUSTRÁLIA, 222, 28 P. PARA CALCULAR A DENSIDADE  %
    % DA ÁGUA ESTUARINA À PRESSÃO HIDROSTÁTICA ZERO    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	dens=sw_dens(salt,temp,0);         % Densidade em (kg/m^3)
	
	sigma_t=(dens-1000);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Cálculo da Prof. Adm. Z=z/h %
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
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%        Velocidade                   %
    %  Usa os dados de velocidade         %
    %  processados pelo decomp.m          %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   u=vel(:,2); % velocidade longitudinal em (m/s)
   prov=-vel(:,1);
   %prov=vel(:,1);

   
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Cálculo da velocidade nas Prof. Adm. Z         %
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
	% Interpolação Pelo Método "Cubic Spline" todas  % 
    % variáveis S, T, densidade e velocidade         % 
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
	%  hidro: Estação  Prof.Adm.  Temp.  Salin. Dens.  Sigma_T   Veloc. %
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


	% Cálculo da diferença entre densidade de superficie e fundo (Delta-ro)
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
	% Cálculo das médias por integração numérica        %
    % ao longo da coluna de água                        %
    % Miranda et al. (2002, 2011 2nd) eq. 5.19 - p. 169 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%    Média no espaço     % 
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
	
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Número de Richardson por Camada - Medida da estratificação vertical da %
% Coluna de água - Mais detalhes em Miranda et al. (2002, 2011 2nd ed.)  %
% eq. 2.35 - p. 85                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	RiL=(9.8.*(abs(h))./vme.^2).*(Dro./rome)	%valores adimensionais
	

	  tempo=5.5:.5:18.; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estratificação da coluna d'água - Número de Richardson por Camada,    %
% Mais detalhes em Miranda et al. (2002, 2011 2nd ed.) eq. 2.35 - p. 85 %
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
   
   title('Número de Richardson, Ri_{L}# sizígia')

      
   ylabel('Ri_{L}','fontsize',14)
   xlabel('Tempo (h)','fontsize',14)
   
	hold off
   

	%print -dbitmap c:\testecuri\riL_si.bmp
   
   print -dbitmap e:\ACADEMICA42\CAP02\richcamada\riL.bmp

   
    diary off
	
