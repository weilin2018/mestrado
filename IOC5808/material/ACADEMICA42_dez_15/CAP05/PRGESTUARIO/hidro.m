%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA hidro.m                                    %
% PROCESSAMENTO E INTERPOLA��O DE DADOS HIDROGR�FICOS %
% ALESSANDRO LUVIZON B�RGAMO - JULHO DE 1998/2001     %
% LUIZ BRUNER DE MIRANDA                              %
% MAR�O DE 2001/2012                                  % 
% TEORIA EM Miranda et al. (2012)                     %
% M�dias espacial e temporal p. 166 e seguintes       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PROCESSAMENTO E INTERPOLA��O DE DADOS HIDROGR�FICOS AMOSTRADOS EM
% PERFIS VERTICAIS ARMAZENADOS NO ARQUIVO "hidro_quad.dat"
% TRATA-SE DE UM PROGRAMA AUXILIAR PARA O PROCESSAMENTO
% DO PROGRAMA "estuario.m", PARA MAIS DETALHES, CONSULTAR
% B�rgamo, Miranda & Corr�a et al. (2002) - 
% Relat.t�c. inst. oceanogr. S�o Paulo, no. 49 pp.1-16.                                    
    
% NOTA��O E LEGENDA:  

%	p  : profundidade|
%	t  : temperatura | dados de
%	s  : salinidade  | entrada	  
%	st : densidade   |
%       hor: tempo (horas)
%	z  : profundidade adimensional
%	Z  :      "            "       em intervalos decimais
%	T  : temperatura interpolada em profundidade adimensional
%	S  : salinidade        "      "      "            "  
%	ST : densidade         "      "      "            "  
%	M  : matriz contendo Z,T,S,ST para cada esta��o
%	MT : matriz final contendo 14 esta��es (ciclo de mar�) de cada
%		dia de coleta de dados
%	curva em preto: dado real
%	curva em vermelho: dado interpolado

  
	clear
    
    % Carrega o arquivo de dados experimentais "hidro_quad.dat" 
	
	load g:\Academica42_maio_14\cap05\prgestuario\hidro_quad.dat; 

	diary g:\Academica42_maio_14\cap05\prgestuario\stp_siz.dat; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Leitura dos dados hidrogr�ficos          %
    % separando da matriz de dados             %
    % as vari�veis: profundidade, temperatura, % 
    % salinidade e densidade                   %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   

   prof=hidro_quad(:,1);     
   temp=hidro_quad(:,2);  
   salt=hidro_quad(:,3);  
   dens=hidro_quad(:,4);

	
	
	 % Encontra a posi��o onde prof.=0
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     x=find(prof==0);  
	
     tam=length(x);      
   
     for i=1:tam        
     if i~=tam       
	I=num2str(i);
	
	 % Separa os dados de p, t, s, st de todas as esta��es           
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
	pi=['p',I];
	profundidade=[pi,'=prof(x(',I,'):x(',I,'+1)-1);'];
	eval(profundidade);
	        
	ti=['t',I];
	temperatura=[ti,'=temp(x(',I,'):x(',I,'+1)-1);'];
	eval(temperatura);  
        	
	si=['s',I];
	salinidade=[si,'=salt(x(',I,'):x(',I,'+1)-1);'];
	eval(salinidade);
	
	sti=['st',I];
	densidade=[sti,'=dens(x(',I,'):x(',I,'+1)-1);'];
	eval(densidade);

    
	% Calcula a diferen�a entre densidade de superficie e fundo (Dst)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	dsti=['dst',I];
	delta=[dsti,'=dens(x(',I,'+1)-1)-dens(x(',I,'));'];
	eval(delta);

	
	% Transforma p em prof. adim.(Z)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	zi=['z',I];
	prof_adim=[zi,'=(prof(x(',I,'):x(',I,'+1)-1))/(prof(x(',I,'+1)-1));'];
	eval(prof_adim);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Interpola os componentes u e v em d�cimos de prof. adim. (Z) % 
    % pelo m�todo%do "cubic spline"                                %
    % Mais detalhes em Miranda et al. (2002, 2011) p. 157          %                                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	Z=0:.1:1;
	
	Ti=['T',I];
	temp_inter=[Ti,'=spline(prof(x(',I,'):x(',I,'+1)-1)/prof(x(',I,'+1)-1),temp(x(',I,'):x(',I,'+1)-1),Z);'];
	eval(temp_inter);
		
	Si=['S',I];
	salt_inter=[Si,'=spline(prof(x(',I,'):x(',I,'+1)-1)/prof(x(',I,'+1)-1),salt(x(',I,'):x(',I,'+1)-1),Z);'];
	eval(salt_inter);

	STi=['ST',I];
	dens_inter=[STi,'=spline(prof(x(',I,'):x(',I,'+1)-1)/prof(x(',I,'+1)-1),dens(x(',I,'):x(',I,'+1)-1),Z);'];
	eval(dens_inter);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compara graficamente os dados de entrada com os interpolados %
    % permitindo o controle de qualidade das interpola��es         %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	B=['red'];
	K=['black'];
	
	nome=['temperatura'];
	graph=['plot(t',I,',-z',I,',K);hold;plot(T',I,',-Z,B);title(nome);hold;pause;'];
	eval(graph);
	
	nome=['salinidade'];
	graph=['plot(s',I,',-z',I,',K);hold;plot(S',I,',-Z,B);title(nome);hold;pause;'];
	eval(graph);
	
	nome=['densidade'];
	graph=['plot(st',I,',-z',I,',K);hold;plot(ST',I,',-Z,B);title(nome);hold;pause;'];
	eval(graph);	

	% Organiza cada esta��o com respectivas prof.adim., sal., temp.,dens. interpoladas
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	Mi=['M',I];
	mat_parc=[Mi,'=[ -Z ; T',I,' ; S',I,' ; ST',I,' ];'];
	eval(mat_parc);
	
      end
	
     end
     
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   % Organiza todas as esta��es para um ciclo de mar� (13 perfis hor�rios)%
   % separando em dados amostrados em intervalos hor�rios (primeira op��o)%
   % ou a cada meia hora (26 perfis), segunda op��o escolhida             %                
   % neste exemplo                                                        % 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %hor=[ones(1,11)*0, ones(1,11)*1, ones(1,11)*2, ones(1,11)*3, ones(1,11)*4, ones(1,11)*5, ...
    %          ones(1,11)*6, ones(1,11)*7, ones(1,11)*8, ones(1,11)*9, ones(1,11)*10, ones(1,11)*11, ...
    %         ones(1,11)*12, ones(1,11)*13, ones(1,11)*14, ones(1,11)*15, ones(1,11)*16, ones(1,11)*17,...
    %         ones(1,11)*18, ones(1,11)*19, ones(1,11)*20, ones(1,11)*21, ones(1,11)*22, ones(1,11)*23, ...
    %          ones(1,11)*24, ones(1,11)*25];
        
        
   hor=[ones(1,11)*0, ones(1,11)*.5, ones(1,11)*1, ones(1,11)*1.5, ones(1,11)*2, ones(1,11)*2.5, ...
              ones(1,11)*3, ones(1,11)*3.5, ones(1,11)*4, ones(1,11)*4.5, ones(1,11)*5, ones(1,11)*5.5, ...
              ones(1,11)*6, ones(1,11)*6.5, ones(1,11)*7, ones(1,11)*7.5, ones(1,11)*8, ones(1,11)*8.5,...
              ones(1,11)*9, ones(1,11)*9.5, ones(1,11)*10, ones(1,11)*10.5, ones(1,11)*11, ones(1,11)*11.5, ...
              ones(1,11)*12, ones(1,11)*12.5];
     

    % Vari�vel (hor) estabelecendo a hora inicial do experimento 
    % que deve ser alterado 
              
	hor=hor+5.5;
    
    for i=1        
 
%   	mat_tt=['mt=[M',num2str(i),' M',num2str(i+1),' M',num2str(i+2),' M',num2str(i+3)...
%,' M',num2str(i+4),' M',num2str(i+5),' M',num2str(i+6),' M',num2str(i+7),' M',num2str(i+8)...
%,' M',num2str(i+9),' M',num2str(i+10),' M',num2str(i+11),' M',num2str(i+12),' M',num2str(i+13),' ];'];
%   	eval(mat_tt);

   	mat_tt=['mt=[M',num2str(i),' M',num2str(i+1),' M',num2str(i+2),' M',num2str(i+3)...
,' M',num2str(i+4),' M',num2str(i+5),' M',num2str(i+6),' M',num2str(i+7),' M',num2str(i+8)...
,' M',num2str(i+9),' M',num2str(i+10),' M',num2str(i+11),' M',num2str(i+12),' M',num2str(i+13)...
,' M',num2str(i+14),' M',num2str(i+15),' M',num2str(i+16),' M',num2str(i+17),' M',num2str(i+18)...
,' M',num2str(i+19),' M',num2str(i+20),' M',num2str(i+21),' M',num2str(i+22),' M',num2str(i+23)...
,' M',num2str(i+24),' M',num2str(i+25),'];'];
   	eval(mat_tt);

     end

	% Organiza uma matriz para delta_ro (diferen�a entre sup. e fundo de sigma-t) %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    for i=1        
 
   	DST=['Dst=[dst',num2str(i),' dst',num2str(i+1),' dst',num2str(i+2),' dst',num2str(i+3)...
,' dst',num2str(i+4),' dst',num2str(i+5),' dst',num2str(i+6),' dst',num2str(i+7),' dst',num2str(i+8)...
,' dst',num2str(i+9),' dst',num2str(i+10),' dst',num2str(i+11),' dst',num2str(i+12),' dst',num2str(i+13)...
,' dst',num2str(i+14),' dst',num2str(i+15),' dst',num2str(i+16),' dst',num2str(i+17),' dst',num2str(i+18)...
,' dst',num2str(i+19),' dst',num2str(i+20),' dst',num2str(i+21),' dst',num2str(i+22),' dst',num2str(i+23)...
,' dst',num2str(i+24),' dst',num2str(i+25),'];'];
   	eval(DST);
  									     
     end



   	MT=[hor;mt];						
	matriz=MT';
 
	
    hor=matriz(:,1);
	 Z=matriz(:,2);
	 T=matriz(:,3);
	 S=matriz(:,4);
	ST=matriz(:,5);

	pstdi22_08=[hor,Z,T,S,ST];
	temp_22_08=[hor,Z,T];                
	salt_22_08=[hor,Z,S];                 
   	dens_22_08=[hor,Z,ST];
   	dst22_08=Dst;
	        

	%  Organiza uma matriz com as prof. adim. n�o interpoladas  %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %	for i=1
        
 
   %	profundidade=['Z=[z',num2str(i),' ;z',num2str(i+1),' ;z',num2str(i+2),' ;z',num2str(i+3)...
%,' ;z',num2str(i+4),' ;z',num2str(i+5),' ;z',num2str(i+6),' ;z',num2str(i+7),' ;z',num2str(i+8)...
%,' ;z',num2str(i+9),' ;z',num2str(i+10),' ;z',num2str(i+11),' ;z',num2str(i+12),';z',num2str(i+13)...
%,' z',num2str(i+14),' z',num2str(i+15),' z',num2str(i+16),' z',num2str(i+17),' z',num2str(i+18)...
%,' z',num2str(i+19),' z',num2str(i+20),' z',num2str(i+21),' z',num2str(i+22),' z',num2str(i+23)...
%,' z',num2str(i+24),' z',num2str(i+25),'];'];
%	eval(profundidade);
  									
%	end


	% Prof. local instant�nea  
    % Observe que s�o em n�mero de 26
    % pois o intervalo de amostragem �
    % � 30 minutos
	
	 
	k=0;     
   
   	for i=1:26        
     	
	k=k+1;
	zl_22_08(k)=prof(x(i+1)-1);
	
    end
   
    figure 
    
    % Intervalo de tempo entre o in�cio do experimento 5.5 h
    % e o t�rmino do experimento (18.0 h) em intervalos e 0.5 h
    
    tempo=5.5:.5:18.;
    
    % Plota o gr�fico da varia��o da espessura da camada de �gua
    
     plot(tempo,zl_22_08,'linewidth',2)
     xlabel('Tempo em horas','fontsize',12)
     ylabel('Espessura da coluna de �gua','fontsize',12) 
   
        
   	%  Sa�da do diary  %
	%%%%%%%%%%%%%%%%%%%%


	pstdi22_08
	zl_22_08
	dst22_08

	ts_22_08=pstdi22_08;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%  Edi��o final dos arquivos e grava��o  %
    %  Quando necess�rio, basta desbloquear  %
    %  a grava��o desses arquivos            %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % Salva e grava os resultados de entrada para o programa gr�fico
    % "Surfer": temperatura, salinidade e densidade
       
    %save e:\Academica42\cap05\temp_22_08.dat temp_22_08 -ascii      
    %save e:\Academica42\cap05\salt_22_08.dat salt_22_08 -ascii   
    %save e:\Academica42\cap05\dens_22_08.dat dens_22_08 -ascii   

    % Salva e grava os resultados de entrada para o programa "estuario.m" 
    % Dados termohalinos (ts), profundidade (zl) e dst (densidade)
    % salinidade e temperatura
		
		
	    %save
	    save g:\Academica42\cap05\prgestuario\ts_22_08.dat ts_22_08 -ascii    
        save g:\Academica42\cap05\prgestuario\zl_22_08.dat zl_22_08 -ascii  
        %save f:\estuario\programas\cap5\dst22_08.dat dst22_08 -ascii  

	
	diary off
