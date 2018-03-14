%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA decomp.m                                   %
% DECOMPOSI��O DE VETORES DE VELOCIDADE               %
% ALESSANDRO LUVIZON B�RGAMO - JULHO DE 2001          %
% LUIZ BRUNER DE MIRANDA                              %
% ALESSANDRO LUVIZON B�RGAMO                          %
% MAR�O DE 2001/2012                                  % 
% TEORIA Miranda et al. (2012 - 2nd ed.)              %
% Equa��es 5.1 e 5.2 p. 154 e seguintes               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PROCESSAMENTO E INTERPOLA��O DE DADOS CORRENTOGR�FICOS AMOSTRADOS EM
% PERFIS VERTICAIS ARMAZENADOS NO ARQUIVO "vel_quad.dat"
% TRATA-SE DE UM PROGRAMA AUXILIAR PARA O PROCESSAMENTO
% DO PROGRAMA "estuario.m", PARA MAIS DETALHES, CONSULTAR
% B�rgamo, Miranda & Corr�a et al. (2002)- 
% Relat.t�c. inst. oceanogr. S�o Paulo, no. 49 pp.1-16.  

% NOTA��O E LEGENDA:  
%	p : profundidade  |
%	u : veloc. trans. | dados de entrada
%	v : veloc. long.  |   
%	Z : profundidade adimensional
%	Z :      "            "       em intervalos decimais
%	U : veloc. trans. interpolada em profundidade adimensional
%       V : veloc. long.        "      "      "            "  
%	MU: contendo Z,U para cada esta��o
%	MV:   "       "     Z,V  "     "     "  
%  hora: hora inteira mais pr�xima da medida" 
%  curva preta :dado de entrada
%  curva vermelha:dado interpolado
%  os dados de entrada est�o na matriz com os dados%
%  hora (opcional), profundidade, velocidade e dire��o%
%  neste programa o nome desse arquivo de dados � vel_quad.dat%	  

		
   clear all

   % Leitura dos dados	
   
   
  load g:\Academica42_maio_14\cap05\prgestuario\vel_quad.dat
	
  diary veloc.txt
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Leitura dos dados correntogr�ficos     %
    % separando da matriz de dados           %
    % as vari�veis: hora, profundidade,      % 
    % intensidade da velocidade (vel) e      %
    % dire��o (dir)                          %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     hora=vel_quad(:,4); 
     prof=-vel_quad(:,5);
     vel=vel_quad(:,6);
     dir=vel_quad(:,7);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Decomposi��o do vetor velocidade nas componentes "u" e "v"                %
	% levando em conta o angulo de declina��o magn�tica (D=23 graus W)          %
    % referente a orienta��o do canal estuarino de Caravelas (BA)no ano de 2008 %
    % e a rota��o do eixo x (B) em 23.0 graus no sentido anti-hor�rio.          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
     D=23.0*pi/180;         
     A=pi/2;
     rad= dir*pi/180;
     B=10*pi/180;
     
     % vel_u e  vel_v denotam os componentes u e v de velocidade

     vel_u=vel.*cos(A-(rad-D)+B);
     vel_v=vel.*sin(A-(rad-D)+B);
     
    %Prepara��o da matriz de armazenamento dos dados 

	vel_22=[hora prof vel dir vel_u vel_v];
	
   
	 % Encontra a posi��o onde prof.=0
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     x=find(prof==0);  
	
     tam=length(x);      
   
     for i=1:tam        
     if i~=tam       
	I=num2str(i);
	
	% Separa de todas as esta��es os dados de prof. (p), u e v             
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
	pi=['p',I];
	profundidade=[pi,'=prof(x(',I,'):x(',I,'+1)-1);'];
	eval(profundidade);
	        
	ui=['u',I];
	vel_trans=[ui,'=vel_u(x(',I,'):x(',I,'+1)-1);'];
	eval(vel_trans);  
       	
        
	vi=['v',I];
	vel_long=[vi,'=vel_v(x(',I,'):x(',I,'+1)-1);'];
	eval(vel_long);
	         
         
    % Transforma a profunidade (p) em prof. adim.(Z)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	zi=['z',I];
	prof_adim=[zi,'=(prof(x(',I,'):x(',I,'+1)-1))/(prof(x(',I,'+1)-1));'];
	eval(prof_adim);
          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
	% Interpola os componentes u e v em d�cimos da prof. adim. (Z) % 
    % pelo m�todo do "cubic spline"                                %
    % Mais detalhes em Miranda et al. (2002, 2011) p. 157          %                                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Z=0:.1:1;

    Ui=['U',I];
	u_interp=[Ui,'=interp1(prof(x(',I,'):x(',I,'+1)-1)/prof(x(',I,'+1)-1),vel_u(x(',I,'):x(',I,'+1)-1),Z);'];
	eval(u_interp);
	  
	Vi=['V',I];
	v_interp=[Vi,'=interp1(prof(x(',I,'):x(',I,'+1)-1)/prof(x(',I,'+1)-1),vel_v(x(',I,'):x(',I,'+1)-1),Z);'];
	eval(v_interp);
        
     Ui=['U',I];
     u_interp=[Ui,'=spline(prof(x(',I,'):x(',I,'+1)-1)/prof(x(',I,'+1)-1),vel_u(x(',I,'):x(',I,'+1)-1),Z);'];
     eval(u_interp);
         
      Vi=['V',I];
      v_interp=[Vi,'=spline(prof(x(',I,'):x(',I,'+1)-1)/prof(x(',I,'+1)-1),vel_v(x(',I,'):x(',I,'+1)-1),Z);'];
      eval(v_interp);
         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
	% Compara graficamente os dados de entrada com os interpolados %
    % Permitindo o controle de qualidade dos resultado obtidos     %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	R=['red'];
	B=['black'];
    
    xxx=[0,0];
    yyy=[0,-1];
   
	
	nome=['componente u (m/s)'];
   graph=['plot(u',I,',-z',I,',B);hold;plot(U',I,',-Z,R);title(nome);line(xxx,yyy);hold;pause(1);pause'];
   %graph=['plot(u',I,',-z',I,',B);hold;plot(U',I,',-Z,R);title(nome);hold;pause(1);pause'];
     
	eval(graph);
	
	nome=['componente v (m/s)'];
   graph=['plot(v',I,',-z',I,',B);hold;plot(V',I,',-Z,R);title(nome);line(xxx,yyy);hold;pause(1);pause'];
   %graph=['plot(v',I,',-z',I,',B);hold;plot(V',I,',-Z,R);title(nome);hold;pause(1);pause'];
     
	eval(graph);
        
        
        
	% Organiza cada esta��o com respectivas prof. adim., e veloc. long. u e veloc. trans.v %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%

	MUi=['MU',I];
	mat_u_parc=[MUi,'=[-Z;U',I,'];'];
	eval(mat_u_parc);

	MVi=['MV',I];
	mat_v_parc=[MVi,'=[-Z;V',I,'];'];
	eval(mat_v_parc);

	MUVi=['MUV',I];
	mat_tt=[MUVi,'=[-Z;U',I,';V',I,'];'];
	eval(mat_tt);

	
      end
	
    end


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Organiza todas as esta��es para um ciclo de mar� ( 13 horas )        %
	% ou intervalos hor�rios separando em dados amostrados                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
	hor=[ones(1,11)*0, ones(1,11)*.5, ones(1,11)*1, ones(1,11)*1.5, ones(1,11)*2, ones(1,11)*2.5, ...
              ones(1,11)*3, ones(1,11)*3.5, ones(1,11)*4, ones(1,11)*4.5, ones(1,11)*5, ones(1,11)*5.5, ...
              ones(1,11)*6, ones(1,11)*6.5, ones(1,11)*7, ones(1,11)*7.5, ones(1,11)*8, ones(1,11)*8.5,...
              ones(1,11)*9, ones(1,11)*9.5, ones(1,11)*10, ones(1,11)*10.5, ones(1,11)*11, ones(1,11)*11.5, ...
              ones(1,11)*12, ones(1,11)*12.5];

   % Vari�vel (hor) estabelecendo a hora inicial do experimento 
   % que deve ser alterado 
              
	hor=hor+5.5;

    for i=1
        
 
   	mat_u_tt=['mut=[MU',num2str(i),' MU',num2str(i+1),' MU',num2str(i+2),' MU',num2str(i+3)...
,' MU',num2str(i+4),' MU',num2str(i+5),' MU',num2str(i+6),' MU',num2str(i+7),' MU',num2str(i+8)...
,' MU',num2str(i+9),' MU',num2str(i+10),' MU',num2str(i+11),' MU',num2str(i+12),' MU',num2str(i+13)...
,' MU',num2str(i+14),' MU',num2str(i+15),' MU',num2str(i+16),' MU',num2str(i+17),' MU',num2str(i+18)...
,' MU',num2str(i+19),' MU',num2str(i+20),' MU',num2str(i+21),' MU',num2str(i+22),' MU',num2str(i+23)...
,' MU',num2str(i+24),' MU',num2str(i+25),'];'];
   	eval(mat_u_tt);

	 mat_v_tt=['mvt=[MV',num2str(i),' MV',num2str(i+1),' MV',num2str(i+2),' MV',num2str(i+3)...
,' MV',num2str(i+4),' MV',num2str(i+5),' MV',num2str(i+6),' MV',num2str(i+7),' MV',num2str(i+8)...
,' MV',num2str(i+9),' MV',num2str(i+10),' MV',num2str(i+11),' MV',num2str(i+12),' MV',num2str(i+13)...
,' MV',num2str(i+14),' MV',num2str(i+15),' MV',num2str(i+16),' MV',num2str(i+17),' MV',num2str(i+18)...
,' MV',num2str(i+19),' MV',num2str(i+20),' MV',num2str(i+21),' MV',num2str(i+22),' MV',num2str(i+23)...
,' MV',num2str(i+24),' MV',num2str(i+25),'];'];
   	eval(mat_v_tt);
                      
           
	 mat_uv_tt=['muvt=[MUV',num2str(i),' MUV',num2str(i+1),' MUV',num2str(i+2),' MUV',num2str(i+3)...
,' MUV',num2str(i+4),' MUV',num2str(i+5),' MUV',num2str(i+6),' MUV',num2str(i+7),' MUV',num2str(i+8)...
,' MUV',num2str(i+9),' MUV',num2str(i+10),' MUV',num2str(i+11),' MUV',num2str(i+12),' MUV',num2str(i+13)...
,' MUV',num2str(i+14),' MUV',num2str(i+15),' MUV',num2str(i+16),' MUV',num2str(i+17),' MUV',num2str(i+18)...
,' MUV',num2str(i+19),' MUV',num2str(i+20),' MUV',num2str(i+21),' MUV',num2str(i+22),' MUV',num2str(i+23)...
,' MUV',num2str(i+24),' MUV',num2str(i+25),'];'];
   	eval(mat_uv_tt);

   	    
    end

   	MUT=[hor;mut];
	u_22_08=MUT';
        
    vl_22_08=mut';
        
	MVT=[hor;mvt];
	v_22_08=MVT';

	
	MUVT=[hor;muvt];
	vint_22_08=MUVT';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%  Organiza uma matriz com as prof. adim. n�o interpoladas  %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    	for i=1
        
 
   	profundidade=['Z=[z',num2str(i),' ;z',num2str(i+1),' ;z',num2str(i+2),' ;z',num2str(i+3)...
,' ;z',num2str(i+4),' ;z',num2str(i+5),' ;z',num2str(i+6),' ;z',num2str(i+7),' ;z',num2str(i+8)...
,' ;z',num2str(i+9),' ;z',num2str(i+10),' ;z',num2str(i+11),' ;z',num2str(i+12),' ;z',num2str(i+13)...
,' ;z',num2str(i+14),' ;z',num2str(i+15),' ;z',num2str(i+16),' ;z',num2str(i+17),' ;z',num2str(i+18)...
,' ;z',num2str(i+19),' ;z',num2str(i+20),' ;z',num2str(i+21),' ;z',num2str(i+22),' ;z',num2str(i+23)...
,' ;z',num2str(i+24),' ;z',num2str(i+25),'];'];
   	eval(profundidade);
  
     	end

	y=[0];
	Z=[Z;y];	

	vad_22_08=[ Z vel dir vel_u vel_v];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Varia��o da espessura da coluna de �gua no ciclo de mar�%
    % Em intervalos de tempo de medida iguais a 30 minutos    %
    % S�o 26 medidas no intrervalo de 5.5 h  a 18.0 horas     % 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	k=0;     
   
   	for i=1:26        
   	k=k+1;
	zl_22_08(k)=-prof(x(i+1)-1);
	
        end

	%  Sa�da do diary  %
	%%%%%%%%%%%%%%%%%%%%



    vel_22
	vad_22_08
	vint_22_08
   u_22_08
	v_22_08
   vl_22_08
	zl_22_08

    
    % Intervalo de tempo entre o in�cio do experimento 5.5 h
    % e o t�rmino do experimento (18.0 h) em intervalos e 0.5 h
    
    tempo=5.5:.5:18.;
    
    figure
    
    % Plota o gr�fico da varia��o da espessura da camada de �gua
    
        tempo=5.5:.5:18.0;
        plot(tempo,zl_22_08,'linewidth',2)
        xlabel('Tempo em horas','fontsize',12)
        ylabel('Espessura da coluna de �gua','fontsize',12) 
        
        
        print -dbitmap e:\Academica42\cap05\prgestuario\h_t_quad.bmp
     
        
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%  Edi��o final dos arquivos e grava��o  %
    %  Quando necess�rio, basta desbloquear  %
    %  a grava��o desses arquivos            %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 	
     % Salva e grava os resultados de entrada para o programa gr�fico
     % "Surfer": componentes u e v de velocidade e profundidade

    %save f:\estuario\programas\cap5\u_22_08.dat u_22_08 -ascii		    
	%save f:\estuario\programas\cap5\v_22_08.dat v_22_08 -ascii		   
    %save f:\estuario\programas\cap5\zl_22_08.dat zl_22_08 -ascii	
    
    % Salva e grava os resultados de entrada para o programa "estuario.m" 
    % resultados do componente longitudinal de velocidade u %
    % e a profundidade (zl)
    		

	%save e:\Academica42\cap05\prgestuario\zl_22_08.dat zl_22_08 -ascii    % dados de entrada para
	 save e:\Academica42\cap05\prgestuario\vl_22_08.dat vl_22_08 -ascii      %	estuario.m


    diary off
