  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % PLOTA OS PAR�METROS NUMERO DE ESTRATIFICA��O(st)x ESTRATIFICA��0 (pe) %
  % NO DIAGRAMA CORRESPONDENTE:  IPPEN & HARLEMAN (1961) e PRANDLE (1985),%
  % ESTABELECENDO UM CRIT�RIO PARA CLASSIFICA��O DE ESTU�RIOS             %
  % TEORIA em Miranda et al. (2002, 2012 2nd ed.)                         %
  % Equa��es 3.4 e 3.5 p. 107 e 108                                       %
  %%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	function diagr_Prandle(st,pe)
    
    clear
    
    close all

  %COEFICIENTES, PAR�METROS (geometria, descarga fluvial)

   k=2.5e-003; % Par�metro adimensional de mistura
   uo=0.8; % Amplitude da velocidade barotr�pica (m/s)
   uf=0.1; % Velocidade gerada pela descarga fluvial
   L=1.0e+006; % Dimens�o longitudinal do estu�rio (m)
   ho=10.0; % Profundidade m�dia do estu�rio
   rob=1030.0; % Densidade na boca do estu�rio (kg/m3)
   roc=1000.0; % Densidade na cabeceira do estu�rio (kg/m3)  
   g=9.8; % acelera��o da gravidade (m/s2)
   ro=1000; % Densidade de refer�ncia (kg/m3) 
	
   % C�lculo do N�mero de estratifica��o "st" 
   % iniciando pela gravidade reduzida "gr"
   
   gr=g*((rob-roc)/ro)
   
   st=0.85*((k*uo^3*L)/(gr*ho^2)*uf)
   
   % C�lculo do par�metro estratifica��o (pe)
   
   %pe=0.02
   
   pe=4.0*st^(-0.55)
   
   % Diagrama com eixos em log x log
   
   loglog(st,pe,'*')
   
   % Tra�a as linhas que delitam a transi��o entre
   % os estu�rios estratificados e os bem misturados
   % ou seja, no intervalo 100<st<400 temos os estu�rio
   % parcialmente misturados
   
    x=[100,100]
    y=[200,0.01]
    
    line(x,y)
        
    xx=[400,400]
    yy=[400,0.01]
    
    line(xx,yy)
        
    ylabel('Par�metro de estratifica��o - pe','fontsize',14)
    xlabel('N�mero de estratifica��o - st ','fontsize',14)
   
	
	axis([1 10^5 10^(-2) 10])
    
    hold on
    
    % Com o comando "get" "s�o indicados no diagrama os 
    % tipos de estu�rios " Estratificado" e "Bem misturado
    
        
    gtext('Estratificado','fontsize',14)	
   
    gtext('Bem misturado','fontsize',14)
      
    print -dbitmap e:\Academica42\cap03\Figura-Prandle

    
