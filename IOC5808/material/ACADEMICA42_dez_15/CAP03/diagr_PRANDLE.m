  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % PLOTA OS PARÂMETROS NUMERO DE ESTRATIFICAÇÃO(st)x ESTRATIFICAÇÃ0 (pe) %
  % NO DIAGRAMA CORRESPONDENTE:  IPPEN & HARLEMAN (1961) e PRANDLE (1985),%
  % ESTABELECENDO UM CRITÉRIO PARA CLASSIFICAÇÃO DE ESTUÁRIOS             %
  % TEORIA em Miranda et al. (2002, 2012 2nd ed.)                         %
  % Equações 3.4 e 3.5 p. 107 e 108                                       %
  %%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	function diagr_Prandle(st,pe)
    
    clear
    
    close all

  %COEFICIENTES, PARÂMETROS (geometria, descarga fluvial)

   k=2.5e-003; % Parâmetro adimensional de mistura
   uo=0.8; % Amplitude da velocidade barotrópica (m/s)
   uf=0.1; % Velocidade gerada pela descarga fluvial
   L=1.0e+006; % Dimensão longitudinal do estuário (m)
   ho=10.0; % Profundidade média do estuário
   rob=1030.0; % Densidade na boca do estuário (kg/m3)
   roc=1000.0; % Densidade na cabeceira do estuário (kg/m3)  
   g=9.8; % aceleração da gravidade (m/s2)
   ro=1000; % Densidade de referência (kg/m3) 
	
   % Cálculo do Número de estratificação "st" 
   % iniciando pela gravidade reduzida "gr"
   
   gr=g*((rob-roc)/ro)
   
   st=0.85*((k*uo^3*L)/(gr*ho^2)*uf)
   
   % Cálculo do parâmetro estratificação (pe)
   
   %pe=0.02
   
   pe=4.0*st^(-0.55)
   
   % Diagrama com eixos em log x log
   
   loglog(st,pe,'*')
   
   % Traça as linhas que delitam a transição entre
   % os estuários estratificados e os bem misturados
   % ou seja, no intervalo 100<st<400 temos os estuário
   % parcialmente misturados
   
    x=[100,100]
    y=[200,0.01]
    
    line(x,y)
        
    xx=[400,400]
    yy=[400,0.01]
    
    line(xx,yy)
        
    ylabel('Parâmetro de estratificação - pe','fontsize',14)
    xlabel('Número de estratificação - st ','fontsize',14)
   
	
	axis([1 10^5 10^(-2) 10])
    
    hold on
    
    % Com o comando "get" "são indicados no diagrama os 
    % tipos de estuários " Estratificado" e "Bem misturado
    
        
    gtext('Estratificado','fontsize',14)	
   
    gtext('Bem misturado','fontsize',14)
      
    print -dbitmap e:\Academica42\cap03\Figura-Prandle

    
