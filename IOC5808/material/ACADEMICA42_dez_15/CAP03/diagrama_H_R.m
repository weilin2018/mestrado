  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % PLOTA OS PAR�METROS ESTRATIFICA��O (pe)x ESTRATIFICA��0 (pc) NO %
  % DIAGRAMA ESTRATIFICA��0-CIRCULA��O DE HANSEN & RATTRAY (1966)   %
  % PROGRAMA PREPARADO POR ALESSANDRO LUVISON B�RGAMO (1998)        %
  % TEORIA EM Miranda et al. (2002, 2012 2nd ed.)                   %
  % Equa��es 3.6 p. 111 ou eq. 11.123 p. 407                        %
  % PARA MAIS DETALHES, CONSULTAR                                   %
  % B�rgamo et al. (2002), Relat.t�c. inst. oceanogr. S�o Paulo,    % 
  % no. 49 pp.1-16.                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	function diagrama(pc,pe)
    
   % Entrada dos par�metros "pc" e "pe" calculados com o programa "estuario.m 
   % Esses valores s�o apresentados no final do processamento
   
      pe=0.015;
      pc=5.33;
   
     
    % k � o par�metro "ni" associado ao par ordenado (parcrc,parstr)
    %
	% n � o par�metro "ni"- estabelece a propor��o relativa dos processos
    % difusivo e advectivo para o transporte de sal estu�rio acima
	%
	% C � o eixo das abscissas (par�metro circula��o)
	%
	% E � o eixo das ordenadas (par�metro estratifica��o)
    % 
    % O programa calcula a equa��o do segundo grau incompleta
    %
    % na inc�gnita (ni) 
       
   	C=[1:10^3];
   
   
   	k=-(-pe^(-1)*(210+252*(pc-3/2))+76*(pc-3/2)+152/3*(pc-3/2)^2)/(pe^(-1)*(210+252*(pc-3/2)));
   
   	E=(((210+252*(C-(3/2))).*k.^2)-((210+252*(C-(3/2)))*k))./(-(76*(C-(3/2)).*k+(152/3)*((C-(3/2)).^2).*k));
      
      		loglog(C,E,'b-')
            
                        
            hold on
            
      
      	n=0.01;		

	E=(((210+252*(C-(3/2))).*n.^2)-((210+252*(C-(3/2)))*n))./(-(76*(C-(3/2)).*n+(152/3)*((C-(3/2)).^2).*n));
				
		loglog(C,E,'k:')


	n=0.1;
	
	E=(((210+252*(C-(3/2))).*n.^2)-((210+252*(C-(3/2)))*n))./(-(76*(C-(3/2)).*n+(152/3)*((C-(3/2)).^2).*n));

		loglog(C,E,'k:')


	n=0.5;
	
	E=(((210+252*(C-(3/2))).*n.^2)-((210+252*(C-(3/2)))*n))./(-(76*(C-(3/2)).*n+(152/3)*((C-(3/2)).^2).*n));
				
		loglog(C,E,'k:')


	n=0.9;
	
	E=(((210+252*(C-(3/2))).*n.^2)-((210+252*(C-(3/2)))*n))./(-(76*(C-(3/2)).*n+(152/3)*((C-(3/2)).^2).*n));

		loglog(C,E,'k:')


	n=0.99;
	
	E=(((210+252*(C-(3/2))).*n.^2)-((210+252*(C-(3/2)))*n))./(-(76*(C-(3/2)).*n+(152/3)*((C-(3/2)).^2).*n));

		loglog(C,E,'k:')

	x1(1)=1.5;
	x1(2)=1.5;
	y1(1)=.001;
	y1(2)=100;
   
		plot(x1,y1,'k:')

	x2(1)=1;
	x2(2)=1000;
	y2(1)=.1;
	y2(2)=.1;

	plot(x2,y2,'k-.')		


	plot(pc,pe,'bo')
       
  
	ylabel('Par�metro Estratifica��o-pe','fontsize',14)
    xlabel('Par�metro Circula��o - pc','fontsize',14)
  	title('Diagrama Estratifica��o-Circula��o')

	axis([1 10^3 10^(-2) 2])
    
    % Os eixos podem ser reconfigurados
	% axis([10^(-1) 10^3 10^(-2) 2])
    
    % A inc�gnita ni � plotada no canto superior direito
    % Utilize o comando "get" para plotar os valores 
    % dessa inc�gnita correspondentes a 1.0, 0.01. O valor
    % do "get" igual a 1.5 indica o valor m�nimo te�rico do 
    % par�metro circula��o (pc) que, quando ni=1.0,
    % � independente do par�metro pe,
    
	for i=k        
     	       
	I=num2str(i);
          
    K=[I];
	texto=['legend(','I',')'];
	eval(texto);

	end

	gtext('1.0')
	gtext('0.01')
	gtext('1.5')

   	 
   print -dbitmap e:\Academica42\cap03\diagrama

   
   hold 