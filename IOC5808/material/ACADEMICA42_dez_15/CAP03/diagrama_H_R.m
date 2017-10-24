  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % PLOTA OS PARÂMETROS ESTRATIFICAÇÃO (pe)x ESTRATIFICAÇÃ0 (pc) NO %
  % DIAGRAMA ESTRATIFICAÇÃ0-CIRCULAÇÃO DE HANSEN & RATTRAY (1966)   %
  % PROGRAMA PREPARADO POR ALESSANDRO LUVISON BÉRGAMO (1998)        %
  % TEORIA EM Miranda et al. (2002, 2012 2nd ed.)                   %
  % Equações 3.6 p. 111 ou eq. 11.123 p. 407                        %
  % PARA MAIS DETALHES, CONSULTAR                                   %
  % Bérgamo et al. (2002), Relat.téc. inst. oceanogr. São Paulo,    % 
  % no. 49 pp.1-16.                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	function diagrama(pc,pe)
    
   % Entrada dos parâmetros "pc" e "pe" calculados com o programa "estuario.m 
   % Esses valores são apresentados no final do processamento
   
      pe=0.015;
      pc=5.33;
   
     
    % k é o parâmetro "ni" associado ao par ordenado (parcrc,parstr)
    %
	% n é o parâmetro "ni"- estabelece a proporção relativa dos processos
    % difusivo e advectivo para o transporte de sal estuário acima
	%
	% C é o eixo das abscissas (parâmetro circulação)
	%
	% E é o eixo das ordenadas (parâmetro estratificação)
    % 
    % O programa calcula a equação do segundo grau incompleta
    %
    % na incógnita (ni) 
       
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
       
  
	ylabel('Parâmetro Estratificação-pe','fontsize',14)
    xlabel('Parâmetro Circulação - pc','fontsize',14)
  	title('Diagrama Estratificação-Circulação')

	axis([1 10^3 10^(-2) 2])
    
    % Os eixos podem ser reconfigurados
	% axis([10^(-1) 10^3 10^(-2) 2])
    
    % A incógnita ni é plotada no canto superior direito
    % Utilize o comando "get" para plotar os valores 
    % dessa incógnita correspondentes a 1.0, 0.01. O valor
    % do "get" igual a 1.5 indica o valor mínimo teórico do 
    % parâmetro circulação (pc) que, quando ni=1.0,
    % é independente do parâmetro pe,
    
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