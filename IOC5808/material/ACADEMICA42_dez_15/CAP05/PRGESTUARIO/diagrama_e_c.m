  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % PLOTA OS PARÂMETROS ESTRATIFICAÇÃO (pe)x ESTRATIFICAÇÃ0 (pc) NO %
  % DIAGRAMA ESTRATIFICAÇÃ0-CIRCULAÇÃO DE HANSEN & RATTRAY (1966)   %
  % PROGRAMA PREPARADO POR ALESSANDRO LUVISON BÉRGAMO (1998)        %
  % TEORIA EM Miranda et al. (2002, 2012 2nd ed.)                   %
  % Equações 5.20 p. 170, e 11.123 p. 399                           %
  % PARA MAIS DETALHES, CONSULTAR                                   %
  % Bérgamo, Miranda & Corrêa et al. (2002) -                       % 
  % Relat.téc. inst. oceanogr. São Paulo, no. 49 pp.1-16.           %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	function diagrama_e_c(pc,pe)
    
    % São dados de entrada os parâmetros "pe"e "pc" calculados pelo
    % programa "estuario.m" como exemplificado a seguir
   
    pe=0.036;
    pc=5.2527;
   
     
   % k é o parâmetro "ni" associado ao par ordenado (parcrc,parstr)
   %
	% n é o parâmetro "ni"
	%
	% C é o eixo das abscissas (parâmetro circulação)
	%
	% E é o eixo das ordenadas (parâmetro estratificação)
    
    figure
   
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
    plot(pc,pe,'bx')
   
  
	ylabel('Parâmetro estratificação','fontsize',16)
    xlabel('Parâmetro circulação','fontsize',16)
   
	%title('Stratification-circulation diagram # spring tide')

	%title('Stratification-circulation diagram # neap tide')

    axis([1 10^3 10^(-2) 2])

	%axis([10^(-1) 10^3 10^(-2) 2])

	for i=k        
     	       
	I=num2str(i);

	K=[I];
	texto=['legend(','I',')'];
	eval(texto);

	end

	gtext('1.0')
	gtext('0.01')
	gtext('1.5')

   
	%print -dbitmap c:\testecuri\dgm_si.bmp
   
   print -dbitmap e:\estuario-03-mar\programas-03-mar\cap05\prgestuario\diagrama

   
hold 