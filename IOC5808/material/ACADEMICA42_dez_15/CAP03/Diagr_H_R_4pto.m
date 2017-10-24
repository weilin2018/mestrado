 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 % PLOTA UM CONJUNTO DE QUATRO PARES DE PARÂMETROS ESTRATIFICAÇÃO  (pe)  %
 % x CIRCULAÇÃO (pc) NO DIAGRAMA ESTRATIFICAÇÃ0-CIRCULAÇÃO               %
 % DE HANSEN & RATTRAY (1966)                                            %
 % PROGRAMA PREPARADO POR FERNANDO PINHEIRO ANDUTTA (2010)               %
 % TEORIA EM Miranda et al. (2002, 2012 - 2nd ed.)                         %
 % Equações 3.6 p. 111 ou eq. 11.123 p. 407                              %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all

% Entrada dos quatro pares de parâmetros  

pe1 = 0.073;
pc1 = 11.42;

pe2 = 0.066;
pc2 = 12.49;

pe3 = 0.063;
pc3 = 7.49;

pe4 = 0.059;
pc4 = 7.88;

% Especificação das plotagens e idioma das figuras
% Português (1) ou Inglês (0), e formato dos caracteres

plot_PcPe = 1;
language = 1;
esp_title = 12;
esp_xlabel = 16;
esp_ylabel = 16;
language = 1;
salva = 1;

if plot_PcPe==1;


% Os comandos abaixo traçam as isolinhas do parâmetro "ni"    
figure     
C=[1:10^3];

n=0.01;		
E=(((210+252*(C-(3/2))).*n.^2)-((210+252*(C-(3/2)))*n))./(-(76*(C-(3/2)).*n+(152/3)*((C-(3/2)).^2).*n));				
loglog(C,E,'k:','linewidth',2)

hold on

% A isolina n=0.1 abaixo pode ser ativada

%n=0.1;	
%E=(((210+252*(C-(3/2))).*n.^2)-((210+252*(C-(3/2)))*n))./(-(76*(C-(3/2)).*n+(152/3)*((C-(3/2)).^2).*n));
%loglog(C,E,'k:')

hold on

n=0.4;	
E=(((210+252*(C-(3/2))).*n.^2)-((210+252*(C-(3/2)))*n))./(-(76*(C-(3/2)).*n+(152/3)*((C-(3/2)).^2).*n));				
loglog(C,E,'k:','linewidth',2)

hold on

n=0.9;	
E=(((210+252*(C-(3/2))).*n.^2)-((210+252*(C-(3/2)))*n))./(-(76*(C-(3/2)).*n+(152/3)*((C-(3/2)).^2).*n));
loglog(C,E,'k:','linewidth',2)

hold on

n=0.99;	
E=(((210+252*(C-(3/2))).*n.^2)-((210+252*(C-(3/2)))*n))./(-(76*(C-(3/2)).*n+(152/3)*((C-(3/2)).^2).*n));
%loglog(C,E,'k:')
loglog(C,E,'k:','linewidth',2)
hold on

%%%%% 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
k1=-(-pe1^(-1)*(210+252*(pc1-3/2))+76*(pc1-3/2)+152/3*(pc1-3/2)^2)/(pe1^(-1)*(210+252*(pc1-3/2)));   
E1=(((210+252*(C-(3/2))).*k1.^2)-((210+252*(C-(3/2)))*k1))./(-(76*(C-(3/2)).*k1+(152/3)*((C-(3/2)).^2).*k1));      
loglog(C,E1,'k-','linewidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                
hold on
%%%% 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k2=-(-pe2^(-1)*(210+252*(pc2-3/2))+76*(pc2-3/2)+152/3*(pc2-3/2)^2)/(pe2^(-1)*(210+252*(pc2-3/2)));   
E2=(((210+252*(C-(3/2))).*k2.^2)-((210+252*(C-(3/2)))*k2))./(-(76*(C-(3/2)).*k2+(152/3)*((C-(3/2)).^2).*k2));      
%loglog(C,E2,'b-')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
%%%% 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k3=-(-pe3^(-1)*(210+252*(pc3-3/2))+76*(pc3-3/2)+152/3*(pc3-3/2)^2)/(pe3^(-1)*(210+252*(pc3-3/2)));   
E3=(((210+252*(C-(3/2))).*k3.^2)-((210+252*(C-(3/2)))*k3))./(-(76*(C-(3/2)).*k3+(152/3)*((C-(3/2)).^2).*k3));      
%loglog(C,E3,'r-')


x1(1)=1.5;
x1(2)=1.5;
y1(1)=.001;
y1(2)=100;   
plot(x1,y1,'k:','linewidth',2)

x2(1)=1;
x2(2)=1000;
y2(1)=.1;
y2(2)=.1;
plot(x2,y2,'k-.','linewidth',2)		

% arredondamento de k para a segunda casa decimal    
                    % Multilicação feita de acordo com o 
k1=k1*100;        	% número significativo desejado;
k2=k2*100;
k3=k3*100;
                    % transformando o decimal em inteiro                   
k1=round(k1);		% Arredondamento do 
k2=round(k2);      	% número inteiro    
k3=round(k3);      	       

                    % Divisão retorna para decimal
k1=k1/100;   		% com os números significativos
k2=k2/100;         	% esperados   
k3=k3/100;         	   


plot(pc1,pe1,'ko')
plot(pc2,pe2,'ko')
plot(pc3,pe3,'ks')
plot(pc4,pe4,'ks')

plot(pc1,pe1,'kx')
plot(pc2,pe2,'kx')
plot(pc3,pe3,'ks')
plot(pc4,pe4,'ks')

if language==1;
   xlabel('Parâmetro circulação','fontsize',16)
   axe_Y = ['Parâmetro estratificação'];
   ylabel(axe_Y,'fontsize',16); 
   Ht_text2=get(gca,'ylabel');
   set(Ht_text2,'fontsize',esp_ylabel);
  % title('Diagrama Estratificação-Circulação')
   else
   xlabel('Circulation parameter','fontsize',16)
   axe_Y = ['Stratification parameter'];
   %ylabel(axe_Y);
   %ylabel('Stratification parameter','fontsize',16)
   Ht_text2=get(gca,'ylabel');
   set(Ht_text2,'fontsize',esp_ylabel);
  % title('Diagram Stratification-Circulation')
   end
end

axis([10^0 10^3 10^(-2) 2])

text(7,0.2,'2b')
text(7,0.03,'2a')
    
text(1.6,0.2,'1.0')
text(1.3,0.008,'1.5')
text(29,0.2,'0.01')

k1_str=num2str(k1);
%k2_str=num2str(k2);
k3_str=num2str(k3);
%k4_str=num2str(k4);
pri = ['gtext(','k1_str',')'];
eval(pri);
%segun = ['gtext(','k2_str',')'];
%eval(segun);    
terc = ['gtext(','k3_str',')'];
eval(terc);    
%quar = ['gtext(','k4_str',')'];
%eval(quar);    

%if salva==1
%print -dbitmap c:\Pc_Pe.bmp

if salva==1
%print -dbitmap c:\Pc_Pe.bmp

print -dbitmap e:\academica42\cap03\Diagr_4pto.bmp

end %diag.m
