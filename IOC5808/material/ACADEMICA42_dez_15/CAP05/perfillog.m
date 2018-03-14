
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROGRAMA perfillog.m                       %
% AJUSTE LORARÍTMICO DE VELOCIDADE           %
% LUIZ bRUNER DE MIRANDA                     %
% ALESSANDRO LUVIZON BÉRGAMO                 %
% MARÇO DE 2001                              %                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

clear

% carrega o arquivo de dados - apenas um perfil adimensional (0-1)

% o nome do arquivo seguido do numero


cr=input('Digite o numero do arquivo:');

     YY=num2str(cr);
     
     ldarq=['load perfil',YY,'.dat'];
     nomearq=['perfil',YY];
    
     eval(ldarq);   
     X=eval(nomearq); 
      
   z=X(:,1);
   vu=X(:,2);
   
   % Converte para modulo da velocidade e de cm/s para m/s
   
  %vu=abs(vu)/100;
  vu=abs(vu);

   % Muda o sentido do eixo 
      
   u=flipud(vu);
   
   % Traça os perfis com ambos os eixos
   
   plotyy(vu,-z,u,z)
      
   xlabel('|u| # (m s^{-1})','fontsize',14)
   ylabel('z # (m)','fontsize',14)
   
   legend ('Origem Superficie', 'Origem Fundo',4) 
   
   print -dbitmap f:\cap5\perfillog\figura1.bmp
   
   pause   
         
   figure
   
   % Calculo da derivada du/dz
   
   du=diff(u)./diff(z);
   zd=z(2:end);
   
   hd=1./zd;
   
   % Ajuste polinomial du/dz x 1/z
   
   p=polyfit(du,hd,1)
   
   hi=polyval(p,du);
   
   subplot(1,2,1)
   
   plot(u,z)
   xlabel('Velocidade (m s^{-1})','fontsize',14)
   ylabel('Espessura da Coluna de Agua -z (m)','fontsize',14)
   
   hold

   subplot(1,2,2)
         
   plot(hd,du,'b')
   
   hold
   
   plot(hi,du,'g')
   
   hold
   
   legend('du/dz','Linear',4)
   ylabel('(du/dz) - (s^{-1}) e Ajuste Linear','fontsize',14)
   xlabel('(1/z) - (m^{-1})','fontsize',14)
   
   print -dbitmap f:\CAP5\perfillog\figura2.bmp

   
   pause
   
   figure
   
   plot(hi,du,'g')
   ylabel('du/dz # (s^{-1})','fontsize',14)
   xlabel('1/z # (m^{-1})','fontsize',14)
   
    
 % Calculo da velocidade de atrito u*
   
   a=p(1);
   
   k=0.4;
   
   % Observe o fator (1/10) de escala entre os eixos coordenados da correlação linear
   % entre du/dz e 1/z. Atençao esse fator poderá mudar, por ex. para 1/1000 !!!
   
   ux=(a*k)/10;
   
   % ou
   %ux=(a*k)/1000
   
   I=num2str(ux);
   
   gtext('u_*(m s^{-1})=')
   texto=['gtext(','I',')'];
   eval(texto);
   
   print -dbitmap f:\CAP5\perfillog\figura3.bmp

   
   pause
   
    
   % Cálculo do coeficiente de arrasto, CD, com diferentes velocidades de referência:
   
   % a-velocidade de referência calculada pela media integrada do perfil de velocidade
   
   
   l=length(vu);
   
   up=vu(2:end-1);
   
   urm=(((vu(1)+vu(end))/2)+sum(up))/(l-1)
   
   cdm=(ux/urm)^2
   
   % Atençao: Falta incluir comando para saida de CDM
   
   
   % b-velocidade de referência calculada pela velocidade a 1 m do fundo
   
      
   ur100=vu(end)

   cd100=(ux/ur100)^2
   
   
   % Atençao: Falta incluir comando para saida de CD100
   
   
   % Calculo do perfil logaritmico u(z)
   % Quando solicitado deve-se entrar com o valor de z0,
   % o programa roda até com 10 valores diferentes
   % caso não haja necessidade interrompe-se o processamento
   % teclando ctrl+C
   
   zr1=input('Digite o valor de z0:');

   Y=num2str(zr1);
   eq=['uz=(ux/k)*log(zd/',Y,');'];
   eval(eq);
   
   figure
   
   plot(u,z,'b*')
   hold
   plot(uz,zd,'g')
   hold
   legend ('u','u(z)',2)
   xlabel('velocidade (m s^{-1})','fontsize',14)
   ylabel('Espessura da coluna de água -z (m)','fontsize',14)
   
   print -dbitmap f:\CAP5\perfillog\figura4.bmp

   
   I=num2str(Y);
   gtext('z_{o}(m)=')
   texto=['gtext(','I',')'];
   eval(texto);
   
   U=num2str(ux);
   gtext('u_*(m s^{-1})=')

   texto=['gtext(','U',')'];
   eval(texto);       
      
   zr2=input('Digite um novo valor de z0:');

   Y=num2str(zr2);
   eq=['uz=(ux/k)*log(zd/',Y,');'];
   eval(eq);
     
   plot(u,z,'b*')
   hold
   plot(uz,zd,'g')
   hold
   legend ('u','u(z)',2)
   xlabel('u : u(z)')
   ylabel('z')
   
   gtext('z_{o}(m)=')
   I=num2str(Y);
   texto=['gtext(','I',')'];
   eval(texto);
   
   U=num2str(ux);
   gtext('u_*(m s^{-1})=')
   texto=['gtext(','U',')'];
   eval(texto);       
   
    
   zr3=input('Digite um novo valor de z0:');

   Y=num2str(zr3);
   eq=['uz=ux/k*log(zd/',Y,');'];
   eval(eq);
     
   plot(u,z,'b*')
   hold
   plot(uz,zd,'g')
   hold
   legend ('u','u(z)',2)
   xlabel('u : u(z)')
   ylabel('z')
   
   gtext('z_{o}(m)=')
   I=num2str(Y);
   texto=['gtext(','I',')'];
   eval(texto);

   U=num2str(ux);
   gtext('u_*(m s^{-1})=')
   texto=['gtext(','U',')'];
   eval(texto);       
   
   zr4=input('Digite um novo valor de z0:');

   Y=num2str(zr4);
   eq=['uz=ux/k*log(zd/',Y,');'];
   eval(eq);
     
   plot(u,z,'b*')
   hold
   plot(uz,zd,'g')
   hold
   legend ('u','u(z)',2)
   xlabel('u : u(z)')
   ylabel('z')
    
   gtext('z_{o}(m)=')
   I=num2str(Y);
   texto=['gtext(','I',')'];
   eval(texto);
   
   U=num2str(ux);
   gtext('u_*(m s^{-1})=')
   texto=['gtext(','U',')'];
   eval(texto);       
   
   zr5=input('Digite um novo valor de z0:');

   Y=num2str(zr5);
   eq=['uz=ux/k*log(zd/',Y,');'];
   eval(eq);
     
   plot(u,z,'b*')
   hold
   plot(uz,zd,'g')
   hold
   legend ('u','u(z)',2)
   xlabel('u : u(z)')
   ylabel('Z')
   
   gtext('z_{o}(m)=')
   I=num2str(Y);
   texto=['gtext(','I',')'];
   eval(texto);
   
   U=num2str(ux);
   gtext('u_*(m s^{-1})=')
   texto=['gtext(','U',')'];
   eval(texto);       
   
   zr6=input('Digite um novo valor de z0:');

   Y=num2str(zr6);
   eq=['uz=ux/k*log(zd/',Y,');'];
   eval(eq);
     
   plot(vu,z,'b*')
   hold
   plot(uz,zd,'g')
   hold
   legend ('u','u(z)',2)
   xlabel('u : u(z)')
   ylabel('z')
   
   gtext('z_{o}(m)=')
   I=num2str(Y);
   texto=['gtext(','I',')'];
   eval(texto);
   
   U=num2str(ux);
   gtext('u_*(m s^{-1})=')
   texto=['gtext(','U',')'];
   eval(texto);       
   
   zr7=input('Digite um novo valor de z0:');

   Y=num2str(zr7);
   eq=['uz=ux/k*log(zd/',Y,');'];
   eval(eq);
     
   plot(u,z,'b*')
   hold
   plot(uz,zd,'g')
   hold
   legend ('u','u(z)',2)
   xlabel('u : u(z)')
   ylabel('z')
   
   gtext('z_{o}=')
   I=num2str(Y);
   texto=['gtext(','I',')'];
   eval(texto);
   
   U=num2str(ux);
   gtext('u_*(m s^{-1})=')
   texto=['gtext(','U',')'];
   eval(texto);       
   
   
   zr8=input('Digite um novo valor de z0:');

   Y=num2str(zr8);
   eq=['uz=ux/k*log(zd/',Y,');'];
   eval(eq);
     
   plot(u,z,'b*')
   hold
   plot(uz,zd,'g')
   hold
   legend ('u','u(z)',2)
   xlabel('u : u(z)')
   ylabel('z')
   
   gtext('z_{o}=')
   I=num2str(Y);
   texto=['gtext(','I',')'];
   eval(texto);

   U=num2str(ux);
   gtext('u_*(m s^{-1})=')
   texto=['gtext(','U',')'];
   eval(texto);       
   
   zr9=input('Digite um novo valor de z0:');

   Y=num2str(zr9);
   eq=['uz=ux/k*log(zd/',Y,');'];
   eval(eq);
     
   plot(u,z,'b*')
   hold
   plot(uz,zd,'g')
   hold
   legend ('u','u(z)',2)
   xlabel('u : u(z)')
   ylabel('z')
   
   gtext('z_{o}=')
   I=num2str(Y);
   texto=['gtext(','I',')'];
   eval(texto);
   
   U=num2str(ux);
   gtext('u_*(m s^{-1})=')
   texto=['gtext(','U',')'];
   eval(texto);       
   
   zr10=input('Digite um novo valor de z0:');

   Y=num2str(zr10);
   eq=['uz=ux/k*log(zd/',Y,');'];
   eval(eq);
     
   plot(u,z,'b*')
   hold
   plot(uz,zd,'g')
   hold
   legend ('u','u(z)',2)
   xlabel('u : u(z)')
   ylabel('z')
   
   gtext('z_{o}=')
   I=num2str(Y);
   texto=['gtext(','I',')'];
   eval(texto);
   
   U=num2str(ux);
   gtext('u_*(m s^{-1})=')
   texto=['gtext(','U',')'];
   eval(texto);  
   
   



   
   