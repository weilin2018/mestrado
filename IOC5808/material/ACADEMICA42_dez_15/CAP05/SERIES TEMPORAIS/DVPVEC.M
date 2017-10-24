function Y=dvpvec(u,v,dt)
%dt = delta T, intervalo de amostragem.
%dt deve ser transformado em segundos:
dt = dt*3600;          %para dt originalmente em minutos.

clc
np=max(size(u)); 
plot(0,0)
clf
hold

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dessoma=0;

%k=1;

x(2)=0;
y(2)=0;

cntemp=0;
for i=1:np ,
 x(1)=x(2);
 y(1)=y(2);
 x(2)=u(i)*dt/10^5 + x(1);
 y(2)=v(i)*dt/10^5 + y(1);
 plot(x,y,'b')
 
 % Faz uma marca a cada 24horas 
 cntemp=cntemp+1;
 if cntemp == (24*3600/dt)
    %plot(x(2),y(2),'or')
    %plot(x(2),y(2),'.r')
    cntemp=0;
 end
 
 %Calcula a distância total percorrida  
 desloca(i)=sqrt((u(i)*dt/10^5)^2 + (v(i)*dt/10^5)^2)+ dessoma;
 dessoma=desloca(i);
end

axis('equal')
xlabel('Distância (km)','fontsize',14)
ylabel('Distância (km)','fontsize',14)
%title('Diagrama vetorial progressivo da corrente, (+) início e (o) a cada 24h')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORIGEM %
tam=.01*max(abs(desloca));  
x(1)=0;
y(1)=0-tam;
x(2)=0;
y(2)=0+tam;
plot(x,y,'r')

x(1)=0-tam;
y(1)=0;
x(2)=0+tam;
y(2)=0;
plot(x,y,'r')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Y=desloca;

