%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% PROGRAMA vel_eq_9_41.m                           %
% AJUSTE PERFIL DE VELOCIDADE                      %
% DA CAMADA DE FUNDO DO ESTU�RIO CUNHA SALINA      %
% LUIZ BRUNER DE MIRANDA                           %
% TEORIA Miranda et al. (2011) - (eq. 9.41 p. 328) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  Simulacao de perfis de velocidade estacion�rios 
%  Foi utilizada a equa��o de Miranda et al. (2011) (eq. 9.41 p. 328)
%  Programa��o de Fernando Pinheiro Andutta & Luiz Bruner de Miranda
%  Estu�rios tipo cunha salina 
%  Equil�brio dos componentes barotr�pico, descarga fluvial, tens�o
%  interfacial de atrito e o atrito interno na cunha salina e atrito m�ximo no fundo"
%  For�antes e parcelas que dependem da profundidade (z) tiveram o sinal trocado
%  em rela��o � orienta��o utilizada no livro (Miranda et al., 2011), logo,
%  a profundidade z est� orientada positivamente para cima com origem na superf�cie -1<Z<0
%  Tens�o de cisalhamento do vento foi considerada desprez�vel

clear
clc
close all

%COEFICIENTES, PAR�METROS (geometria, descarga fluvial)

nr=1; % input('Numero da simulacao:');
Qf=1040.; % Descarga fluvial em m3/s;
Ho=10.0; % espessura da  camada na cabeceira do estu�rio (m) - Ho:');
h1=4.0; % espessura da  camada sobrejacente � cunha salina (m) - h1:');
h2=6.0; % espessura da  camada da cunha salina (m) - h2:');
B=1000.0 % Largura do estu�rio na posi��o x=7,8 km da boca do estu�rio


Y=num2str(nr);
diario=['diary simulacao',Y,'.txt'];
eval(diario)    

% C�lculo da velocidade gerada pela descarga fluvial "u1" na camada
% sobrejacente � cunha salina

u1=Qf/(h1*B)

% Formula��o matem�tica - perfil de velocidade - (u=u(x,z)) Equa��o 9.41 p.328

z=-4.0:-0.1:-10.0

coef1=-(2*u1/h2)

coef2=(3.0*u1/h2^2)

u2z=coef1*(Ho+z)+coef2*[(Ho+z).^2]

vmedia=mean(u2z)

% Resultados

x=[0,0]
X=[0,-10]
 
xx=[0.26,0.26]
yy=[0.26,-4.0]

xxx=[-0.1,0.3]
yyy=[0.0,0.0]

figure

f1 = plot(u2z,z,'k');
set(f1,'linewidth',3);
hold

line(x,X)
line(xx,yy)
line(xxx,yyy)

legend('Barotr�pico+Tens�o interfacial',4)
xlabel('Componente, u(m s^{-1})','fontsize',14);
ylabel('Profundidade, z','fontsize',14);

gtext('Vazante','fontsize',14)
gtext('Enchente','fontsize',14)

salva=['print -dbitmap simulacao',Y,'.bmp'];
eval(salva)   

diary off




 


