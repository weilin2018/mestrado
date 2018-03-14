% Exercicio 01.01 - calcular primeiras derivadas de f(x)=sin(x) através de fórmula analíticas
% e de diferencas finitas de 2a e 3a ordem.


clear all; close all;

% configurando os parametros gerais
dx = 0.1;
x = 0:dx:pi;
f = sin(x);

% calculo usando a formula analitica
deranalit = cos(x);
deranalit = deranalit';
% '

% definindo as matrizes que armazenarao as soluções discretas
jmax = size(x,2);
fava = zeros(jmax,1);            % avancada
fret = zeros(jmax,1);            % retardada
fcen = zeros(jmax,1);            % centrada

% calculos
fava(1:jmax-2) = ( -3*f(1:jmax-2) + 4*f(2:jmax-1) - f(3:jmax) ) / (2*dx);
fret(3:jmax)   = (  3*f(3:jmax) - 4*f(2:jmax-1) + f(1:jmax-2) ) / (2*dx);
fcen(3:jmax-2) = 2*( (f(4:jmax-1) - f(2:jmax-3)) / (2*dx) ) - (f(5:jmax) - f(1:jmax-4)) / (4*dx);

% plotar os dados calculados
figure(1);
plot(x,f,'LineWidth',2)
grid on
hold
plot(x(1:jmax-2),fava(1:jmax-2),'r','LineWidth',2)
plot(x(3:jmax),fret(3:jmax),'g','LineWidth',2)
plot(x(3:jmax-2),fcen(3:jmax-2),'k','LineWidth',2)
axis([x(1) x(jmax) -inf inf])
title('Funcao [em azul] e diferencas finitas [ avancada em verm, retardada em verde e centrada em preto]', 'fontsize', 12)
ylabel('f, df/dx alta ordem','fontsize',12)
xlabel('x','fontsize',12)

%% calcular as estatísticas de erro
% de forma simples, podemos calcular o erro apenas ao subtrair a solucao analitica
% da solucao discreta:

% definir matrizes de erro
erro_ava = zeros(jmax,1);            % avancada
erro_ret = zeros(jmax,1);            % retardada
erro_cen = zeros(jmax,1);            % centrada

% calcular erro de cada metodo
erro_ava(1:jmax-2) = deranalit(1:jmax-2) - fava(1:jmax-2);
erro_ret(3:jmax)   = deranalit(3:jmax) - fret(3:jmax);
erro_cen(3:jmax-2) = deranalit(3:jmax-2) - fcen(3:jmax-2);

% plotar o erro
figure(2)
plot(x(1:jmax-2),erro_ava(1:jmax-2),'r','LineWidth',2)
grid on
hold
plot(x(3:jmax),erro_ret(3:jmax),'g','LineWidth',2)
plot(x(3:jmax-2),erro_cen(3:jmax-2),'k','LineWidth',2)
axis([x(1) x(jmax) -0.004 0.004])
title('Erros p/ deriv analit: av(verm), ret(verde), centr(preto)','fontsize',12)
ylabel('Erro - alta ordem','fontsize',12)
xlabel('x','fontsize',12)
