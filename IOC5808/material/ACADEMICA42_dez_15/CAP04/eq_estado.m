% Equação de estado não linear ro=ro(S,T)
% Mamayev (1975) e Miranda (2002, 2012)
% Mais detalhes em Miranda et al. (2012) p. 147-161 es. 4.15 e 4.16
% Vide também Morgan (1994) - Sub-rotina "sw_dens0"


% Leitura dos dados termohalinos %(temperatura, salinidade) 
% Na forma de matrizes coluna

load dados_TS.dat

temp=dados_TS(:,3)
sal=dados_TS(:,4)

s35=linspace(35,35,99);

% Calculo da matriz transposta

s35=s35';
ss=(sal-s35)

%Cálculo da equação de estado de Mamayev (1975)

sigma=28.152-temp.*7.35e-002-temp.^2*4.69e-003+(0.802-temp.*2.0e-003).*ss
roma=sigma+1000.0

%Cálculo da equação de estado de Miranda et al. (2012) 


romi=1028.106-temp.*8.68243e-002-temp.^2*3.819566e-003+(8.036307e-001-temp.*1.7117e-003).*ss
sigmi=romi-1000.0

close all
