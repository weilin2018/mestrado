Queridos alunos,

A tarefa referente ao mapeamento de função de corrente a partir dos dados sinóticos é a seguinte. O cruzeiro a ser utilizado é o WESTRAX 2 e tenham em mãos o meu velho artigo de 2000. Adotem para a AOE e AOV os parâmetros que forneço. Cuidado com as unidades das saídas! Sugiro a conversão da velocidade do perfilador Pegasus para graus/s quando entrarem na AOV. Depois da interpolação, reconvertam para m/s.

1)  Mapeamento de Função de Corrente Geostrófica​

​a) ​ Calculem a função de corrente geostrófica nos pontos dos dados, ou seja, peguem os pontos de obs. e calculem a anomalia do geopotencial para 10 dbar relativamente a 1200 dbar. A posteriori, dividam por f0 com f avaliado a 5N (tal qual no artigo da minha tese). Removam a média do vetor formado pelos valores de psi_g.

b) Precisamos satisfazer as condições de contorno. Comecemos com a mais fácil de implementar: no slip. Obtenham de arquivo a isóbata de 200 m, seja pelo ETOPO ou GEBCO. Suavizem-na com uma média corrida, janela móvel ou reinterpolem usando um esquema cúbico. Vcs escolhem.

c) Selecionem lats e lons dessa isóbata e atribuam o valor de zero na borda.

d) Componha o vetor juntando os dados de verdade com os dados falsos da c.c.

e) Rode a rotina scaloa.m ou a versão python que o Iurizinho.

​f) Plote e verifique se o padrão é semelhante ao do artigo.

g) Agora repita o mesmo caso, mas use imagens. Faça a aproximação linear da isóbata de 200 m e a use como eixo de simetria.​


​2) Mapeamento de Função de Corrente Observada

a) o procedimento é semelhante ao da AOE, mas não removam a média dos vetores. Mapeiem a psi_obs em 10 dbar. A rotina tem um parâmetro especial que estabelece o nível médio de psi_obs=0.

b) Usem tanto a no-slip como as imagens (free slip)​ como cc.

c) Comparemos os padrões.

d) se quiserem, podem fazer também a velocidade relativa a 1200 dbar para direta comparação entre psi_g e psi_obs. Para tanto, basta subtrairem o vetor de 1200 dbar do vetor de 10 m.

​NÃO EMPAQUEM. EM CASO DE DÚVIDA, PROCUREM-ME. ESTOU SEM AULAS PELA TARDE, A MENOS DAS 2 QUARTAS DE REPOSIÇÃO DE DFG2.​
