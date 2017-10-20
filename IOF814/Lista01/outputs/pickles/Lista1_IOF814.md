
# Lista 1 - IOF814: Modelos Numéricos Aplicados a Processos Costeiros e Estuarinos

##### Q.01 - Demonstre que a solução da equação da difusão unidimensional linear que utiliza esquema centrado no tempo e espaço é incondicionalmente instável para uma solução explícita e incondicionalmente estável para uma solução implícita.

A equação da difusão unidimensional linear é dada por:

\begin{equation}
\frac{\partial{f}}{\partial{x}} = D\frac{\partial^2{f}}{\partial{x^2}}
\end{equation}

onde $D$ ($D > 0$) é o coeficiente de difusão.

Tomando a solução de (1) de forma explícita, com esquema centrado no espaço e no tempo, obtemos:

\begin{equation}
\frac{f^{n+1}_{j} - f^{n-1}_{j}}{2\Delta{t}} = \frac{D}{\Delta{x^2}}(f^{n}_{j+1} - 2f^{n}_{j} + f^{n}_{j-1})
\end{equation}

Ainda com (1), podemos tomar a solução implícita com equema centrado no espaço e no tempo, obtendo:

\begin{equation}
\frac{f^{n+1}_{j} - f^{n-1}_{j}}{2\Delta{t}} = \frac{D}{\Delta{x^2}}(f^{n+1}_{j+1} - 2f^{n+1}_{j} + f^{n+1}_{j-1})
\end{equation}

Para analisar a estabilidade das equações (2) e (3), precisamos aplicar o Método de Von Neumann para estabilidade.


**Para a Eq. (2)**

Admitindo-se uma solução harmônica no espaço e polinomial no tempo, da forma: 
\begin{equation}
    f^{n}_{j} = B^{n}_{j}e^{iKj\Delta{x}}
\end{equation}
onde $K$ é o número de onda.

Aplicamos (4) em (2), obtendo:

\begin{equation}
    \frac{B^{n+1}e^{iKj\Delta{x}} - B^{n-1}e^{iKj\Delta{x}} }{2\Delta{x}} = \frac{D}{\Delta{x^2}}(B^{n}e^{iK(j+1)\Delta{x}} - 2B^{n}e^{iKj\Delta{x}} + B^{n}e^{iK(j-1)\Delta{x}})
\end{equation}

A seguir, assumimos uma relação entre as amplitudes em instantes de tempo sucessivos, onde:

\begin{equation}
    B^{n+1} = \lambda B^{n} e B^{n} = \lambda B^{n-1} \Rightarrow B^{n-1} = \lambda^{-1}B^n
\end{equation}

Aplicando (6) em (5) e dividindo o resultado por (4), teremos:

\begin{equation}
    \frac{\lambda - \lambda^{-1}}{2\Delta{x}} = \frac{D}{\Delta{x^2}}(e^{iK\Delta{x}} - 2 + e^{-iK\Delta{x}})
\end{equation}

Lembrando que $e^{i\theta} = cos(\theta) - isin(\theta)$, então (7) é reescrita da seguinte forma:

\begin{equation}
    \frac{\lambda - \lambda^{-1}}{2\Delta{x}} = \frac{D}{\Delta{x^2}}(2cos(K\Delta{x}) - 2)
\end{equation}

Desenvolvendo (8) e considerando $q = 4\frac{\Delta{t}D}{\Delta{x^2}}$, temos:

\begin{equation}
    \lambda - \lambda^{-1} - q(cos(K\Delta{x} -1)) = 0
\end{equation}

Manipulando (9), chegamos em uma equação de 2º, da seguinte forma:

\begin{equation}
    \lambda^{2} - \lambda q(cos(K\Delta{x} -1)) -1 = 0
\end{equation}

A soluçao de (10) é dada por:

\begin{equation}
    \lambda = \frac{q[cos(K\Delta{x} -1)]^2 \pm \sqrt{q^2[cos(K\Delta{x})-1]^2 + 4} }{2}
\end{equation}

Analisando (11), nota-se que, para qualquer $\Delta{x}$ e $\Delta{t}$, a solução de $|\lambda| > 1$, ou seja, para qualquer passo de tempo e de espaço a solução implícita com esquema centrado no espaço e tempo será incondicionalmente instável.


**Para a Eq. (3)**

De forma análoga ao realizado para (2), fazemos para (3).

Aplicando (4) em (3):

\begin{equation}
        \frac{B^{n+1}e^{iKj\Delta{x}} - B^{n-1}e^{iKj\Delta{x}} }{2\Delta{x}} = \frac{D}{\Delta{x^2}}(B^{n+1}e^{iK(j+1)\Delta{x}} - 2B^{n+1}e^{iKj\Delta{x}} + B^{n+1}e^{iK(j-1)\Delta{x}})
\end{equation}

Assumindo a relação (6) e dividindo o resultado por (4), finalmente obtemos:

\begin{equation}
    \frac{\lambda - \lambda^{-1}}{2\Delta{x}} = \frac{D}{\Delta{x^2}}(\lambda e^{iK\Delta{x}} - 2\lambda + \lambda e^{-iK\Delta{x}})
\end{equation}

Manipulando (13), também chegamos em uma equação 2º:

\begin{equation}
    \lambda^{2}[1 - q(cos(K\Delta{x}) - 1)] = 1
\end{equation}

Portando, a solução de (14) será:

\begin{equation}
    \lambda = \frac{1}{\sqrt{1 - q(cos(K\Delta{x}) - 1)}}
\end{equation}

Analisando a equação (15), nota-se que, para qualquer passo de tempo e espaço, $|\lambda| \le 1$ e, portanto, a solução implícita com esquema centrado no tempo e espaço será incondionalmente estável.

__________________________________________________________________________________________

##### Q.02 - Implemente um modelo para a solução da equação da advecção bi-dimensional por esquema semi-implícito centrado no tempo e no espaço.

\begin{equation}
    \frac{\partial{f}}{\partial{t}} + u\frac{\partial{f}}{\partial{x}} + v\frac{\partial{f}}{\partial{y}} = 0
\end{equation}

O esquema semi-implícito centrado no tempo e espaço é dado através da média ponderada entre o esquema explícito e implícito, onde a soma do peso deve ser 1. Desta forma, precisamos desenvolver ambos os esquemas, para termos o requisitado no exercício.

**Esquema Explícito**

\begin{equation}
    \frac{ f^{n+1}_{j,k} - f^{n-1}_{j,k} }{2\Delta{x}} + u\frac{ f^{n}_{j+1,k} - f^{n}_{j-1,k} }{2\Delta{x}} + v\frac{f^{n}_{j,k+1} - f^{n}_{j,k-1}}{2\Delta{y}} = 0
\end{equation}


**Esquema Implícito**

\begin{equation}
    \frac{ f^{n+1}_{j,k} - f^{n-1}_{j,k} }{2\Delta{x}} + u\frac{ f^{n+1}_{j+1,k} - f^{n+1}_{j-1,k} }{2\Delta{x}} + v\frac{f^{n+1}_{j,k+1} - f^{n+1}_{j,k-1}}{2\Delta{y}} = 0
\end{equation}

A diferença entre (2) e (3) está no nível de tempo utilizado para calcular as derivadas espaciais, respectivamente, n (atual) e n+1 (renovado).

**Esquema Semi-Implícito**

\begin{equation}
    \frac{f^{n+1}_{j,k} - f^{n-1}_{j,k}}{2\Delta{t}} + \frac{u}{4\Delta{x}}(f^{n+1}_{j+1,k} - f^{n-1}_{j+1,k} + f^{n}_{j+1,k} - f^{n}_{j-1,k}) + \frac{v}{4\Delta{y}}( f^{n}_{j,k+1} - f^{n}_{j,k-1} + f^{n+1}_{j,k+1} -f^{n+1}_{j,k-1} ) = 0
\end{equation}

Isolando o termo $f^{n+1}_{j,k}$, obteremos:

\begin{equation}
    f^{n+1}_{j,k} = f^{n-1}_{,k} - \frac{\Delta{t}}{2\Delta{x}}(f^{n+1}_{j+1,k} - f^{n+1}_{j-1,k} + f^{n}_{j+1,k} - f^{n}_{j-1,k}) - \frac{\Delta{t}}{2\Delta{y}}( f^{n}_{j,k+1} - f^{n}_{j,k-1} + f^{n+1}_{j,k+1} - f^{n+1}_{j,k-1})
\end{equation}

Para resolver a equação (5), devemos reescrevê-la da seguinte forma:

\begin{equation}
   a_jf^{n+1}_{j-1} + b_jf^{n+1}_{j} + c_jf^{n+1}_{j+1} = d_j
\end{equation}

\begin{equation}
    a_kf^{n+1}_{k-1} + b_kf^{n+1}_{k} + c_kf^{n+1}_{k+1} = d_k
\end{equation}

Nesta etapa do desenvolvimento, convém trabalharmos com cada equação de forma independente. Sendo assim, reescrevemos (5) somente para __as variações em $x$__, de forma que obteremos:

\begin{equation}
    \frac{1}{2}\frac{\Delta{t}}{\Delta{x}}uf^{n+1}_{j-1} + f^{n+1}_{j} + \frac{1}{2}\frac{\Delta{t}}{\Delta{x}}uf^{n+1}_{j+1}  = f^{n-1}_{j} - \frac{1}{2}\frac{\Delta{t}}{\Delta{x}}u(f^{n}_{j+1} - f^{n}_{j-1})
\end{equation}

Onde, comparando (8) e (6), temos os seguintes coeficientes:

\begin{equation}
a_j = \frac{1}{2}\frac{\Delta{t}}{\Delta{x}}u; b_j = 1; c_j = \frac{1}{2}\frac{\Delta{t}}{\Delta{x}}u; d_j = f^{n-1}_{j} - \frac{1}{2}\frac{\Delta{t}}{\Delta{x}}u(f^{n}_{j+1} - f^{n}_{j-1})
\end{equation}

Como (8) é uma equação de três pontos, utilizamos o __método de inversão de linha__ (Lindzen & Kuo, 1969) para transformá-la em uma equação de dois pontos.

Inicialmente, definimos: 
\begin{equation}
    f^{n+1}_{j} = s_{j}f^{n+1}_{j+1} + p_j \Rightarrow f^{n+1}_{j-1} = s_{j-1}f^{n+1}_{j} + p_{j-1}
\end{equation}

Aplicando (9) e (10) em (8), teremos:

\begin{equation}
    a_j(s_{j-1}f^{n+1}_{j} + p_{j-1}) + b_jf^{n+1}_{j} + c_jf^{n+1}_{j+1} = d_j
\end{equation}

E, por fim, obtemos:

\begin{equation}
    s_j = \frac{-c_j}{b_j + a_js_{j-1}}; p_j = \frac{d_j - a_jp_{j-1}}{b_j + a_js_{j-1}}
\end{equation}

De modo análogo, podemos reescrever (5) para as __variações em $y$__, obtendo $s_k$ e $p_k$:

\begin{equation}
    \underbrace{-\frac{1}{2}\frac{\Delta{t}}{\Delta{y}}v}_{a_k}f^{n+1}_{k-1} + 
    \underbrace{\frac{1}{2}\frac{\Delta{t}}{\Delta{y}}v}_{c_k}f^{n+1}_{k+1}  = 
    \underbrace{-\frac{1}{2}\frac{\Delta{t}}{\Delta{y}}v(f^{n}_{k+1} - f^{n}_{k-1})}_{d_k}
\end{equation}

Sendo assim, obtemos:

\begin{equation}
    s_k = \frac{-c_k}{b_k + a_ks_{k-1}}; p_k = \frac{d_k - a_jp_{k-1}}{b_k + a_ks_{k-1}}
\end{equation}

Por fim, com as equações (11)-(14) podemos obter a solução geral do problema ao reescrever (10) na forma:

\begin{equation}
    f^{n+1}_{j,k} = s_{j}f^{n+1}_{j+1,k} + s_{k}f^{n+1}_{j,k+1}- p_{j} - p_{k}
\end{equation}

A solução numérica do exercício será realizada para a equação (16) e, todo o desenvolvimento numérico pode ser conferido no arquivo __lista1_Q02.m__. A seguir, analisamos os resultados obtidos no modelo elaborado.

<center>
    <h1><font color='red'>falta finalizar o código para colocar as imagens aqui e discuti-las</font></h1>
</center>

__________________________________________________________________________________________

##### Q.03 - O que são condições de contorno nas formas não gradiente, extrapolação linear e radiacional ? Dê exemplos para a equação da difusão uni-dimensional linear. Que consequências e que restrições deve ser consideradas ao usar estas condições de contorno ?



__________________________________________________________________________________________

##### Q.04 - Implemente dois modelos de advecção de um sinal retangular numa grande uni-dimensional, através de esquemas avançados no tempo, que utilizam diferenças finitas no espaço de 1ª ordem e de 4ª ordem. Compare os resultados dos dois modelos.


A equação da advecção é dada por:

\begin{equation}
    \frac{\partial{f}}{\partial{t}} + c\frac{\partial{f}}{\partial{x}} = 0
\end{equation}

Discretizando (1) para o esquema avançado no tempo e espaço de __1ª ordem__, teremos:

_____

\begin{equation}
    \frac{f^{n+1}_{j} - f^{n-1}_{j}}{\Delta{t}} + c\Biggl( \frac{f^{n}_{j+1} - f^{n}_{j}}{\Delta{x}} \Biggl) = 0
\end{equation}

A fórmula de recorrência será então:

\begin{equation}
    f^{n+1}_{j} = f^{n}_{j} - c\frac{\Delta{t}}{\Delta{t}}(f^{n}_{j+1} - f^{n}_{j})
\end{equation}

E, na forma matricial, (3) ficará:

\begin{equation}
fren(2:jmax-1) = fatu(2:jmax-1) - c\frac{dt}{dx}\biggr[ fatu(3:jmax) - fatu(2:jmax-1) \biggr]
\end{equation}

_________

Discretizando (1), para o mesmo esquema, mas desta vez de __4ª ordem__, é necessário trabalhar de forma diferente. Sendo assim, desenvolvemos inicialmente a discretização no esquema avançado no tempo:

\begin{equation}
    \frac{\partial{f}}{\partial{t}} = \frac{f^{n+1}_{j} - f^{n}_{j}}{\Delta{t}}
\end{equation}

Já para o esquema de 4ª ordem no espaço, teremos um conjunto de 5 equações:

Uma geral, em que os pontos de grade calculados são de 3:jmax-2, pensando em um aplicação no matlab, onde jmax é o a quantidade máxima de pontos em j:
\begin{equation}
    \frac{\partial{f}}{\partial{x}} = \frac{\biggl( f^{n}_{j-2} - 8f^{n}_{j-1} + 8f^{n}_{j+1} - f^{n}_{j+2} \bigg)}{12\Delta{x}}
\end{equation}

Para, respectivamente, o 1º, 2º, penúltimo e último ponto de grade, teremos as seguintes equações:

\begin{equation}
\begin{cases}
    \frac{\partial{f}}{\partial{x}} = \frac{\biggl( - 25f^{n}_{j} + 48f^{n}_{j+1} - 36f^{n}_{j+2} + 16f^{n}_{j+3} - 3f^{n}_{j+4} \bigg)}{12\Delta{x}} & \text{(A)}\\
    \frac{\partial{f}}{\partial{x}} = \frac{\biggl( - 3f^{n}_{j-1} - 10f^{n}_{j} + 18f^{n}_{j+1} - 6f^{n}_{j+2} + f^{n}_{j+3} \bigg)}{12\Delta{x}} & \text{(B)}\\
    \frac{\partial{f}}{\partial{x}} = \frac{\biggl( - f^{n}_{j-3} + 6f^{n}_{j-2} - 18f^{n}_{j-1} + 10f^{n}_{j} + 3f^{n}_{j+1} \bigg)}{12\Delta{x}} & \text{(C)}\\
    \frac{\partial{f}}{\partial{x}} = \frac{\biggl( 3f^{n}_{j-4} - 16f^{n}_{j-3} + 36f^{n}_{j-2} - 48f^{n}_{j-1} + 25f^{n}_{j} \bigg)}{12\Delta{x}} & \text{(D)}
\end{cases}
\end{equation}


Sendo assim, ao programarmos um modelo com a equação de advecção discretiza avançada no tempo e de 4ª ordem no espaço, teremos o seguinte conjunto de equações a ser implementado.

Iniciamos o desenvolvimento da discretização de (1) pelo caso mais geral, equações (5) e (6):

\begin{equation}
    \frac{f^{n+1}_{j} - f^{n}_{j}}{\Delta{t}} + \frac{c}{12\Delta{x}}\biggl( f^{n}_{j-2} - 8f^{n}_{j-1} + 8f^{n}_{j+1} - f^{n}_{j+2} \bigg) = 0
\end{equation}

Que, isolando $f^{n+1}_{j}$ obtemos a função de renovação:

\begin{equation}
    f^{n+1}_{j} = f^{n}_{j}-c\frac{\Delta{t}}{12\Delta{x}}\biggl( f^{n}_{j-2} - 8f^{n}_{j-1} + 8f^{n}_{j+1} - f^{n}_{j+2} \bigg)    
\end{equation}

De forma análoga, podemos reescrever o conjunto de equações (7), isolando o nível renovado de tempo:

\begin{equation}
    \begin{cases}
            f^{n+1}_{j} = f^{n}_{j}-c\frac{\Delta{t}}{12\Delta{x}}\biggl( - 25f^{n}_{j} + 48f^{n}_{j+1} - 36f^{n}_{j+2} + 16f^{n}_{j+3} - 3f^{n}_{j+4} \bigg) & (A) \\
            f^{n+1}_{j} = f^{n}_{j}-c\frac{\Delta{t}}{12\Delta{x}}\biggl( - 3f^{n}_{j-1} - 10f^{n}_{j} + 18f^{n}_{j+1} - 6f^{n}_{j+2} + f^{n}_{j+3}  \bigg) & (B) \\
            f^{n+1}_{j} = f^{n}_{j}-c\frac{\Delta{t}}{12\Delta{x}}\biggl( - f^{n}_{j-3} + 6f^{n}_{j-2} - 18f^{n}_{j-1} + 10f^{n}_{j} + 3f^{n}_{j+1} \bigg) & (C) \\
            f^{n+1}_{j} = f^{n}_{j}-c\frac{\Delta{t}}{12\Delta{x}}\biggl( 3f^{n}_{j-4} - 16f^{n}_{j-3} + 36f^{n}_{j-2} - 48f^{n}_{j-1} + 25f^{n}_{j} \bigg) & (D)
    \end{cases}
\end{equation}

A solução numérica de ambos os esquemas de discretização é realizada no arquivo lista1_Q04.m e, a seguir, seguem os resultados obtidos nos modelos elaborados.

<font color='red'>PROGRAMAR</font>

__________________________________________________________________________________________

##### Q.05 - Determine os coeficientes da forma discretiza da equação da difusão

\begin{equation}
    \frac{\partial{f}}{\partial{t}} - D\frac{\partial^2{f}}{\partial{x^2}} = \alpha f^{n+1}_{j} + \beta f^{n}_{j+1} + \gamma f^{n}_{j-1} + \delta f^{n-1}_{j}
\end{equation}

##### a partir de expressões da Série de Taylor. Qual é a principal utilidade do "Método dos coeficientes indefinidos"?

Inicialmente expandimos os termos do lado direito da equação (1), truncando a expansão no terceiro termo:

\begin{equation}
\begin{aligned}
    f^{n+1} = f^{n} + \Delta{t}\frac{\partial{f}}{\partial{t}} + \frac{\Delta{t}^2}{2}\frac{\partial^2{f}}{\partial{t^2}} + ...
    \\
    f^{n-1} = f^{n} - \Delta{t}\frac{\partial{f}}{\partial{t}} + \frac{\Delta{t}^2}{2}\frac{\partial^2{f}}{\partial{t^2}} - ...
    \\
    f_{j+1} = f_{j} + \Delta{x}\frac{\partial{f}}{\partial{x}} + \frac{\Delta{t}^2}{2}\frac{\partial^2{f}}{\partial{x^2}} + ...
    \\
    f_{j-1} = f_{j} - \Delta{x}\frac{\partial{f}}{\partial{x}} + \frac{\Delta{t}^2}{2}\frac{\partial^2{f}}{\partial{x^2}} - ...
\end{aligned}
\end{equation}

Aplicando o conjunto de equações (2) em (1), obteremos:

\begin{equation}
\begin{aligned}
    \frac{\partial{f}}{\partial{t}} - D\frac{\partial^2{f}}{\partial{x^2}} = 
    \alpha (f^{n} + \Delta{t}\frac{\partial{f}}{\partial{t}} + \frac{\Delta{t}^2}{2}\frac{\partial^2{f}}{\partial{t^2}}) + 
    \beta (f_{j} + \Delta{x}\frac{\partial{f}}{\partial{x}} + \frac{\Delta{t}^2}{2}\frac{\partial^2{f}}{\partial{x^2}}) + 
    \gamma (f^{n} - \Delta{t}\frac{\partial{f}}{\partial{t}} + \\ + \frac{\Delta{t}^2}{2}\frac{\partial^2{f}}{\partial{t^2}}) + 
    \delta (f_{j} - \Delta{x}\frac{\partial{f}}{\partial{x}} + \frac{\Delta{t}^2}{2}\frac{\partial^2{f}}{\partial{x^2}})
\end{aligned}
\end{equation}

Associando os termos do lado direito aos termos do lado esquerdo de (3), obteremos o seguinte sistema de equações:

\begin{equation}
  \begin{cases}
    \alpha + \beta + \gamma + \delta = 0 & \text{(a)}\\
    \alpha \Delta{t} - \delta \Delta{t} = 0  & \text{(b)} \\
    \beta \Delta{x} - \gamma \Delta{x} = 0 \to \beta = \gamma & \text{(c)}\\
    \beta \frac{\partial{\Delta{x^2}}}{2} + \gamma \frac{\partial{\Delta{x^2}}}{2} = -D \to \beta + \gamma = \frac{-2D}{\Delta{x^2}} & \text{(d)}
  \end{cases}
\end{equation}

Para determinar os coeficientes $\alpha$, $\beta$, $\gamma$ e $\delta$, precisamos resolver o sistema de equações (4).

__Determinando $\alpha$__

Multiplica-se (4.a) por $\Delta{t}$ e soma-se (4.b):

\begin{equation}
    (\alpha + \beta + \gamma + \delta)\Delta{t} + \alpha \Delta{t} - \gamma \Delta{t} = 0 \to 2\alpha{\Delta{t}} + \beta{\Delta{t}} + \gamma{\Delta{t}} = 0 \\ 
    \therefore 2\alpha + \beta + \gamma = \frac{1}{\Delta{t}}
\end{equation}

Multiplicando (4.d) por (-1) e somando (5), obtemos:

\begin{equation}
    \alpha = \frac{1}{2\Delta{t}} + \frac{D}{\Delta{x^2}}
\end{equation}

__Determinando $\delta$__

Aplica-se (6) em (4.b):

\begin{equation}
    \begin{aligned}
        (\frac{1}{2\Delta{t}} + \frac{D}{\Delta{x^2}})\Delta{t} - \delta \Delta{t} = 0 \to  \frac{1}{2} + \frac{\Delta{t}D}{\Delta{x^2}} - \delta{\Delta{t}} = 1 \to \delta = -\frac{1}{\Delta{t}} + \frac{D}{\Delta{x^2}}
    \end{aligned}
\end{equation}

__Determinando $\beta$__

Aplicando (4.c) em (4.d):

\begin{equation}
    \begin{aligned}
        \beta{\frac{\Delta{x^2}}{2}} + \beta{\frac{\Delta{x^2}}{2}} = -D \to \beta{\Delta{x^2}} = -D \\ 
        \therefore \beta = -\frac{D}{\Delta{x^2}}
    \end{aligned}
\end{equation}

__Determinando $\gamma$__

Como de (4.c) $\beta = \gamma$, então:

\begin{equation}
    \gamma = -\frac{D}{\Delta{x^2}}
\end{equation}

Finalmente, aplicamos as equações (6) a (9) em (1) e obteremos:

\begin{equation}
\frac{\partial{f}}{\partial{t}} - D\frac{\partial^2{f}}{\partial{x^2}} = 
    (\frac{1}{2\Delta{t}} + \frac{D}{\Delta{x^2}})f^{n+1}_{j} + 
    (-\frac{D}{\Delta{x^2}})f^{n}_{j+1} + 
    (-\frac{D}{\Delta{x^2}})f^{n}_{j-1} + 
    (-\frac{1}{2\Delta{t}} + \frac{D}{\Delta{x^2}})f^{n-1}_{j}
\end{equation}

Desenvolvendo a equação (10), obteremos:

\begin{equation}
    \frac{\partial{f}}{\partial{t}} - D\frac{\partial^2{f}}{\partial{x^2}} = 
    \frac{f^{n+1}_{j}}{2\Delta{t}} + \frac{Df^{n+1}_{j}}{\Delta{x^2}} -
    \frac{Df^{n}_{j+1}}{\Delta{x^2}} - \frac{Df^{n}_{j-1}}{\Delta{x^2}} -
    \frac{f^{n-1}_{j}}{2\Delta{t}} + \frac{Df^{n-1}_{j}}{\Delta{x^2}} 
\end{equation}

Reorganizando a equação (11), finalmente obtemos a equação da difusão discretizada com esquema centrado no tempo e no espaço pseudo-implícito:

\begin{equation}
    \frac{f^{n+1}_{j} - f^{n-1}_{j}}{2\Delta{t}} = 
    \frac{D}{\Delta{x^2}}(f^{n}_{j+1} + f^{n}_{j-1} - f^{n+1}_{j} - f^{n-1}_{j})
\end{equation}

A principal utilidade de se utilizar o __"Método dos coeficientes indefinidos"__ é que podemos usar mais pontos no esquema de discretização. Neste caso, foram utilizados 4 pontos e, portanto, 4 coeficientes foram definidos. Em caso de de usar 6, 8, 10 pontos, poderíamos ter, respectivamente, 6, 8 ou 10 coeficientes a serem determinados.

__________________________________________________________________________________________

##### Q.06 - Implemente um modelo numérico de advecção-difusão-decaimento 2D para a área costeira ao largo de Santos (SP), para os limites (46.5ºW - 46.2ºW; 23.95ºS - 24.15ºS). Processe o modelo para a dispersão de um contaminante despejado de forma contínua no ponto 46.35ºW 24.01ºS (simulando a operação do emissário submarino) e a dispersão de um contaminante despejado instantaneamente em 46.35ºW 24.10ºS (simulando acidente com embarcação tem trânsito). Processe o modelo com correntes de 1m/s para Norte, Nordeste e Noroeste. Forneça mapas das plumas de contaminantes e detecte que áreas costeiras podem ser atingidas em cada caso.


