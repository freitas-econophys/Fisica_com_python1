# Tarefa 1

- **Disciplina:** Física com Python 1 (2026.2) — CBPF
- **Autor:** Isaque Porto de Freitas
- **Obs.:** Código desenvolvido para fins de avaliação acadêmica. Disponibilizado publicamente sob a Licença MIT para consulta educacional.

# Título: Aplicação do ridge plot à capacidade térmica isocórica de líquidos simples

# 1) Objetivo:
Este notebook implementa numericamente a capacidade térmica isocórica
$C_V(T;\sigma)$ obtida a partir do modelo de líquidos simples como sistemas
desordenados. Além disso, investiga como as curvas $C_V(T)$ evoluem com o strength
da desordem $\sigma$ por meio de um ridge plot.

# 2) Descrição:
A capacidade térmica isocórica de um líquido de composição molecular simples
(baixo peso e número atômico) e com viscosidade desprezível é calculada por meio
de métodos da teoria estatística de campos aplicada a sistemas com desordem.
Neste contexto, foi obtida uma fórmula para a capacidade térmica isocórica
em termos da temperatura normalizada do líquido. Sendo o líquido simples
modelado como um sistema desordenado, o strength da desordem também aparece
como um parâmetro livre na fórmula da capacidade térmica.

Portanto, através do método do ridge plot, este projeto irá calcular
diferentes curvas teóricas da capacidade térmica, explorando sua evolução
com o strength da desordem.

# 3) Ambientação teórica

# 3.a) O líquido como sistema desordenado

Consideramos o seguinte funcional de partição para modelar o líquido simples
sem viscosidade:

$$Z(h) = \int [d\varphi] \exp \left(-\int_0^{\hbar \beta} d\tau \int_V d\mathbf{x} [ \frac{1}{2}\varphi(\tau,\mathbf{x}) ( -\frac{1}{u^2}\frac{\partial^2}{\partial \tau^2} -\Delta ) \varphi(\tau,\mathbf{x}) - h(\tau,\mathbf{x})\varphi(\tau,\mathbf{x}) ] \right)$$

sendo $h(\tau, \x)$ um campo de desordem que modela estatisticamente os processos 
interatômicos do líquido. Tais processos estão em uma escala de energia para além 
da energia característica dos processos estudados na hidrodinâmica clássica.
Consideramos a seguinte correlação para o campo de desordem $h(\tau, x)$:

$$\mathbb{E}\left[h(\tau, \x) h(\tau', x') \right] =  \sigma^2\delta (\tau-\tau')\delta(\x - \x')$$

no qual $\sigma$ é uma quantidade real estatística que modela a intensidade da 
desordem no sistema. 

# 3.b) A capacidade térmica e o método da zeta-distribucional

A capacidade térmica a volume constante do líquido é dada por:

$$C_V(h)  = \beta^2 \frac{\partial^2}{\partial \beta^2} \ln{Z(h)}.$$

Para um descrição efetiva, devemos calcular o $C_v(h)$ medido sobre todas 
as possíveis configurações da desordem, ou seja, 

$$ \mathbb{E}\left[C_V(h) \right] =  \beta^2 \frac{\partial^2}{\partial \beta^2} \mathbb{E}\left[ \ln{Z(h)} \right]$$

Para o cálculo da média $\mathbb{E}[\ln{Z(h)}]$ vamos recorrer ao método da função 
zeta distribucional. A derivação deste método pode ser vista com mais detalhes em
https://doi.org/10.1142/S0217751X1650144X . Entretanto, por simplicidade, vamos 
recorrer a um modo mais simples para verificar a fórmula final de $\mathbb{E}[\ln{Z(h)}]$.

A função integral exponecial $\mathrm{Ei}(-ax)$, onde $a>0$ é dada por:

$$ \mathrm{Ei}(-ax) = \int_{ax}^\infty \frac{e^{-t}}{t} dt$$

e sua expansão em série é:

$$\mathrm{Ei}(-ax) = \gamma + \ln{a|x|} + \sum_{k=1}^\infty \frac{(-1)^{k}}{k!k} (ax)^k $$

onde $\gamma$ é a constante de Euler-Mascheroni. Considerando $ax \rightarrow \frac{b}{Z_0} Z(h)$, 
onde $b>0$ e $Z_0$ é a função de partição calculada quando a desordem é nula, teremos a seguinte 
representação

$$ \mathrm{Ei}\left(-b\frac{Z(h)}{Z_0} \right) = \int_{b\frac{Z(h)}{Z_0}}^\infty \frac{e^{-t}}{t} dt $$

e

$$  \mathrm{Ei}\left(-b\frac{Z(h)}{Z_0} \right) = +\gamma - \ln{Z_0} +\ln{b} +\ln{Z(h)} + \sum_{k=1}^\infty \frac{(-b)^{k}}{k!k} \frac{[Z(h)]^k}{Z_0^k} $$

onde 

$$ \mathrm{Ei}\left(-b\frac{Z(h)}{Z_0} \right) \leq \frac{e^{-b}}{b} $$

Aplicando-se a média na equação de $\mathrm{Ei}\left(-b\frac{Z(h)}{Z_0} \right)$, teremos:

$$\mathbb{E}[\ln Z(h)] = \ln Z_0 + \sum_{k=1}^{\infty} \frac{(-1)^{k+1}b^k}{k!k}\mathbb{E}\left[\frac{Z(h)^k}{Z_0^k}\right] - \ln b - \gamma - R(b).$$

sendo

$$R(b) =  -\mathbb{E}\left[\int_{\frac{b Z(h)}{Z_0}}^\infty \frac{e^{-t}}{t} dt \right].$$

Para valores suficientes grandes de $b>1$, a contribuição de $R(b)$ se torna desprezível. 
Por meio da regularização de determinantes, podemos representar $Z_0$ e $\mathbb{E}[Z^k(h)]$
em termos das seguintes exponenciais:
$$Z(0) = \exp\{ V\alpha\int_0^1 dy y^{2}\ln[2\sinh(\frac{T_D}{2T}y)]\}$$

e
$$\mathbb{E}[\frac{Z(h)^k}{Z_0^k}] = \exp\{V\alpha \int_0^1 dy y [\sqrt{y^2 + k(\frac{\sigma}{(2\pi^2 \alpha)^{1/3}})^2} - y]\ln[2\sinh(\frac{T_D}{2T}y)]\}$$

onde $V$ é o volume total do sistema, $\alpha$ é o inverso do volume característico (comprimento de onda de Debye ao cubo) 
e $T_D$ é a temperatura de Debye. 

# 4) Metodologia

## Metodologia

A partir da expressão teórica para a capacidade térmica isocórica $C_V(T)$, serão calculadas 
numericamente diferentes curvas $C_V(T; \sigma)$ para valores distintos da intensidade 
(strength) da desordem $\sigma$. Assim, as curvas obtidas serão organizadas e 
representadas por meio de um ridge plot, permitindo visualizar o do comportamento da
capacidade térmica de acordo com a mudança da intensidade da desordem.
