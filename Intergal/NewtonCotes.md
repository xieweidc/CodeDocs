---
title: Newton-Cotes Quadrature formula
---

# Introduction

The Newton-Leibniz formulation can used to compute definite integral.

$$
I = \int_a^b f(x) \mathrm{d}x = F(b) - F(a)
$$

where $F$ is the primitive function of $f$.

But for practice, the above is hard to achive, because of 
we can't express the $F$ explicit.
For this case, using another quadrature method instead, 
numerical integration.

- Trapezoid formula, 

$$
I \approx \frac{b-a}{2} (f(a) + f(b)).
$$

- Mid-rectangle formula,

$$
I \approx (b-a) f(\frac{b-a}{2}).
$$

- Simpson formula

$$
I \approx \frac{b-a}{6} (f(a) + 4f(\frac){b-a}{2}+f(b)).
$$

# Appendix

​	Divide the interval $[a, b]$ into $n$ part.

$$
x_k = a+kh, \quad k=0,1,2, \cdots, n.
$$

where $h = \frac{b-a}{n}.$ The n-order Newton-Cotes quadrature formula is 

$$
\int_a^b f(x) ~ \mathrm{d}x = 
(b-a) \sum_{k=0}^n c_k^{(n)} f(x_k),
$$

where

$$
c_k^{(n)} = \frac{(-1)^{n-k}}{n\cdot k!(n-k)!}
\int_0^n \prod_{i=0 \atop i\ne k}^n (t-i) 
~\mathrm{d}t.
$$

​	Obviously,  $c_k^{(n)}$ is not related to integral interval.

<center class="half">
<img src="../pics/NewtonCotesCoefficient.png"  width="80%"/>
</center>