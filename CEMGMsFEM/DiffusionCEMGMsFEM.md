---
title:	Constraint Energy Minimizing Generalized Multiscale Finite Element Method to solve Diffusion Equation
---

# PDE Model

The diffusion equation as follows:

$$
\begin{cases}
\begin{aligned}
-\mathrm{div}(\kappa{\nabla}u) &= f,
&\text{in} ~\Omega,\\
u &= g, &\text{on} ~\partial\Omega.
\end{aligned}
\end{cases}
$$

Construct CEM-GMsFEM space are mainly in two steps:

- auxilary mutliscale space, $V_{aux}$.
- multiscale space, $V_{ms}$.

# Simulation

**To be consistent with the notation of other documents,** 
**we still using $\phi$ to denote the multiscale basis,** 
**rather than $\psi$ in the CEM-GMsFEM paper.**

### Auxiliary multiscale space

Using $K_i$ denote $ith$ coarse cell,

$$
a_i(\psi_j^{(i)}, w) = 
\lambda_j^{(i)} s_i(\psi_j^{(i)}, w),
\quad \forall w \in V(K_i)
$$

where

$$
a_i(v, w)=\int_{K_i} \kappa
{\nabla}v{\cdot}{\nabla}w, \quad
s_i(v,w)=\int_{K_i} \tilde{\kappa} vw, 
\quad
\tilde{\kappa}=\sum_{j=1}^{N_c} \kappa
|\nabla \chi_j^{ms}|^2.
$$

We will use the first $l_i$ eigenfunctions to construct our local auxiliary 
multiscale space $V_{aux}^{(i)}$. where 
$V_{aux}^{(i)}=\mathrm{span}\{\psi_j^{(i)}|j\le l_i\}$.
The global auxiliary multiscale space $V_{aux}$ is the sum of these local 
auxiliary multiscale space, namely 
$V_{aux}=\oplus_{i=1}^N V_{aux}^{(i)}$.

### Multiscale space

Next define the multiscale basis function $\phi_{j,ms}^{(i)}$ by：

$$
\phi_{j,ms} =  \mathrm{argmin}
\{a(\phi,\phi) | \psi \in V_0(K_{i,m}),
\phi ~\text{is}~ \psi_j^i-\text{orthogonal} \}
$$

By using Lagrange Multiplier, we can get 

$$
(\phi_{j,ms}, \lambda) =  \mathrm{argmin}
\{a(\phi,\phi) + 2s(\phi-\psi_j^{(i)}, \lambda) 
| \psi \in V_0(K_{i,m}),
\phi ~\text{is}~ \psi_j^i-\text{orthogonal}, 
\lambda \in V_{aux} \}
$$


The matrix form as follows:

$$
\begin{bmatrix}
A_h^i & M_h^i P^i \\
(M_h^i P^i)^T & 0
\end{bmatrix}
\begin{bmatrix}
\psi_h^i \\ \lambda_h^i
\end{bmatrix} = 
\begin{bmatrix}
 0 \\ I_i
\end{bmatrix}
$$

其中 $A_h^i, M_h^i$ 为细网格线性基函数在 $K_{i,m}$ 上的刚度，质量矩阵。
$P^i$ 为辅助空间的离散形式。

# Appendix

Defined $V=H_0^1.$

## Orthogonal of two Spaces

The operator $\pi_i$ is given by

$$
\pi_i(u) = \sum_{j=1}^{l_i} 
\frac{s_i(u,\psi_j^{(i)})}{s_i(\psi_j^{(i)},\psi_j^{(i)})}
\psi_j^{(i)}, \quad \forall u \in V.
$$

The null space of the projection $\pi$, namely, $\tilde{V}=\{v\in V| \pi(v)=0\}$

Then for any $\phi_j^{(i)}\in V_{glo}$, we have 

$$
a(\phi_j^{(i)}, v)=0, \quad \forall v\in \tilde{V}.
$$

## Rewrite of Minimizing problem

$$
\mathcal{L}(\phi, \lambda) = 
a(\phi, \phi)+2s(\phi-\psi_k^{(i)}, \lambda),
$$

原问题就转化为找一 
$(\phi, \lambda) \in V_0(K_{i,m}) \times V_{aux}(K_{i,m}),$ 使得
$\mathcal{L}(\phi, \lambda)$ 函数值最小，
那么 $\mathcal{L}(\phi, \lambda)$ 分别对 $\phi, \lambda$ 求导。
$$
\begin{aligned}
\mathcal{L}_{\phi}(\phi, \lambda) &= 
\lim_{t\to 0} \frac{a(\phi+tp,\phi+tp)+
2s(\phi+tp-\psi_k^{(i)},\lambda)-a(\phi,\phi)-
2s(\phi-\psi_k^{(i)},\lambda)}{t} \\
&= 2a(\phi,p) + 2(p,\lambda) \\
\mathcal{L}_{\lambda}(\phi,\lambda) &= 
\lim_{t\to 0}
\frac{a(\phi,\phi)+
2s(\phi-\psi_k^{(i)},\lambda+tq)
-a(\phi,\phi)-
2s(\phi-\psi_k^{(i)},\lambda)}{t} \\
&= 2s(\phi-\psi_k^{(i)}, q)
\end{aligned}
$$

这样问题就变成了

$$
\begin{cases}
\begin{aligned}
a(\phi_{j,ms}^i, p) + s(p,\lambda) &= 0,
&\forall p\in V_0(K_{i,m}), \\
s(\phi_{j,ms}^i-\psi_j^i, q) &= 0,
&\forall q\in V_{aux}(K_{i,m}),
\end{aligned}
\end{cases}
$$