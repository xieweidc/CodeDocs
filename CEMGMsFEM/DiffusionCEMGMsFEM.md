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

| Notation |           Meaning            |
| :----:   | :----------------------------- |
|  NXC     | Number of segments in the x direction for coarse mesh|
|  NYC     | Number of segments in the y direction for coarse mesh|
|  NXF     |  Number of segments in the x direction for fine mesh |
|  NYF     |  Number of segments in the y direction for fine mesh |
|  DNN     | Number of fine node in $D_i$ |
|  DNA     | Number of auxiliary multiscale basis in $D_i$ |

## Multiscale space

**To be consistent with the notation of other documents,** 
**we still using $\phi$ to denote the multiscale basis,** 
**rather than $\psi$ in the CEM-GMsFEM paper.**
Besides, $\varphi$ is stand for the basis in fine grid.


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

Using Lagrange multiplied, it comes to

$$
\begin{cases}
\begin{aligned}
a(\phi_{j,ms}^i, p) + s(p,\lambda) &= 0,
&\forall p\in V_0(K_{i,m}), \\
s(\phi_{j,ms}^i-\psi_j^i, q) &= 0,
&\forall q\in V_{aux}(K_{i,m}).
\end{aligned}
\end{cases}
$$

Considering $\phi_{j,ms}$ can be a linear combination of $\varphi$ in 
$D_i$, while $\lambda$ is the linear combination of $\psi_j^i$. 
Denote $M_h$ be the matrix such that $M_{h,ij}=s(q_j,q_i)$, 
$A_h$ is the stiffmatrix, i.e. $A_{h,ij}=a_(q_j,q_i)$. 
$A_h^i$ and $M_h^i$ be the restriction of $A_h$ and $M_h$ on $K_{i,m}$ 
respectively. $P^i$ is the matrix that includes all the discrete 
auxiliary basis in space $V_{aux}(K_{i,m})$.
The matrix form as follows:

$$
\begin{bmatrix}
A_h^i & M_h^i P^i \\
(M_h^i P^i)^T & 0
\end{bmatrix}
\begin{bmatrix}
\phi_h^i \\ \lambda_h^i
\end{bmatrix} = 
\begin{bmatrix}
 0 \\ I
\end{bmatrix}
$$

where $I$ is a unit matrix, which dimension is (DNN, DNN). 
In this case, we have to point out that the order of 
$\psi_{j,h}^i$ is consistent with 
the order of $\phi_{j,ms}^i$.

In addition, using the previously mentioned marks, let's talk about 
the dimension of the above matrix.

$$
\begin{aligned}
A_h^i &: \text{(DNN, DNN)} \\
M_h^i &: \text{(DNN, DNN)} \\
P^i &: \text{(DNN, DNA)} \\
\phi_h^i &: \text{(DNN, DNA)} \\
\lambda_h^i &: \text{(DNA, DNA)}
\end{aligned}
$$

# Appendix

Defined $V=H_0^1.$

## Rewrite of Minimizing problem

The original form as follows:

$$
\phi_{j,ms} =  \mathrm{argmin}
\{a(\phi,\phi) | \psi \in V_0(K_{i,m}),
\phi ~\text{is}~ \psi_j^i\text{-orthogonal} \}
$$

By using Lagrange Multiplier, we can get 

$$
(\phi_{j,ms}, \lambda) =  \mathrm{argmin}
\{a(\phi,\phi) + 2s(\phi-\psi_j^{(i)}, \lambda) 
| \psi \in V_0(K_{i,m}),
\phi ~\text{is}~ \psi_j^i\text{-orthogonal}, 
\lambda \in V_{aux} \}
$$

Define

$$
\mathcal{L}(\phi, \lambda) = 
a(\phi, \phi)+2s(\phi-\psi_k^{(i)}, \lambda),
$$

Then the original problem has transferd to, looking for a 
$(\phi, \lambda) \in V_0(K_{i,m}) \times V_{aux}(K_{i,m})$, 
such that $\mathcal{L}(\phi, \lambda)$ is minimizing.

The partial of $\mathcal{L}(\phi, \lambda)$ with respect 
to $\phi, \lambda$.

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

Then it comes to

$$
\begin{cases}
\begin{aligned}
a(\phi_{j,ms}^i, p) + s(p,\lambda) &= 0,
&\forall p\in V_0(K_{i,m}), \\
s(\phi_{j,ms}^i-\psi_j^i, q) &= 0,
&\forall q\in V_{aux}(K_{i,m}).
\end{aligned}
\end{cases}
$$

## Orthogonal of two Spaces

In the following, we will discuss why did the two space 
$V_{glo}$ and $\tilde{V}$ are orthogonal.
The operator $\pi_i$ is given by

$$
\pi_i(u) = \sum_{j=1}^{l_i} 
\frac{s_i(u,\psi_j^{(i)})}{s_i(\psi_j^{(i)},\psi_j^{(i)})}
\psi_j^{(i)}, \quad \forall u \in V.
$$

The null space of the projection $\pi$, namely, $\tilde{V}=\{v\in V| \pi(v)=0\}$. 
Then for any $\phi_j^{(i)}\in V_{glo}$, we have 

$$
a(\phi_j^{(i)}, v)=0, \quad \forall v\in \tilde{V}.
\tag{1.1}
$$

As discussed before, 

$$
a(\phi_j^{(i)}, v) = -s(v, \lambda)
$$

where $\lambda \in V_{aux}$, while $v \in \tilde{V}$, so 
we can conclude that (1.1).
