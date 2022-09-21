---
title:  MsFEM solve Linear Poisson
---

# Introduction

​	Multiscale finite element method(MsFEM) is proposed by Hou and Wu in 1997. The original parper could be found in, A Multiscale Finite Element Method for Elliptic Problems in Composite Materials and Porous Media, Vol 134, 169-189.

​	In this document, we will use a 2D problem to illustrate the basic ideas of the method. 1D and 3D also be concerned. The 2D problem is:

$$
\begin{cases}
\begin{aligned}
L_{\epsilon}u =
-\nabla \cdot 
\left(
\kappa
(\frac{\boldsymbol{x}}{\epsilon})
\nabla u
\right)
&= f \quad in \quad \Omega \\
\quad u &= g \quad on\ \ \partial\Omega
\end{aligned}
\end{cases}
$$

​	When we need to solve a heterogeneity problem by Finite Element Method on a fine scale, probability, we get a very large algebraic equations system. For example, if we assume the area that we need to solve is $[0, 1] \times [0, 1]$, and using $1024 \times 1024$ rectangle to divide that area. In that case, we store it in a csr matrix form, the memory is 

$$
2^3 \times 2^2 \times 2^2
\times 2^{10} \times 2^{10} Byte
= 2^{27} Byte
= 2^7 MB.
$$

​	Although it's not too big for today's computers, but you should note that it's just a csr matrix and not too big a partition in engineering practice. If we want to compute it at the same time, or other complicated operations,  it will take a long time.

​	To reduce the matrix size, they put forward a generalized fem with a new name, multiscale finite element method. 

# PDE model

Follow the previous example, 

$$
\begin{cases}
\begin{aligned}
L_{\epsilon}u =
-\nabla \cdot 
\left(
\kappa
(\frac{\boldsymbol{x}}{\epsilon})
\nabla u
\right)
&= f \quad in \quad \Omega \\
\quad u &= g \quad on\ \ \partial\Omega
\end{aligned}
\end{cases}
$$

where $\Omega = [0, 1] \times [0, 1].$

# Calculation

The computation step is

- Variation
- Discretization
- Space
- Assemble
- Boundary condition
- Solve

	The computation step is same as the FEM, the only different thing is space.

## 1 Variation

Like the linear FEM, given a subspace of $H^1(\Omega),$  named $V_H.$ Then

$$
-\int_{\Omega}
\left(
\nabla \cdot 
\left(
a(\frac{\boldsymbol{x}}{\epsilon})
\nabla u 
\right) 
\right) 
\cdot v ~\mathrm{d} \boldsymbol{x} 
= 
\int_{\Omega} f \cdot v 
~\mathrm{d} \boldsymbol{x}, 
\quad \forall v \in V_{H,0},
$$

where $V_{H,0} = V_H \cap H_0^1.$  Using green formulation, 

$$
\int_{\Omega} 
a(\frac{\boldsymbol{x}}{\epsilon}) 
\nabla u \cdot \nabla v
~\mathrm{d} \boldsymbol{x} 
= 
\int_{\Omega} f \cdot v
~\mathrm{d} \boldsymbol{x}
$$

​	We will talk about the space $V_H,$ which composed by multiscale finite element basis, at subsection 4.

## 2 Discretization

​	For MsFEM, we need coarse mesh and fine mesh. Coarse mesh is used to approximate the numerical solution, while fine mesh used to construct the MsFEM basis. Let's give some notation to illustrate this method better.
| Notation |           Meaning            |
| :----:   | :----------------------------- |
|  NXC     | Number of segments in the x direction for coarse mesh|
|  NYC     | Number of segments in the y direction for coarse mesh|
|  NCC     |       Number of elements in coarse mesh       |
|  NNC     |       Number of nodes in coarse mesh       |
| LDC      | Number of local basis functions of cell in coarse mesh |
|  NXF     |  Number of segments in the x direction for fine mesh |
|  NYF     |  Number of segments in the y direction for fine mesh |
|  NCF     |       Number of elements in fine mesh       |
|  NNF     |       Number of nodes in fine mesh       |
| LDF 	   | Number of local basis functions of cell in fine mesh |
| $\mathcal{T}_H$ | coarse mesh |
| $K_H$ | element in coarse mesh|
| $\mathcal{T}_h$ | fine mesh in a coarse element |
| $K_h$ | element in fine mesh |

​	Note that the mesh is not restricted to rectangular elements, we can use triangle, tetrahedron, hexahedron either. We gave the local degree of freedom number, the edge of the serial number.

<div align=center>
<img src="../pics/fealpy_mesh_display.png" alt="参考单元网格节点编号" width="40%"/>
<img src="../pics/RectangleMesh.png" alt="笛卡尔坐标下的单元网格" width="40%"/>
</div>

## 3 Space

### 3.1 Multiscale finite element basis

​	Maybe you are interested in how the MsFEM basis looks like, I give three pictures in the next. From left to right, is linear FEM basis, linear boundary MsFEM basis, oscillating boundary MsFEM basis.

<center class="half">
<img src="../pics/P1_FEM_basis.png"  width="30%"/>
<img src="../pics/LinearBC_MsFEM_basis.png"  width="30%"/>
<img src="../pics/OscillationBC_MsFEM_basis.png"  width="30%"/>
</center>

​	MsFEM can capture the local scale because of the basis, but unfortunately, the basis function doesn't have a explicit expression. So it's the key to construct the basis. 

​	The basis satisfies 

$$
\begin{cases}
L_{\epsilon} \phi_i = 0 ~\quad\ in\quad K_H\\
\quad \phi_i = \mu_i \quad on\ \partial K_H
\end{cases}
\qquad i=0,1,2,3
$$

where $\mu_i$ is artificial boundary condition.

### 4.2 子问题构造

在粗网格单元 $K_H$ 上，我们通过微分算子 $L_{\epsilon}$ 构造基函数来抓住细尺度的信息。
为了确保基函数的全局连续性，需给它设定一个边界条件 $\mu_i$，
由此，我们可以得到原问题的子问题：



$\mu_i$ 为粗网格节点 $x_i$ 处对应的多尺度基函数在 $\partial K_H$ 处的定义。
值得指出的是，一个粗网格单元内子问题个数为 $K_H$ 节点个数
（四边形网格内就有四个子问题）。
对于子问题的边界条件，我们有两种选择，

+ 线性边界条件
+ 振荡边界条件

$\mu_i$ 将在下一小节具体讨论，

### 4.3 子问题边界选择

+ **线性边界条件**

对于粗网格单元的第 $j = 0,1,2,3$ 个节点，再结合

$$
\phi_i(\boldsymbol{x}_j) = 
\delta_{ij}
$$

其多尺度基函数在边界处的函数表达式如下：

$$
\begin{cases}
\mu_0(x,y) = 
\cfrac{(x-x_3)(y-y_3)}{(x_0-x_3)(y_0-y_3)}
\\
\mu_1(x,y) = 
\cfrac{(x-x_3)(y-y_0)}{(x_0-x_3)(y_3-y_0)} 
\\
\mu_2(x,y) = 
\cfrac{(x-x_0)(y-y_3)}{(x_3-x_0)(y_0-y_3)} 
\\
\mu_3(x,y) = 
\cfrac{(x-x_0)(y-y_0)}{(x_3-x_0)(y_3-y_0)}
\end{cases}
$$

其中 $\mu_i$ 的编号如下图所示：

<div align=center>
<img src="../pics/RectangleMesh.png" alt="笛卡尔坐标下的单元网格" width="40%"/>
</div>

+ **振荡边界条件**

对于基函数的振荡边界条件，为了较好的保留尺度信息。
我们将原问题去除掉正交于 $\partial K$ 方向的偏导数，
只保留正切于 $\partial K$ 方向的偏导数，得到振荡边界条件。
对于上文中给出的特殊单元，表达式如下：

$$
\begin{aligned}
\frac{\mathrm{d}}{\mathrm{d}x} 
\left( 
a(\frac{\boldsymbol{x}}{\epsilon}) 
\frac{\mathrm{d} \mu_i}{x}
\right) 
= 0 \quad 
\mu_i ~on~\Gamma_0~or~\Gamma_3  
\\
\frac{\mathrm{d}}{\mathrm{d}y} 
\left( 
a(\frac{\boldsymbol{x}}{\epsilon}) 
\frac{\mathrm{d} \mu_i}{y}
\right) 
= 0 \quad 
\mu_i ~on~\Gamma_1~or~\Gamma_2
\end{aligned}
$$

再结合

$$
\phi_i(\boldsymbol{x}_j) = 
\delta_{ij}
$$

这样就能求出基函数在边界处的表达式:

$$
\begin{cases}
\mu_0(x,y) =
\cfrac{\int_{x}^{x_3} \cfrac{1}{a} ~\mathrm{d}x}
{\int_{x_0}^{x_3} \cfrac{1}{a} ~\mathrm{d}x}
\cfrac{\int_{y}^{y_3} \cfrac{1}{a} ~\mathrm{d}y}
{\int_{y_0}^{y_3} \cfrac{1}{a} ~\mathrm{d}y}
\\
\mu_1(x,y) =
\cfrac{\int_{x}^{x_3} \cfrac{1}{a} ~\mathrm{d}x}
{\int_{x_0}^{x_3} \cfrac{1}{a} ~\mathrm{d}x}
\cfrac{\int_{y_0}^{y} \cfrac{1}{a} ~\mathrm{d}y}
{\int_{y_0}^{y_3} \cfrac{1}{a} ~\mathrm{d}y}
\\
\mu_2(x,y) =
\cfrac{\int_{x_0}^{x} \cfrac{1}{a} ~\mathrm{d}x}
{\int_{x_0}^{x_3} \cfrac{1}{a} ~\mathrm{d}x}
\cfrac{\int_{y}^{y_3} \cfrac{1}{a} ~\mathrm{d}y}
{\int_{y_0}^{y_3} \cfrac{1}{a} ~\mathrm{d}y}
\\
\mu_3(x,y) =
\cfrac{\int_{x_0}^{x} \cfrac{1}{a} ~\mathrm{d}x}
{\int_{x_0}^{x_3} \cfrac{1}{a} ~\mathrm{d}x}
\cfrac{\int_{y_0}^{y} \cfrac{1}{a} ~\mathrm{d}y}
{\int_{y_0}^{y_3} \cfrac{1}{a} ~\mathrm{d}y}
\end{cases}
$$

+ **两种边界条件的优缺点**

|   边界   |         优点         |         缺点         |
| :------: | :------------------: | :------------------: |
| 线性边界 |       计算简单       | 放弃边界处细尺度信息 |
| 振荡边界 | 保留边界处细尺度信息 |       计算复杂       |


### 4.4 子问题变分

由于子问题中的多尺度基函数构造的特殊性，无法获得确切的数学表达式。
对此，我们采用线性元来数值逼近多尺度基函数。

给定一个试探函数空间 $V_h$，

$$
-\int_{K_H}
\left(
\nabla \cdot 
\left(
a(\frac{\boldsymbol{x}}{\epsilon}) 
\nabla \phi
\right) 
\right) 
\cdot v ~\mathrm{d} \boldsymbol{x} = 
\int_{K_H} 0\cdot v
~\mathrm{d}\boldsymbol{x}, 
\quad \forall v \in V_{h,0},
$$

分部积分可得

$$
\int_{K_H} 
a(\frac{\boldsymbol{x}}{\epsilon}) 
\nabla \phi \cdot \nabla v
~\mathrm{d} \boldsymbol x 
= 0
$$

## 5 参数说明

对于这两套网格，我们进行必要的参数说明



那么，在求解域 $\Omega$ 内，所有粗网格节点组成的有限元空间为

$$
V_H = 
\begin{bmatrix}
\phi_0, \phi_1, \cdots, \phi_{nnc-1}
\end{bmatrix}
$$

一个粗网格单元 $K_H$ 内，其有限元空间为：

$$
V_h = 
\begin{bmatrix}
\varphi_0, \varphi_1, \cdots, \varphi_{nnf-1}
\end{bmatrix}
$$

## 6 刚度矩阵以及载荷向量的组装

### 6.1子问题求解

用线性元求解子问题，我们可以得到：

$$
\phi_i = 
\sum\limits_{l=0}^{nnf-1} q_{i,l} \varphi_l\ 
, \quad  i=0,1,2,3
$$

其中的 $d_{il}$ 表示第 $i$ 个多尺度基函数由细网格内的 $l$ 个线性元线性表达的系数，
写成矩阵形式为：

$$
\begin{bmatrix}
\phi_0, & \phi_1, & \phi_2, & \phi_3
\end{bmatrix} = 
\begin{bmatrix} 
\varphi_0, & \varphi_1, & \cdots, & \varphi_{ncf-1}
\end{bmatrix}
\begin{bmatrix}
q_{0,0} & q_{1,0} & q_{2,0} & q_{3,0} \\
q_{0,1} & q_{1,1} & q_{2,1} & q_{3,1} \\
q_{0,2} & q_{1,2} & q_{2,2} & q_{3,2} \\
\vdots  & \vdots  & \vdots  & \vdots  \\
q_{0,nnf-1} & q_{1,nnf-1} & q_{2,nnf-1} & q_{3,nnf-1} \\
\end{bmatrix}
$$

其中 $\varphi_i$ 为细网格单元的基函数，我们记

$$
\boldsymbol{\Phi}_H = 
\begin{bmatrix}
\phi_0, & \phi_1, & \phi_2, & \phi_3
\end{bmatrix}, \quad 
\boldsymbol{\varPhi}_h = 
\begin{bmatrix} 
\varphi_0, & \varphi_1, & \cdots, & \varphi_{nnf-1}
\end{bmatrix}, \quad 
D_h = 
\begin{bmatrix}
q_{0,0} & q_{1,0} & q_{2,0} & q_{3,0} \\
q_{0,1} & q_{1,1} & q_{2,1} & q_{3,1} \\
q_{0,2} & q_{1,2} & q_{2,2} & q_{3,2} \\
\vdots  & \vdots  & \vdots  & \vdots  \\
q_{0,nnf-1} & q_{1,nnf-1} & q_{2,nnf-1} & q_{3,nnf-1} \\
\end{bmatrix}
$$

$\Phi_H, \varPhi_h$ 都是行向量，则有：

$$
\boldsymbol{\Phi}_H =  
\boldsymbol{\varPhi}_h Q
$$

### 6.2 单元刚度矩阵与载荷向量组装

$\Phi_H$ 的梯度记为

$$
\nabla \Phi_H = 
\begin{bmatrix}
\nabla \phi_0, & \nabla \phi_1, & 
\nabla \phi_2, & \nabla \phi_3
\end{bmatrix}
$$

这里标量函数的梯度默认是列向量的形式，那么 $\nabla \Phi_H$ 
实际上是一个 $d \times 4$ 的矩阵函数（d为求解域维数，这里 $d=2$ ) . 

类似的，$\varPhi_h$ 的梯度记为

$$
\nabla \varPhi_h = 
\begin{bmatrix}
\nabla \phi_0, & \nabla \phi_1, & 
\cdots, & \nabla \phi_{nnf-1}
\end{bmatrix}
$$

$\nabla \varPhi_h$ 是一个 $d \times nnf$ 的矩阵函数。

我们记子问题在粗网格单元 $K_H$ 上的总刚度矩阵 $A_h$，
总载荷向量 $b_h$ 分别为

$$
A_h = 
\int_{K_H} a \cdot
\left(
\nabla \varPhi_h
\right)^T
\nabla \varPhi_h
~\mathrm{d}\boldsymbol{x},
\qquad
b_h = 
\int_{K_H} f \cdot
\varPhi_h^T
~\mathrm{d}\boldsymbol{x}
$$


##### 多尺度单元刚度矩阵

$$
\begin{aligned}
A_{K_H} =&~
\int_{K_H} a \cdot
\left( 
\nabla \Phi_H
\right)^T
\nabla \Phi_H
~\mathrm{d}\boldsymbol{x}  \\
=&~ \int_{K_H} a \cdot
\left( 
\nabla \varPhi_h Q
\right)^T
\nabla \varPhi_h
Q
~\mathrm{d}\boldsymbol{x} \\
=&~ Q^T A_h Q
\end{aligned}
$$

##### 多尺度单元载荷向量

$$
\begin{aligned}
b_{K_H} &= 
\int_{K_H} f \cdot \Phi_H^T
~\mathrm{d}\boldsymbol{x} \\
&= 
\int_{K_H} f \cdot 
\left(
\varPhi_h Q
\right)^T
~\mathrm{d}\boldsymbol{x} \\
&= Q^T b_h
\end{aligned}
$$

再将单元刚度矩阵与单元载荷向量按照索引存储到总刚度矩阵与总载荷向量中去，即可求解。

## 6 数值实验结果

这里我们分别给出两个案例，检验在真解是否振荡的两种情形下，
高次元 、线性边界多尺度有限元 、振荡边界多尺度有限元 这三种方法的性态。
值得注意的是，在两种边界条件下的多尺度有限元法，我们都选择细网格剖分段数为 
$ nxf = nyf = 8$.

$$
\begin{cases}
L_{\epsilon}u = f \quad in \quad \Omega \\
\quad u = g \quad on\ \ \partial\Omega
\end{cases} 
$$

其中 

$$
L_{\epsilon}u = 
-\nabla \cdot 
\left(
a(\frac{\boldsymbol{x}}{\epsilon_1})\nabla u
\right), \quad 
a(\frac{\boldsymbol{x}}{\epsilon_1}) = 
\frac{1}
{(2+P\sin(\frac{x}{\epsilon_1}))
(2+P\sin(\frac{y}{\epsilon_1}))}, \quad
u(\frac{\boldsymbol{x}}{\epsilon_2}) =
\sin(\frac{x}{\epsilon_2})
\sin(\frac{y}{\epsilon_2})
$$

#### 6.1 真解非振荡情形

$P = 1.5, \epsilon_1= \epsilon_2 = 1$，结果如下：

<table>
	<tr align="center">
		<td> NX = NY </td>
		<td colspan="2"> P4 元 </td>
        <td colspan="2"> 线性边界多尺度有限元 </td>
		<td colspan="2"> 振荡边界多尺度有限元 </td>
	</tr>
	<tr align='center'>
        <td> 2 </td>
        <td> 1.09252023e-04 </td>
        <td> </td>
        <td> 0.123393 </td>
        <td> </td>
        <td> 0.120459 </td>
        <td> </td>
    </tr>
    <tr align='center'>
        <td> 4 </td>
        <td>  3.52844059e-06 </td>
        <td> 4.95 </td>
        <td> 0.037659 </td>
        <td> 1.71 </td>
        <td> 0.064074 </td>
        <td> 0.91 </td>
    </tr>
    <tr align='center'>
        <td> 8 </td>
        <td> 1.06368504e-07 </td>
        <td>  5.05 </td>
        <td> 0.009496 </td>
        <td> 1.99 </td>
        <td> 0.015543 </td>
        <td> 2.04 </td>
    </tr>
    <tr align='center'>
        <td> 16 </td>
        <td> 3.30508628e-9 </td>
        <td> 5.01 </td>
        <td> 0.002391 </td>
        <td> 1.99 </td>
        <td> 0.003893 </td>
        <td> 2.00 </td>
    </tr>
    <tr align='center'>
        <td> 32 </td>
        <td> 1.03151488e-10 </td>
        <td> 5.00 </td>
        <td> 0.000600 </td>
        <td> 2.00 </td>
        <td> 0.000977 </td>
        <td> 1.99 </td>
    </tr>
    <tr align='center'>
        <td> 64 </td>
        <td> 3.22538390e-12 </td>
        <td> 5.00 </td>
        <td> 0.000150 </td>
        <td> 2.00 </td>
        <td> 0.000244 </td>
        <td> 2.00 </td>
    </tr>
</table>

#### 6.2 真解强振荡情形

$P = 1.8,\epsilon_1=0.008, \epsilon_2=0.08$，结果如下：

<table>
	<tr align="center">
		<td> NX = NY </td>
		<td colspan="2"> P8 元 </td>
        <td colspan="2"> 线性边界多尺度有限元 </td>
		<td colspan="2"> 振荡边界多尺度有限元 </td>
	</tr>
        <tr align='center'>
        <td> 2 </td>
        <td> 178.23570728 </td>
        <td> </td>
        <td> 3.40821887e+01 </td>
        <td> </td>
        <td> 3.18629635e+01  </td>
        <td> </td>
    </tr>
    <tr align='center'>
        <td> 4 </td>
        <td> 76.15654157 </td>
        <td>  1.23  </td>
        <td> 5.65486603e-01 </td>
        <td> 5.91 </td>
        <td> 6.35698253e-01 </td>
        <td> 5.65  </td>
    </tr>
    <tr align='center'>
        <td> 8 </td>
        <td>  35.76403924 </td>
        <td> 1.09 </td>
        <td> 5.36255946e-01  </td>
        <td>  0.08   </td>
        <td> 5.41033265e-01 </td>
        <td> 0.23  </td>
    </tr>
    <tr align='center'>
        <td> 16 </td>
        <td> 61.42495307 </td>
        <td> -0.78  </td>
        <td> 3.05249901e-01 </td>
        <td> 0.81 </td>
        <td>  3.06064551e-01 </td>
        <td> 0.82 </td>
    </tr>
    <tr align='center'>
        <td> 32 </td>
        <td>  3.54661894</td>
        <td> 4.11 </td>
        <td> 8.22935495e-02  </td>
        <td> 1.89  </td>
        <td> 7.86093119e-02 </td>
        <td> 1.96  </td>
    </tr>
    <tr align='center'>
        <td> 64 </td>
        <td> 0.90121616 </td>
        <td> 1.98 </td>
        <td> 3.07820424e-02 </td>
        <td> 1.42 </td>
        <td> 2.97493282e-02 </td>
        <td> 1.40 </td>
    </tr>
</table>

#### 6.3 结论

可以看出，对于非振荡的情形，用高阶有限元就能得到很好的结果；
而对于强振荡的情形，多尺度有限元更加具有优势。