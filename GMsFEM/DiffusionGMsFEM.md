---
title:	GMsFEM solve Diffusion Equation
---

# Introduction

​	Different to the multiscale finite element method(MsFEM), 
generalized multiscale finite element method(GMsFEM) 
can have more than one basis in each coarse node. 

# PDE model

  In this document, like the MsFEM, we still use a 2D problem to 
illustrate the basic ideas of the method. 
1D and 3D also be concerned. The 2D problem is:

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

# Calculation

The computation step is

- Variation
- Discretization
- Space
- Assemble
- Boundary condition
- Solve

  The only different thing is still how to construct multiscale space.
	
## 2 Discretization
   
  For GMsFEM, we need coarse mesh and fine mesh. Coarse mesh is used to 
approximate the numerical solution, while fine mesh used to construct 
the MsFEM basis. Let's give some notation to illustrate this method better.

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
| DNC | Number of coarse cell in the support of one multiscale basis |
| $\mathcal{T}_H$ | coarse mesh |
| $K_H$ | element in coarse mesh|
| $\mathcal{T}_h$ | fine mesh in a coarse element |
| $K_h$ | element in fine mesh |

  Note that the mesh is not restricted to rectangular elements, 
we can use triangle, tetrahedron, hexahedron either. 
We gave the local degree of freedom number, the edge of the serial number.

<div align=center>
<img src="../pics/fealpy_mesh_display.png" alt="参考单元网格节点编号" width="40%"/>
<img src="../pics/RectangleMesh.png" alt="笛卡尔坐标下的单元网格" width="40%"/>
</div>

## 3 Space

  In constructing the GMsFEM basis, divides it into two stages: 
offline and online.

1. Offline computation:
   - 1.1 Coarse grid generation;
   - 1.2 Construction of snapshot space;
   - 1.3 Performing dimension reduction of the snapshot space.
2. Online computation:
   - 2.1 Performing dimension reduction of the offline space;
   - 2.2 Times the partition unity;
   - 2.3 Iterative solvers, if needed.

  As for the actual computation, we need to compute 
partition unity $\chi$, snapshot basis function $\psi^{snap}$, 
offline basis function $\psi^{off}$, 
online basis function $\psi^{on}$ if needed.


### 3.1 Partition Unity

  We can take the original MsFEM as partition unity, the computation of 
MsFEM basis already examplified in this folder. Using $\chi$ as the partition unity notation, to different with the multiscale basis.

### 3.2 Offline computation

#### 3.2.1 Snapshot Space

  First, we using the same grid to compute the snapshot space,
when compute the partition unity. 
There are two options to construct the snapshot space. 
The first choice is to solve a set of local problems with 
boundary conditions in each coarse neighborhood.
The second choice is to use fine-grid nodal bases 
in each coarse region $K_H$.

  In here, we use $P_{K_{H,i}}^{snap}$ to denote the multiscale basis value 
in each fine node. We talk the multiscale basis for one cell,
so abbreviate the subscript $K_{H,i}$ is a convient choice.

**Choice 1**

  Solve following 

**Choice 2**

  In this case, $P^{snap}$ is very simple.

$$
P^{snap} = I
$$

#### 3.2.2 Offline Space

  Then, we need a spectrum decomposition to reduce the dimension
of snapshot space.