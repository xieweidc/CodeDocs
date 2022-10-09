---
title: Newton-Cotes Quadrature formula
---

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
