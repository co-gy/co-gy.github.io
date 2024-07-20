# 矩阵求解迭代法

# 材料

华东师范课件：[https://math.ecnu.edu.cn/~jypan/Teaching/NA/slides_04A_Jacobi_GS_SOR.pdf](https://math.ecnu.edu.cn/~jypan/Teaching/NA/slides_04A_Jacobi_GS_SOR.pdf)

博客：[https://leungyukshing.cn/archives/求解矩阵的三种迭代法.html](https://leungyukshing.cn/archives/%E6%B1%82%E8%A7%A3%E7%9F%A9%E9%98%B5%E7%9A%84%E4%B8%89%E7%A7%8D%E8%BF%AD%E4%BB%A3%E6%B3%95.html)

东北大学网课

# 概要

矩阵求解的迭代法是一种用于求解线性方程组 $Ax^{*} = b$ 的数值方法。通过迭代公式逐步逼近解，直到满足精度要求。用于大型稀疏矩阵。常见的迭代方法包括**雅可比迭代法（Jacobi）**、**高斯-赛德尔迭代法（G-S）**、**SOR迭代法**、**共轭梯度法**。迭代通式如下：

$$
x^{(k+1)} \leftarrow Bx^{(k)}+f
$$

其中：k表示迭代次数，$x^{(0)}$任取

迭代收敛定理：收敛的充要条件是迭代矩阵B的谱半径小于1，即：$\rho(B) < 1$

矩阵分裂：给定$A = M - N（M非奇异）$，称为矩阵$A$的一个分裂，原方程就有：$Mx = Nx + b$，迭代式就有：$x^{(k+1)} = M^{-1}Nx^{(k)} + M^{-1}b$，此时$B = M^{-1}N$，$f = M^{-1}b$。根据不同的分裂，就有不同的方法

给定：$A = D - L - U$

$$
D=\left(
	\begin{matrix}
	a_{11} & 0 & \cdots & 0\\
	0 & a_{22} & \cdots & 0\\
	\vdots & \vdots & \ddots & \vdots\\
	0 & 0 & \cdots & a_{nn} \\
    \end{matrix}
    \right)
$$

$$
L = -\left(
	\begin{matrix}
	0 & 0 & \cdots & 0\\
	a_{21} & 0 & \cdots & 0\\
	\vdots & \vdots & \ddots & \vdots\\
	a_{n1} & \cdots & a_{nn-1} & 0 \\
    \end{matrix}
    \right)
$$

$$
U = -\left(
	\begin{matrix}
	0 & a_{12} & \cdots & a_{1n}\\
	0 & 0 & \cdots & a_{2n}\\
	\vdots & \vdots & \ddots & \vdots\\
	0 & \cdots & 0 & 0 \\
    \end{matrix}
    \right)
$$

# 雅可比

令：

$$
M = D, N = L+U
$$

则：

$$
\begin{gather*}B = M^{-1}N = D^{-1}(L+U)\\ f = M^{-1}b = D^{-1}b\end{gather*}
$$

迭代式为：

$$
x^{(k+1)} \leftarrow D^{-1}(L+U)x^{(k)} + D^{-1}b
$$

$$
x_i^{(k+1)} \leftarrow \frac{1}{a_{ii}}(b_i - \sum_{j=1}^{i-1}a_{ij}x_j^{(k)}-\sum_{j=i+1}^n a_{ij}x_j^{(k)})
$$

将矩阵 $B$ 进行变换，得到如此形式的迭代式：

$$
B = D^{-1}(L+U) = D^{-1}[D-(D-L-U)] = I-D^{-1}A
$$

$$
x^{(k+1)} \leftarrow x^{(k)} + D^{-1}(b-Ax^{(k)})
$$

可以看到的是，$D^{-1}(b-Ax^{(k)})$是对$x^{(k)}$的一种所谓的“修正”，使得$x^{(k+1)}$比$x^{(k)}$更靠近真实解$$x^{*}$$。
观察迭代式可以发现
$$
D^{-1}(b-Ax^{(k)}) = D^{-1}(Ax^{*}-Ax^{(k)}) = D^{-1}A(x^{*}-x^{(k)})
$$
直观去理解的话，有点类似于是在真实误差$$x^{*}-x^{(k)}$$上乘了一个范数较小的矩阵$$D^{-1}A$$。然后通过迭代步步逼近$$x^{*}$$

但是其中有个逻辑还没有太理解的是，为什么矩阵$D^{-1}A$是如何作为修正项的系数矩阵的。

# 高斯-赛德尔

$$
M = D-L, N = U
$$

则：

$$
\begin{gather*}B = M^{-1}N = (D-L)^{-1}U\\ f = M^{-1}b = (D-L)^{-1}b\end{gather*}
$$

迭代式的分量形式如下：

$$
Dx^{(k+1)} \leftarrow Lx^{(k+1)} + Ux^k + b
$$

$$
x_i^{(k+1)} \leftarrow \frac{1}{a_{ii}}(b_i - \sum_{j=1}^{i-1}a_{ij}x_j^{(k+1)}-\sum_{j=i+1}^n a_{ij}x_j^{(k)})
$$

从迭代式的元素形式能够发现，G-S迭代充分利用了已经算出的新数据，从而获得更快的收敛速度

# SOR

基本思想：将 G-S 迭代法的 $x^{(k+1)}$ 与 $x^{(k)}$ 加权平均，以获得更好的近似解

$$
M = D-L, N = U
$$

迭代式：

$$
x^{(k+1)} = (D-\omega L^{-1})((1-\omega)D+\omega U)x^{(k)}+\omega(D-\omega L)^{-1}b
$$

分量形式：

$$
x_i^{(k+1)} \leftarrow (1-\omega)x_i^{(k)} + \frac{\omega}{a_{ii}}(b_i-\sum_{j=1}^{i-1}a_{ij}x_j^{(k+1)}-\sum_{j=i+1}^n a_{ij}x_j^{(k)})
$$

- $\omega=1$时，就是G-S
- $\omega < 1$，称为低松弛方法
- $\omega > 1$，称为超松弛方法，一般会有比较好的收敛效果，不同的系数矩阵$A$会有不同的最优的$\omega$
