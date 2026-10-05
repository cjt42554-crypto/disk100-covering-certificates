# 两实际单位锚点的非零 tau 引理与强范数域 two-face fee

## 0. 定理、适用范围、没有使用的预算

本文件固定

\[
 R=\frac{5661}{50000},\qquad \lambda=\frac7{50},\qquad
 t=\frac{1071883}{20000000}.
\]

两张实际正向非退化三角形为 \((p,c_0,c_1)\)、\((p,c_1,c_2)\)，其外接圆半径均不超过 \(R\)。对每张面的底边，存在实际单位锚点 \(q_i\)，满足

\[
 |q_i|=1,\quad |q_i-c_i|,|q_i-c_{i+1}|\le R.
\]

不假定两个端点距离相等；不假定锚点在 apex 的另一侧；不假定边界链凸、镜像对称、反射最优、圆心在 assigned-arc lens 内端，或最近-owner 条件。

使用 actual signed A：

\[
 A_i=-\tfrac12\det(c_{i+1}-c_i,q_i-c_i),\qquad
 S_2=\sum_{i=0}^1\left[\frac{|T_i|+A_i}{R^2}-4\lambda x_i^2+t\gamma_i\right],
\quad D=2C-S_2.
\]

其中 \(x_i=|c_{i+1}-c_i|/(2R)\)，\(\gamma_i\) 是实际 apex 角，\(C\) 是第 1 节的完整 apex/skew 支持。

**两面条件定理。** 如果另外

\[
 |c_1|^2\le1-R^2,\qquad |c_0|^2,|c_2|^2\ge\frac{83}{100},
\tag{0.1}
\]

则

\[
 \boxed{D>\frac9{400}=.0225.}
\tag{0.2}
\]

该结论已经处理非零 actual-anchor slack。它不使用中间弧 variance，也不把 J-A / tau 放入其他全局面积或 Q 预算。

用户给出的外 BC 范数下界 \(.911088875\) 的平方严格大于 \(83/100\)。完整 assigned arc 长度严格大于 \(.147\) 还直接给出更强的范数下界（第 9 节）。因此 (0.1) 是所给 highM 链的一个较弱几何接口。

定义目标 fee

\[
 f_{\rm req}=\frac{178523206195158271774904671901}
 {8011730250000000000000000000000}.
\]

精确有

\[
 \frac9{400}-f_{\rm req}
 =\frac{1740724429841728225095328099}
 {8011730250000000000000000000000}>\frac1{5000}.
\tag{0.3}
\]

所以在 (0.1) 下，任意实际中间 HALFarc 都满足

\[
 D+\frac{1221}{640R^2}\left(a_{\rm mid}-\frac\pi{33}\right)^2
 >f_{\rm req}+\frac1{5000}.
\tag{0.4}
\]

(0.4) 只是由非负 variance 得出的推论；variance 在证明中未支付任何项。

这不是全局 \(r_D(100)\) 认证：其他行、其他拓扑及 highM 全链依旧须独立完成。验证程序仅检查这里新引理的标量常数，没有运行旧三面或旧全行 checker。

---

## 1. 完整支持与单份 deficit 预算

令 \(\theta=\arcsin x\)。外接圆半径限制给出

\[
 0<\theta\le\pi/2,\qquad \theta\le\gamma\le\pi-\theta.
\]

定义

\[
 z=\gamma-\theta,\quad
 \sigma=x^2\cot(\gamma/2)-|T|/R^2\ge0,\quad
 \tau=x\cos\theta-A/R^2\ge0.
\]

面积上界来自固定底边和 apex 角的等腰最大面积；\(A\) 上界来自两个端点圆的 lens，不要求 \(A\ge0\)。单面 score 恒等于

\[
 x^2\cot(\gamma/2)+x\cos\theta-4\lambda x^2+t\gamma-\sigma-\tau.
\tag{1.1}
\]

固定 \(\theta\) 后，去掉 \(-\sigma-\tau\) 的部分关于 \(\gamma\) 凸。于是完整支持为两个 apex 端点分支之最大值：

\[
 F(\theta)=\sin\theta(1+2\cos\theta)-4\lambda\sin^2\theta+t\theta,
\]
\[
 G(\theta)=\sin\theta-4\lambda\sin^2\theta+t(\pi-\theta).
\]

\(G\le1+t\pi<6/5\)，而 \(F(5/6)>7/5\)。因此

\[
 C=\max_{0\le\theta\le\pi/2}F(\theta).
\tag{1.2}
\]

不使用全区间凹性：事实上 \(F''(\pi/2)=8\lambda-1>0\)。正确的一维归约是

\[
 F'=\cos\theta+2\cos2\theta-4\lambda\sin2\theta+t,
\]
\[
 F''=-\sin\theta-4\sin2\theta-8\lambda\cos2\theta.
\]

在 \([0,\pi/4]\) 上 \(F''<0\)；在 \([\pi/4,\pi/3]\) 上
\(F''\le-\sqrt2/2-2\sqrt3+4\lambda<0\)；在 \([\pi/3,\pi/2]\) 上
\(F'\le-1/2+t<0\)。精确 guards 给出

\[
 F'(.828)>0>F'(.829).
\]

故唯一峰 \(\theta_*\in(.828,.829)\)，\(C=F(\theta_*)\)。

以下反设 \(D\le d_*=9/400\)。各单面 deficit 非负，所以各自不超过 \(d_*\)。由

\[
 F(5/6)-F(.7)>.035,\qquad F(5/6)-F(.96)>.035
\]

及上述单峰性、G 分支分离，得到 \(.7<\theta_i<.96\)。

定义 apex 损失

\[
 L_t(\theta,z)=\sin^2\theta\left[\cot\frac\theta2-
 \cot\frac{\theta+z}{2}\right]-tz.
\]

在整个允许的 z 区间，

\[
 \partial_zL_t\ge1-\cos\theta-t>0.
\]

当 \(0\le z\le.1\) 时，

\[
 \frac{\sin\theta}{\sin((\theta+.1)/2)}
 =\frac{2\cos(\theta/2)}{\cos(.05)+\cot(\theta/2)\sin(.05)}
 \ge\frac{3871}{2500},
\]

其中使用 \(\cos(\theta/2)\ge1-.48^2/2=553/625\)、
\(\cot(\theta/2)\le20/7\)、\(\sin(.05)\le.05\)、\(\cos(.05)\le1\)。于是

\[
 \partial_zL_t\ge\tfrac12(3871/2500)^2-t>1.
\]

若某面 \(z\ge.1\)，其 deficit 至少为 .1，矛盾。因此 \(z_i<.1\) 且 \(L_t\ge z_i\)。

在 \([.7,.96]\) 上，

\[
 -F''(\theta)\ge\sin(.7)+4\sin(1.92)+(28/25)\cos(1.92)>4.
\]

这里 \(\sin(2\theta)\) 的最小值是两端较小的 \(\sin(1.92)\)。故 Taylor 公式给出
\(C-F(\theta)\ge2(\theta-\theta_*)^2\)。得到唯一一份预算

\[
 \boxed{2\sum_i(\theta_i-\theta_*)^2+
 \sum_i z_i+\sum_i\sigma_i+\sum_i\tau_i\le D.}
\tag{1.3}
\]

本节只给两面的支持/局部化，未使用任何三面几何 gap。

---

## 2. 精确 skew 约定及预算推出的紧区域

设 \(k_i=x_i/\sin\gamma_i\le1\)。定义 skew \(d_i\) 为三角形两个底角之差的一半，具体约定为

\[
 |c_i-p|/R=2k_i\cos(\gamma_i/2+d_i),\quad
 |c_{i+1}-p|/R=2k_i\cos(\gamma_i/2-d_i).
\tag{2.1}
\]

因此 \(|d_i|<(\pi-\gamma_i)/2<\pi/2\)。若在 apex 角平分线框架中写外 tip，这正对应

\[
 e^{im_i}[U_i\cos d_i+iV_i\sin d_i],
\]

其中 \(U_i=k_i(1+\cos\gamma_i)+\cos\theta_i\)、
\(V_i=k_i(1-\cos\gamma_i)-\cos\theta_i\)。

实际面积恒等式给出

\[
 \boxed{\sigma_i=\frac{2\sin^2\theta_i}{\sin\gamma_i}\sin^2d_i.}
\tag{2.2}
\]

记

\[
 s=(\theta_0+\theta_1)/2,\quad \Delta=\theta_0-\theta_1,\quad
 Z=(z_0+z_1)/2,\quad\eta=(d_0-d_1)/2,\quad T=\tau_0+\tau_1.
\]

从 (1.3)、\(D\le9/400\)、\(.828<\theta_*<.829\) 得到

\[
 .72<\theta_i<.94,\quad .753<s<.904,\quad
 |\Delta|\le.15,\quad 0\le Z\le9/800,\quad T\le9/400.
\tag{2.3}
\]

因为 \(0\le z_i\le9/400\)，且 \(\theta_i+z_i<\pi/2\)，函数
\(2\sin^2\theta/\sin(\theta+z)\) 随 \(\theta\) 增加、随 z 减少。因此

\[
 \frac{2\sin^2\theta_i}{\sin\gamma_i}
 \ge\frac{2\sin^2(.72)}{\sin(.7425)}>\frac54.
\]

于是 \(\sin^2d_i\le9/500\)，由 \(\sin^2(.14)>9/500\) 得到 \(|d_i|<.14\)。进一步

\[
 |\sin d_i|\ge(1-.14^2/6)|d_i|=(14951/15000)|d_i|.
\]

故

\[
 \boxed{\sigma_0+\sigma_1\ge c_\eta\eta^2,
 \qquad c_\eta=\frac52(14951/15000)^2
 =\frac{223532401}{90000000}.}
\tag{2.4}
\]

同时 \(|\eta|<12/125=.096\)。

定义

\[
 \boxed{\alpha=\frac{3s}{2}+\frac Z2+\eta,\qquad
 h=\frac s2-\frac Z2-\eta=2s-\alpha.}
\tag{2.5}
\]

令 \(A_0=\alpha-\theta_0\)、\(A_1=\alpha-\theta_1\)。可直接得到下列有理区间：

\[
 \alpha\in[2067/2000,11661/8000],\quad
 h\in[2199/8000,137/250],
\]
\[
 A_i\in[411/2000,5029/8000].
\tag{2.6}
\]

全都在第一象限。另有

\[
 x_i^2>13/30,\quad \sin\alpha>6/7,\quad\cos\alpha>1/10,
 \quad\sin A_i>81/400.
\tag{2.7}
\]

**符号审查结果。** 与 \(U\cos d+iV\sin d\) 及 (2.1) 配套的是 (2.5) 中的 \(+\eta\)、\(-\eta\)。用户写的 \(\alpha=3s/2+Z/2-\eta\)、\(h=s/2-Z/2+\eta\) 对应把两个 \(d_i\) 同时反号。两种约定可用，但不能将一种的 tip 公式与另一种的相位公式混用。

---

## 3. 零 tau 公式的精确作用域

将共同边 \(pc_1\) 的方向暂取为正实轴，则 apex 中相位为
\(m_0=-\gamma_0/2\)、\(m_1=\gamma_1/2\)。底边向量和右法向为

\[
 (c_{i+1}-c_i)/R=2ix_i e^{i(m_i-d_i)},\quad
 N_i=e^{i(m_i-d_i)}.
\]

外 tip 相对中间圆心满足

\[
 (q_0^*-c_1)/R=e^{-i(\theta_0+\gamma_0/2+d_0)},\quad
 (q_1^*-c_1)/R=e^{i(\theta_1+\gamma_1/2-d_1)}.
\]

再作一个共同旋转后，这两个单位向量恰为

\[
 a_0^*=e^{-i\alpha},\qquad a_1^*=e^{i\alpha},
\]

且

\[
 N_0=e^{-iA_0},\quad N_1=e^{iA_1},\quad T_i=iN_i.
\tag{3.1}
\]

这里未对两个三角形作镜像化；\(\theta_0,\theta_1,d_0,d_1,z_0,z_1\) 均可不同。

设 \(B_0=\sin\alpha\)、\(u_0=\sqrt{1-R^2B_0^2}\)。如果 tips 本身是单位锚点，则其平均值的两个平移分支为 \(\pm u_0\)；中间圆心为

\[
 c_1=\pm u_0-R\cos\alpha.
\]

在近支持区域 \(\cos\alpha>0\)，负分支使 \(|c_1|^2>1-R^2\)，所以 (0.1) 强制正分支。

令

\[
 a_* =\arcsin(R\sin\alpha),\qquad\chi=h+a_*.
\]

在正分支上精确有

\[
 \boxed{|c_0|^2=1+R^2-2R\cos(\chi+\Delta),\qquad
 |c_2|^2=1+R^2-2R\cos(\chi-\Delta).}
\tag{3.2}
\]

所以平均值为

\[
 \boxed{M_*:=\tfrac12(|c_0|^2+|c_2|^2)
 =1+R^2-2R\cos\chi\cos\Delta.}
\tag{3.3}
\]

最小值则是

\[
 1+R^2-2R\{\cos\chi\cos\Delta+|\sin\chi\sin\Delta|\}.
\]

在本近支持域，\(0<\chi<\pi/2\) 且 \(|\Delta|\le2\chi\)，故确可进一步上界为
\(1+R^2-2R\cos\chi\)。本证明不把该最后一步当成无条件代数推论。

固定 \(s,\eta\) 时，

\[
 \boxed{\partial_Z\chi=-\tfrac12+
 \frac{R\cos\alpha}{2\sqrt{1-R^2\sin^2\alpha}}<0.}
\tag{3.4}
\]

严格号来自 \(u_0^2-R^2\cos^2\alpha=1-R^2>0\)，不是镜像极值假设。

---

## 4. 非零 actual-anchor slack 的精确 lens 形式

仍在 (3.1) 的框架中，实际锚点误差写成

\[
 e_i=(q_i-q_i^*)/R=-h_iN_i+\ell_iT_i,\qquad
 h_i=\tau_i/x_i\ge0.
\tag{4.1}
\]

底边中点到锚点的归一化法向高度为 \(y_i-h_i\)，其中 \(y_i=\cos\theta_i\)。两个实际距离限制**等价于**

\[
 \boxed{(x_i+|\ell_i|)^2+(y_i-h_i)^2\le1.}
\tag{4.2}
\]

也就是

\[
 \ell_i^2+2x_i|\ell_i|+h_i^2\le2y_ih_i.
\]

因此保留方向符号的精确约束是

\[
 |\ell_i|\le w_i(h_i):=
 \sqrt{x_i^2+2y_ih_i-h_i^2}-x_i,
\tag{4.3}
\]

而线性放宽为

\[
 |\ell_i|\le\frac{y_i}{x_i}h_i,\qquad |e_i|\le\frac{\tau_i}{x_i^2}.
\tag{4.4}
\]

关键不是单独使用最后的模长估计，而是在第 6 节先用 (4.1) 的法向符号和切向范围求线性支持。模长界只用于二阶余项。

这里即使 \(A_i<0\) 也没有任何变化。近支持预算随后自动限制 \(h_i\) 很小。

---

## 5. 两单位锚点消除平移后的精确公式

记

\[
 g=(e_0+e_1)/2,\quad d_b=(e_1-e_0)/2,\quad
 b=(0,B_0)+d_b,\quad B=|b|,
\]
\[
 n=(b_y,-b_x)/B,\quad u=\sqrt{1-R^2B^2},\quad
 k=(\cos\alpha,0)+g.
\]

两实际单位锚点之差为 \(2Rb\)，平均值必为 \(\pm un\)，故

\[
 c_1=\pm un-Rk.
\tag{5.1}
\]

这就是无平移自由参数的精确解。并非把实际圆心挪到径向内端。

由 (4.4) 及 (2.7)，

\[
 |g|,|d_b|\le P:=\frac{15}{13}T\le\frac{27}{1040}.
\tag{5.2}
\]

沿 \(e_i(v)=ve_i\)、\(0\le v\le1\) 的同伦，两个 lens 都由凸性保留。尤其两个 anchor-minus-middle 向量仍各自长度不超过 1，因此 \(|b(v)|\le1\)。同时

\[
 |b(v)|\ge6/7-27/1040>83/100,\qquad u(v)>99/100.
\tag{5.3}
\]

还有 \(n_x\ge83/100\)，所以

\[
 k\cdot n\ge(1/10)(83/100)-27/1040>0.
\]

故 (0.1) 同样排除了实际非零 tau 情形的负分支，实际必须使用 (5.1) 的正号。

定义

\[
 E_*=-\cos h(\cos\Delta,\sin\Delta),
\]
\[
 F=(-\sin h\sin\Delta,\ B_0+\sin h\cos\Delta),\quad E=E_*-g.
\tag{5.4}
\]

由实际底边向量直接相加得到

\[
 (c_0+c_2)/2=un+RE,\qquad (c_2-c_0)/2=RF.
\]

从而

\[
 \boxed{M:=\tfrac12(|c_0|^2+|c_2|^2)
 =u^2+2Ru\,n\cdot E+R^2(|E|^2+|F|^2).}
\tag{5.5}
\]

并且还有完整的最小值公式

\[
 \boxed{\min(|c_0|^2,|c_2|^2)
 =M-2R|u\,n\cdot F+R E\cdot F|.}
\tag{5.6}
\]

(5.5) 与用户给出的通用平移恒等式一致，但此处 E、F、b 已明确连接到 actual tau 的法向/切向误差。

---

## 6. 单份 tau 的较小代价：线性支持与二阶余项

定义 \(f(b)=\sqrt{1-R^2|b|^2}\,J_{\rm cw}b/|b|=un\)。将 (5.5) 在 \(e_0=e_1=0\) 展开，记 \(c_h=\cos h\)、\(c_\Delta=\cos\Delta\)、\(s_\Delta=\sin\Delta\)。线性部分为

\[
 L_0\cdot e_0+L_1\cdot e_1,
\]

其中

\[
 L_{0x}=-Ru_0+R^2c_hc_\Delta-Ru_0c_hs_\Delta/B_0,
\]
\[
 L_{1x}=-Ru_0+R^2c_hc_\Delta+Ru_0c_hs_\Delta/B_0,
\]
\[
 L_{0y}=R^2c_hs_\Delta+R^2B_0-R^3B_0c_hc_\Delta/u_0,
\]
\[
 L_{1y}=R^2c_hs_\Delta-R^2B_0+R^3B_0c_hc_\Delta/u_0.
\tag{6.1}
\]

由 (2.6)—(2.7)，

\[
 L_0\cdot T_0<0,\qquad L_1\cdot T_1>0.
\tag{6.2}
\]

例如设 \(a=R(99/100)-R^2-7R/40\)，则 \(L_{ix}\le-a\)，
\(L_{0y}\le23R^2/20\)、\(L_{1y}\ge-23R^2/20\)，而
\(a/5-23R^2/20>0\)。这证明 (6.2) 在整个区域成立。

用 (4.1)—(4.4) **带符号**地取支持后，得到

\[
 L_0\cdot e_0+L_1\cdot e_1\le\kappa_0\tau_0+\kappa_1\tau_1,
\]

其中不是粗范数系数，而是

\[
 \boxed{\kappa_0=
 \frac{R\rho_m}{\sin^2\theta_0}
 \left[B_0\left(1-\frac{R\cos h\cos\Delta}{u_0}\right)
 +\cos h\sin\Delta\right],}
\]
\[
 \boxed{\kappa_1=
 \frac{R\rho_m}{\sin^2\theta_1}
 \left[B_0\left(1-\frac{R\cos h\cos\Delta}{u_0}\right)
 -\cos h\sin\Delta\right],\quad
 \rho_m=u_0-R\cos\alpha.}
\tag{6.3}
\]

推导中的关键约分是

\[
 -N_0-\cot\theta_0\,T_0=(-\sin\alpha,-\cos\alpha)/x_0,
\]
\[
 -N_1+\cot\theta_1\,T_1=(-\sin\alpha,\cos\alpha)/x_1.
\]

### 6.1 仅八个完整区间 slab 的精确线性 guard

本文件证明

\[
 \boxed{\kappa_0,\kappa_1<\frac{41}{200}.}
\tag{6.4}
\]

令 \(\omega=Z/2+\eta\in[-12/125,813/8000]\)，则
\(\alpha=3s/2+\omega\)、\(h=s/2-\omega\)。把
\([753/1000,113/125]\) 等分成 8 个 s 区间；omega 保留整个区间，不选取点。

设 \(k=41/200\)。不等式 \(\kappa_0\le k\) 等价于

\[
 A+Q\cos\Delta+P\sin\Delta\ge0,
\]

其中

\[
 A=k/2-R\rho_m B_0,
\quad Q=-k\cos2s/2+R^2\rho_mB_0\cos h/u_0,
\]
\[
 P=k\sin2s/2-R\rho_m\cos h.
\tag{6.5}
\]

\(\kappa_1\) 对应 \(\Delta\mapsto-\Delta\)。在每个完整 s slab 上，
\(\rho_m(\alpha)\)、\(\rho_m(\alpha)\sin\alpha\) 随 alpha 严格递增：这由
\(u_0^2-R^2\cos^2\alpha=1-R^2>0\)、\(\cos2\alpha<0\)、\(1-2R^2\sin^2\alpha>0\) 直接微分得到。

所以可用 alpha 的两个端点包围 rho 和 rho B；cos h 用 h 端点包围；sin 2s 的下界取 s 两端较小者（sin 在这里凹），上界取 1。于是获得严格区间 \(A\ge A_-\)、\(Q\ge Q_-\)、\(P\in[P_-,P_+]\)。整个 \(|\Delta|\le3/20\) 上的下界为

\[
 A_-+\min\{Q_-,Q_-\cos(3/20)\}
 -\max\{|P_-|,|P_+|\}\sin(3/20).
\tag{6.6}
\]

八个 slab 全部严格大于 \(1/4000\)；最小区间下界约 .0003061330347，精确有理值在 verification.json 内。
这不是端点采样：omega 的完整范围及 Delta 的完整范围已通过 (6.5)—(6.6) 包围。

### 6.2 二阶余项，只用剩余的二次误差

对 b 的径向、切向分量记为 a、v，可直接微分得到

\[
 \|Df\|\le\frac{100}{83},
\]

以及

\[
 \|D^2f(b)[w,w]\|
 \le\left(\frac{R^2}{u^3}+\frac{2}{\sqrt3\,u|b|^2}\right)|w|^2
 <\frac74|w|^2.
\tag{6.7}
\]

这里的 \(2/\sqrt3\) 来自
\(\sqrt{v^4+4a^2v^2}\le(2/\sqrt3)(a^2+v^2)\)，而不是模糊的整体平移 Lipschitz 常数。用 \(u>.99\)、\(|b|>.83\)、\(2/\sqrt3<7/6\) 即可严格验证 (6.7)。

将 (5.5) 完整展开，线性部分之外的余项上界为

\[
 R^2|g|^2+R(7/4)|d_b|^2+2R(100/83)|d_b||g|.
\]

因此由 (5.2)，

\[
 \text{余项}\le
 \left[R^2+\frac74R+\frac{200}{83}R\right](15/13)^2T^2
 <\frac{13}{20}T^2.
\tag{6.8}
\]

最终得到可复用的 refined tau 引理：

\[
 \boxed{M\le M_*+\kappa_0\tau_0+\kappa_1\tau_1+\frac{13}{20}T^2.}
\tag{6.9}
\]

以及统一版本

\[
 \boxed{M\le M_*+\frac{41}{200}T+\frac{13}{20}T^2
 \le M_*+\frac{11}{50}T.}
\tag{6.10}
\]

最后一步使用 \(T\le9/400\)，且
\(41/200+(13/20)(9/400)=1757/8000<11/50\)。

只有一份 \(T=\tau_0+\tau_1\)。它没有另进全局 Q、BI、独立 J 费或弧面积账本。

---

## 7. 从 tau 引理闭合强范数域 fee，不使用 variance

(1.3)、(2.4) 给出

\[
 4(s-\theta_*)^2+c_\eta\eta^2+\Delta^2+2Z+T\le D.
\tag{7.1}
\]

设 \(E=4(s-\theta_*)^2+c_\eta\eta^2\)。由 (3.4)，丢掉非负 Z 只会增大 chi。
对 \(g(s,\eta)=s/2-\eta+\arcsin(R\sin(3s/2+\eta))\)，整个区域内

\[
 0<\partial_sg<3/5,\qquad -1<\partial_\eta g<0.
\]

第一个上界由 alpha 最小端点的精确 guard 得到。又因为

\[
 R\sin(2487/2000)<\sin(43/400),
\]

并且 \(\theta_*<.829\)，有 \(g(\theta_*,0)<261/500=.522\)。因此

\[
 \chi\le\frac{261}{500}+rac35|s-\theta_*|+|\eta|
 \le\frac{261}{500}+\frac{71}{100}\sqrt E.
\tag{7.2}
\]

这里使用精确不等式

\[
 \frac9{100}+\frac1{c_\eta}<\left(\frac{71}{100}\right)^2.
\]

在近支持区 cos chi 为正。由 (3.3)、(6.10)、\(1-\cos\Delta\le\Delta^2/2\)，

\[
 M\le1+R^2-2R\cos\chi+R\Delta^2+\frac{11}{50}T.
\]

注意 \(R<11/50\)，故 (7.1)—(7.2) 给出

\[
 M\le\Phi(E):=1+R^2-
 2R\cos\left(\frac{261}{500}+\frac{71}{100}\sqrt E\right)
 +\frac{11}{50}(d_*-E),\quad 0\le E\le d_*.
\tag{7.3}
\]

这一步在同一预算内合并了 \(\Delta^2+T\)，没有同时把它们各自上界为一整份 D。

对 \(0<E\le d_*\)，

\[
 \Phi'(E)\ge
 \frac{R(71/100)\sin(261/500)}{3/20}-\frac{11}{50}>\frac1{25}>0.
\]

所以无需优化器，最大值就是右端点 E=d_*。精确 guard 给出

\[
 \boxed{M\le1+R^2-2R\cos(1257/2000)
 <\frac{83}{100}-\frac3{10000}.}
\tag{7.4}
\]

较窄的实际包围是

\[
 1+R^2-2R\cos(1257/2000)
 =.8296491164252835007624107\ldots.
\]

但 (0.1) 蕴含 \(M\ge83/100\)，矛盾。故反设 \(D\le9/400\) 不成立，定理 (0.2) 得证。

---

## 8. 可复用的实际 HALFarc 与 signed tau 接口

这一节不是 fee 证明的前提，也不额外支付 variance。若 q0、q1 是中间完整 assigned arc 的两个端点，那么

\[
 \boxed{a_{\rm mid}=\arcsin(R|b|),\qquad 0<a_{\rm mid}\le\arcsin R.}
\tag{8.1}
\]

这里是 HALFarc；与底边 \(x_i=\sin\theta_i\) 无关。

写 \(d_b=(r,w)\)。由 (4.1)、(3.1)，

\[
 w=-\tfrac12(h_0\sin A_0+h_1\sin A_1)
 +\tfrac12(\ell_1\cos A_1-\ell_0\cos A_0).
\]

因此有带方向的界

\[
 -\frac12\left[
 \frac{\tau_0\cos(h+\Delta)}{x_0^2}+
 \frac{\tau_1\cos(h-\Delta)}{x_1^2}\right]
 \le w\le
 \frac{\cos\alpha}{2}
 \left(\frac{\tau_0}{x_0^2}+\frac{\tau_1}{x_1^2}\right).
\tag{8.2}
\]

结合 \(B_0+w>83/100\) 和
\(\sqrt{a^2+r^2}\le a+r^2/(2a)\)，得到

\[
 B_-:=B_0-\frac12\left[
 \frac{\tau_0\cos(h+\Delta)}{x_0^2}+
 \frac{\tau_1\cos(h-\Delta)}{x_1^2}\right]
 \le |b|,
\]
\[
 |b|\le B_+:=\min\left\{1,\ B_0+
 \frac{\cos\alpha}{2}\left(\frac{\tau_0}{x_0^2}+\frac{\tau_1}{x_1^2}\right)
 +\frac{81}{100}T^2\right\}.
\tag{8.3}
\]

\(B_->0\) 由同一已证区域保证。所以

\[
 \boxed{\arcsin(RB_-)\le a_{\rm mid}\le\arcsin(RB_+).}
\tag{8.4}
\]

更精细的无需线性放宽版本，直接使用 (4.3) 的 \(w_i(h_i)\)：

\[
 w_-^{\rm exact}=-\tfrac12\sum_i[h_i\sin A_i+w_i(h_i)\cos A_i],
\]
\[
 w_+^{\rm exact}=-\tfrac12\sum_i[h_i\sin A_i-w_i(h_i)\cos A_i],
\]
\[
 |r|\le\tfrac12\left[|h_0\cos A_0-h_1\cos A_1|
 +w_0(h_0)\sin A_0+w_1(h_1)\sin A_1\right].
\]

将这些代入 \(|b|=\sqrt{(B_0+w)^2+r^2}\) 即得更窄区间。

如果以后使用 variance，可以只添加
\(K\operatorname{dist}(\pi/33,[\arcsin(RB_-),\arcsin(RB_+)])^2\) 一次；不能将同一个 T 同时当作另一个独立面积亏损。

**真正等距 switch 的可选强化。** 若父线程确有
\(|q_i-c_i|=|q_i-c_{i+1}|\)，则由两距离平方相减严格得到
\(\ell_i=0\)。这不是 (0.2) 的前提。此时线性系数可直接改为
\(-L_i\cdot N_i/x_i\)，余项里也可使用 \(|e_i|=\tau_i/x_i\)。不要在仅给出“两端都被覆盖”的一般 anchor 域中擅自使用这个强化。

---

## 9. 完整弧长度到范数下界，不强制 radial 内端

若 c 的 R 圆覆盖一条半角为 a 的完整单位圆周弧，取其两个端点的平均距离平方，得到

\[
 1+|c|^2-2|c|\cos a\le R^2.
\]

这里仅使用圆心方向与弧中方向的余弦不超过 1，不要求两个方向相同。因此

\[
 |c|\ge\cos a-\sqrt{R^2-\sin^2a}.
\]

右端在 \(0<a<\arcsin R\) 上严格递增。完整弧长 >147/1000 时，

\[
 |c|>
 \cos(147/2000)-\sqrt{R^2-\sin^2(147/2000)}
 >.911124459.
\]

这当然强于使用 \(.911088875\) 的较弱输入，也强于本定理只要求的平方下界 83/100。该计算没有从 assigned-cell lens 扣取任何面积或亏损。

---

## 10. 闭包、代码证明范围与独立复核接口

x=0 的参数化退化闭包满足
\(0\le x^2\cot(\gamma/2)\le2x\to0\)、\(|A|/R^2\le x\sqrt{1-x^2}\to0\)，故其 score 至多为 t pi，远离完整支持。x=1 时两 apex 端点合并，已经包含在 F、G 内。

反设 D<=9/400 后，证明自动获得 theta、gamma、d 的严格紧区间、b>.83、u>.99。所有低-deficit 相关三角形和单位锚点对都因此远离所用分母的零点；退化序列不产生未检查的低-deficit 分支。非严格外接圆半径、非严格端点覆盖和非严格范数输入均已允许。

verify.py 的接受计算只有：精确有理四则运算、整数向外舍入、整数平方根，以及附 Lagrange 余项的 Taylor 多项式。48 个 guards 包含 8 个覆盖完整连续参数范围的 s-slabs。程序没有通过“残差区间含零”宣称几何恒等式；这些恒等式均在正文代数推导。代码也不声称其常数 guards 本身代替了正文的全域几何归约。

复跑命令：

```bash
python verify.py
```

verification.json 记录精确有理上下端点；verification.log 为实际运行输出。结论仍只是一条可以嵌入 highM 链的两面定理，不能独自宣布全局 .11322 或精确最优已完成。
