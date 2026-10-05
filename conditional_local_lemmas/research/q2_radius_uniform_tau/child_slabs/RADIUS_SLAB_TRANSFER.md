# 冻结八个 slab 的半径传输：纯解析推导

本文只处理 two-anchor tau lemma 第 6.1 节的线性支持，作用域仍是两张实际相邻三角形及其条件参数域。它不是全局 disk-radius 认证。没有执行或导入数学程序、重跑 checker、搜索样点，亦未修改冻结来源。

来源是 `research/pro_sources_20261005/sources/two_anchor_tau_lemma/two_anchor_tau_lemma/PROOF.md` 第 6.1 节、`verify.py` 的 `rho_product` 与八个 slab 循环，以及 `verification.json` 的已审有理下界。

## 1. 固定角域及精确连续函数

固定

\[
R_0=\frac{5661}{50000},\qquad k_0=\frac{41}{200},\qquad
\lambda=\frac7{50},\quad t=\frac{1071883}{20000000},\quad D_{\max}=\frac9{400}.
\]

保留原来的八个完整 slab，并保留

\[
s\in[753/1000,113/125],\quad
\omega\in[-12/125,813/8000],\quad |\Delta|\le d:=3/20,
\]

\[
\alpha=3s/2+\omega,
\qquad h=s/2-\omega.
\]

所有半径传输估计均在整个

\[
\boxed{R_0\le R\le3/25=.12}
\]

上成立。角域及已审的角度不等式不依赖 R，特别是

\[
\sin\alpha>6/7,\qquad \cos\alpha>1/10,
\qquad 0<h<\pi/2,
\qquad \pi/4<\alpha<\pi/2.
\]

对固定角 a，定义精确实函数

\[
b=\sin a,\quad c=\cos a,\quad
u(R,a)=\sqrt{1-R^2b^2},\quad
\rho(R,a)=u(R,a)-Rc.
\]

在 slab \([s_-,s_+]\) 内写

\[
a_-={3s_-\over2}-{12\over125},\quad
a_+={3s_+\over2}+{813\over8000},
\]

\[
h_-={s_-\over2}-{813\over8000},\quad
h_+={s_+\over2}+{12\over125}.
\]

定义下列精确端点函数，对应冻结循环中各区间表达式的理想实数值：

\[
\mathcal A(R)=k_0/2-R\rho(R,a_+)\sin a_+,
\]

\[
\mathcal Q(R)=-{k_0\over2}\cos(2s_-)
+R^2\rho(R,a_-)\sin a_-\cos h_+,
\]

\[
\mathcal P_-(R)={k_0\over2}
\min\{\sin(2s_-),\sin(2s_+)\}
-R\rho(R,a_+)\cos h_-,
\]

\[
\mathcal P_+(R)=k_0/2-R\rho(R,a_-)\cos h_+.
\]

最后令

\[
\boxed{\mathcal L(R)=\mathcal A(R)
+\min\{\mathcal Q(R),\mathcal Q(R)\cos d\}
-\max\{|\mathcal P_-(R)|,|\mathcal P_+(R)|\}\sin d.}
\tag{1}
\]

这里 \(\mathcal Q\) 的正项没有除以 u；这是 `verify.py` 的保守下界，而非正文精确 Q。正项为正且 u<1，故省略 1/u 的方向正确。

## 2. 端点包围在半径变化后仍然有效

\[
\partial_a\rho=Rb(1-Rc/u)>0,
\]

因为 \(u^2-R^2c^2=1-R^2>0\)。同时

\[
\partial_a(\rho\sin a)
=\frac{c(1-2R^2b^2)}u-R\cos(2a)>0.
\]

后一个严格号使用 a 的第二倍角在 \((\pi/2,\pi)\)、c>0 及 R<=3/25。因此，原来的 rho 与 rho*sin(alpha) 端点包围在整个半径区间保留。cos h、cos 2s、sin 2s 的端点包围没有半径依赖。

令 \(G_i(k)=\sin^2\theta_i(k-\kappa_i)\)，其中 \(\theta_0=s+\Delta/2\)、\(\theta_1=s-\Delta/2\)。正文的精确 A,Q,P 恒等式及上述保守端点包围给出

\[
\boxed{G_0(k_0)\ge\mathcal L(R),\qquad
G_1(k_0)\ge\mathcal L(R).}
\tag{2}
\]

第二个式子只是 Delta 反号；(1) 同时覆盖两个符号及完整 Delta 区间。

## 3. 半径导数：A、P 控制损失，Q 增加

在 R<=3/25 时有

\[
u>99/100,
\qquad {R^2b^2\over u}<\frac4{275}<\frac3{200}.
\]

直接微分得到

\[
H(R,a):=\partial_R(R\rho)
=u-2Rc-{R^2b^2\over u},
\]

\[
J(R,a):={1\over R}\partial_R(R^2\rho)
=2u-3Rc-{R^2b^2\over u}.
\]

两者满足

\[
\boxed{\frac{147}{200}<H<m:=1-\frac{R_0}{5}
=\frac{244339}{250000}<1,}
\tag{3}
\]

\[
\boxed{\frac{321}{200}<J<2.}
\tag{4}
\]

H 的下界来自 .99-.24-.015=.735；J 的下界来自 1.98-.36-.015=1.605。H 的上界使用 c>1/10、R>=R0 及 u<=1。

于是

\[
-m<\mathcal A'(R)<-63/100,
\]

因为 \(\sin a_+>6/7\)，并且 \((147/200)(6/7)=63/100\)。同时

\[
0<\mathcal Q'(R)=R J(R,a_-)\sin a_-\cos h_+<6/25.
\]

每个 P 端点的导数为 \(-H\cos h_{\pm}\)，故绝对值小于 m。绝对值与最大值运算保留 m-Lipschitz 界。

函数 \(q\mapsto\min(q,q\cos d)\) 非降、1-Lipschitz；这也覆盖 Q 穿过零点的情形。因此 Q 部分不会削弱半径增加后的下界，也不必确认其符号。

对任意 R2>R1>=R0、R2<=3/25，(1) 给出

\[
\mathcal L(R_2)-\mathcal L(R_1)
>-m(1+\sin d)(R_2-R_1).
\]

因为 sin d<d=3/20，

\[
m(1+\sin d)
<\frac{244339}{250000}\frac{23}{20}
=\frac{5619797}{5000000}<\frac98<\frac{23}{20}.
\]

故可直接使用

\[
\boxed{\mathcal L(R)\ge\mathcal L(R_0)
-\frac98(R-R_0)
\ge\mathcal L(R_0)-\frac{23}{20}(R-R_0).}
\tag{5}
\]

还可从同一导数获得上界：A 的降幅至少 (63/100)*(R2-R1)，Q 的增幅至多 (6/25)*(R2-R1)，P 绝对值项给 L 的增幅至多 (3/20)*(R2-R1)。因此

\[
\boxed{-\frac98(R_2-R_1)
<\mathcal L(R_2)-\mathcal L(R_1)
<-\frac6{25}(R_2-R_1).}
\tag{6}
\]

每个 slab 的精确保守函数因此严格递减且双侧 9/8-Lipschitz。结论不要求 A,Q,P 的绝对值选择分支保持固定。

## 4. 继承冻结有理下界

设 g_j 是 `verification.json` 对第 j 个 slab 输出的已审 gap lower。原区间算术包含上述精确实函数的值，故只需继承

\[
\mathcal L_j(R_0)\ge g_j>1/4000.
\]

不需要重新计算任何三角函数、平方根或新 guard；对于所有 R0<=R<=3/25，(5) 直接给出

\[
\boxed{G_i(k_0)\ge g_j-\frac98(R-R_0)
\ge g_j-\frac{23}{20}(R-R_0).}
\tag{7}
\]

特别在 R<=1133/10000=.1133 时，R-R0<=1/12500，故

\[
G_i(k_0)>\frac1{4000}-\frac98\frac1{12500}
=\frac1{6250}>0.
\]

所以原 k0=.205 支持系数至少连续推广至 .1133。仅继承 g_j>1/4000，亦可推广至 .11344，残差 >1/400000。

若进一步直接比较冻结 json 中已存在的有理数，八个 g_j 都大于 30613/100000000=.00030613：第一项下界是 .0003061330347...；后七项依次大于 .003、.005、.008、.011、.014、.017、.020。这些是现有有理输出的比较，并非新浮点求值。在 R<=.11349 时，R-R0<=27/100000，故

\[
G_i(k_0)>\frac{30613}{100000000}
-\frac98\frac{27}{100000}
=\frac{119}{50000000}>0.
\]

## 5. 增加 k 的 pointwise 余量：主线程所用 .116 范围

对任意 epsilon>=0，精确恒等式为

\[
\boxed{G_i(k_0+\epsilon)
=G_i(k_0)+\epsilon\sin^2\theta_i.}
\tag{8}
\]

加价应施于精确 G_i；无需声称重新定义后的保守 slab 函数会自动增加相同的量。在原条件参数域内已审有 \(\sin^2\theta_i>13/30\)。取

\[
k_1=213/1000=k_0+1/125,
\qquad R\le29/250=.116,
\]

则 R-R0<=139/50000。仅用较松的 23/20 版本及所有已审 g_j>1/4000，得到

\[
\begin{aligned}
G_i(k_1)
&>\frac1{4000}-\frac{23}{20}\frac{139}{50000}
+\frac1{125}\frac{13}{30}\\
&=\boxed{\frac{1559}{3000000}>0.}
\end{aligned}
\tag{9}
\]

所以同一条件角域内，原已审八个 slab 支持 \(\kappa_i<213/1000\) 在完整 R 范围 [.11322,.116] 成立，而不用重跑任何 slab guard。

## 6. 带符号切向支持的半径保留

正文 (6.2) 的充分量为

\[
a=R(99/100)-R^2-7R/40.
\]

在 R<=3/25 上，

\[
\frac a5-\frac{23R^2}{20}
=R\left(\frac{163}{1000}-\frac{27R}{20}\right)
\ge\frac R{1000}>0.
\]

因此原始带符号切向支持的方向也在此半径范围保留。所有其他几何前提、tau 二阶余项及最终范数闭合由主证明分别处理；本文没有将本线性传输扩张成全局结论。
