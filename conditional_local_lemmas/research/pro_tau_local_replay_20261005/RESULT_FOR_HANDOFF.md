# 父线程接口：非零 actual-anchor tau 已处理，强范数域两面 fee > .0225

参数固定 R=5661/50000、lambda=7/50、t=1071883/20000000；C 是整个 apex/skew 域的同参数支持。

## 可直接复核的新定理

两实际连续正向三角形 (p,c0,c1),(p,c1,c2)，各 circumradius<=R；每个实际 unit anchor qi 被对应底边两端的 R 圆覆盖。允许 signed A<0，不要求 anchors 等距，不要求镜像、radial-inner、凸边界链或 nearest-owner。

如果 |c1|²<=1-R² 且 |c0|²,|c2|²>=83/100，则

    D=2C-S2 > 9/400 = .0225.

目标 fee 与 .0225 的精确差为

    1740724429841728225095328099 /
    8011730250000000000000000000000 > 1/5000.

所以题设 joint D+K(a_mid-pi/33)^2 直接超过目标，且没有动用 a_mid variance。

highM 已给范数 .911088875 的平方 >83/100。完整 assigned arc>.147 直接给更强的 norm>.911124459；证明不强制圆心处于内端。

## 相位修正

若 d 的约定与 tip e^{im}[U cos d+iV sin d] 一致，则

    alpha=3s/2+Z/2+eta,
    h=s/2-Z/2-eta,
    eta=(d0-d1)/2.

先前 alpha 中的 -eta 对应将两个 d 同时反号。正 unit-anchor 平移分支零 tau 精确式：

    rho0²=1+R²-2R cos(chi+Delta),
    rho2²=1+R²-2R cos(chi-Delta),
    chi=h+asin(R sin alpha).

chi 对 Z 严格递减。负平移分支由中间圆心 norm<=sqrt(1-R²) 排除。未假设一般构型有镜像极值。

## 处理 actual tau 的新引理

取 Ni 为底边右法向，Ti=iNi，ei=(qi-tip_i)/R=-hi Ni+li Ti，hi=tau_i/xi_base。
精确 lens：

    (x_i+|li|)²+(cos(theta_i)-hi)²<=1.

其法向符号与切向界用于线性项，而不是先取粗模长。令 T=tau0+tau1，D<=9/400 时：

    M=(rho0²+rho2²)/2
      <= M_star+kappa0*tau0+kappa1*tau1+(13/20)T²
      <= M_star+(11/50)T,

    M_star=1+R²-2R cos(chi)cos(Delta),
    kappa_i<41/200.

kappa 的完整表达式、精确消平移公式、8 个连续 s-slabs 的守卫见 PROOF.md 第 5、6 节。

单份预算给出

    4(s-theta*)²+c_eta eta²+Delta²+2Z+T<=D,
    c_eta=223532401/90000000.

由此无需优化器，严格推得在 D<=9/400 时

    M<=1+R²-2R cos(1257/2000)
      <83/100-3/10000,

与两侧范数下界矛盾。

## 计费与范围

只有一份 T。未把 J-A/tau 再放入 Q、BI、独立J费用或 assigned lens 面积；未支付 variance。还给出实际 HALFarc 的 signed-tau 区间，供未来联合预算使用。若是严格真实等距 switch，可额外令切向 li=0；主定理不需要这一前提。

新 verify.py 已实际通过 48 个精确 guards（含 8 个全区间 slabs）。没有跑旧三面、402 或其他旧全行 checker，没有做构造搜索或优化。

不是全局100圆证明。highM 外部桥接、其余行及七面区域仍应由原执行线独立复核；本轮不改变已认证全局区间。
