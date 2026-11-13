# 哥德巴赫猜想的圆法框架及次弧困难分析

## 1. 题目、定义与本稿结论

题目为：对每一个大于 2 的偶数 n，寻找素数 p、q，使 n=p+q。本稿沿用主弧、次弧和高次素数幂的分解，逐项修正此前引理 3、引理 6 及最终结论。已证结论包括：圆法恒等式；带完整误差项的主弧下界；次弧核的精确公式；Type II 双线性式的有效初等上界；高次素数幂的统一上界；以及全部偶数 4≤n≤10000 的显式素数证书。无限区间的充分条件列在第 3 卷末，作为继续研究的准确接口。

以下 log 表示自然对数，e(x)=exp(2πix)。素数为恰有两个正因子的正整数。记 Λ(p^k)=log p（p 为素数，k≥1），其余正整数的 Λ 值为零，并约定 Λ(0)=0。

对给定偶数 n，定义

$$S_n(\alpha)=\sum_{k=1}^{n}\Lambda(k)e(k\alpha),\qquad \psi(n)=\sum_{k=1}^{n}\Lambda(k).$$

$$R(n)=\sum_{k=1}^{n-1}\Lambda(k)\Lambda(n-k),\quad R_{\rm prime}(n)=\sum_{p+q=n\atop p,q\in\mathbb P}\log p\log q,\quad R_{\rm pow}=R-R_{\rm prime}.$$

和式按有序对计数。所有权重非负，因此 R_prime(n)>0 等价于 n 有两个素数之和的表示。

### 定义 1.1：主弧和次弧

在圆周 T=R/Z 上取

$$P=\tfrac12\log n,\quad Q=\frac{n}{(\log n)^2},\quad \eta_q=\frac1{qQ}.$$

对既约中心 a/q（1≤q≤P），定义 M(q,a)={α∈T:dist_T(α,a/q)≤η_q}，M 为这些弧的并，m=T\M。圆周距离自动处理 0 与 1 的接合。

对可测集合 A⊂T，定义

$$W_A(h)=\int_A e(h\alpha)d\alpha,\qquad a_A(m)=\int_A S_n(\alpha)e(-m\alpha)d\alpha.$$

### 引理 1：正交性、互不相交性与精确分解

对 n≥400，上述主弧两两不交，并且

$$R(n)=I_M(n)+I_m(n),\qquad I_A(n)=\int_A S_n(\alpha)^2e(-n\alpha)d\alpha.$$

**证明。** 整数 h 的单位周期积分为 δ(h,0)。展开 S_n² 并逐项积分，恰好保留 k+l=n 的项，得到 R(n)。

函数 f(x)=(log x)^3/x 的导数为 (log x)^2(3-log x)/x²，在 x≥400 为负；f(400)<1，因此 Q>log n=2P。不同既约中心的圆周距离至少为 1/(qq')，而

$$\eta_q+\eta_{q'}=\frac{q+q'}{qq'Q}\le\frac{2P}{qq'Q}<\frac1{qq'}.$$

故主弧互不相交。其测度为

$$u:=|M|=\sum_{q\le P}\frac{2\varphi(q)}{qQ}\le\frac{2P}{Q}=\frac{(\log n)^3}{n}<1.$$

由积分可加性即得分解。圆周关于 α→−α 对称，而 S_n(−α)=overline(S_n(α))，所以 I_M、I_m 均为实数。证毕。

### 引理 2：奇异级数及其显式上下界

令 μ 为莫比乌斯函数，φ 为欧拉函数，c_q(h)=Σ_(a mod q,(a,q)=1)e(ah/q)。对固定偶数 n，

$$\mathfrak S(n):=\sum_{q\ge1}\frac{\mu(q)^2c_q(n)}{\varphi(q)^2}=2C_2\prod_{p\mid n,\ p>2}\frac{p-1}{p-2},\qquad C_2=\prod_{p>2}\left(1-\frac1{(p-1)^2}\right).$$

该级数绝对收敛，而且 0.6557<C₂≤3465/5120。

**证明。** 对素数 p，c_p(n)=p−1（p|n）或 −1（p∤n），中国剩余定理给出关于 q 的积性。p|n 仅有有限项，其余绝对局部项 1/(p−1)² 可求和。展开绝对收敛的欧拉乘积得到公式。p=2 贡献因子 2。

前四个奇素数的乘积为

$$P_4=\frac34\frac{15}{16}\frac{35}{36}\frac{99}{100}=\frac{3465}{5120}=0.6767578125.$$

其余因子在 (0,1)，给出上界。对 p≥13，设 x_p=1/(p−1)²≤1/144。由幂级数，log(1−x_p)≥−x_p−x_p²。13、17、19、23、29、31、37、41、43、47 的 x_p 之和小于 0.020831。p≥53 时 p−1 是至少为 52 的偶数，故

$$\sum_{p\ge53}x_p\le\frac14\sum_{k=26}^{\infty}\frac1{k^2}<\frac14\sum_{k=26}^{\infty}\frac1{k(k-1)}=0.01.$$

所以 Σx_p<0.030831，Σx_p²<0.000215。于是 C₂≥P₄ exp(−0.031046)≥P₄(1−0.031046)>0.6557。特别地 S(n)>1.3114。证毕。

### 引理 3（修订）：保留全部误差项的主弧下界

定义有限三角多项式 T_n(β)=Σ_(k=1)^n e(kβ)，并在弧 M(q,a) 上写成精确等式

$$S_n(a/q+\beta)=V_q(\beta)+E_{q,a}(\beta),\quad V_q(\beta)=\frac{\mu(q)}{\varphi(q)}T_n(\beta).$$

误差 E 由此等式定义。记

$$D_{q,a}=\int_{-\eta_q}^{\eta_q}|E_{q,a}(\beta)|^2d\beta,\qquad \mathfrak S_P(n)=\sum_{q\le P}\frac{\mu(q)^2c_q(n)}{\varphi(q)^2},$$

$$A_{\rm arc}(n)=\sum_{q\le P}\frac{\mu(q)^2|c_q(n)|}{\varphi(q)^2}\left(\frac1{2\eta_q}-1\right),$$

$$A_E(n)=\sum_{q\le P}\sum_{(a,q)=1}\left(\frac{2|\mu(q)|\sqrt{nD_{q,a}}}{\varphi(q)}+D_{q,a}\right).$$

这些量均有明确的有限积分或有限和定义；有

$$I_M(n)\ge (n-1)\mathfrak S_P(n)-A_{\rm arc}(n)-A_E(n).\tag{3.1}$$

**证明。** 对每一条弧，S²−V²=2VE+E²。Cauchy–Schwarz 给出

$$\left|\int_{-\eta_q}^{\eta_q}(S^2-V^2)e(-n\beta)d\beta\right|\le2\|V\|_2\|E\|_2+\|E\|_2^2\le\frac{2|\mu(q)|\sqrt{nD_{q,a}}}{\varphi(q)}+D_{q,a},$$

因为完整周期上 ∫|T_n|²=n，而子区间上的非负积分至多为 n。

令 J_q=∫_(-η_q)^(η_q)T_n(β)²e(−nβ)dβ。对全部既约中心求和后，V² 项成为 Σ_(q≤P) μ(q)²c_q(−n)J_q/φ(q)²；c_q(−n)=c_q(n) 为实数。完整周期上 ∫T_n²e(−nβ)=n−1。

由几何级数公式和 sin(π|β|)≥2|β|（|β|≤1/2），有 |T_n(β)|≤min(n,1/(2|β|))。由于 η_q<1/2，

$$|J_q-(n-1)|\le2\int_{\eta_q}^{1/2}\frac{d\beta}{4\beta^2}=\frac1{2\eta_q}-1.$$

逐弧求和便得到 (3.1)。证毕。

#### 3.1 全部误差均可还原为有限和

设 b_k=Λ(k)e(ka/q)−μ(q)/φ(q)，K_η(h)=2η（h=0），K_η(h)=sin(2πhη)/(πh)（h≠0）。有限和交换积分得到

$$D_{q,a}=\sum_{k,l=1}^n b_k\overline{b_l}K_{\eta_q}(k-l).\tag{3.2}$$

同理，令 v_n(s)=s−1（2≤s≤n+1），v_n(s)=2n+1−s（n+1<s≤2n），则

$$J_q=\sum_{s=2}^{2n}v_n(s)K_{\eta_q}(s-n).\tag{3.3}$$

因此，主弧估计中的每一个量都能按原定义复算。式 (3.1) 的有效性直接来自代数恒等式、正交性及 Cauchy–Schwarz。

#### 3.2 对零点密度包络与平滑归一化的修正

此前的 0.001882、0.0030 与 0.65 改由 (3.1) 中的显式误差量取代。把高度划分为十段只给出各段误差之和；各段误差的大小需要各自的已证估计。每一个数值系数都应由这些误差量的界推出。

平滑算子的尺度也需完整保留。固定基准 n，令

$$H(t)=\frac1{24(\log n)^5}\sum_{k\le t}(t-k)^4\Lambda(k).$$

在非整数 t 处，逐段四次求导得到 H^(4)(t)=ψ(t)/(log n)^5。将 H 的上界转成 ψ 的上界需要控制四阶导数；归一化因子 (log n)^−5 随复原过程抵消。式 (3.1) 直接处理原始 S_n，避免该转换中的尺度缺项。
