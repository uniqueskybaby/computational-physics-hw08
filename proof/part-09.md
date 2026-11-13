# 任意大于 2 的偶数的二素数表示：主弧估计的完整修订

## 题目与量词

目标命题为

$$\forall n\in\mathbb Z,\quad n>2\ \land\ 2\mid n\quad\Longrightarrow\quad\exists p,q\in\mathbb P,\ n=p+q.$$

本卷处理该证明链中的主弧实部下界。第 1–3 卷的定义、核恒等式及有限证书继续使用。本卷给出已经证明的充分大 n 主弧下界，并精确保留二素数全称命题需要的剩余次弧条件。

### 定理 10：新的主弧参数与无条件实部下界

令 L=log n，改取

$$P_*=L^{100},\qquad Q_*=\frac{n}{L^{202}},\qquad \eta_{*,q}=\frac1{qQ_*}.$$

主弧 M_* 为圆周上所有既约 a/q（q≤P_*）的半宽 η_(＊,q) 邻域，m_* 为其补集。存在固定阈值 N_*，对所有偶数 n≥N_*，有

$$I_{M_*}(n)\ge0.8\,\mathfrak S(n)n.\tag{10.1}$$

证明采用经典 Siegel–Walfisz 定理的存在型版本；本卷逐项推导它与所选 P_*、Q_* 的兼容性及全部余项。N_* 保持该定理所给出的存在量词。

#### 第一步：使用的标准定理

Siegel–Walfisz 定理：对任意固定正数 A、B，存在仅依赖 A、B 的常数 C_(A,B)，使对充分大的 x、q≤(log x)^B、(r,q)=1，

$$\left|\psi(x;q,r)-\frac{x}{\varphi(q)}\right|\le C_{A,B}\frac{x}{(\log x)^A}.\tag{SW}$$

这里 ψ(x;q,r)=Σ_(k≤x,k≡r mod q)Λ(k)。采用的定理版本可见 Terence Tao，254A, Notes 2: Complex-analytic multiplicative number theory，2014 年 12 月 9 日，Exercise 64。其指数误差形式推出 (SW)，因为 exp(−c sqrt(log x)) 乘以任意固定的 log x 幂次都趋于零。

#### 第二步：统一控制全部部分和

取 A=800、B=101。对 sqrt(n)≤t≤n，充分大 n 时有 L^100≤(log t)^101，且 log t≥L/2。因此 (SW) 给出

$$\left|\psi(t;q,r)-\frac{\lfloor t\rfloor}{\varphi(q)}\right|\le C\frac n{L^{800}}\qquad(q\le L^{100}).\tag{10.2}$$

常数吸收 2^800 及取整误差。对 0≤t<sqrt(n)，用初等界 ψ(t;q,r)≤t log n，误差至多 sqrt(n)(L+1)。由于 L^801/sqrt(n)→0，同一个常数在充分大 n 时仍给出 (10.2)。因此它对所有 0≤t≤n 同时成立。

令 F_(q,r)(t)=ψ(t;q,r)−floor(t)/φ(q)。对 e(βt) 做 Stieltjes 分部积分，得到

$$\left|\sum_{k\le n\atop k\equiv r\ (q)}\Lambda(k)e(k\beta)-\frac{T_n(\beta)}{\varphi(q)}\right|\le \sup_{t\le n}|F_{q,r}(t)|(1+2\pi|\beta|n).$$

对全部 φ(q) 个互素剩余类求和，并使用 Σ_(r mod q,(r,q)=1)e(ar/q)=μ(q)，得到主项 μ(q)T_n(β)/φ(q)。非互素类的 Λ 支持位于 p|q 的素数幂上，其总权重至多 ω(q)L≤100L log L/log 2。由于 |β|n≤L^202、φ(q)≤L^100，最终得到统一近似

$$S_n(a/q+\beta)=\frac{\mu(q)}{\varphi(q)}T_n(\beta)+E_{q,a}(\beta),\qquad |E_{q,a}(\beta)|\le C\frac n{L^{498}}.\tag{10.3}$$

这一步的幂次为 800−100−202=498；非互素类的多项式对数项也被右端吸收。

#### 第三步：证明奇异级数截断误差足够小

记 a_q=μ(q)²|c_q(n)|/φ(q)²，A(n)=Σ_(q≥1)a_q。由第 1 卷的欧拉乘积，

$$\frac{A(n)}{\mathfrak S(n)}=\prod_{p\nmid n}\frac{1+(p-1)^{-2}}{1-(p-1)^{-2}}\le C_2^{-2}<2.33.\tag{10.4}$$

使用带 log q 的绝对一阶矩，可得

$$\frac{\sum_q a_q\log q}{A(n)}=\sum_p\frac{a_p\log p}{1+a_p}\le\sum_{p\mid n}\frac{\log p}{p}+\sum_p\frac{\log p}{(p-1)^2}.\tag{10.5}$$

等式来自绝对收敛欧拉乘积逐个素因子展开；p|n 时 a_p/(1+a_p)=1/p。

下面给出这两个和的完整初等控制。首先 θ(x)≤4x 对所有 x≥1 成立：对整数 m，超过 ceil(m/2) 的每个素数都整除 binom(m,floor(m/2))，而该二项式系数≤2^m。因此

$$\vartheta(m)\le\vartheta(\lceil m/2\rceil)+m\log2.$$

从 m=1 开始归纳，偶数情形用 2m+m log 2≤4m，奇数 m≥3 用 2m+2+m log 2≤4m；实数情形取整即可。于是分部求和给出

$$\sum_{p\le L}\frac{\log p}{p}=\frac{\vartheta(L)}L+\int_2^L\frac{\vartheta(t)}{t^2}dt\le4+4\log L.$$

对 p|n、p>L，Σ log p≤log n=L，所以对应和≤1。另一方面，函数 f(x)=log(x+1)/x² 在 x≥1 递减，故

$$\sum_p\frac{\log p}{(p-1)^2}\le\sum_{j\ge1}\frac{\log(j+1)}{j^2}\le\log2+\int_1^\infty\frac{\log(x+1)}{x^2}dx=3\log2<3.$$

结合 (10.5)，得到 Σa_q log q≤A(n)(4log L+8)。对 q>P_*，log q>100log L，因而

$$\left|\mathfrak S(n)-\mathfrak S_{P_*}(n)\right|\le\sum_{q>P_*}a_q\le\frac{A(n)(4\log L+8)}{100\log L}.$$

当 log L≥8 时，最后一项≤0.05A(n)<0.1165S(n)。因此

$$\mathfrak S_{P_*}(n)\ge0.8835\mathfrak S(n).\tag{10.6}$$

此处完整保留了截断误差对 n 的一致性，避免将固定 n 的级数收敛直接用于变化中的 n。

#### 第四步：弧截断误差与近似误差逐项合并

充分大 n 时 Q_*>2P_*，第 1 卷引理 1 的互不相交性证明适用，主弧测度至多 2P_*/Q_*=2L^302/n。

对所有正整数 q，有 φ(q)≥sqrt(q/2)。证明为逐素数幂相乘：奇素数 p 的 p−1≥sqrt(p)；p^a 中的因子 p^(a−1)(p−1)≥p^(a/2)；素数 2 仅在 a=1 时损失 sqrt(2) 因子。于是由 |c_q(n)|≤φ(q)，

$$A_{\rm arc,*}(n)\le\frac{Q_*}{2}\sum_{q\le P_*}q a_q\le\frac{Q_*}{2}\sum_{q\le P_*}\sqrt{2q}\le\frac n{\sqrt2 L^{52}}=o(n).\tag{10.7}$$

用 (10.3) 及 |μ(q)T_n/φ(q)|≤n，主弧总近似误差至多

$$|M_*|\left(2n\sup|E|+(\sup|E|)^2\right)\le Cn\left(L^{-196}+L^{-694}\right)=o(n).\tag{10.8}$$

应用第 1 卷引理 3 的同一分解，再结合 (10.6)–(10.8)，得到

$$I_{M_*}(n)\ge0.8835(n-1)\mathfrak S(n)-o(n).$$

因为 S(n)>1.3114，存在统一阈值 N_* 使右端至少为 0.8S(n)n。这证明 (10.1)。证毕。

### 定理 11：将全称命题归结为精确剩余任务

对 n≥max(N_*,10^12)，第 3 卷已证 R_pow(n)<0.001n。如果进一步建立

$$I_{m_*}(n)\ge-0.79\,\mathfrak S(n)n,\tag{11.1}$$

则由 (10.1) 得到

$$R_{\rm prime}(n)\ge0.01\mathfrak S(n)n-0.001n>0.012114n>0.$$

完整全称证明的剩余任务是：(i) 对全部充分大偶数证明 (11.1) 或更强的次弧下界；(ii) 给出可使用的阈值并核对其以下偶数。本稿的有限证书已覆盖 4≤n≤10000；(11.1) 保持待证状态。

### 引理 12：旧窄扇区零点声明的几何校正

旧稿使用 η=(log n)²/(q²sqrt(H))，q>sqrt(n)，H≥n(log n)^7。由这些条件直接得到

$$0<n\eta<\frac1{\sqrt n(\log n)^{3/2}}\le\frac1{20(\log400)^{3/2}}<\frac1{293}.$$

所以 sin(πnη)>0。该半宽处的边界项按完整数值保留。几何级数核的真实零点要求 nη 为正整数，与这一区间不相交。

### 引理 13：截断投影矩阵的范数下界

对任意可测 A⊂T，定义 G_(j,k)=W_A(k−j)，1≤j,k≤n。则 G 是 Hermitian 半正定矩阵，全部特征值位于 [0,1]，并且

$$\|G\|_{\rm op}\ge\frac{\operatorname{tr}G}{n}=|A|.$$

**证明。** 对 x∈C^n，x*Gx=∫_A|Σx_k e(kα)|²dα，介于 0 与完整周期积分 Σ|x_k|² 之间，给出特征值区间。每个对角元均为 |A|，因此最大特征值至少为平均特征值 |A|。证毕。

对 A=m，得到范数下界 1−|M|。扇区算子的额外小因子需要额外的系数限制或算术抵消，并在 Type II 非对角和中证明。
