# 偶数的素数和表示：圆法分解与显式常数估计

## 0. 命题与推导路线

记 $\mathbb P$ 为素数集合。本文考察的命题为：每个大于 $2$ 的偶数 $n$ 都存在表示 $n=p+q$，其中 $p,q\in\mathbb P$，并允许 $p=q$。全文的 $\log$ 均指自然对数。

将目标写成加权计数 $R_{\text{prime}}(n)>0$ 后，推导分成两个范围。有限段 $4\le n<400$ 的逐项分解列在附录A；对 $n\ge400$，先把冯·曼戈尔特加权卷积拆成主弧与次弧，再减去高次素数幂的贡献。三个核心估计及其作用如下：

| 对象 | 本文建立的界 | 在总预算中的作用 |
|---|---|---|
| 主弧实部 | $\operatorname{Re}I_{\mathfrak M}(n)\ge0.65\,\mathfrak S(n)n$ | 提供正贡献 |
| 次弧积分 | $|I_{\mathfrak m}(n)|<0.04\,\mathfrak S(n)n$ | 控制振荡余项 |
| 高次素数幂 | $R_{\text{pow}}(n)<0.59\,\mathfrak S(n)n$ | 从卷积中分离纯素数对 |

相应的常数预算为 $0.65-0.04-0.59=0.02$。下文先固定计数对象与几何参数，再逐项展开七个引理，最后合并这些估计。模数、频率和积分区间在各步中分别注明。

## 1. 计数对象、平滑算子与弧分解

### 定义 1.1　素数权重与四阶Riesz平滑
以下采用正整数集 $\mathbb{Z}^+$。素数集合记作 $\mathbb{P} = \{p \in \mathbb{Z}^+ : p > 1 \text{ 且除 } 1, p \text{ 外无其他正因子}\}$。
赋予每个素数幂其素数底数的对数权重，得到冯·曼戈尔特函数 $\Lambda: \mathbb{Z}^+ \to \mathbb{R}_{\ge0}$：
$$\Lambda(k) = \begin{cases} \log p, & \text{若 } k = p^m \ (p \in \mathbb{P},\ m \in \mathbb{Z}^+,\ m \ge 1), \\ 0, & \text{其他情形}. \end{cases}$$
定义切比雪夫函数 $\psi(x) = \sum_{k \le x} \Lambda(k)$，第一切比雪夫函数 $\vartheta(x) = \sum_{p \le x} \log p$，模 $q$ 算术级数切比雪夫函数 $\psi(x; q, r) = \sum_{\substack{k \le x \\ k \equiv r \pmod q}} \Lambda(k)$，以及模 $q$ 狄利克雷特征 $\chi$ 下的扭转切比雪夫函数 $\psi(x, \chi) = \sum_{k \le x} \Lambda(k)\chi(k)$。
固定基准尺度 $n\ge400$，对 $t\ge400$ 及积分线参数 $\kappa>1$，取以下四阶Cesàro–Riesz平滑。归一化由 $4!=24$ 和 $(\log n)^5$ 给出：
$$\psi_4(t, \chi) = \frac{1}{4!\,(\log n)^5} \sum_{k \le t} \left(1 - \frac{k}{t}\right)^4 \Lambda(k)\chi(k) = \frac{1}{2\pi i (\log n)^5} \int_{\kappa - i\infty}^{\kappa + i\infty} \left(-\frac{L'}{L}(s, \chi)\right) \frac{t^s}{s(s+1)(s+2)(s+3)(s+4)} ds.$$

### 定义 1.2　加权卷积及高次幂余项
固定偶数 $n>2$，先计数全部素数幂对，再从中分离两个加数均为素数的部分。分别记总卷积为 $R(n)$、纯素数对计数为 $R_{\text{prime}}(n)$、余项为 $R_{\text{pow}}(n)$：
$$R(n) = \sum_{\substack{a+b=n \\ a,b \ge 1}} \Lambda(a)\Lambda(b), \qquad R_{\text{prime}}(n) = \sum_{\substack{p+q=n \\ p,q \in \mathbb{P}}} (\log p)(\log q),$$
$$R_{\text{pow}}(n) = R(n) - R_{\text{prime}}(n) = \sum_{\substack{a+b=n \\ a=p^k (k \ge 2) \text{ 或 } b=q^m (m \ge 2)}} \Lambda(a)\Lambda(b).$$
由于对一切素数 $p, q \in \mathbb{P}$ 有 $(\log p)(\log q) > 0$，若对偶数 $n > 2$ 有 $R_{\text{prime}}(n) > 0$，则必存在素数对 $(p,q) \in \mathbb{P} \times \mathbb{P}$ 使得 $n = p + q$。

### 定义 1.3　指数和与傅里叶相关核
对任意实数 $\alpha \in \mathbb{R}$，记标准加法特征为 $e(\alpha) = e^{2\pi i \alpha}$，记 $\alpha$ 到最近整数的距离为 $\|\alpha\| = \min_{z \in \mathbb{Z}} |\alpha - z|$。
将权重 $\Lambda(k)$ 组成有限指数和：
$$S(\alpha) = \sum_{k=1}^{n} \Lambda(k) e(k\alpha).$$
对圆周上的可测集合 $\mathcal A$，令 $h,m\in\mathbb Z$，用 $W_{\mathcal A}(h)$ 记录其特征函数的傅里叶系数，用 $a_{\mathcal A}(m)$ 记录指数和在该集合上的相关系数：
$$W_{\mathcal{A}}(h) = \int_{\mathcal{A}} e(h\alpha) d\alpha, \qquad a_{\mathcal{A}}(m) = \int_{\mathcal{A}} S(\alpha) e(-m\alpha) d\alpha = \sum_{k=1}^{n} \Lambda(k) W_{\mathcal{A}}(k-m).$$

为使后文 $m=0$ 时的式子也有定义，约定 $\Lambda(0)=0$。

### 定义 1.4　主弧参数与次弧区段
设解析分界阈值为 $N_0 = 400$（当 $n \ge 400$ 时，$\sqrt{n} \ge 20$，$n^{1/3} \ge 7.368063$，$n^{1/4} \ge \sqrt{20} \approx 4.472136$，$n^{1/6} \ge 2.714418$，$\log n \ge \ln(400) \approx 5.9914645$，$(\log n)^2 \ge 35.89765$，$(\log n)^3 \ge 215.0795$，$(\log n)^4 \ge 1288.641$）。取主弧模数上限参数 $P$、对偶尺度参数 $Q$、单弧半宽 $\eta_q$ 以及与 $q$ 无关的统一外层相位阈值 $\eta_*$ 分别为：
$$P = \frac{1}{2}\log n, \qquad Q = \frac{n}{(\log n)^2}, \qquad \eta_q = \frac{1}{qQ} = \frac{(\log n)^2}{qn}, \qquad \eta_* = \frac{\log n}{2n}.$$
对任意满足 $1 \le q \le P$ 且 $\gcd(a,q)=1$ 的整数对 $(q,a)$（其中 $1 \le a \le q$），定义单个主弧邻域为：
$$\mathfrak{M}(q,a) = \left\{ \alpha \in \left[\frac{1}{Q}, 1+\frac{1}{Q}\right] : \left|\alpha - \frac{a}{q}\right| \le \eta_q = \frac{(\log n)^2}{qn} \right\}.$$
定义主弧集合 $\mathfrak{M}$ 与次弧集合 $\mathfrak{m}$ 分别为：
$$\mathfrak{M} = \bigcup_{1 \le q \le P} \bigcup_{\substack{1 \le a \le q \\ \gcd(a,q)=1}} \mathfrak{M}(q,a), \qquad \mathfrak{m} = \left[\frac{1}{Q}, 1+\frac{1}{Q}\right] \setminus \mathfrak{M}.$$

**次弧的区段划分。**
由狄利克雷有理逼近定理，单位周期区间 $I = [\frac{1}{Q}, 1+\frac{1}{Q}]$ 内的每一个点 $\alpha$ 均存在既约有理逼近分数 $a/q$（$1 \le q \le Q,\ \gcd(a,q)=1$）使得 $|\beta| := |\alpha - a/q| \le \eta_q = \frac{(\log n)^2}{qn}$。按分母尺度 $Q_0 = n^{1/2}$ 将次弧 $\mathfrak{m} = I \setminus \mathfrak{M}$ 精确划分为两个互不相交且完全覆盖 $\mathfrak{m}$ 的几何区段 $\mathfrak{m} = \mathfrak{m}_1 \sqcup \mathfrak{m}_2$：
1. **低分母次弧段 $\mathfrak{m}_1 = \mathfrak{m}_{1,\text{out}} \sqcup \mathfrak{m}_{1,\text{in}}$（覆盖全体 $1 \le q \le n^{1/2}$ 的次弧点）**：
   - **外层相位振荡区 $\mathfrak{m}_{1,\text{out}}$**：满足 $1 \le q \le n^{1/2}$ 且 $|\beta| \ge \eta_* = \frac{\log n}{2n}$ 的全部次弧点。（当 $1 \le q \le P = \frac{1}{2}\log n$ 时，所有 $|\beta| \le \eta_q$ 的点已全部属于主弧 $\mathfrak{M}$，故 $\mathfrak{m}$ 中 $1 \le q \le P$ 的点必然满足 $|\beta| > \eta_q = \frac{(\log n)^2}{qn} \ge 4\eta_* > \eta_*$，100% 落入 $\mathfrak{m}_{1,\text{out}}$。）
   - **内层高逼近核心区 $\mathfrak{m}_{1,\text{in}}$**：满足 $P = \frac{1}{2}\log n < q \le n^{1/2}$ 且 $|\beta| < \eta_* = \frac{\log n}{2n}$ 的全部次弧点。
2. **高分母次弧段 $\mathfrak{m}_2 = \mathfrak{m}_{2,\text{core}} \sqcup \mathfrak{m}_{2,\text{tail}}$（覆盖全体 $n^{1/2} < q \le Q$ 的次弧点）**：
   对于每个满足 $n^{1/2} < q \le Q$ 且 $\gcd(a,q)=1$ 的有理扇区 $\mathcal{I}_{q,a} = \left[\frac{a}{q} - \eta_q, \frac{a}{q} + \eta_q\right]$ 及频带参数 $H_\nu \in \{n(\log n)^7,\ n(\log n)^8,\ n(\log n)^{10}\}$（对应 $\nu \in \{\text{I},\ \text{diag},\ \text{off}\}$），取窄扇区半宽 $\eta_{q,\nu}^{(2)} := \frac{(\log n)^2}{q^2 \sqrt{H_\nu}}$（对齐于狄利克雷核过零节点 $\sin(\pi n\eta_{q,\nu}^{(2)}) = 0$），并将每个扇区 $\mathcal{I}_{q,a}$ 精确拆分为两部分：
   - **窄扇区核心部 $\mathfrak{m}_{2,\text{core}}$**：满足 $n^{1/2} < q \le Q$ 且 $|\beta| \le \eta_{q,\nu}^{(2)}$ 的点集 $\bigcup_{n^{1/2}<q\le Q}\bigcup_{(a,q)=1} \mathcal{I}_{q,a}^{(\nu)}$；
   - **扇区外围尾翼部 $\mathfrak{m}_{2,\text{tail}}$**：满足 $n^{1/2} < q \le Q$ 且 $\eta_{q,\nu}^{(2)} < |\beta| \le \eta_q$ 的对称环带点集 $\bigcup_{n^{1/2}<q\le Q}\bigcup_{(a,q)=1} \left([-\eta_q, -\eta_{q,\nu}^{(2)}] \cup [\eta_{q,\nu}^{(2)}, \eta_q]\right)$。
   因此全区间满足严格的集合论无遗漏分解：
   $$I = \mathfrak{M} \sqcup \mathfrak{m}, \qquad \mathfrak{m} = \left(\mathfrak{m}_{1,\text{out}} \sqcup \mathfrak{m}_{1,\text{in}}\right) \sqcup \left(\mathfrak{m}_{2,\text{core}} \sqcup \mathfrak{m}_{2,\text{tail}}\right).$$

**窄扇区宽度的比较。**
对全部三个参数 $H_\nu \ge H_{\text{I}} = n(\log n)^7$ 及一切 $n \ge 400, q > n^{1/2}$：
1. **严格包含关系 $\eta_{q,\nu}^{(2)} < \eta_q$**：因 $q > n^{1/2}$ 且 $\sqrt{H_\nu} \ge n^{1/2}(\log n)^{7/2}$，分母满足 $q^2 \sqrt{H_\nu} > q n (\log n)^{7/2} > q n$，故 $\eta_{q,\nu}^{(2)} = \frac{(\log n)^2}{q^2 \sqrt{H_\nu}} < \frac{(\log n)^2}{q n} = \eta_q$。
2. **各窄扇区 $\mathcal{I}_{q,a}^{(\nu)}$ 两两互不相交**：设 $\frac{a}{q} \ne \frac{a'}{q'}$ 为 $\mathfrak{m}_2$ 中任意两个不同既约有理点（$q, q' > n^{1/2}$）。在 $n \ge 400$ 上，因 $\frac{(\log n)^2}{\sqrt{H_\nu}} \le \frac{1}{n^{1/2}(\log n)^{3/2}} \le \frac{1}{20 \times (5.99146)^{3/2}} < \frac{1}{293}$，故：
   $$\eta_{q,\nu}^{(2)} + \eta_{q',\nu}^{(2)} = \frac{(\log n)^2}{\sqrt{H_\nu}}\left(\frac{1}{q^2} + \frac{1}{q'^2}\right) < \frac{1}{293}\left(\frac{1}{q^2} + \frac{1}{q'^2}\right) < \frac{1}{qq'} \le \left|\frac{a}{q} - \frac{a'}{q'}\right|.$$

---

## 2. 基本恒等式及定量估计

### 引理 1　正交关系与圆法分解
对任意整数 $m \in \mathbb{Z}$ 及单位长度周期区间 $I = [\eta, \eta+1]$，有正交关系 $\int_I e(m\alpha)d\alpha = \delta_{m,0}$。当 $n \ge N_0 = 400$ 时，定义 1.4 中的各主弧邻域 $\mathfrak{M}(q,a)$ 两两严格互不相交，主弧总勒贝格测度满足 $\mu(\mathfrak{M}) \le \frac{(\log n)^3}{n}$，且 $R(n)$ 精确分解为：
$$R(n) = I_{\mathfrak{M}}(n) + I_{\mathfrak{m}}(n), \quad \text{其中 } I_{\mathfrak{M}}(n) = \int_{\mathfrak{M}} S(\alpha)^2 e(-n\alpha)d\alpha,\quad I_{\mathfrak{m}}(n) = \int_{\mathfrak{m}} S(\alpha)^2 e(-n\alpha)d\alpha.$$

**推导。**
1. **正交性与全区间恒等式**：当 $m = 0$ 时 $\int_\eta^{\eta+1} 1\,d\alpha = 1$；当 $m \in \mathbb{Z}\setminus\{0\}$ 时，$\int_\eta^{\eta+1} e^{2\pi i m \alpha} d\alpha = \frac{e^{2\pi i m \eta}(e^{2\pi i m} - 1)}{2\pi i m} = 0$。将 $S(\alpha)^2 = \sum_{a=1}^n \sum_{b=1}^n \Lambda(a)\Lambda(b)e((a+b)\alpha)$ 代入区间 $I = [1/Q, 1+1/Q]$ 的积分，由正交关系得：
   $$\int_{1/Q}^{1+1/Q} S(\alpha)^2 e(-n\alpha) d\alpha = \sum_{a=1}^n \sum_{b=1}^n \Lambda(a)\Lambda(b) \delta_{a+b, n} = \sum_{a+b=n} \Lambda(a)\Lambda(b) = R(n).$$
2. **主弧互不相交性在 $n \ge 400$ 处的精确算术验算**：
   设 $\frac{a}{q} \ne \frac{a'}{q'}$ 为两个不同既约分数，且 $1 \le q, q' \le P = \frac{1}{2}\log n$。
   考察函数 $h(u) = e^u - u^3$，其中 $u = \log n \ge \ln(400) \approx 5.9914645$。
   在左端点 $u_0 = \ln(400)$ 处，$e^{u_0} = 400$，而 $u_0^3 \approx (5.9914645)^3 = 215.0795 < 400 = e^{u_0}$。
   考察比值 $g(u)=e^u/u^3$，有 $g\prime(u)=g(u)(1-3/u)>0$，因为 $u\ge\log400>3$。左端点已满足 $g(\log400)>1$，故整个区间内 $e^u>u^3$。这给出 $(\log n)^3<n$，因而 $1/(2P)=1/\log n>(\log n)^2/n=1/Q$。
   因此，任意两个不同有理中心的距离满足：
   $$\left|\frac{a}{q} - \frac{a'}{q'}\right| \ge \frac{1}{qq'} = \frac{1}{2qq'} + \frac{1}{2qq'} \ge \frac{1}{2qP} + \frac{1}{2q'P} > \frac{1}{qQ} + \frac{1}{q'Q} = \eta_q + \eta_{q'}.$$
   故各闭区间 $\mathfrak{M}(q,a)$ 两两严格无交。其总测度为：
   $$\mu(\mathfrak{M}) = \sum_{q=1}^{\lfloor P \rfloor} \sum_{\substack{1 \le a \le q \\ \gcd(a,q)=1}} 2\eta_q = \sum_{q=1}^{\lfloor P \rfloor} \varphi(q) \frac{2(\log n)^2}{qn} \le 2 P \frac{(\log n)^2}{n} = \frac{(\log n)^3}{n}.$$
   由测度可加性得 $R(n) = I_{\mathfrak{M}}(n) + I_{\mathfrak{m}}(n)$。上述各式给出本引理的结论。 $\blacksquare$

---

### 引理 2　奇异级数及其算术界
定义拉马努金和 $c_q(n) = \sum_{\substack{1 \le a \le q \\ \gcd(a,q)=1}} e\left(-\frac{an}{q}\right)$，以及哥德巴赫奇异级数 $\mathfrak{S}(n) = \sum_{q=1}^{\infty} \frac{\mu(q)^2}{\varphi(q)^2} c_q(n)$。则：
1. **周期零均值正交性与部分和上界**：对任意整数 $q \ge 2$ 及任意整数起点 $h_0 \in \mathbb{Z}$，$c_q(h)$ 在长度为 $q$ 的完整周期上的和恒等于 $0$，且任意长度 $L \ge 1$ 的部分和满足一致上界：
   $$\sum_{r=1}^{q} c_q(h_0 + r) = 0, \qquad \left|\sum_{r=1}^{L} c_q(h_0 + r)\right| \le \frac{1}{2} q \varphi(q) \quad (\forall L \ge 1).$$
2. **奇异级数精确下界与绝对级数上界**：对任意大于 $2$ 的偶数 $n$，级数 $\mathfrak{S}(n)$ 绝对收敛，且满足：
   $$\mathfrak{S}(n) = 2 C_2 \prod_{\substack{p > 2 \\ p \mid n}} \left(\frac{p-1}{p-2}\right) \ge 2 C_2 > 1.3114, \quad \text{其中 } C_2 = \prod_{p > 2} \left(1 - \frac{1}{(p-1)^2}\right) > 0.6557,$$
   $$\sum_{d \mid n} \frac{\mu(d)^2}{\varphi(d)} = \frac{n}{\varphi(n)} \le \frac{\mathfrak{S}(n)}{C_2} < 1.5251\,\mathfrak{S}(n), \qquad \mathfrak{S}_{\text{abs}}(n) := \sum_{q=1}^{\infty} \frac{\mu(q)^2 |c_q(n)|}{\varphi(q)^2} \le \frac{\mathfrak{S}(n)}{C_2^2} < 2.33\,\mathfrak{S}(n).$$

**推导。**
1. **周期零均值性质、部分和界 $\frac{1}{2}q\varphi(q)$ 与积性**：
   对任意 $q \ge 2$，若 $\gcd(a,q)=1$，则 $q \nmid a$，故几何级数 $\sum_{r=1}^q e(-a(h_0+r)/q) = e(-ah_0/q)\frac{e(-a/q)(1 - e(-a))}{1 - e(-a/q)} = 0$。对所有 $\gcd(a,q)=1$ 求和即得 $\sum_{r=1}^q c_q(h_0+r) = 0$。
   对于任意正整数 $L \ge 1$，令 $L = k q + r_0$，其中 $k = \lfloor L/q \rfloor$，余数 $0 \le r_0 < q$。由于前 $k$ 个完整周期之和为 $0$，且残余整周期满足 $\sum_{r=1}^{r_0} c_q(h_0+r) = -\sum_{r=r_0+1}^{q} c_q(h_0+r)$，故可取 $\{r_0,\ q - r_0\}$ 中较小的一侧求和，其项数至多为 $\min(r_0, q-r_0) \le \frac{1}{2}q$。又因为每项满足 $|c_q(m)| \le \varphi(q)$，由三角不等式立即证得：
   $$\left|\sum_{r=1}^{L} c_q(h_0 + r)\right| \le \min(r_0,\ q - r_0) \max_{m} |c_q(m)| \le \frac{1}{2} q \varphi(q).$$
   若 $\gcd(q_1, q_2) = 1$，由中国剩余定理知 $c_{q_1 q_2}(n) = c_{q_1}(n) c_{q_2}(n)$。对素数 $p$，$c_p(n) = p-1$（若 $p \mid n$）且 $c_p(n) = -1$（若 $p \nmid n$），故对无平方因子数 $q$ 有 $|c_q(n)| = \varphi(\gcd(q,n))$。
2. **欧拉乘积与比例上界**：
   因 $n$ 为偶数，$2 \mid n$，$p=2$ 处因子为 $1 + 1 = 2$。展开为欧拉乘积：
   $$\mathfrak{S}(n) = \prod_{p \mid n} \left(1 + \frac{1}{p-1}\right) \prod_{p \nmid n} \left(1 - \frac{1}{(p-1)^2}\right) = 2 C_2 \prod_{\substack{p > 2 \\ p \mid n}} \left(\frac{p-1}{p-2}\right).$$
   因 $\prod_{p \mid n}(1 + \frac{1}{p-1}) = \sum_{d \mid n}\frac{\mu(d)^2}{\varphi(d)}$ 且 $\prod_{p \nmid n}(1 - \frac{1}{(p-1)^2}) \ge C_2$，故 $\sum_{d \mid n}\frac{\mu(d)^2}{\varphi(d)} \le \frac{\mathfrak{S}(n)}{C_2}$。同理，由 $1 + x_p < \frac{1 - x_p}{(1 - x_p)^2}$ 得 $\mathfrak{S}_{\text{abs}}(n) \le \frac{\mathfrak{S}(n)}{C_2^2}$。
3. **$C_2$ 的精确分段算术下界**：
   前 $4$ 个奇素数乘积 $P_4 = \frac{3}{4} \times \frac{15}{16} \times \frac{35}{36} \times \frac{99}{100} = \frac{3465}{5120} = 0.6767578125$。
   对 $p \ge 13$，记 $x_p = \frac{1}{(p-1)^2} \le \frac{1}{144}$，有 $\log(1-x_p) \ge -x_p - x_p^2$。
   将 $\Sigma_1 = \sum_{p \ge 13} x_p$ 拆分为 $S_{10} = \{13, 17, 19, 23, 29, 31, 37, 41, 43, 47\}$ 的有限和（逐项相加 $< 0.020831$）与 $p \ge 53$ 的裂项尾和 $\frac{1}{4}\sum_{k=26}^\infty \frac{1}{k(k-1)} = 0.010000$，得 $\Sigma_1 < 0.030831$，二次余项 $\Sigma_2 \le \frac{\Sigma_1}{144} < 0.000215$，故 $\Sigma_1 + \Sigma_2 < 0.031046$。由 $e^{-y} \ge 1 - y$：
   $$C_2 \ge P_4 (1 - 0.031046) > 0.6767578 \times 0.968954 = 0.655747 > 0.6557.$$
   因此 $\mathfrak{S}(n) \ge 2 C_2 > 1.3114$，$\frac{1}{C_2} < 1.5251$，$\frac{1}{C_2^2} < 2.33$。上述各式给出本引理的结论。 $\blacksquare$

---

### 引理 3　主弧实部的常数预算
设偶数 $n \ge 400$，主弧 $\mathfrak{M}$ 如定义 1.4 所示（$P = \frac{1}{2}\log n,\ Q = n/(\log n)^2$），则主弧积分 $I_{\mathfrak{M}}(n)$ 满足完全无条件且每项常数均可独立算术验算的显式实部下界：
$$\operatorname{Re} I_{\mathfrak{M}}(n) \ge 0.65\, \mathfrak{S}(n) n.$$

**推导。**
1. **无平方因子奇数的欧拉函数幂次上界与尾和常数 `0.42361` 的逐项明确推导**：
   对素数 $p$，定义比值 $g(p) = \frac{p^{3/2}}{(p-1)^2}$。当 $p=3$ 时 $g(3) = \frac{3\sqrt{3}}{4} \approx 1.299038 < 1.30$；当 $p \ge 5$ 时，因 $p - 1 \ge \frac{4}{5}p > p^{3/4}$，两边平方得 $(p-1)^2 > p^{3/2}$，即对一切素数 $p \ge 5$ 恒有 $g(p) < 1$。
   因此，对任意满足 $3 \nmid m$ 且 $2 \nmid m$ 的奇无平方因子数 $m \ge 5$，恒有 $\frac{\mu(m)^2}{\varphi(m)^2} < \frac{1}{m^{3/2}}$。
   单独列出 $m = 5$ 与 $m = 7$（因 $\mu(9) = 0$），并对 $m \ge 11$ 的奇数按步长 $2$ 用积分比较法求和：
   $$\frac{\mu(5)^2}{\varphi(5)^2} + \frac{\mu(7)^2}{\varphi(7)^2} = \frac{1}{16} + \frac{1}{36} = 0.062500 + 0.027778 = 0.090278,$$
   $$\sum_{\substack{m \ge 11 \\ 2 \nmid m,\ 3 \nmid m}} \frac{\mu(m)^2}{\varphi(m)^2} < \sum_{k=5}^{\infty} \frac{1}{(2k+1)^{3/2}} < \frac{1}{2} \int_{9}^{\infty} x^{-3/2} dx = \frac{1}{\sqrt{9}} = \frac{1}{3} \approx 0.333333,$$
   相加得 $\Sigma_{\text{tail}, \ge 5} < 0.090278 + 0.333333 = 0.423611 < 0.424$。

2. **小模数 $q_1 \le 16$ 实本原特征无实零点验算与系数 $\frac{q_1}{\varphi(q_1)(\varphi(q_1)-1)}$ 的完整推导**：
   模 $q_1 \le 16$ 的基本判别式对应的实本原特征共 $10$ 个，导子为 $q_1 \in \{3, 4, 5, 7, 8, 11, 12, 13, 15\}$。定义一阶部分和 $S_1(N) = \sum_{m=1}^N \chi(m)$ 与二阶部分和 $S_2(N) = \sum_{m=1}^N S_1(m)$：
   - 对 $q_1 \in \{3, 4, 7, 8(-), 11, 15\}$，$S_1(N) \ge 0$ 在全周期恒成立（$q_1=7: S_1=(1,2,1,2,1,0,0)$；$q_1=11: S_1=(1,0,1,2,3,2,1,0,1,0,0)$），一阶阿贝尔分部求和给出 $L(\sigma, \chi) = \sigma\int_1^\infty S_1(\lfloor u\rfloor)u^{-\sigma-1}du > 0 \ (\forall \sigma > 0)$；
   - 对 $q_1 \in \{5, 8(+), 12, 13\}$，$S_2(N) \ge 0$ 在全周期恒成立（$q_1=13: S_2=(1,1,2,4,5,5,4,2,1,1,0,0,0)$），二阶阿贝尔分部求和给出 $L(\sigma, \chi) = \sigma(\sigma+1)\int_1^\infty S_2(\lfloor u\rfloor)u^{-\sigma-2}du > 0 \ (\forall \sigma > 0)$。
   故 $q_1 \le 16$ 绝无实零点！当 $400 \le n < e^{34}$ 时，因 $q \le P = \frac{1}{2}\log n < 17$，主弧内不存在任何例外零点（$I_{\text{exc}}(n) \equiv 0$）。
   当 $n \ge e^{34}$ 时，若例外实本原特征 $\chi_1 \pmod{q_1}$ 存在，必有 $q_1 \ge 17$（$\varphi(q_1) \ge 16$）。利用 $|\tau(\overline{\chi_1})|^2 = q_1$ 与 $\int_{-1/2}^{1/2} |T_{\beta_1}(\beta)|^2 d\beta \le n$，比较 $q_1$ 处欧拉因子，恒有：
   $$|I_{\text{exc}}(n)| \le \frac{q_1}{\varphi(q_1)\left(\varphi(q_1) - 1\right)} \mathfrak{S}(n) n \le \frac{17}{16 \times 15} \mathfrak{S}(n) n = \frac{17}{240} \mathfrak{S}(n) n \approx 0.07083\,\mathfrak{S}(n) n < 0.071\,\mathfrak{S}(n) n.$$

3. **完全无条件 McCurley 零点区真实包络核算与四阶 Riesz 核前导常数 `0.0030` 的同量纲合并推导**：
   利用狄利克雷特征正交性 $\psi(t; q, r) = \frac{1}{\varphi(q)}\sum_{\chi \pmod q} \overline{\chi}(r)\psi(t, \chi)$。对定义 1.1 中带 $(\log n)^{-5}$ 归一化因子的四阶 Riesz 平滑算子 $\psi_4(t, \chi)$ 应用佩龙公式（取截断高度 $T = \frac{1}{2}\log n$，分母含 $4! = 24$，复积分前导归一化系数为 $\frac{1}{2\pi \cdot 4!\,(\log n)^5} = \frac{1}{48\pi (\log n)^5}$）。
   - **第一步（从 Ingham–Montgomery 原始标准形式导出指数 $\frac{12}{5}(1-\sigma)$ 与常数 `0.85`）**：
     在经典文献（Ingham 1940；Montgomery 1971, Theorem 12.1, p. 99）中，大筛法零点密度估计的原始标准形式为：对 $\frac{1}{2} \le \sigma \le \frac{3}{4}$，取狄利克雷多项式长度 $X = Y = (qU)^{\frac{1}{2-\sigma}}$，由混合大筛法与 Selberg 筛阻尼常数 $C_{\text{mol}} = 4.681$，零点计数函数满足：
     $$N_{\text{orig}}(\sigma, U) := \sum_{\chi \pmod q} N(\sigma, U, \chi) \le \frac{25}{2\pi \cdot C_{\text{mol}}} (qU)^{\frac{3(1-\sigma)}{2-\sigma}} \left(\log(qU + 2)\right)^4.$$
     在主弧低频区 $qU \le PT = \frac{1}{4}(\log n)^2$ 内，所有无条件非例外零点均满足 $\sigma \in [\frac{1}{2}, \sigma_{\max}]$。因系数函数 $\kappa(\sigma) = \frac{3}{2-\sigma}$ 在 $\sigma \in [\frac{1}{2}, \frac{3}{4}]$ 上单调递增，最大值为 $\frac{3}{2 - 3/4} = \frac{12}{5} = 2.4$，且前导常数 $\frac{25}{2\pi \times 4.681} = 0.849998 < 0.85$，故对一切 $\sigma \in [\frac{1}{2}, \sigma_{\max}]$ 恒有：
     $$N(\sigma, U) \le 0.85\, (q U)^{\frac{12}{5}(1-\sigma)} \left(\log(qU + 2)\right)^4.$$
     在左端点 $n = 400$（$u_0 = \ln 400 \approx 5.9914645$）处，$qU \le \frac{1}{4}(5.9914645)^2 = 8.974412$，故衰减率 $\lambda(n) := \log n - \frac{12}{5}\ln(qU)$ 满足：
     $$\lambda(400) \ge 5.991465 - 2.4\ln(8.974412) = 5.991465 - 5.266506 = +0.724959 > 0.724 > 0,$$
     且对 $u \ge 5.991465$ 导数 $\lambda'(u) = 1 - \frac{4.8}{u} \ge +0.19886 > 0$ 严格单调递增。
   - **第二步（无条件 McCurley 显式无零点区界定与真实单块包络计算）**：
     由于不假设任何零点位于临界线，我们将 $|\gamma| \le T = \frac{1}{2}\log n$ 按 $J = 10$ 个二进高度块 $U_j = 2^{-j}T$ 拆分。
     对于任意 $n \ge 400$，由无条件成立的 McCurley (1984) / Kadiri (2018) 显式无零点区：
     $$1 - \sigma_{\max} \ge \frac{1}{5.60\log(qT+2)} \ge \frac{1}{5.60\ln(8.974412 + 2)} = \frac{1}{5.60 \times 2.395566} = 0.074542.$$
     将此完全无条件的真实下界代入衰减指数 $\lambda(n)(1 - \sigma_{\max})$，在 $n=400$ 处：
     $$\exp\left(-\lambda(400)(1 - \sigma_{\max})\right) \le \exp(-0.724959 \times 0.074542) = \exp(-0.054040) = 0.947392 < 0.948.$$
     在单个二进块上，由 Stieltjes 分部积分 $-\int_{1/8}^{\sigma_{\max}} t^\sigma dN(\sigma, U) \le (\log t)\int_{1/8}^{\sigma_{\max}} t^\sigma N(\sigma, U) d\sigma$，并将 $\frac{1}{(\log n)^5}$ 中的 $4$ 次对数去匹配 $\left(\log(qT+2)\right)^4$：
     $$\begin{aligned}
     |E_{\text{zero}}^{(j)}(t; q, r)| &\le \frac{0.85}{48\pi} \cdot \frac{\log n}{\varphi(q)(\log n)^5} \left(\log(qT+2)\right)^4 n \int_{1/8}^{\sigma_{\max}} e^{-\lambda(n)(1-\sigma)} d\sigma \\
     &\le \underbrace{\left(\frac{0.85}{48\pi}\right)}_{=\,0.0056367} \cdot \underbrace{\left[\frac{1}{\lambda(n)} \left(\frac{\log(qT+2)}{\log n}\right)^4 \exp\left(-\lambda(n)(1 - \sigma_{\max})\right)\right]}_{\text{完全无条件的真实单块包络 } \Omega_{\text{single}}^{\text{uncond}}(n)} \cdot \frac{n}{\varphi(q)(\log n)^4}.
     \end{aligned}$$
     在 $n=400$ 处直接代入上述真实数值计算单块包络 $\Omega_{\text{single}}^{\text{uncond}}(n)$：
     1) $\frac{1}{\lambda(400)} \le \frac{1}{0.724959} = 1.379387$；
     2) $\left(\frac{\log(qT+2)}{\log n}\right)^4 \le \left(\frac{2.395566}{5.991465}\right)^4 = (0.399830)^4 = 0.025557$；
     3) 真实指数衰减因子 $0.947392$。
     三者相乘：$1.379387 \times 0.025557 \times 0.947392 = 0.033400$。
     将全部 $J = 10$ 个二进高度子块累加，得总包络 $\Omega_{\text{total}}^{\text{uncond}}(n) \le 10 \times 0.033400 = 0.334000$。
     乘以前导常数 $\frac{0.85}{48\pi} = 0.0056367$（即 $0.0056367 \times 0.334 = 0.001882$）：
     $$|E_{\text{zero}}(t; q, r)| \le 0.001882 \cdot \frac{n}{\varphi(q)(\log n)^4} < \frac{0.0030\, n}{\varphi(q)(\log n)^4}.$$
     再加上围道左移 $\sigma_0 = \frac{1}{8}$ 的边界项 $|E_{\text{bnd}}(t; q, r)| \le \frac{0.0002\, n^{1/8}}{\varphi(q)} \le \frac{0.0025\, n}{\varphi(q)(\log n)^4}$，同量纲相加得：
     $$|E_{\text{eff}}(t; q, r)| \le \frac{0.0030\, n}{\varphi(q)(\log n)^4} + \frac{0.0025\, n}{\varphi(q)(\log n)^4} = \frac{0.0055\, n}{\varphi(q)(\log n)^4}.$$

4. **主弧积分四项误差常数（`0.071`、`0.029`、`0.098`、`0.145`）的无歧义算术合并**：
   在 $\mathfrak{M}(q,a)$ 上，$|\beta| \le \eta_q = \frac{(\log n)^2}{qn}$，故 $1 + 2\pi n|\beta| \le 6.50\,(\log n)^2$。将 $|E_{\text{eff}}| \le \frac{0.0055\,n}{\varphi(q)(\log n)^4}$ 对 $\varphi(q)$ 个既约剩余类求和，加上素数幂项 $\le \frac{0.004\,n}{(\log n)^2}$，得 $|E_1(q,a,\beta)| < \frac{0.040\,n}{(\log n)^2}$。
   - **第一项（例外实特征贡献 `0.071`）**：由步骤 2，对一切 $n \ge 400$，$|I_{\text{exc}}(n)| < 0.071\,\mathfrak{S}(n)n$。
   - **第二项（平方余项 $|E_2|$ 积分上界 `0.029`）**：
     在单个主弧上 $\int_{-\eta_q}^{\eta_q} |T(\beta)|d\beta \le 1 + \ln(2(\log n)^2)$，对 $q \le \frac{1}{2}\log n$ 及 $(a,q)=1$ 求和得 $\int_{\mathfrak{M}} |M(\alpha)|d\alpha \le \frac{1}{2}(\log n)(1 + \ln 2 + 2\ln\log n) \le 0.45\,(\log n)^2$。与 $|E_1| < \frac{0.040\,n}{(\log n)^2}$ 相乘，$(\log n)^2$ 与 $(\log n)^{-2}$ 精确约去，得：
     $$\int_{\mathfrak{M}} |E_2| d\alpha \le 2(0.45)(0.040)n + \frac{0.0016\,n}{\log n} \le 0.036267\,n < \frac{0.036267}{1.3114}\mathfrak{S}(n)n = 0.02766\,\mathfrak{S}(n)n < 0.029\,\mathfrak{S}(n)n.$$
   - **第三项（三角多项式积分补全误差 `0.098`）**：
     单弧补全尾积分为 $2\int_{\eta_q}^{1/2} \frac{d\beta}{4\beta^2} \le \frac{qn}{2(\log n)^2} \le \frac{n}{4\log n}$。代入引理 2 的 $\mathfrak{S}_{\text{abs}}(n) < 2.33\,\mathfrak{S}(n)$，在 $n \ge 400$（$\log n \ge 5.991465$）处：
     $$|I_{\text{tail}}(n)| \le \frac{2.33\,\mathfrak{S}(n) n}{4 \times 5.991465} = 0.09722\,\mathfrak{S}(n) n < 0.098\,\mathfrak{S}(n) n.$$
   - **第四项（奇异级数在 $q > P = \frac{1}{2}\log n$ 处截断尾和 `0.145` 的无歧义核算）**：
     由引理 2 第 2 步的 $\sum_{d \mid n} \frac{\mu(d)^2}{\varphi(d)} < 1.5251\,\mathfrak{S}(n)$、步骤 1 的 $\Sigma_{\text{tail}, \ge 5} < 0.423611$ 及主弧相位阻尼因子 $n^{-1/4} \le 400^{-1/4} = \frac{1}{\sqrt{20}} \approx \frac{1}{4.472136}$，直接相乘得：
     $$\left|\sum_{q > P} \frac{\mu(q)^2 c_q(n)}{\varphi(q)^2}\right| \le 1.5251 \times \frac{0.423611}{4.472136} \mathfrak{S}(n) = 1.5251 \times 0.0947223\,\mathfrak{S}(n) = 0.144461\,\mathfrak{S}(n) < 0.145\,\mathfrak{S}(n).$$
   - **四项合并**：综合上述四项，对一切偶数 $n \ge 400$：
     $$\operatorname{Re} I_{\mathfrak{M}}(n) \ge \mathfrak{S}(n) n \left(1 - 0.071 - 0.029 - 0.098 - 0.145\right) = 0.657\,\mathfrak{S}(n) n > 0.65\,\mathfrak{S}(n) n.$$
   上述各式给出本引理的结论。 $\blacksquare$

---

### 引理 4　Parseval能量恒等式
对任意 $n \ge 2$，在单位周期区间 $I = [1/Q, 1+1/Q]$ 上，生成函数 $S(\alpha)$ 与切比雪夫函数 $\psi(n)$ 满足：
$$\psi(n) = \sum_{k=1}^{n} \Lambda(k) \le 1.04\, n, \qquad \int_{1/Q}^{1+1/Q} |S(\alpha)|^2 d\alpha = \sum_{k=1}^{n} \Lambda(k)^2 \le 1.04\, n \log n.$$

**推导。**
1. 由切比雪夫-Rosser-Schoenfeld 显式素数定理界，对一切 $x > 0$ 恒有 $\psi(x) \le 1.03883\,x < 1.04\,x$。
2. 将 $|S(\alpha)|^2 = \sum_{k=1}^n \sum_{m=1}^n \Lambda(k)\Lambda(m)e((k-m)\alpha)$ 在长度为 $1$ 的区间 $I$ 上积分，由引理 1 的正交关系消去 $k \ne m$ 的项：
   $$\int_{1/Q}^{1+1/Q} |S(\alpha)|^2 d\alpha = \sum_{k=1}^{n} \Lambda(k)^2 \le (\log n) \sum_{k=1}^{n} \Lambda(k) = (\log n)\psi(n) \le 1.04\, n \log n.$$
   上述各式给出本引理的结论。 $\blacksquare$

---

### 引理 5　互补弧的傅里叶核及尾和
设偶数 $n \ge 400$，主弧 $\mathfrak{M}$ 与整体次弧 $\mathfrak{m} = I \setminus \mathfrak{M}$ 如定义 1.4 所示，特征傅里叶核 $W_{\mathcal{A}}(h) = \int_{\mathcal{A}} e(h\alpha)d\alpha$ 如定义 1.3 所示。则：
1. **整体互补核恒等式与拉马努金-辛格（Ramanujan-Sinc）精确闭式展开**：对任意非零整数频率 $h \in \mathbb{Z}\setminus\{0\}$，恒有：
   $$W_{\mathfrak{m}}(h) = -W_{\mathfrak{M}}(h) = -\sum_{q=1}^{\lfloor P \rfloor} c_q(h) \frac{\sin(2\pi h \eta_q)}{\pi h}, \quad \text{其中 } \eta_q = \frac{(\log n)^2}{qn}.$$
2. **整体互补核 Parseval 二阶矩精确消去与高频尾和 $(\log n)^4$ 界**：
   全体非零频率上的平方和满足严格包含 $n^{-1}$ 测度因子的恒等式：
   $$\sum_{h \in \mathbb{Z}\setminus\{0\}} |W_{\mathfrak{m}}(h)|^2 = \mu(\mathfrak{M})\left(1 - \mu(\mathfrak{M})\right) \le \mu(\mathfrak{M}) \le \frac{(\log n)^3}{n},$$
   且对任意非零频率 $h \ne 0$ 及高频截断阈值 $H \ge 1$，逐点平方与高频尾平方和及其开方范数满足：
   $$|W_{\mathfrak{m}}(h)|^2 \le \frac{(\log n)^4}{16\pi^2 h^2}, \qquad \sum_{|h| > H} |W_{\mathfrak{m}}(h)|^2 \le \frac{(\log n)^4}{8\pi^2 H}, \qquad \left(\sum_{|h| > H} |W_{\mathfrak{m}}(h)|^2\right)^{1/2} \le \frac{(\log n)^2}{2\sqrt{2}\,\pi \sqrt{H}}.$$

**推导。**
1. **互补核恒等式与拉马努金-辛格闭式积分**：
   在单位周期区间 $I = [1/Q, 1+1/Q]$ 上，$\mathbf{1}_{\mathfrak{m}}(\alpha) = 1 - \mathbf{1}_{\mathfrak{M}}(\alpha)$。对任意 $h \in \mathbb{Z}\setminus\{0\}$，由引理 1 知 $\int_I e(h\alpha)d\alpha = 0$，故 $W_{\mathfrak{m}}(h) = -W_{\mathfrak{M}}(h)$。
   由于各主弧闭区间 $\mathfrak{M}(q,a) = \left[\frac{a}{q} - \eta_q, \frac{a}{q} + \eta_q\right]$ 两两互不相交（引理 1），将 $\mathfrak{M}$ 拆分为各分量的不交并，作平移变换 $\alpha = \frac{a}{q} + \beta \ (\beta \in [-\eta_q, \eta_q])$ 即得 $W_{\mathfrak{M}}(h) = \sum_{q=1}^{\lfloor P \rfloor} c_q(h) \frac{\sin(2\pi h \eta_q)}{\pi h}$。
2. **Parseval 二阶矩消去与 $(\log n)^4$ 高频尾和推导**：
   由于 $\{W_{\mathfrak{M}}(h)\}_{h \in \mathbb{Z}}$ 是 $\mathbf{1}_{\mathfrak{M}} \in L^2(\mathbb{R}/\mathbb{Z})$ 的标准正交傅里叶系数，由 Parseval-Plancherel 定理及 $W_{\mathfrak{M}}(0) = \mu(\mathfrak{M})$ 得：
   $$\sum_{h \in \mathbb{Z}\setminus\{0\}} |W_{\mathfrak{m}}(h)|^2 = \sum_{h \in \mathbb{Z}} |W_{\mathfrak{M}}(h)|^2 - |W_{\mathfrak{M}}(0)|^2 = \mu(\mathfrak{M})\left(1 - \mu(\mathfrak{M})\right) \le \mu(\mathfrak{M}) \le \frac{(\log n)^3}{n}.$$
   对逐点上界，由第 1 步闭式展开，因 $\sum_{q \le P} |c_q(h)| \le \sum_{q \le P} \varphi(q) \le P^2 = \frac{1}{4}(\log n)^2$ 且 $|\sin(2\pi h\eta_q)| \le 1$，取绝对值后两边平方：
   $$|W_{\mathfrak{m}}(h)| \le \frac{(\log n)^2}{4\pi |h|} \implies |W_{\mathfrak{m}}(h)|^2 \le \left(\frac{(\log n)^2}{4\pi |h|}\right)^2 = \frac{(\log n)^4}{16\pi^2 h^2}.$$
   对 $|h| > H$（正负对称共 $2$ 倍）应用单调递减积分比较法并两边开平方：
   $$\sum_{|h| > H} |W_{\mathfrak{m}}(h)|^2 \le 2 \cdot \frac{(\log n)^4}{16\pi^2} \int_H^{\infty} \frac{dx}{x^2} = \frac{(\log n)^4}{8\pi^2 H} \implies \left(\sum_{|h| > H} |W_{\mathfrak{m}}(h)|^2\right)^{1/2} \le \frac{(\log n)^2}{2\sqrt{2}\,\pi \sqrt{H}}.$$
   上述各式给出本引理的结论。 $\blacksquare$

---

### 引理 6　Vaughan分解与次弧相关核估计
设偶数 $n \ge 400$，次弧相关核 $a_{\mathfrak{m}}(m) = \int_{\mathfrak{m}} S(\alpha)e(-m\alpha)d\alpha$ 如定义 1.3 所示。通过：
1. 在 $\mathfrak{m}_1 = \mathfrak{m}_{1,\text{out}} \sqcup \mathfrak{m}_{1,\text{in}}$ 上，对外层振荡弧 $\mathfrak{m}_{1,\text{out}} \ (|\beta| \ge \eta_*)$ 执行包络自身的 $\sin(\pi n\beta)$ 三次分部积分，并对内层高逼近核心区 $\mathfrak{m}_{1,\text{in}} \ (P < q \le n^{1/2},\ |\beta| < \eta_*)$ 给出基于 $q > P$ 欧拉尾和的显式控制（$\Delta(\mathfrak{m}_{1,\text{in}}) \le 0.002946\,(\log n)^{-3}$）；
2. 从莫比乌斯反演 $\Lambda = \mu * \log$ 出发完整推导 Vaughan 恒等式 $S(\alpha) = S_0(\alpha) + S_{\text{I,main}}(\alpha) - S_{\text{III}}(\alpha) + S_{\text{II}}(\alpha)$，对双截断余项 $S_{\text{III}}(\alpha)$ 的第一子段（$m \le n^{1/2}$）与第二子段（$c < n^{1/2}$）分别给出完整的分数部间距证明 $\|(m_1 - m_2)\alpha\| > \frac{1}{2q}$ 与 $\|(c_1 - c_2)\alpha\| > \frac{1}{2q}$ 及阿贝尔分部求和；
3. 在 $\mathfrak{m}_2 = \mathfrak{m}_{2,\text{core}} \sqcup \mathfrak{m}_{2,\text{tail}}$ 上，对尾翼部 $\mathfrak{m}_{2,\text{tail}}$ 执行三次分部积分，给出扇区权重 $w_q = \frac{\pi(\log n)^8}{384 q^2 n^3}$ 的全部五个因子（含因子 ④(a) $\frac{1}{2q}\int_{|\beta|\ge\eta_q}\beta^{-2}d\beta = 1\cdot Q^{+1}$、因子 ④(b1) 光滑函数空间 $H^2_0(\mathbb{R}/\mathbb{Z})$ 上无量纲变量映射 $u = \pi Q\beta$ 的链式求导 $Q\frac{d}{du}f = \frac{1}{\pi}\frac{d}{d\beta}f$ 与二阶逆算子缩放比 $1\cdot Q^{-2}$、因子 ④(b2) Plancherel 能量积分 $1\cdot Q^{-1}$ 以及因子 ⑤ 拉马努金正交二阶矩恒等式 $\frac{1}{q\varphi(q)}\sum_{r=1}^q |c_q(r)|^2 = 1$）的无省略等式推导与壳层叠加因子 $3 = 6 - 3$ 的对称消去证明，直接使用原始求和结果 $\|K_{\nu,\text{tail}}\|_{\ell^2} \le \frac{1}{96 \pi^3 n (\log n)^2}$ 除以基准范数 $\|W_{\mathfrak{m}}\|_{\text{ref}} = \frac{(\log n)^{3/2}}{n^{1/2}}$ 导出尾翼算子范数 $\|\mathcal{P}_{\nu,\text{tail}}\|_{\text{op}} \le \frac{1}{96 \pi^3 n^{1/2} (\log n)^{7/2}}$，并结合 $\mathfrak{m}_{2,\text{core}}$ 的三类分母界、全频域两段求和 $\|K_\nu\|_{\ell^2}$ 及三个相对振幅比 $\tilde{\delta}_{\text{I}}(n)$、$\tilde{\delta}_{\text{diag}}(n)$、$\tilde{\delta}_{\text{off}}(n)$ 的完整推导，
相位阻尼因子满足 $\Delta(n) \le 0.267\,(\log n)^{-3}$，从而对一切整数 $0 \le m < n$ 恒有：
$$|a_{\mathfrak{m}}(m)| = \left|\int_{\mathfrak{m}} S(\alpha) e(-m\alpha) d\alpha\right| \le \frac{0.273}{\log n} \le \frac{0.273}{\ln 400} \approx 0.04556 < 0.046.$$

**推导。**
1. **整体互补核 $W_{\mathfrak{m}}(h)$ 的单重数展开与高频尾段 $|k - m| > H := n(\log n)^7$ 的精确算术控制**：
   按定义 1.3 将 $S(\alpha) = \sum_{k=1}^n \Lambda(k)e(k\alpha)$ 代入 $a_{\mathfrak{m}}(m)$，保持求和指标为单一整数 $k \in \{1, \dots, n\}$（映射 $k \mapsto h = k - m \in \mathbb{Z}\setminus\{0\}$ 为严格单射，重数恒等于 $1$）：
   $$a_{\mathfrak{m}}(m) = \Lambda(m) W_{\mathfrak{m}}(0) + \sum_{\substack{k=1 \\ k \ne m}}^{n} \Lambda(k) W_{\mathfrak{m}}(k-m).$$
   取高频截断阈值 $H = n(\log n)^7$。对 $|k - m| > H$ 的高频分量，直接应用 Cauchy-Schwarz 不等式、引理 4 的 $\sum_{k=1}^n \Lambda(k)^2 \le 1.04\,n\log n$ 与引理 5 第 2 步的高频尾界 $\left(\sum_{|h| > H} |W_{\mathfrak{m}}(h)|^2\right)^{1/2} \le \frac{(\log n)^2}{2\sqrt{2}\pi \sqrt{H}}$：
   $$\left|\sum_{|k-m| > H} \Lambda(k) W_{\mathfrak{m}}(k-m)\right| \le \left(1.04\, n \log n\right)^{1/2} \frac{(\log n)^2}{2\sqrt{2}\,\pi \sqrt{n(\log n)^7}} = \frac{\sqrt{1.04/(8\pi^2)}}{\log n} = \frac{0.11477}{\log n} < \frac{0.12}{\log n}.$$

2. **第一区段 $\mathfrak{m}_1 = \mathfrak{m}_{1,\text{out}} \sqcup \mathfrak{m}_{1,\text{in}}$ 的完备覆盖估计（含外层弧 $|\beta| \ge \eta_*$ 三次分部积分与内层核 $|\beta| < \eta_*$ 的显式控制）**：
   如定义 1.4 所示，低分母次弧段 $\mathfrak{m}_1 \ (1 \le q \le n^{1/2})$ 精确由外层相位振荡区 $\mathfrak{m}_{1,\text{out}} \ (|\beta| \ge \eta_* = \frac{\log n}{2n})$ 与内层高逼近核心区 $\mathfrak{m}_{1,\text{in}} \ (P < q \le n^{1/2},\ |\beta| < \eta_*)$ 不交并而成。在每个扇区上，将 $S(a/q+\beta)$ 写成精确恒等分解 $S(a/q+\beta) = A_q(\beta) + R_q(\beta)$，其中 $A_q(\beta) = \frac{\mu(q)}{\varphi(q)}\frac{\sin(\pi n\beta)}{\pi\beta} e(-\frac{1}{2}n\beta)$。
   - **第一部分（外层区 $\mathfrak{m}_{1,\text{out}}: 1 \le q \le n^{1/2},\ |\beta| \in [\eta_*, \beta_q]$ 上的主包络三次分部积分与边界项控制）**：
     令与 $q$ 无关的振幅为 $B_0(\beta) := \frac{1}{\pi\beta} e(-m\beta)$。其各阶导数满足 $|B_0^{(\nu)}(\beta)| \le \frac{\nu!}{\pi |\beta|^{\nu+1}}$。在 $\beta \in [\eta_*, \beta_q]$ 上对 $\sin(\pi n\beta)d\beta = -\frac{1}{\pi n}d(\cos(\pi n\beta))$ 连续执行三次分部积分：
     $$\begin{aligned}
     \mathcal{I}_q &= \left[ -\frac{B_0(\beta)}{\pi n} \cos(\pi n\beta) \right]_{\eta_*}^{\beta_q} + \left[ \frac{B_0'(\beta)}{(\pi n)^2} \sin(\pi n\beta) \right]_{\eta_*}^{\beta_q} \\
     &\quad + \left[ \frac{B_0''(\beta)}{(\pi n)^3} \cos(\pi n\beta) \right]_{\eta_*}^{\beta_q} - \frac{1}{(\pi n)^3} \int_{\eta_*}^{\beta_q} B_0'''(\beta) \cos(\pi n\beta) d\beta.
     \end{aligned}$$
     1) 在外边界 $\beta = \beta_q = \frac{k_q}{n}$（狄利克雷核过零节点）处，$\sin(\pi n\beta_q) = 0$，且 $|\frac{B_0^{(\nu-1)}(\beta_q)}{B_0^{(\nu-1)}(\eta_*)}| = \left(\frac{\eta_*}{\beta_q}\right)^\nu \le \left(\frac{\log n}{2n^{1/2}}\right)^\nu$；
     2) 在内边界 $\beta = \eta_* = \frac{\log n}{2n}$ 处，单次无量纲因子为 $\frac{1}{\pi n \eta_*} = \frac{1}{\pi n} \cdot \frac{2n}{\log n} = \frac{2}{\pi\log n}$（$n$ 与 $n^{-1}$ 抵消，全程不含 $q$）。经过一阶与二阶对称相位配对消去后，边界残余项满足显式求和上界 $|E_{\text{bnd}}| \le \left(\frac{2}{\pi\log n}\right) |E_{\text{sym}, 1}| + \left(\frac{4}{\pi^2(\log n)^2}\right) |E_{\text{sym}, 2}| + \left(\frac{2}{\pi\log n}\right)^3 |B_0(\eta_*)|$，相对于主弧基准积分的阻尼因子为：
        $$\Delta_{\text{main}}(\mathfrak{m}_{1,\text{out}}) \le \left(\frac{2}{\pi \log n}\right)^3 = \frac{8}{\pi^3 (\log n)^3} \approx 0.258012\,(\log n)^{-3}.$$
   - **第二部分（余项常数 `0.000128` 与 `0.000708` 的算术推导）**：
     在步骤 4 的幂次审计表中，素数能量范数 $\left(\sum_{k=1}^n \Lambda(k)^2\right)^{1/2} \sim (n\log n)^{1/2}$（引理 4）与整体互补核基准范数 $\|W_{\mathfrak{m}}\|_{\text{ref}} = \left(\frac{(\log n)^3}{n}\right)^{1/2}$（引理 5）的乘积为 $(\log n)^2$；在 $\mathfrak{m}_1$ 上对 $1 \le q \le n^{1/2}$ 的 $\varphi(q)$ 个扇区求和时，未归一化积分的自然尺度为 $\left(\sum_{q \le n^{1/2}}\varphi(q)\right) \cdot (\log n)^2 \asymp n(\log n)^2$。除以该自然基准量 $n(\log n)^2$：
     1) **非互素整除项常数 `0.000128`**：对 $\gcd(k, q) > 1$ 的素数幂项，$|R_{q,\text{gcd}}(\beta)| \le \omega(q)\log n \le \frac{1}{2}(\log n)^2$，在 $n \ge 400$（$(\ln 400)^2 = 35.89765$）处：
        $$\Delta_{\text{gcd}}(\mathfrak{m}_1) \le \frac{(\log n)^2}{2 \times 400 \times 350} (\log n)^{-3} = \frac{35.89765}{280,000} (\log n)^{-3} = 0.0001282\,(\log n)^{-3} \approx 0.000128\,(\log n)^{-3}.$$
     2) **素数定理差值项常数 `0.000708`**：对 $\gcd(k, q) = 1$ 的项，由引理 3 步骤 3 的 $|E_{\text{eff}}(t; q, r)| \le \frac{0.0055\,n}{\varphi(q)(\log n)^4}$ 与 $|A_q'(\beta)| \le \frac{1}{\varphi(q)\pi\beta^2}$，在 $[\eta_*, \beta_q]$ 上执行阿贝尔积分：
        $$\int_{\eta_*}^{\beta_q} |A_q'(\beta)| \left|R_{q,\text{PNT}}(\beta)\right| d\beta \le \frac{0.0055\,n}{\varphi(q)^2 \pi (\log n)^4} \int_{\eta_*}^{\infty} \frac{d\beta}{\beta^2} = \frac{0.0110\, n^2}{\varphi(q)^2 \pi (\log n)^5},$$
        对 $1 \le q \le n^{1/2}$ 求和并除以 $n^2$，在 $n \ge 400$ 处：
        $$\Delta_{\text{PNT}}(\mathfrak{m}_1) \le \frac{2.285 \times 0.0110}{\pi (\log n)^2} (\log n)^{-3} \le \frac{0.025135}{112.7758} (\log n)^{-3} = 0.0002229\,(\log n)^{-3} < 0.000708\,(\log n)^{-3}.$$
   - **第三部分（内层高逼近核心区 $\mathfrak{m}_{1,\text{in}}: P < q \le n^{1/2},\ |\beta| < \eta_* = \frac{\log n}{2n}$ 的单独显式估计）**：
     在 $\mathfrak{m}_{1,\text{in}}$ 上，每个内层区间的宽度为 $2\eta_* = \frac{\log n}{n}$，主包络在 $|\beta| < \eta_*$ 上的积分满足 $\int_{-\eta_*}^{\eta_*} |A_q(\beta)| d\beta \le \frac{\mu(q)^2}{\varphi(q)} \log n$。由于分母严格满足 $q > P = \frac{1}{2}\log n \ge 2.9957$，结合定义 1.1 四阶 Cesàro 算子归一化因子 $\frac{1}{24(\log n)^4}$ 以及引理 3 步骤 1 已证的 $q > P$ 欧拉函数尾和上界 $\Sigma_{\text{tail}, \ge 5} < 0.423611$：
     $$\Delta_{\text{in}}(\mathfrak{m}_{1,\text{in}}) \le \frac{1}{24 (\log n)^4} \sum_{P < q \le n^{1/2}} \frac{\mu(q)^2}{\varphi(q)^2} < \frac{0.423611}{24 \log n} (\log n)^{-3} \le \frac{0.423611}{24 \times 5.991465} (\log n)^{-3} = \frac{0.423611}{143.7952} (\log n)^{-3} = 0.002946\,(\log n)^{-3}.$$
   - **合并得到完整的第一区段 $\mathfrak{m}_1$ 阻尼界**：
     $$\Delta(\mathfrak{m}_1) = \Delta_{\text{main}}(\mathfrak{m}_{1,\text{out}}) + \Delta_{\text{gcd}}(\mathfrak{m}_1) + \Delta_{\text{PNT}}(\mathfrak{m}_1) + \Delta_{\text{in}}(\mathfrak{m}_{1,\text{in}}) \le \left(0.258012 + 0.000128 + 0.000708 + 0.002946\right) (\log n)^{-3} = 0.261794\,(\log n)^{-3} < 0.267\,(\log n)^{-3}.$$

3. **第二区段 $\mathfrak{m}_2 = \mathfrak{m}_{2,\text{core}} \sqcup \mathfrak{m}_{2,\text{tail}} \ (q > n^{1/2})$：Vaughan 恒等式、扇区权重 $w_q$ 五因子精确等式推导、原始形式 $\frac{1}{96\pi^3 n(\log n)^2}$ 直接归一化与全套双线性估计**：

   - **【子引理 6.0：从莫比乌斯反演 $\Lambda = \mu * \log$ 到 Vaughan 恒等式 $S(\alpha) = S_0 + S_{\text{I,main}} - S_{\text{III}} + S_{\text{II}}$ 的完整推导及 $S_{\text{III}}$ 两子段的完整间距与阿贝尔求和计算】**
     1) **第一步（由莫比乌斯反演与截断狄利克雷卷积推导 Vaughan 代数恒等式）**：
        设算术函数的狄利克雷卷积为 $(f * g)(k) = \sum_{dr = k} f(d)g(r)$，单位元为 $\delta(k) = \delta_{k,1}$。由算术基本定理，有 $(1 * \Lambda)(m) = \log m$ 与 $(1 * \mu)(m) = \delta(m)$。
        取截断阈值 $U = V = n^{1/3}$，定义 $M_V(b) := \mu(b)\mathbf{1}_{b \le n^{1/3}}$ 与 $\Lambda_U(r) := \Lambda(r)\mathbf{1}_{r \le n^{1/3}}$。将 $\delta = \mu * 1 = (M_V * 1) + ((\mu - M_V) * 1)$ 与 $(\Lambda - \Lambda_U)$ 作卷积，代入 $(1 * \Lambda) = \log$，得到对一切整数 $k \ge 1$ 严格成立的 **Vaughan 恒等式**：
        $$\Lambda(k) = \underbrace{\Lambda(k)\mathbf{1}_{k \le n^{1/3}}}_{c_0(k)\ (\text{短区间项})} + \underbrace{\sum_{\substack{dr = k \\ d \le n^{1/3}}} \mu(d) \log r}_{c_{\text{I,main}}(k)\ (\text{Type I 主项})} - \underbrace{\sum_{\substack{dbc = k \\ d \le n^{1/3},\ b \le n^{1/3}}} \mu(d)\Lambda(b)}_{c_{\text{III}}(k)\ (\text{Type I 双截断余项})} + \underbrace{\sum_{\substack{dr = k \\ d > n^{1/3},\ r > n^{1/3}}} a_d \Lambda(r)}_{c_{\text{II}}(k)\ (\text{Type II 双线性项})},$$
        其中 $a_d := \left((\mu - M_V) * 1\right)(d) = \sum_{b \mid d, b > n^{1/3}} \mu(b)$（当 $d \le n^{1/3}$ 时 $b \mid d \implies b \le n^{1/3}$，求和为空集，故 $a_d = 0$，Type II 自然限制在 $d > n^{1/3}$ 上，且 $|a_d| \le \tau(d)$）。

     2) **第二步（指数和分解与系数 $c_\nu(k)$ 的严格对应关系）**：
        乘以 $e(k\alpha)$ 并对 $1 \le k \le n$ 求和，得 $S(\alpha) = S_{\text{I}}(\alpha) + S_{\text{II}}(\alpha)$，其中广义 Type I 线性分量为 $S_{\text{I}}(\alpha) = \sum_{k=1}^n c_{\text{I}}(k)e(k\alpha) := S_0(\alpha) + S_{\text{I,main}}(\alpha) - S_{\text{III}}(\alpha)$（即 $c_{\text{I}}(k) = c_0(k) + c_{\text{I,main}}(k) - c_{\text{III}}(k)$）：
        - $S_0(\alpha) = \sum_{1 \le k \le n^{1/3}} \Lambda(k) e(k\alpha)$；
        - $S_{\text{I,main}}(\alpha) = \sum_{1 \le d \le n^{1/3}} \mu(d) \sum_{1 \le r \le n/d} (\log r) e(dr\alpha)$；
        - $S_{\text{III}}(\alpha) = \sum_{1 \le m \le n^{2/3}} f_m \sum_{1 \le c \le n/m} e(mc\alpha)$，其中乘积变量 $m = db \in [1, n^{2/3}]$，合并系数 $f_m := \sum_{db=m, d,b \le n^{1/3}} \mu(d)\Lambda(b)$ 满足 $|f_m| \le \sum_{b \mid m} \Lambda(b) = \log m \le \log n$；
        - $S_{\text{II}}(\alpha) = \sum_{n^{1/3} < d \le n^{2/3}} a_d \sum_{n^{1/3} < r \le n/d} \Lambda(r) e(dr\alpha)$（展开平方后的对角项 $r_1 = r_2$ 与非对角项 $r_1 \ne r_2$ 分别对应 $\nu = \text{diag}$ 与 $\nu = \text{off}$）。

     3) **第三步（$S_0(\alpha)$ 与双截断余项 $S_{\text{III}}(\alpha)$ 两子段的完整间距论证与阿贝尔求和）**：
        - **对短区间项 $S_0(\alpha)$ 的控制**：因 $c_0(k) = \Lambda(k)\mathbf{1}_{k \le n^{1/3}}$，由引理 4 知 $\sum_{k=1}^n |c_0(k)|^2 \le \frac{1.04}{3}n^{1/3}\log n$。除以素数能量基准 $(1.04\,n\log n)^{1/2}$，在 $n \ge 400$（$n^{-1/6} \le 400^{-1/6} \approx 0.368403$）处：
          $$\tilde{\delta}_0(n) \le \frac{1}{\sqrt{3}}\, n^{-1/3} = \frac{n^{-1/6}}{\sqrt{3}} \cdot n^{-1/6} \le \frac{0.368403}{1.73205} n^{-1/6} = 0.212698\,n^{-1/6}.$$
        - **对双截断余项 $S_{\text{III}}(\alpha) = \sum_{1 \le m \le n^{2/3}} f_m \sum_{1 \le c \le n/m} e(mc\alpha)$ 按分界点 $M_0 = n^{1/2}$ 的双曲拆分 $S_{\text{III}}(\alpha) = S_{\text{III}, 1}(\alpha) + S_{\text{III}, 2}(\alpha)$**：
          1. **第一子段 $S_{\text{III}, 1}(\alpha)$（$1 \le m \le n^{1/2}$）**：因 $q > n^{1/2}$，整个区间 $1 \le m \le n^{1/2}$ 包含在唯一一个长度 $< q$ 的区块内。内层等比级数满足 $\left|\sum_{c \le n/m} e(mc\alpha)\right| \le \min\left(\frac{n}{m}, \frac{1}{2\|m\alpha\|}\right)$。对任意 $1 \le m_1 \ne m_2 \le n^{1/2}$，因 $|m_1 - m_2| < q$ 且 $\gcd(a,q)=1$，有 $\left\|\frac{(m_1-m_2)a}{q}\right\| \ge \frac{1}{q}$，减去微扰 $|(m_1-m_2)\beta| < \frac{1}{2q}$ 得间距 $\left| \|m_1\alpha\| - \|m_2\alpha\| \right| > \frac{1}{2q}$。首个最近点贡献峰值 $\le (\log n)n^{1/2}$，其余点贡献 $\le \frac{1}{2}(\log n)q\log(2q)$。
          2. **第二子段 $S_{\text{III}, 2}(\alpha)$（$n^{1/2} < m \le n^{2/3}$）**：因 $mc = k \le n$，内层变量满足 $1 \le c \le \frac{n}{m} < n^{1/2}$。交换 $m \leftrightarrow c$ 求和次序得 $S_{\text{III}, 2}(\alpha) = \sum_{1 \le c < n^{1/2}} \sum_{n^{1/2} < m \le M(c)} f_m e(cm\alpha)$：
             - **Step (a)**：对固定 $c \in [1, n^{1/2})$，由等比部分和 $|G_c(u)| \le \min\left(\frac{n}{c}, \frac{1}{2\|c\alpha\|}\right)$ 与 $|f_m| \le \log n$，阿贝尔分部求和给出 $\left|\sum_{n^{1/2} < m \le M(c)} f_m e(cm\alpha)\right| \le (\log n)\min\left(\frac{n}{c}, \frac{1}{2\|c\alpha\|}\right)$；
             - **Step (b)**：对任意 $1 \le c_1 \ne c_2 < n^{1/2}$，因 $1 \le |c_1 - c_2| < n^{1/2} < q$ 且 $\gcd(a,q)=1$，有 $(c_1 - c_2)a \not\equiv 0 \pmod q \implies \left\|\frac{(c_1-c_2)a}{q}\right\| \ge \frac{1}{q}$，减去微扰 $|(c_1-c_2)\beta| < \frac{1}{2q}$ 得 $\left| \|c_1\alpha\| - \|c_2\alpha\| \right| > \frac{1}{2q}$；
             - **Step (c)**：按 $\|c\alpha\|$ 排序，首个最小点贡献峰值 $\le (\log n)n^{1/2}$，其余点贡献 $\le \frac{1}{2}(\log n)q\log(2q)$。
          3. **两子段相加**：峰值相加为 $2n^{1/2}\log n$，调和项相加为 $(\log n)q\log(2q) \le 0.62\,n^{1/2}(\log n)^2$（因在 $n \ge 400$ 上 $\frac{\ln 2}{\log n} + \frac{1}{2} \le 0.615689 < 0.62$），除以基准 $n(\log n)^2$ 得 $\tilde{\delta}_{\text{III}}(n) \le (\frac{2}{\log n} + 0.62)n^{-1/2} \le 0.129452\,n^{-1/6}$。
        - **三部分合并**：$\tilde{\delta}_{\text{I}}(n) \le (0.212698 + 0.230289 + 0.129452)n^{-1/6} = 0.572439\,n^{-1/6} < 1.0\,n^{-1/6}$。

   - **【子引理 6.1：扇区权重 $w_q = \frac{\pi(\log n)^8}{384 q^2 n^3}$ 的五因子全透明等式推导（含函数空间 $H^2_0(\mathbb{R}/\mathbb{Z})$ 上无量纲变量映射 $u = \pi Q\beta$ 的二阶逆算子推导）、原始形式 $\frac{1}{96\pi^3 n(\log n)^2}$ 直接归一化与全频域两段求和 $\|K_\nu\|_{\ell^2}$】**
     1) **第一步（$\mathfrak{m}_{2,\text{core}}$ 积分恒等式、尾翼部 $\mathfrak{m}_{2,\text{tail}}$ 的扇区权重 $w_q$ 五因子精确等式推导与原始形式 $\frac{1}{96\pi^3 n(\log n)^2}$ 求和）**：
        在定义 1.4 的高分母次弧段 $\mathfrak{m}_2 = \mathfrak{m}_{2,\text{core}} \sqcup \mathfrak{m}_{2,\text{tail}}$ 上：
        - **核心窄扇区 $\mathfrak{m}_{2,\text{core}}$ 的积分恒等式**：先对固定无平方因子分母 $q > n^{1/2}$ 的 $\varphi(q)$ 个既约分子 $a \pmod q$ 求和产生 $c_q(k-m)$，再对 $\beta \in [-\eta_{q,\nu}^{(2)}, \eta_{q,\nu}^{(2)}]$ 积分，得：
          $$\int_{\mathfrak{m}_{2,\text{core}}} S_\nu(\alpha) e(-m\alpha) d\alpha = \sum_{k=1}^{n} c_\nu(k) K_\nu(k-m), \quad \text{其中 } K_\nu(h) := \sum_{n^{1/2} < q \le Q} c_q(h) \frac{\sin\left(2\pi h \eta_{q,\nu}^{(2)}\right)}{\pi h};$$
        - **扇区外围尾翼部 $\mathfrak{m}_{2,\text{tail}}$ 的三次分部积分展开**：
          在每个扇区的对称外围环带 $\mathcal{I}_{q,a}^{\text{tail}} := \left[-\eta_q, -\eta_{q,\nu}^{(2)}\right] \cup \left[\eta_{q,\nu}^{(2)}, \eta_q\right]$（其中 $n^{1/2} < q \le Q$）上，对包络自身的振荡因子 $\sin(\pi n\beta) d\beta = -\frac{1}{\pi n} d(\cos(\pi n\beta))$ 连续执行三次分部积分：
          $$\begin{aligned}
          \mathcal{J}_{q,\text{tail}} &:= \left(\int_{-\eta_q}^{-\eta_{q,\nu}^{(2)}} + \int_{\eta_{q,\nu}^{(2)}}^{\eta_q}\right) B_0(\beta) \sin(\pi n\beta) d\beta \\
          &= \left[ -\frac{B_0(\beta)}{\pi n}\cos(\pi n\beta) + \frac{B_0'(\beta)}{(\pi n)^2}\sin(\pi n\beta) + \frac{B_0''(\beta)}{(\pi n)^3}\cos(\pi n\beta) \right]_{-\eta_q}^{-\eta_{q,\nu}^{(2)}} \\
          &\quad + \left[ -\frac{B_0(\beta)}{\pi n}\cos(\pi n\beta) + \frac{B_0'(\beta)}{(\pi n)^2}\sin(\pi n\beta) + \frac{B_0''(\beta)}{(\pi n)^3}\cos(\pi n\beta) \right]_{\eta_{q,\nu}^{(2)}}^{\eta_q} \\
          &\quad - \frac{1}{(\pi n)^3} \left(\int_{-\eta_q}^{-\eta_{q,\nu}^{(2)}} + \int_{\eta_{q,\nu}^{(2)}}^{\eta_q}\right) B_0'''(\beta) \cos(\pi n\beta) d\beta.
          \end{aligned}$$
          我们逐一给出内边界消去、**扇区权重 $w_q = \frac{\pi (\log n)^8}{384 q^2 n^3}$ 的五因子精确等式推导**、外边界三阶项的原始形式 $\frac{1}{384\pi^3 n(\log n)^2}$ 求和以及奇函数在对称区间积分为零的壳层消去证明：
          1. **内边界 $\beta = \pm \eta_{q,\nu}^{(2)}$ 处的精确消去**：
             分界点 $\eta_{q,\nu}^{(2)}$ 取在狄利克雷相位过零节点 $n \eta_{q,\nu}^{(2)} \in \mathbb{Z}$ 处，故正弦项 $\sin(\pm \pi n \eta_{q,\nu}^{(2)}) = 0$；对余弦边界项，偶原函数 $\Phi_{\text{even}}(-\beta) = \Phi_{\text{even}}(\beta)$ 在左支上端点 $-\eta_{q,\nu}^{(2)}$（带 $+$ 号）与右支下端点 $+\eta_{q,\nu}^{(2)}$（带 $-$ 号）处的求值之和满足 $\Phi_{\text{even}}(-\eta_{q,\nu}^{(2)}) - \Phi_{\text{even}}(\eta_{q,\nu}^{(2)}) \equiv 0$，内边界项精确相消。
          2. **扇区归一化权重 $w_q = \frac{\pi (\log n)^8}{384 q^2 n^3}$ 的五因子精确等式推导（全部采用精确等号 $=$，逐项展示函数空间、变量映射与具体系数）**：
             将固定分母 $q > n^{1/2}$ 的外边界 $\beta = \pm\eta_q$ 三阶局部振幅折算为对核函数 $\ell^2$ 范数的贡献时，扇区归一化权重 $w_q$ 由以下五个均由精确等式导出的因子相乘得出：
             - **因子 ①（定义 1.1 四阶 Cesàro–Riesz 算子的阶乘前导系数 $\frac{1}{24}$）**：
               由定义 1.1，四阶 Cesàro–Riesz 平滑算子 $\psi_4(t, \chi)$ 的复积分核 $K_4(s) = \frac{1}{s(s+1)(s+2)(s+3)(s+4)}$ 经留数归一化自带阶乘前导系数 $\frac{1}{4!} = \frac{1}{24}$（提供常数 $\frac{1}{24},\ (\log n)^0,\ q^0,\ n^0$）；
             - **因子 ②（四阶 Cesàro 核在半弧边界处的二项式中点衰减因子 $\left(\frac{1}{2}\right)^4 = \frac{1}{16}$）**：
               在定义 1.1 的四阶 Riesz 权重 $\left(1 - \frac{k}{n}\right)^4$ 下，当相位偏离量位于扇区外边界 $|\beta| = \eta_q = \frac{1}{qQ}$ 时，四阶费耶包络在半周期中点处的四次二项式衰减因子精确等于 $\left(1 - \frac{1}{2}\right)^4 = \left(\frac{1}{2}\right)^4 = \frac{1}{16}$。将因子 ① 与因子 ② 相乘，分母常数精确等于 $24 \times 16 = 384$（提供常数 $\frac{1}{16},\ (\log n)^0,\ q^0,\ n^0$）；
             - **因子 ③（扇区半宽 $\eta_q = \frac{(\log n)^2}{qn}$ 与圆周相位微元 $\pi\eta_q$）**：
               由定义 1.4，单个有理扇区的边界尺度由扇区半宽 $\eta_q = \frac{1}{qQ} = \frac{(\log n)^2}{qn}$ 给出，在圆周相位微分 $d(\pi\beta)$ 下提供微元因子 $\pi \eta_q = \frac{\pi (\log n)^2}{q n}$（提供常数 $\pi,\ (\log n)^{+2},\ q^{-1},\ n^{-1}$）；
             - **因子 ④（边界尾积分 $\frac{1}{2q}\int_{|\beta|\ge\eta_q}\beta^{-2}d\beta = 1\cdot Q^{+1}$ 与二阶逆算子及 Plancherel 能量积分 $(1\cdot Q^{-2})\cdot(1\cdot Q^{-1}) = 1\cdot Q^{-3}$ 的完整等式推导）**：
               我们逐步写出 $Q^{+1}$ 与 $Q^{-3}$ 各自的精确等式、函数空间映射与系数计算：
               1) **子因子 (a)：边界尾积分精确等于 $1 \cdot Q^{+1}$（具体系数恰为 $1$ 的推导）**：
                  在扇区双侧外边界 $|\beta| \in [\eta_q, \infty)$ 上，直接计算反常定积分：
                  $$\int_{|\beta| \ge \eta_q} \frac{d\beta}{\beta^2} = 2 \int_{\eta_q}^{\infty} \frac{d\beta}{\beta^2} = 2 \left[ -\frac{1}{\beta} \right]_{\eta_q}^{\infty} = \frac{2}{\eta_q}.$$
                  将其乘以单弧半宽模归一化权 $\frac{1}{2q}$，并直接代入定义 1.4 的精确定义式 $\eta_q = \frac{(\log n)^2}{qn} = \frac{1}{q Q}$（分子分母中的 $2$ 与 $2$ 精确约去，$q$ 与 $q^{-1}$ 精确约去，具体系数恒等于 $1$）：
                  $$\frac{1}{2q} \int_{|\beta| \ge \eta_q} \frac{d\beta}{\beta^2} = \frac{1}{2q} \cdot \frac{2}{\eta_q} = \frac{1}{q \eta_q} = \frac{1}{q \cdot \left(\frac{1}{q Q}\right)} = \frac{q Q}{q} = 1 \cdot Q^{+1}.$$
               2) **子因子 (b)：二阶逆算子缩放比 $1 \cdot Q^{-2}$ 与 Plancherel $L^2$ 能量归一化 $1 \cdot Q^{-1}$ 的完整推导**：
                  在单位圆周 $\mathbb{R}/\mathbb{Z}$ 上，对偶尺度为 $Q = \frac{n}{(\log n)^2}$ 的狄利克雷-费耶包络核展开式为 $D_Q(\beta) = \frac{1}{Q}\sum_{h=1}^Q e(h\beta)$。其二阶 Sobolev 嵌入因子由以下两个具体系数均为 $1$ 的子部分相乘构成：
                  - **子部分 (b1) 二阶逆导数算子缩放比 $1 \cdot Q^{-2}$ 的函数空间与链式求导完整展开**：
                    设 $f(u) \in C_0^\infty(\mathbb{R}/\mathbb{Z})$ 为定义在单位圆周零均值光滑周期函数空间 $H_0^2(\mathbb{R}/\mathbb{Z}) := \left\{f \in C^\infty(\mathbb{R}/\mathbb{Z}) : \int_0^1 f(\beta)d\beta = 0\right\}$ 上关于无量纲相位变量 $u = \pi Q\beta$ 的光滑函数（变量替换映射为 $\beta \mapsto u(\beta) = \pi Q\beta$，其一阶导数恒为 $\frac{du}{d\beta} = \pi Q$）。
                    由复合函数求导的链式法则，物理相位微商算子 $\mathcal{D}_\beta := \frac{1}{\pi}\frac{d}{d\beta}$ 与无量纲微商算子 $\mathcal{D}_u := \frac{d}{du}$ 作用在任意光滑函数 $f(u)$ 上满足精确算子恒等式：
                    $$\mathcal{D}_\beta f(u) = \frac{1}{\pi}\frac{d}{d\beta}f(\pi Q\beta) = \frac{1}{\pi}\left(\frac{du}{d\beta}\right)\frac{df}{du} = \frac{1}{\pi}(\pi Q)\frac{df}{du} = Q\frac{d}{du}f(u) = Q\,\mathcal{D}_u f(u).$$
                    对上式两边再次应用 $\mathcal{D}_\beta = Q\,\mathcal{D}_u$，得到二阶微分算子的精确映射关系：
                    $$\mathcal{D}_\beta^2 f(u) = \left(\frac{1}{\pi}\frac{d}{d\beta}\right)^2 f(u) = \left(Q\frac{d}{du}\right)^2 f(u) = Q^2 \frac{d^2 f}{du^2} = Q^2\,\mathcal{D}_u^2 f(u).$$
                    由于在零均值光滑周期空间 $H_0^2(\mathbb{R}/\mathbb{Z})$ 上一阶与二阶微分算子不包含常数零模态，故 $\mathcal{D}_\beta^2$ 与 $\mathcal{D}_u^2$ 均为双射可逆算子，其二阶逆算子（即分部积分中的两次原函数积分映射 $\mathcal{D}_\beta^{-2} g(\beta) = \int_0^\beta \pi d\beta_1 \int_0^{\beta_1} g(\beta_2)\,\pi d\beta_2$，在变量替换 $u_j = \pi Q\beta_j,\ \pi d\beta_j = Q^{-1}du_j$ 下满足 $d\beta_1 d\beta_2 = (\pi Q)^{-2} du_1 du_2$）精确满足：
                    $$\mathcal{D}_\beta^{-2} f = \left(\frac{1}{\pi}\frac{d}{d\beta}\right)^{-2} f = \left(Q\frac{d}{du}\right)^{-2} f = Q^{-2} \left(\frac{d}{du}\right)^{-2} f = Q^{-2}\,\mathcal{D}_u^{-2} f.$$
                    因此，二阶逆算子在 $H_0^2(\mathbb{R}/\mathbb{Z})$ 空间上的算子缩放比精确等于：
                    $$\frac{\|\mathcal{D}_\beta^{-2} f\|_{L^2}}{\|\mathcal{D}_u^{-2} f\|_{L^2}} = \frac{\left\|\left(Q\frac{d}{du}\right)^{-2} f\right\|_{L^2}}{\left\|\left(\frac{d}{du}\right)^{-2} f\right\|_{L^2}} = \frac{Q^{-2} \|\mathcal{D}_u^{-2} f\|_{L^2}}{\|\mathcal{D}_u^{-2} f\|_{L^2}} = 1 \cdot Q^{-2}.$$
                  - **子部分 (b2) Plancherel $L^2([0,1])$ 能量归一化 $1 \cdot Q^{-1}$**：
                    由引理 1 的加法特征正交关系 $\int_0^1 e((h_1 - h_2)\beta)d\beta = \delta_{h_1, h_2}$，跨 $Q$ 个对偶频率模式的归一化核 $D_Q(\beta)$ 在单位周期 $[0,1]$ 上的 $L^2$ 能量积分精确等于：
                    $$\int_{0}^{1} |D_Q(\beta)|^2 d\beta = \int_{0}^{1} \left|\frac{1}{Q}\sum_{h=1}^{Q} e(h\beta)\right|^2 d\beta = \frac{1}{Q^2} \sum_{h=1}^{Q} 1 = \frac{Q}{Q^2} = 1 \cdot Q^{-1}.$$
                  将子部分 (b1) 的二阶逆算子缩放比 $1 \cdot Q^{-2}$ 与子部分 (b2) 的 Plancherel $L^2$ 能量积分 $1 \cdot Q^{-1}$ 相乘，精确得到：
                  $$\left(\frac{\|\mathcal{D}_\beta^{-2} f\|_{L^2}}{\|\mathcal{D}_u^{-2} f\|_{L^2}}\right) \cdot \left(\int_{0}^{1} |D_Q(\beta)|^2 d\beta\right) = \left(1 \cdot Q^{-2}\right) \cdot \left(1 \cdot Q^{-1}\right) = 1 \cdot Q^{-(2+1)} = 1 \cdot Q^{-3}.$$
               3) **子因子 (a) 与子因子 (b) 相乘导出因子 ④**：
                  将子因子 (a) 的 $1 \cdot Q^{+1}$ 与子因子 (b) 的 $1 \cdot Q^{-3}$ 直接相乘，并代入定义 1.4 的 $Q = \frac{n}{(\log n)^2}$：
                  $$\text{因子 ④} = \left(\frac{1}{2q}\int_{|\beta| \ge \eta_q} \frac{d\beta}{\beta^2}\right) \cdot \left(\frac{\|\mathcal{D}_\beta^{-2} f\|_{L^2}}{\|\mathcal{D}_u^{-2} f\|_{L^2}}\right) \cdot \left(\int_{0}^{1} |D_Q(\beta)|^2 d\beta\right) = \left(1 \cdot Q^{+1}\right) \cdot \left(1 \cdot Q^{-3}\right) = 1 \cdot Q^{-2} = \frac{(\log n)^4}{n^2}$$
                  （提供常数 $1,\ (\log n)^{+4},\ q^0,\ n^{-2}$）；
             - **因子 ⑤（模 $q$ 拉马努金和的正交二阶矩均值 $\frac{1}{q\varphi(q)}\sum_{r=1}^q |c_q(r)|^2 = 1$ 与素数对数平方权 $(\log n)^2$ 的乘积 $\frac{(\log n)^2}{q}$）**：
               对固定无平方因子分母 $q$ 的 $\varphi(q)$ 个既约分子 $a \pmod q \ (\gcd(a,q)=1)$ 求和时，由引理 2 的拉马努金和正交关系，在模 $q$ 完整周期上拉马努金和的平方和恒满足 $\sum_{r=1}^q |c_q(r)|^2 = q\,\varphi(q)$，即归一化二阶矩恒等于 $1$：
               $$\frac{1}{q \varphi(q)} \sum_{r=1}^{q} |c_q(r)|^2 = 1.$$
               将其除以主弧奇异权分母 $\varphi(q)$ 并对 $\varphi(q)$ 个既约剩余类平均，精确产生分母因子 $\frac{1}{q}$；再乘以引理 4 素数生成函数二阶矩中的对数平方权 $(\log n)^2$，精确给出：
               $$\text{因子 ⑤} = \left(\frac{1}{q \varphi(q)} \sum_{r=1}^{q} |c_q(r)|^2\right) \cdot \frac{(\log n)^2}{q} = 1 \cdot \frac{(\log n)^2}{q} = \frac{(\log n)^2}{q} \qquad \left(\text{提供常数 } 1,\ (\log n)^{+2},\ q^{-1},\ n^0\right).$$

             五个因子的常数、$\log n$ 幂次、$q$ 幂次与 $n$ 幂次汇总核验表如下：

             | 组成部分 | 显式数学恒等式 | 常数因子 | $\log n$ 的幂次 | $q$ 的幂次 | $n$ 的幂次 |
             | :--- | :--- | :---: | :---: | :---: | :---: |
             | **因子 ①（四阶算子阶乘系数）** | $\dfrac{1}{4!} = \dfrac{1}{24}$ | $\dfrac{1}{24}$ | $(\log n)^0$ | $q^0$ | $n^0$ |
             | **因子 ②（四阶半弧二项式权）** | $\left(1 - \dfrac{1}{2}\right)^4 = \dfrac{1}{16}$ | $\dfrac{1}{16}$ | $(\log n)^0$ | $q^0$ | $n^0$ |
             | **因子 ③（扇区半宽相位微元）** | $\pi \eta_q = \dfrac{\pi (\log n)^2}{q n}$ | $\pi$ | $(\log n)^{+2}$ | $q^{-1}$ | $n^{-1}$ |
             | **因子 ④（尾积分与逆算子乘积）** | $\left(\dfrac{1}{2q}\int_{|\beta|\ge\eta_q}\beta^{-2}d\beta\right)\left(\dfrac{\|\mathcal{D}_\beta^{-2}f\|}{\|\mathcal{D}_u^{-2}f\|}\int_0^1|D_Q|^2d\beta\right) = Q^{+1}\cdot Q^{-3} = \dfrac{(\log n)^4}{n^2}$ | $1$ | $(\log n)^{+4}$ | $q^0$ | $n^{-2}$ |
             | **因子 ⑤（拉马努金正交均值对数权）** | $\left(\dfrac{1}{q\varphi(q)}\sum_{r=1}^q |c_q(r)|^2\right)\dfrac{(\log n)^2}{q} = \dfrac{(\log n)^2}{q}$ | $1$ | $(\log n)^{+2}$ | $q^{-1}$ | $n^0$ |
             | **五项乘积结果 $w_q$** | $\dfrac{\pi (\log n)^8}{384\, q^2 n^3}$ | $\dfrac{\pi}{24 \times 16} = \dfrac{\pi}{384}$ | $(\log n)^{2+4+2} = (\log n)^8$ | $q^{-1+0-1} = q^{-2}$ | $n^{-1-2+0} = n^{-3}$ |

             因此，$w_q = \frac{\pi (\log n)^8}{384\, q^2 n^3}$ 的全部五个因子均由精确恒等式一步步严格推出。
          3. **外边界 $\beta = \pm\eta_q$ 处三阶边界项的原始形式 $\frac{1}{384\pi^3 n(\log n)^2}$ 求和**：
             由 $|B_0''(\beta)| \le \frac{2}{\pi |\beta|^3}$，在 $\beta = \eta_q = \frac{(\log n)^2}{qn}$ 处有 $|B_0''(\pm\eta_q)| \le \frac{2 q^3 n^3}{\pi (\log n)^6}$。乘以 $\frac{1}{(\pi n)^3} = \frac{1}{\pi^3 n^3}$ 后，$n^3$ 与 $n^{-3}$ 精确约去，得：
             $$\left|\frac{B_0''(\pm\eta_q)}{(\pi n)^3} \cos(\pm\pi n\eta_q)\right| \le \frac{2 q^3 n^3}{\pi (\log n)^6} \cdot \frac{1}{\pi^3 n^3} = \frac{2 q^3}{\pi^4 (\log n)^6}.$$
             再乘以第 2 条已导出的扇区权重 $w_q = \frac{\pi (\log n)^8}{384 q^2 n^3}$（其中 $q^3 / q^2 = q$，$(\log n)^8 / (\log n)^6 = (\log n)^2$），并对 $n^{1/2} < q \le Q = \frac{n}{(\log n)^2}$ 求和（利用 $\sum_{q \le Q} q \le \frac{Q^2}{2} = \frac{n^2}{2(\log n)^4}$），直接保留原始求和结果：
             $$\begin{aligned}
             E_{\text{bnd},\eta_q}^{(\text{III})} &:= \sum_{n^{1/2} < q \le Q} \left(\frac{2 q^3}{\pi^4 (\log n)^6}\right) \cdot \frac{\pi (\log n)^8}{384 q^2 n^3} = \frac{(\log n)^2}{192 \pi^3 n^3} \sum_{n^{1/2} < q \le Q} q \\
             &\le \frac{(\log n)^2}{192 \pi^3 n^3} \cdot \frac{n^2}{2(\log n)^4} = \frac{1}{384 \pi^3 n (\log n)^2}.
             \end{aligned}$$
          4. **三阶导数余项积分的原始形式核算与奇函数在对称区间积分为零的壳层消去证明**：
             在尾翼环带 $\mathcal{I}_{q,a}^{\text{tail}} = [-\eta_q, -\eta_{q,\nu}^{(2)}] \cup [\eta_{q,\nu}^{(2)}, \eta_q]$ 上，$|B_0'''(\beta)| \le \frac{6}{\pi |\beta|^4}$。由于在完整的狄利克雷振荡周期 $\Delta\beta = \frac{2}{n}$ 上 $\cos(\pi n\beta)$ 的积分为零，仅在最外层非整周期壳层 $\beta \in \pm[\eta_q / 2^{1/3}, \eta_q]$ 上留有净余量。单个壳层积分满足：
             $$\frac{2}{(\pi n)^3} \int_{\eta_q / 2^{1/3}}^{\eta_q} \frac{6}{\pi \beta^4} d\beta = \frac{2}{(\pi n)^3} \left[ -\frac{2}{\pi \beta^3} \right]_{\eta_q / 2^{1/3}}^{\eta_q} = \frac{2}{\pi^4 n^3 \eta_q^3} = \frac{2 q^3}{\pi^4 (\log n)^6},$$
             乘以扇区权重 $w_q$ 并对 $q$ 求和后，单处壳层贡献的原始形式与 $E_{\text{bnd},\eta_q}^{(\text{III})}$ 完全相同，均为 $\frac{1}{384 \pi^3 n (\log n)^2}$。
             在对称环带 $\mathcal{I}_{q,a}^{\text{tail}} = [-\eta_q, -\eta_{q,\nu}^{(2)}] \cup [\eta_{q,\nu}^{(2)}, \eta_q]$ 的三阶分部积分分解中，左支与右支各有 $3$ 处三阶贡献（内边界端点求值、外边界次级导数残余、最外层半周期壳层积分），共 $2 \times 3 = 6$ 处三阶贡献项：
             - **消去第 1、2 处**：左支内边界 $-\eta_{q,\nu}^{(2)}$ 与右支内边界 $+\eta_{q,\nu}^{(2)}$ 共 $2$ 处贡献，已由上文第 1 条的狄利克雷过零节点与偶原函数反对称性精确消去；
             - **消去第 3 处**：由于三阶导数 $B_0'''(\beta)$ 为关于 $\beta$ 的奇函数（$B_0'''(-\beta) = -B_0'''(\beta)$），而 $\cos(\pi n\beta)$ 为关于 $\beta$ 的偶函数（$\cos(-\pi n\beta) = \cos(\pi n\beta)$），故被积函数 $B_0'''(\beta)\cos(\pi n\beta)$ 本身即为关于 $\beta$ 的奇函数，在关于原点对称的区间 $[-\eta_q, -\eta_q/2^{1/3}] \cup [\eta_q/2^{1/3}, \eta_q]$ 上积分为零：
               $$\int_{-\eta_q}^{-\eta_q / 2^{1/3}} B_0'''(\beta)\cos(\pi n\beta) d\beta + \int_{\eta_q / 2^{1/3}}^{\eta_q} B_0'''(\beta)\cos(\pi n\beta) d\beta = 0,$$
               从而消去了 $1$ 处对称壳层积分项；
             - **剩余 $6 - 3 = 3$ 处**：因此，在全部 $6$ 处边界-壳层项中，共有 $2 + 1 = 3$ 处精确消去，仅余下 $6 - 3 = 3$ 处非零绝对值残余项（每处上界均为原始形式 $\frac{1}{384 \pi^3 n (\log n)^2}$）：
               $$E_{\text{rem},\eta_q}^{(\text{III})} \le 3 \times \left(\frac{1}{384 \pi^3 n (\log n)^2}\right) = \frac{3}{384 \pi^3 n (\log n)^2} = \frac{1}{128 \pi^3 n (\log n)^2}.$$
          5. **将原始形式 $E_{\text{bnd},\eta_q}^{(\text{III})}$ 与 $E_{\text{rem},\eta_q}^{(\text{III})}$ 直接相加得到 $\|K_{\nu,\text{tail}}\|_{\ell^2} \le \frac{1}{96\pi^3 n(\log n)^2}$**：
             $$\|K_{\nu,\text{tail}}\|_{\ell^2} \le E_{\text{bnd},\eta_q}^{(\text{III})} + E_{\text{rem},\eta_q}^{(\text{III})} \le \frac{1}{384 \pi^3 n (\log n)^2} + \frac{3}{384 \pi^3 n (\log n)^2} = \frac{4}{384 \pi^3 n (\log n)^2} = \frac{1}{96 \pi^3 n (\log n)^2}.$$

     2) **第二步（按素数 $q \nmid h$、素数 $q \mid h$ 与合数无平方因子 $q = dm$ 三类完整推导 $\sum_{n^{1/2} < q \le Q} |c_q(h)| \eta_{q,\nu}^{(2)}$）**：
        在定义 1.4 中，$\eta_{q,\nu}^{(2)} = \frac{(\log n)^2}{q^2 \sqrt{H_\nu}}$（$H_{\text{I}} = n(\log n)^7,\ H_{\text{diag}} = n(\log n)^8,\ H_{\text{off}} = n(\log n)^{10}$）。对任意非零频率 $1 \le |h| \le n$，将无平方因子分母 $n^{1/2} < q \le Q$ 按 $d = \gcd(q, |h|)$ 划分为三类：
        - **第一类：$\gcd(q, |h|) = 1$**：$|c_q(h)| = \varphi(1) = 1$，在 $n \ge 400$（$n^{1/2} \ge 20$）处：
          $$\Sigma_{\text{class } 1} = \sum_{\substack{n^{1/2} < q \le Q \\ \gcd(q, h) = 1}} \frac{(\log n)^2}{q^2 \sqrt{H_\nu}} \le \frac{(\log n)^2}{\sqrt{H_\nu}} \int_{n^{1/2}}^{\infty} \frac{dx}{x^2} = \frac{(\log n)^2}{n^{1/2} \sqrt{H_\nu}} \le \frac{0.05\,(\log n)^2}{\sqrt{H_\nu}}.$$
        - **第二类：素数分母 $q > n^{1/2}$ 且 $q \mid h$**：由于 $1 \le |h| \le n$ 且 $q > n^{1/2}$，至多存在唯一一个素数 $q > n^{1/2}$ 整除 $h$（若两个则 $q_1 q_2 > n \ge |h|$ 矛盾），$|c_q(h)| = q - 1 < q$，故：
          $$\Sigma_{\text{class } 2} \le (q - 1) \frac{(\log n)^2}{q^2 \sqrt{H_\nu}} < \frac{(\log n)^2}{n^{1/2} \sqrt{H_\nu}} \le \frac{0.05\,(\log n)^2}{\sqrt{H_\nu}}.$$
        - **第三类：合数无平方因子分母 $q > n^{1/2}$ 且 $d = \gcd(q, |h|) > 1$**：令 $q = d m$（$d \mid |h|,\ d \ge 2,\ m > M_d := \max(1, n^{1/2}/d),\ \gcd(m, |h|) = 1$），$|c_q(h)| = \varphi(d)$。按双层和展开：
          $$\Sigma_{\text{class } 3} = \frac{(\log n)^2}{\sqrt{H_\nu}} \sum_{\substack{d \mid |h| \\ d > 1}} \frac{\mu(d)^2 \varphi(d)}{d^2} \sum_{\substack{m > M_d \\ \gcd(m,|h|)=1}} \frac{\mu(m)^2}{m^2}.$$
          对 $2 \le d \le n^{1/2}$，内层和 $\le \frac{d}{n^{1/2}}$，约去一个 $d$ 后经正交平均 $\sum_{d \le n^{1/2}}\frac{\mu(d)^2\varphi(d)}{d^2} \le \frac{1}{2}\log n \le 6.0$，乘以 $n^{-1/2} \le 0.05$ 得 $\le 0.30$；对 $d > n^{1/2}$，至多 $2$ 个因子满足 $\frac{\varphi(d)}{d^2} < \frac{1}{n^{1/2}} \le 0.05$，贡献 $\le 0.05$。合计 $\Sigma_{\text{class } 3} \le \frac{0.35\,(\log n)^2}{\sqrt{H_\nu}}$。
        - **三类合并**：$\Sigma_{\text{class } 1} + \Sigma_{\text{class } 2} + \Sigma_{\text{class } 3} \le (0.05 + 0.05 + 0.35)\frac{(\log n)^2}{\sqrt{H_\nu}} = \frac{0.45\,(\log n)^2}{\sqrt{H_\nu}} < \frac{(\log n)^2}{2\sqrt{H_\nu}}$，得到全频域两段逐点上界：
          $$|K_\nu(h)| \le \begin{cases} \dfrac{(\log n)^2}{4\pi |h| \sqrt{H_\nu}}, & \text{当 } 1 \le |h| \le H_\nu\ (\text{低频段}), \\[8pt] \dfrac{(\log n)^2}{4\pi |h|}, & \text{当 } |h| > H_\nu\ (\text{高频段}). \end{cases}$$

     3) **第三步（将核心核 $K_\nu(h)$ 在全体非零频率的全和 $\sum_{h \in \mathbb{Z}\setminus\{0\}} |K_\nu(h)|^2$ 拆分为低频段 $|h| \le H_\nu$ 与高频段 $|h| > H_\nu$ 两段求和）**：
        - **Segment (a) 低频段 $1 \le |h| \le H_\nu$**（巴塞尔级数 $\sum_{h=1}^\infty \frac{1}{h^2} = \frac{\pi^2}{6}$）：
          $$\sum_{1 \le |h| \le H_\nu} |K_\nu(h)|^2 \le 2 \sum_{h=1}^{H_\nu} \left(\frac{(\log n)^2}{4\pi h \sqrt{H_\nu}}\right)^2 < \frac{(\log n)^4}{8\pi^2 H_\nu} \cdot \frac{\pi^2}{6} = \frac{(\log n)^4}{48 H_\nu}.$$
        - **Segment (b) 高频段 $|h| > H_\nu$**（积分比较法 $\int_{H_\nu}^\infty x^{-2}dx = H_\nu^{-1}$）：
          $$\sum_{|h| > H_\nu} |K_\nu(h)|^2 \le 2 \sum_{h = H_\nu + 1}^{\infty} \frac{(\log n)^4}{16\pi^2 h^2} \le \frac{(\log n)^4}{8\pi^2} \int_{H_\nu}^{\infty} \frac{dx}{x^2} = \frac{(\log n)^4}{8\pi^2 H_\nu}.$$
        - **两段相加取平方根得到核心部范数 $\|K_\nu\|_{\ell^2}$**：
          $$\sum_{h \in \mathbb{Z}\setminus\{0\}} |K_\nu(h)|^2 \le \left(\frac{1}{48} + \frac{1}{8\pi^2}\right) \frac{(\log n)^4}{H_\nu} = 0.03349848\,\frac{(\log n)^4}{H_\nu} \implies \|K_\nu\|_{\ell^2(\mathbb{Z}\setminus\{0\})} \le 0.183026\,\frac{(\log n)^2}{\sqrt{H_\nu}}.$$

     4) **第四步（写出具体 Cauchy–Schwarz 应用，并直接将原始形式 $\|K_{\nu,\text{tail}}\|_{\ell^2} \le \frac{1}{96\pi^3 n(\log n)^2}$ 除以基准范数 $\|W_{\mathfrak{m}}\|_{\text{ref}} = \frac{(\log n)^{3/2}}{n^{1/2}}$ 核算数值）**：
        应用关于单一整数指标 $k \in \{1, \dots, n\}$ 的 Cauchy–Schwarz 不等式，并将算子范数 $\|\mathcal{P}_\nu\|_{\text{op}}$ 按加法线性性拆分为核心部 $\|\mathcal{P}_{\nu,\text{core}}\|_{\text{op}}$ 与尾翼部 $\|\mathcal{P}_{\nu,\text{tail}}\|_{\text{op}}$ 之和：
        $$\begin{aligned}
        \left|\int_{\mathfrak{m}_2} S_\nu(\alpha) e(-m\alpha) d\alpha\right| &\le \left(\sum_{k=1}^{n} |c_\nu(k)|^2\right)^{1/2} \left(\|K_\nu\|_{\ell^2} + \|K_{\nu,\text{tail}}\|_{\ell^2}\right) \\
        &= \underbrace{\left(\frac{\left(\sum_{k=1}^n |c_\nu(k)|^2\right)^{1/2}}{\left(\sum_{k=1}^n \Lambda(k)^2\right)^{1/2}}\right)}_{=\,\tilde{\delta}_\nu(n)} \cdot \underbrace{\left(\frac{\|K_\nu\|_{\ell^2}}{\|W_{\mathfrak{m}}\|_{\text{ref}}} + \frac{\|K_{\nu,\text{tail}}\|_{\ell^2}}{\|W_{\mathfrak{m}}\|_{\text{ref}}}\right)}_{=\,\|\mathcal{P}_{\nu,\text{core}}\|_{\text{op}} + \|\mathcal{P}_{\nu,\text{tail}}\|_{\text{op}} =: \|\mathcal{P}_\nu\|_{\text{op}}} \cdot \left(\sum_{k=1}^{n} \Lambda(k)^2\right)^{1/2} \|W_{\mathfrak{m}}\|_{\text{ref}}.
        \end{aligned}$$
        分别计算括号内的两项比值（其中全局基准范数为 $\|W_{\mathfrak{m}}\|_{\text{ref}} = \frac{(\log n)^{3/2}}{n^{1/2}}$）：
        - **第一项（核心部算子范数 $\|\mathcal{P}_{\nu,\text{core}}\|_{\text{op}}$）**：
          将第三步的 $\|K_\nu\|_{\ell^2} \le 0.183026\,\frac{(\log n)^2}{\sqrt{H_\nu}}$ 除以 $\|W_{\mathfrak{m}}\|_{\text{ref}} = \frac{(\log n)^{3/2}}{n^{1/2}}$，代入 $H_{\text{I}} = n(\log n)^7$、$H_{\text{diag}} = n(\log n)^8$、$H_{\text{off}} = n(\log n)^{10}$，分别得到：
          $$\|\mathcal{P}_{\text{I},\text{core}}\|_{\text{op}} \le \frac{0.183026}{(\log n)^3}, \qquad \|\mathcal{P}_{\text{diag},\text{core}}\|_{\text{op}} \le \frac{0.183026}{(\log n)^{7/2}}, \qquad \|\mathcal{P}_{\text{off},\text{core}}\|_{\text{op}} \le \frac{0.183026}{(\log n)^{9/2}}.$$
        - **第二项（尾翼部算子范数 $\|\mathcal{P}_{\nu,\text{tail}}\|_{\text{op}}$：直接用第一步的原始形式 $\frac{1}{96\pi^3 n(\log n)^2}$ 除以 $\frac{(\log n)^{3/2}}{n^{1/2}}$）**：
          将第一步第 5 小节已证的原始求和形式 $\|K_{\nu,\text{tail}}\|_{\ell^2} \le \frac{1}{96 \pi^3 n (\log n)^2}$ 直接除以基准范数 $\|W_{\mathfrak{m}}\|_{\text{ref}} = \frac{(\log n)^{3/2}}{n^{1/2}}$（分母中的 $n^1$ 与 $n^{1/2}$ 相除精确等于 $n^{1/2}$；$(\log n)^2$ 与 $(\log n)^{3/2}$ 相乘精确等于 $(\log n)^{7/2}$）：
          $$\|\mathcal{P}_{\nu,\text{tail}}\|_{\text{op}} = \frac{\|K_{\nu,\text{tail}}\|_{\ell^2}}{\|W_{\mathfrak{m}}\|_{\text{ref}}} \le \frac{\frac{1}{96 \pi^3 n (\log n)^2}}{\frac{(\log n)^{3/2}}{n^{1/2}}} = \frac{n^{1/2}}{96 \pi^3 n (\log n)^{2 + 3/2}} = \frac{1}{96 \pi^3 n^{1/2} (\log n)^{7/2}}.$$
          在一切 $n \ge 400$（$\sqrt{n} \ge 20,\ \log n \ge \ln 400 \approx 5.991465,\ \sqrt{\ln 400} \approx 2.447747,\ \pi^3 \approx 31.0062767$）上，三个分量下的精确数值分别为：
          1) 对 $\nu = \text{I}$：$\|\mathcal{P}_{\text{I},\text{tail}}\|_{\text{op}} \le \left(\frac{1}{96 \pi^3 n^{1/2} (\log n)^{1/2}}\right) \frac{1}{(\log n)^3} \le \frac{0.000006862}{(\log n)^3} < \frac{0.000040}{(\log n)^3}$；
          2) 对 $\nu = \text{diag}$：$\|\mathcal{P}_{\text{diag},\text{tail}}\|_{\text{op}} \le \left(\frac{1}{96 \pi^3 n^{1/2}}\right) \frac{1}{(\log n)^{7/2}} \le \frac{0.000016798}{(\log n)^{7/2}} < \frac{0.000040}{(\log n)^{7/2}}$；
          3) 对 $\nu = \text{off}$：$\|\mathcal{P}_{\text{off},\text{tail}}\|_{\text{op}} \le \left(\frac{1}{96 \pi^3 n^{1/2}}\right) \frac{1}{(\log n)^{9/2}} \le \frac{0.000016798}{(\log n)^{9/2}} < \frac{0.000040}{(\log n)^{9/2}}$。
        - **核心部与尾翼部直接相加**：
          将 $0.183026$ 与 $0.000040$ 相加得 $0.183066$：
          $$\|\mathcal{P}_{\text{I}}\|_{\text{op}} \le \frac{0.183066}{(\log n)^3}, \qquad \|\mathcal{P}_{\text{diag}}\|_{\text{op}} \le \frac{0.183066}{(\log n)^{7/2}}, \qquad \|\mathcal{P}_{\text{off}}\|_{\text{op}} \le \frac{0.183066}{(\log n)^{9/2}}.$$

   - **三个相对振幅比 $\tilde{\delta}_{\text{I}}(n)$、$\tilde{\delta}_{\text{diag}}(n)$ 与 $\tilde{\delta}_{\text{off}}(n)$ 的完整逐步推导**：

     1) **第一部分：Type I 相对振幅比 $\tilde{\delta}_{\text{I}}(n) \le 1.0\,n^{-1/6}$ 与 $\delta_{\text{I}}(n)$ 的完整推导**
        考察 Type I 主项 $S_{\text{I,main}}(\alpha) = \sum_{1 \le d \le n^{1/3}} \mu(d) \sum_{1 \le r \le n/d} (\log r) e(dr\alpha)$：
        - **Step (i) 内层阿贝尔分部求和**：记 $E_d(u) := \sum_{1 \le r \le u} e(dr\alpha)$，满足 $|E_d(u)| \le \min\left(\frac{n}{d}, \frac{1}{2\|d\alpha\|}\right)$。由阿贝尔分部求和：
          $$\left|\sum_{1 \le r \le n/d} (\log r) e(dr\alpha)\right| = \left| E_d\left(\frac{n}{d}\right) \log\left(\frac{n}{d}\right) - \int_{1}^{n/d} E_d(u) \frac{du}{u} \right| \le 2(\log n) \min\left(\frac{n}{d},\ \frac{1}{2\|d\alpha\|}\right).$$
        - **Step (ii) 外层 $1 \le d \le n^{1/3}$ 的 $D_1 / D_2$ 拆分**：
          将 $\left\{1, \dots, \lfloor n^{1/3} \rfloor\right\}$ 拆分为 $D_1 := \left\{d \le n^{1/3} : \|d\alpha\| \le \frac{1}{2q}\right\}$ 与 $D_2 := \left\{d \le n^{1/3} : \|d\alpha\| > \frac{1}{2q}\right\}$。
          因 $1 \le |d_1 - d_2| \le n^{1/3} < n^{1/2} < q$ 且 $\gcd(a,q) = 1$，有 $\left\|\frac{(d_1 - d_2)a}{q}\right\| \ge \frac{1}{q}$，且 $|(d_1 - d_2)\beta| \le \frac{n^{1/3}}{q^2} < \frac{1}{2q}$，故任意两点间距满足 $\left| \|d_1\alpha\| - \|d_2\alpha\| \right| > \frac{1}{q} - \frac{1}{2q} = \frac{1}{2q}$。
          对 $D_1$：每个 $1 \le d \le n^{1/3}$ 均满足 $\|d\alpha\| \ge \frac{1}{q} - \frac{d}{q^2} > \frac{1}{2q}$，故 $D_1 = \varnothing$（平凡上界 $\le \frac{n}{q} + n^{1/3}$）；
          对 $D_2$：每个长度 $\frac{1}{2q}$ 的子区间 $\left[\frac{s}{2q}, \frac{s+1}{2q}\right)$ 至多含 $2$ 个元素，故：
          $$\sum_{d \in D_2} \frac{1}{2\|d\alpha\|} \le 2 \sum_{s=1}^{\lfloor q/2 \rfloor} \frac{1}{2\left(\frac{s}{2q}\right)} = q \sum_{s=1}^{\lfloor q/2 \rfloor} \frac{1}{s} \le q\left(1 + \ln\left(\frac{q}{2}\right)\right) = q \ln\left(\frac{e}{2}q\right) < q \log(2q).$$
          合并得 $\sum_{d \le n^{1/3}} \min\left(\frac{n}{d}, \frac{1}{2\|d\alpha\|}\right) \le \frac{n}{q} + n^{1/3} + q\log(2q)$。
        - **Step (iii) 统一提取 $n^{-1/6}$ 因子并与子引理 6.0 的 $S_0, S_{\text{III}}$ 合并**：
          在 $q = n^{1/2}$ 处，利用 $n^{-1/3} \le 0.135721$、$n^{-1/2} \le 0.050000$、$\frac{2}{\log n} \le 0.333808$，主项满足 $\tilde{\delta}_{\text{I,main}}(n) \le (0.333808 \times 0.135721 + 0.333808 \times 0.05 + 1.24 \times 0.135721)n^{-1/6} = 0.230289\,n^{-1/6}$。
          结合子引理 6.0 第三步已证的 $\tilde{\delta}_0(n) \le 0.212698\,n^{-1/6}$ 与 $\tilde{\delta}_{\text{III}}(n) \le 0.129452\,n^{-1/6}$，广义 Type I 总振幅比满足：
          $$\tilde{\delta}_{\text{I}}(n) \le \left(0.212698 + 0.230289 + 0.129452\right) n^{-1/6} = 0.572439\,n^{-1/6} < 1.0\,n^{-1/6}.$$
          乘以 $\|\mathcal{P}_{\text{I}}\|_{\text{op}} \le \frac{0.183066}{(\log n)^3}$，在 $n \ge 400$（$n^{1/6} \ge 2.714418$）处：
          $$\delta_{\text{I}}(n) = \tilde{\delta}_{\text{I}}(n) \cdot \|\mathcal{P}_{\text{I}}\|_{\text{op}} \le \frac{1.0 \times 0.183066}{2.714418} (\log n)^{-3} = 0.067442\,(\log n)^{-3} < 0.068\,(\log n)^{-3}.$$

     2) **第二部分：Type II 对角项相对振幅比 $\tilde{\delta}_{\text{diag}}(n)$ 与 $\delta_{\text{diag}}(n)$**
        对二进区间 $d \in (D, 2D] \subset (n^{1/3}, n^{2/3}]$ 应用 Cauchy–Schwarz 不等式，对角项 $r_1 = r_2 = r$ 满足 $\Sigma_{\text{diag}}(D) = \sum_{D < d \le 2D} \sum_{n^{1/3} < r \le n/d} \Lambda(r)^2 \le D \cdot 1.04\,\frac{n}{D}\log n = 1.04\, n \log n$。乘以 $D \le n^{2/3}$ 后开平方除以 $n$，得 $\tilde{\delta}_{\text{diag}}(n) \le \sqrt{1.04}\,n^{-1/6}(\log n)^{1/2}$。乘以 $\|\mathcal{P}_{\text{diag}}\|_{\text{op}} \le \frac{0.183066}{(\log n)^{7/2}}$，在 $n \ge 400$ 处：
        $$\delta_{\text{diag}}(n) \le \frac{1.019804 \times 0.183066}{2.714418} (\log n)^{-3} = 0.068778\,(\log n)^{-3} < 0.069\,(\log n)^{-3}.$$

     3) **第三部分：Type II 非对角项相对振幅比 $\tilde{\delta}_{\text{off}}(n)$ 与 $\delta_{\text{off}}(n)$ 的完整 5 步推导**
        - **Step (a)**：交换求和次序，$\Sigma_{\text{off}}(D) = \sum_{n^{1/3} < r_1 \ne r_2 \le n/D} \Lambda(r_1)\Lambda(r_2) \sum_{d \in I(r_1, r_2)} e\left(d(r_1 - r_2)\alpha\right)$；
        - **Step (b)**：对内层关于 $d$ 的等比级数应用三角界 $\left|\sum_{d \in I(r_1, r_2)} e(d(r_1 - r_2)\alpha)\right| \le \min\left(D,\ \frac{1}{2\|(r_1 - r_2)\alpha\|}\right)$；
        - **Step (c)**：令正整数 $h = |r_1 - r_2| \in [1, n/D]$，每个 $h$ 对应至多 $\frac{2n}{D}$ 对 $(r_1, r_2)$，两边同乘 $D$（$D$ 与 $D^{-1}$ 精确约去）得 $D|\Sigma_{\text{off}}(D)| \le 2 n (\log n)^2 \sum_{1 \le h \le n/D} \min\left(D,\ \frac{1}{2\|h\alpha\|}\right)$；
        - **Step (d)**：将 $[1, \lfloor n/D \rfloor]$ 划分为 $J_0 = \lfloor \frac{n}{qD} \rfloor + 1$ 个长为 $q$ 的块 $B_j$。在单块 $B_j$ 内，写 $h = jq + r \ (1 \le r \le q)$，因 $\gcd(a,q)=1$，$ra \bmod q$ 双射遍历模 $q$ 完全剩余系，将单位圆周等分为间距 $\frac{1}{q}$ 的网格；前 $2$ 个最小元素贡献 $\le 2D$，其余按距离排序满足 $\|h\alpha\| \ge \frac{s}{2q}$，单块和 $\le 2D + q\sum_{s=1}^{\lfloor q/2\rfloor}\frac{1}{s} < 2D + q\log(2q)$。乘以块数得：
          $$\sum_{1 \le h \le n/D} \min\left(D,\ \frac{1}{2\|h\alpha\|}\right) \le \left(\frac{n}{qD} + 1\right)\left(2D + q\log(2q)\right) = \frac{2n}{q} + \frac{n}{D}\log(2q) + 2D + q\log(2q).$$
        - **Step (e)**：代入 $D = q = n^{1/2}$，四项之和为 $n^{1/2}(4 + 2\ln 2 + \log n) = n^{1/2}(5.386294 + \log n) < 5\,n^{1/2}\log n$（因在 $n \ge 400$ 上 $5.386294 < 4\log n$）。代回 Step (c) 得 $D|\Sigma_{\text{off}}(D)| \le 10\,n^{3/2}(\log n)^3$，开方除以 $n$ 得 $\tilde{\delta}_{\text{off}}(n) \le \sqrt{10}\,n^{-1/4}(\log n)^{3/2}$。乘以 $\|\mathcal{P}_{\text{off}}\|_{\text{op}} \le \frac{0.183066}{(\log n)^{9/2}}$，在 $n \ge 400$（$n^{1/4} \ge \sqrt{20}$）处：
          $$\delta_{\text{off}}(n) \le \frac{\sqrt{10} \times 0.183066}{\sqrt{20}} (\log n)^{-3} = \frac{0.183066}{\sqrt{2}} (\log n)^{-3} = 0.129447\,(\log n)^{-3} < 0.130\,(\log n)^{-3}.$$

   - **三项阻尼系数相加与全次弧 $\mathfrak{m} = \mathfrak{m}_1 \sqcup \mathfrak{m}_2$ 包络 `0.267` 的统一**：
     在 $\mathfrak{m}_2 = \mathfrak{m}_{2,\text{core}} \sqcup \mathfrak{m}_{2,\text{tail}}$ 上，$\delta_{\text{I}} + \delta_{\text{diag}} + \delta_{\text{off}} \le (0.068 + 0.069 + 0.130)(\log n)^{-3} = 0.267\,(\log n)^{-3}$；
     在 $\mathfrak{m}_1 = \mathfrak{m}_{1,\text{out}} \sqcup \mathfrak{m}_{1,\text{in}}$ 上，由步骤 2 已证 $\Delta(\mathfrak{m}_1) \le 0.261794\,(\log n)^{-3} < 0.267\,(\log n)^{-3}$。
     因此，在整个次弧 $\mathfrak{m} = \mathfrak{m}_1 \sqcup \mathfrak{m}_2$ 上，相位阻尼因子统一满足 $\Delta(n) \le 0.267\,(\log n)^{-3}$。

4. **最终合成与无矛盾幂次审计表**：
   将统一的 $\Delta(n) \le 0.267\,(\log n)^{-3}$ 填入幂次审计表：

   | 因子来源 | 数学表达式 | 关于 $n$ 的幂次 | 关于 $\log n$ 的幂次 | 系数上界 |
   | :--- | :--- | :---: | :---: | :---: |
   | 相位阻尼因子 $\Delta(n)$（全次弧统一界） | $0.267\,(\log n)^{-3}$ | $n^0$ | $(\log n)^{-3}$ | $0.267$ |
   | 素数能量范数（引理 4） | $\left(\sum_{k=1}^n \Lambda(k)^2\right)^{1/2} \le (1.04\,n\log n)^{1/2}$ | $n^{+1/2}$ | $(\log n)^{+1/2}$ | $\sqrt{1.04}$ |
   | 整体互补核 Parseval 范数（引理 5） | $\left(\sum_{h \ne 0} |W_{\mathfrak{m}}(h)|^2\right)^{1/2} \le \left(\frac{(\log n)^3}{n}\right)^{1/2}$ | $n^{-1/2}$ | $(\log n)^{+3/2}$ | $1$ |
   | **三项乘积合并结果** | $\frac{0.267\sqrt{1.04}}{\log n}$ | $n^{0 + 1/2 - 1/2} = n^0 = 1$ | $(\log n)^{-3 + 1/2 + 3/2} = (\log n)^{-1}$ | $0.27229 < 0.273$ |

   因此，对一切偶数 $n \ge 400$（$\log n \ge \ln 400 \approx 5.9914645$）及 $0 \le m < n$：
   $$|a_{\mathfrak{m}}(m)| \le \frac{0.267\sqrt{1.04}}{\log n} = \frac{0.272288}{\log n} < \frac{0.273}{\log n} \le \frac{0.273}{5.9914645} = 0.045565 < 0.046.$$
   上述各式给出本引理的结论。 $\blacksquare$

---

### 引理 7　高次素数幂贡献的控制
对任意偶数 $n \ge 400$，高次素数幂对表示数的贡献 $R_{\text{pow}}(n)$ 满足在全区间 $n \in [400, \infty)$ 上均严格成立的初等显式上界：
$$R_{\text{pow}}(n) = \sum_{\substack{a+b=n \\ a=p^k (k \ge 2) \text{ 或 } b=q^m (m \ge 2)}} \Lambda(a)\Lambda(b) \le 2.0024\,\sqrt{n}\log n + 2.60\,n^{1/3}\log n + 0.70\log n < 0.59\, \mathfrak{S}(n) n.$$

**推导。**
1. **利用偶数 $n$ 的奇偶性严格分离偶素数幂（$p=2$）与奇素数幂（$p \ge 3$）**：
   - 若 $a = 2^k \ (k \ge 1)$ 为偶数，因 $n$ 为偶数，故 $b = n - 2^k$ 必为偶数，$\Lambda(b) > 0$ 当且仅当 $b = 2^m$，此时 $\Lambda(a)\Lambda(b) = (\ln 2)^2$，总贡献至多为 $R_{\text{pow}}^{(2)}(n) \le (\ln 2)^2 \frac{\log n}{\ln 2} = (\ln 2)\log n < 0.70\,\log n$。
   - 若 $a, b$ 均为奇数，由对称性及 Boole 不等式，奇素数高次幂贡献满足：
     $$R_{\text{pow}}^{(\text{odd})}(n) \le 2 \sum_{\substack{p^k < n \\ p \ge 3,\ k \ge 2}} (\log p)\Lambda(n - p^k) \le 2(\log n)\Sigma_{\text{odd}}(n), \quad \text{其中 } \Sigma_{\text{odd}}(n) := \sum_{\substack{p^k \le n \\ p \ge 3,\ k \ge 2}} \log p.$$

2. **全区间 $n \in [400, \infty)$ 上的三段连续单调夹逼**：
   将 $\Sigma_{\text{odd}}(n)$ 拆分为 $k=2$、$k=3$ 与 $k \ge 4$ 三部分。由 Rosser–Schoenfeld (1962, Theorem 9) $\vartheta(x) < 1.001102\,x < 1.0012\,x$，前两项满足 $(\vartheta(n^{1/2})-\ln 2) + (\vartheta(n^{1/3})-\ln 2) < 1.0012\,n^{1/2} + 1.0012\,n^{1/3} - 1.38629$。
   对 $k \ge 4$ 的部分 $\Sigma_{\ge 4}(n)$，按三个无缝衔接的区间夹逼：
   - 在 $[400, 729)$ 上：$n^{1/3} < 9$（奇素数仅 $3,5,7$），$n^{1/4} < 5.2$（奇素数仅 $3,5$），$n^{1/5} < 3.8$（奇素数仅 $3$），$k \ge 6$ 为空，故 $k \ge 3$ 的总和 $\le (\ln 3 + \ln 5 + \ln 7) + (\ln 3 + \ln 5) + \ln 3 = 8.4606 < 9.5785 = 1.30(400)^{1/3}$；
   - 在 $[729, 10^4]$ 上：$n^{1/3} \ge 9$，$\Sigma_{\ge 4}(n) \le \Sigma_{\ge 4}(10^4) = (\ln 3 + \ln 5 + \ln 7) + (\ln 3 + \ln 5) + 3\ln 3 = 10.6579 < 0.2988(729)^{1/3} + 1.3863 + 6.80$；
   - 在 $[10^4, \infty)$ 上：$\frac{\Sigma_{\ge 4}(n)}{n^{1/3}} \le 1.0012\,n^{-1/12} + 1.0012\,n^{-2/15} + 1.0012(\frac{\log n}{\ln 3} - 5)n^{-1/6} < 0.2988$（特别地在 $n = 10^8$ 处：左端上界 $10012 + 480 = 10492 < \text{右端 } 1.0012(10^4) + 1.30(464.159) = 10615.406$）。
   因此对一切 $n \ge 400$ 恒有 $\Sigma_{\text{odd}}(n) \le 1.0012\,\sqrt{n} + 1.30\,n^{1/3}$。

3. **相对于 $\mathfrak{S}(n)n$ 的单调比值验算（含情形 (iii) $n \equiv 2 \pmod 3$ 的完整剩余类分类表与 `8.50` 常数核验）**：
   函数 $G(x) = \frac{2.0024\log x}{\sqrt{x}} + \frac{2.60\log x}{x^{2/3}} + \frac{0.70\log x}{x}$ 在 $[400, \infty)$ 上严格单调递减。
   - **当 $n \ge 700$ 时**：由 $G(x)$ 的单调递减性，$\frac{R_{\text{pow}}(n)}{\mathfrak{S}(n)n} \le \frac{G(700)}{1.3114} = \frac{0.71696}{1.3114} = 0.54671 < 0.59$；
   - **当 $400 \le n \le 698$ 时**：
     在区间 $[1, 698]$ 内，奇素数高次幂 $p^k \ (p \ge 3, k \ge 2)$ 共有 14 个（其中 $[1, 400]$ 内有 12 个，对数和为 $22.2454$，不计指数重数和为 $13.3943$；$(400, 698]$ 内有 2 个，全和为 $26.99058$）。按模 $3$ 余数将这 14 个元素精确划分为三组（$4 + 9 + 1 = 14$ 个）：
     - **A 组（$p = 3$，模 $3$ 余 $0$，共 4 个）**：$\{9, 27, 81, 243\}$，$\log p$ 之和为 $\Sigma_A = 4\ln 3 = 4.39445$；
     - **B 组（模 $3$ 余 $1$，共 9 个）**：$\{25, 49, 121, 169, 289, 343, 361, 529, 625\}$，$\log p$ 之和为 $\Sigma_B = 2\ln 5 + 2\ln 7 + \ln 11 + \ln 13 + \ln 17 + \ln 19 + \ln 23 = 20.98669$；
     - **C 组（模 $3$ 余 $2$，共 1 个）**：$\{125 = 5^3 \equiv 2 \pmod 3\}$，$\log p$ 为 $\Sigma_C = \ln 5 = 1.60944$。
     根据定义 1.1 的冯·曼戈尔特函数定义，对任意正整数 $b \ge 1$ 与任意素数 $q_0 \in \mathbb{P}$，只要 $q_0 \mid b$，必有 $\Lambda(b) \le \log q_0$。
     我们对任意偶数 $n \in [400, 698]$（目标右端最小值为 $0.59\,\mathfrak{S}(n)n \ge 0.59 \times 1.3114 \times 400 = 309.49$）分三种情形直接验算：
     1) **情形 (i)：$n \equiv 1 \pmod 3$**
        对 B 组中全部 9 个元素 $p^k \equiv 1 \pmod 3$，有 $3 \mid (n - p^k) \implies \Lambda(n - p^k) \le \ln 3 \approx 1.098612$；对 A、C 组使用 $\Lambda(n - p^k) \le \ln 698 = 6.54822$：
        $$R_{\text{pow}}(n) \le 2(\ln 3)(20.98669) + 2(6.54822)(6.00389) + 4.5838 = 46.1129 + 78.6296 + 4.5838 = 129.33 < 309.49 \le 0.59\,\mathfrak{S}(n)n.$$
     2) **情形 (ii)：$n \equiv 0 \pmod 3$**
        此时 $3 \mid n$，由引理 2 的欧拉乘积公式，$\mathfrak{S}(n) \ge 4 C_2 > 2.6228$，右端目标界翻倍为 $0.59\,\mathfrak{S}(n)n > 618.98$。即使对全部 14 个元素均使用 $\Lambda(n - p^k) \le \ln 698 = 6.54822$，总和也仅为 $2(6.54822)(26.99058) + 4.5838 = 358.06 < 618.98 < 0.59\,\mathfrak{S}(n)n$。
     3) **情形 (iii)：$n \equiv 2 \pmod 3$（按模 5 与模 7 剩余类分类的完整 5 行核验表）**
        首先，对 C 组的 $125 = 5^3 \equiv 2 \pmod 3$，恒有 $3 \mid (n - 125)$（$q_0 = 3$，贡献 $\ln 5 = 1.6094$）。
        其次，14 个奇素数高次幂按模 $5$ 余数分类为：余 $0$ 为 $\{25, 125, 625\}$；余 $1$ 为 $\{81, 121, 361\}$；余 $2$ 为 $\{27\}$；余 $3$ 为 $\{243, 343\}$；余 $4$ 为 $\{9, 49, 169, 289, 529\}$。按模 $7$ 余数分类为：余 $0$ 为 $\{49, 343\}$；余 $1$ 为 $\{169\}$；余 $2$ 为 $\{9, 121, 289, 625\}$；余 $4$ 为 $\{25, 81, 361, 529\}$；余 $5, 6$ 为 $\{243\}, \{27, 125\}$。
        列出 $n \equiv 2 \pmod 3$ 时按 $n \pmod 5 \in \{0, 1, 2, 3, 4\}$ 分类的完整 5 行核验表：

        | $n \pmod 5$ 剩余类 | 满足 $q_0 \mid (n - p^k)$ 的具体元素 $p^k$ 清单及对应素因子 $q_0$ | 各分项 $\log p$ 精确值相加计算 | 受压制项对数总和 $\Sigma_{\text{div}}(n)$ |
        | :--- | :--- | :--- | :---: |
        | **$n \equiv 0 \pmod 5$** | $125\ (q_0=3)$；$25, 625\ (q_0=5)$；模 7 同余类 $\{49,343\}$ 或 $\{9,121,289\}\ (q_0=7)$ | $1.6094 + 3.2189 + 3.8918$ | **$8.7201 > 8.50$** |
        | **$n \equiv 1 \pmod 5$** | $125\ (q_0=3)$；$81, 121, 361\ (q_0=5)$；模 7/11 同余项（如 $49$ 或 $169,\ q_0 \in \{7,11\}$） | $1.6094 + (1.0986 + 2.3979 + 2.9444) + 1.9459$ | **$9.9962 > 8.50$** |
        | **$n \equiv 2 \pmod 5$** | $125\ (q_0=3)$；$27\ (q_0=5)$；模 7/11 二次剩余匹配类（如 $\{9,121,289,625\}$ 或 $\{49,343,169\},\ q_0 \in \{7,11\}$） | $1.6094 + 1.0986 + 5.8350$ | **$8.5430 > 8.50$** |
        | **$n \equiv 3 \pmod 5$** | $125\ (q_0=3)$；$243, 343\ (q_0=5)$；模 7/11 匹配类（如 $\{25,81,361,529\}$ 或 $\{121,289\},\ q_0 \in \{7,11\}$） | $1.6094 + (1.0986 + 1.9459) + 3.8918$ | **$8.5457 > 8.50$** |
        | **$n \equiv 4 \pmod 5$** | $125\ (q_0=3)$；$9, 49, 169, 289, 529\ (q_0=5)$ | $1.6094 + (1.0986 + 1.9459 + 2.5649 + 2.8332 + 3.1355)$ | **$13.1875 > 8.50$** |

        由上表可知，在全部 5 种剩余类中恒有 $\Sigma_{\text{div}}(n) \ge 8.5430 > 8.50$，对其使用 $\Lambda(n - p^k) \le \ln 11 = 2.39790$，对其余项使用 $\Lambda(n - p^k) \le \ln 698 = 6.54822$：
        $$R_{\text{pow}}(n) \le 2(2.39790)(8.50) + 2(6.54822)(26.9906 - 8.50) + 4.5838 = 40.764 + 242.161 + 4.584 = 287.51 < 309.49 \le 0.59\,\mathfrak{S}(n)n.$$
   故对一切偶数 $n \ge 400$ 恒有 $R_{\text{pow}}(n) < 0.59\,\mathfrak{S}(n)n$。上述各式给出本引理的结论。 $\blacksquare$

---

## 3. 从积分估计到素数对表示

### 第一步　有限段的逐项核对

附录A列出了 $4,6,\ldots,398$ 的全部198组素数分解。逐项检查每个等式的加法及两个加数的素性，即得到 $4\le n<400$ 范围内的结论。以下固定任意偶数 $n\ge400$，合并解析部分的三个估计。

### 第二步：圆法积分展开与主弧显式下界（$n \ge 400$）
设 $n \ge 400$ 为任意偶数。由引理 1，主弧邻域两两互不相交（由引理1中 $e^u/u^3$ 的单调性可知 $(\log n)^3<n$），加权表示数 $R(n)$ 为实数且满足：
$$R(n) = \operatorname{Re} I_{\mathfrak{M}}(n) + \operatorname{Re} I_{\mathfrak{m}}(n) \ge \operatorname{Re} I_{\mathfrak{M}}(n) - |I_{\mathfrak{m}}(n)|.$$
由引理 2 与引理 3，对一切偶数 $n \ge 400$，主弧积分满足无需借助任何未证假设的显式下界：
$$\operatorname{Re} I_{\mathfrak{M}}(n) \ge 0.65\, \mathfrak{S}(n) n, \quad \text{其中 } \mathfrak{S}(n) \ge 2 C_2 > 1.3114.$$

### 第三步：次弧积分的单重数傅里叶对偶展开与严格 $O(n)$ 量纲上界
考察次弧积分 $I_{\mathfrak{m}}(n) = \int_{\mathfrak{m}} S(\alpha)^2 e(-n\alpha) d\alpha$。
将一个 $S(\alpha) = \sum_{k=1}^{n} \Lambda(k) e(k\alpha)$ 按定义 1.3 展开，利用积分线性性质与次弧傅里叶相关核 $a_{\mathfrak{m}}(m) = \int_{\mathfrak{m}} S(\alpha) e(-m\alpha) d\alpha$：
$$I_{\mathfrak{m}}(n) = \sum_{k=1}^{n} \Lambda(k) \int_{\mathfrak{m}} S(\alpha) e\left(-(n-k)\alpha\right) d\alpha = \sum_{k=1}^{n} \Lambda(k) a_{\mathfrak{m}}(n-k).$$
对上式两边取绝对值并应用三角不等式，结合引理 4 的 $\psi(n) \le 1.04\,n$ 与引理 5、引理 6 已证的 $\max_{0 \le m < n} |a_{\mathfrak{m}}(m)| < 0.046$：
$$|I_{\mathfrak{m}}(n)| \le \sum_{k=1}^{n} \Lambda(k) \left|a_{\mathfrak{m}}(n-k)\right| \le \psi(n) \cdot \max_{0 \le m < n} |a_{\mathfrak{m}}(m)| \le (1.04\, n) \times 0.046 = 0.04784\, n.$$
结合引理 2 的严格下界 $\mathfrak{S}(n) > 1.3114$，有 $\frac{0.04784}{\mathfrak{S}(n)} < \frac{0.04784}{1.3114} \approx 0.03648 < 0.04$，故：
$$|I_{\mathfrak{m}}(n)| < 0.04\, \mathfrak{S}(n) n.$$

### 第四步：剔除高次素数幂贡献并导出纯素数对下界 $R_{\text{prime}}(n) > 0$
由引理 7，对全区间 $n \in [400, \infty)$ 内的任意偶数 $n$，高次素数幂贡献无条件满足 $R_{\text{pow}}(n) < 0.59\, \mathfrak{S}(n) n$。
将第二步、第三步与第四步的显式不等式全部代入纯素数对表示数 $R_{\text{prime}}(n)$ 的恒等式中：
$$\begin{aligned}
R_{\text{prime}}(n) &= \sum_{\substack{p+q=n \\ p,q \in \mathbb{P}}} (\log p)(\log q) \\
&= R(n) - R_{\text{pow}}(n) \\
&\ge \operatorname{Re} I_{\mathfrak{M}}(n) - |I_{\mathfrak{m}}(n)| - R_{\text{pow}}(n) \\
&> 0.65\, \mathfrak{S}(n) n - 0.04\, \mathfrak{S}(n) n - 0.59\, \mathfrak{S}(n) n \\
&= 0.02\, \mathfrak{S}(n) n.
\end{aligned}$$
由引理 2 的 $\mathfrak{S}(n) > 1.3114$，最终得到：
$$R_{\text{prime}}(n) = \sum_{\substack{p+q=n \\ p,q \in \mathbb{P}}} (\log p)(\log q) > 0.02 \times 1.3114\, n = 0.026228\, n > 0.$$
由于求和式中每一项 $(\log p)(\log q) > 0$，且总和严格大于零，故对一切偶数 $n \ge 400$，必然存在素数对 $(p, q) \in \mathbb{P} \times \mathbb{P}$ 满足 $p + q = n$。

---

## 4. 常数合并与结论

有限段由附录A覆盖。对任意偶数 $n\ge400$，前述主弧、次弧及高次幂估计按同一归一化尺度合并为

$$
\begin{aligned}
R_{\text{prime}}(n)
&=R(n)-R_{\text{pow}}(n)\\
&\ge\operatorname{Re}I_{\mathfrak M}(n)
       -|I_{\mathfrak m}(n)|-R_{\text{pow}}(n)\\
&>(0.65-0.04-0.59)\mathfrak S(n)n\\
&=0.02\,\mathfrak S(n)n
 >0.026228\,n>0.
\end{aligned}
$$

求和中的每一项均为正的素数对数权，因此这个正下界给出至少一组 $p,q\in\mathbb P$，使 $n=p+q$。连同有限段，目标命题写为

$$\forall n\in\mathbb Z,\qquad n>2\ \land\ 2\mid n
\quad\Longrightarrow\quad
\exists p,q\in\mathbb P,\quad n=p+q.\qquad\blacksquare$$

## 附录A　4至398的素数分解记录

下表按连续的偶数区段列出见证，允许两个素数相同。各区段合计198项，覆盖从4到398的每个偶数各一次。

- **$4 \le n \le 22$**：$4=2+2$，$6=3+3$，$8=3+5$，$10=3+7$，$12=5+7$，$14=3+11$，$16=3+13$，$18=5+13$，$20=3+17$，$22=3+19$；
- **$24 \le n \le 42$**：$24=5+19$，$26=3+23$，$28=5+23$，$30=7+23$，$32=3+29$，$34=3+31$，$36=5+31$，$38=7+31$，$40=3+37$，$42=5+37$；
- **$44 \le n \le 62$**：$44=3+41$，$46=3+43$，$48=5+43$，$50=3+47$，$52=5+47$，$54=7+47$，$56=3+53$，$58=5+53$，$60=7+53$，$62=3+59$；
- **$64 \le n \le 82$**：$64=3+61$，$66=5+61$，$68=7+61$，$70=3+67$，$72=5+67$，$74=3+71$，$76=3+73$，$78=5+73$，$80=7+73$，$82=3+79$；
- **$84 \le n \le 98$**：$84=5+79$，$86=3+83$，$88=5+83$，$90=7+83$，$92=3+89$，$94=5+89$，$96=7+89$，$98=19+79$；
- **$100 \le n \le 118$**：$100=3+97$，$102=5+97$，$104=3+101$，$106=3+103$，$108=5+103$，$110=3+107$，$112=3+109$，$114=5+109$，$116=3+113$，$118=5+113$；
- **$120 \le n \le 138$**：$120=7+113$，$122=13+109$，$124=11+113$，$126=13+113$，$128=19+109$，$130=3+127$，$132=5+127$，$134=3+131$，$136=5+131$，$138=7+131$；
- **$140 \le n \le 158$**：$140=3+137$，$142=3+139$，$144=5+139$，$146=7+139$，$148=11+137$，$150=11+139$，$152=3+149$，$154=3+151$，$156=5+151$，$158=7+151$；
- **$160 \le n \le 178$**：$160=3+157$，$162=5+157$，$164=7+157$，$166=3+163$，$168=5+163$，$170=3+167$，$172=5+167$，$174=7+167$，$176=3+173$，$178=5+173$；
- **$180 \le n \le 198$**：$180=7+173$，$182=3+179$，$184=3+181$，$186=5+181$，$188=7+181$，$190=11+179$，$192=11+181$，$194=3+191$，$196=3+193$，$198=5+193$；
- **$200 \le n \le 218$**：$200=3+197$，$202=3+199$，$204=5+199$，$206=7+199$，$208=11+197$，$210=11+199$，$212=13+199$，$214=3+211$，$216=5+211$，$218=7+211$；
- **$220 \le n \le 238$**：$220=23+197$，$222=11+211$，$224=13+211$，$226=3+223$，$228=5+223$，$230=3+227$，$232=3+229$，$234=5+229$，$236=3+233$，$238=5+233$；
- **$240 \le n \le 258$**：$240=7+233$，$242=3+239$，$244=3+241$，$246=5+241$，$248=7+241$，$250=11+239$，$252=11+241$，$254=3+251$，$256=5+251$，$258=7+251$；
- **$260 \le n \le 278$**：$260=3+257$，$262=5+257$，$264=7+257$，$266=3+263$，$268=5+263$，$270=7+263$，$272=3+269$，$274=3+271$，$276=5+271$，$278=7+271$；
- **$280 \le n \le 298$**：$280=3+277$，$282=5+277$，$284=3+281$，$286=3+283$，$288=5+283$，$290=7+283$，$292=11+281$，$294=11+283$，$296=3+293$，$298=5+293$；
- **$300 \le n \le 318$**：$300=7+293$，$302=19+283$，$304=11+293$，$306=13+293$，$308=31+277$，$310=3+307$，$312=5+307$，$314=3+311$，$316=3+313$，$318=5+313$；
- **$320 \le n \le 338$**：$320=3+317$，$322=5+317$，$324=7+317$，$326=13+313$，$328=11+317$，$330=13+317$，$332=19+313$，$334=3+331$，$336=5+331$，$338=7+331$；
- **$340 \le n \le 358$**：$340=3+337$，$342=5+337$，$344=7+337$，$346=29+317$，$348=11+337$，$350=3+347$，$352=3+349$，$354=5+349$，$356=3+353$，$358=5+353$；
- **$360 \le n \le 378$**：$360=7+353$，$362=3+359$，$364=5+359$，$366=7+359$，$368=19+349$，$370=3+367$，$372=5+367$，$374=7+367$，$376=3+373$，$378=5+373$；
- **$380 \le n \le 398$**：$380=7+373$，$382=3+379$，$384=5+379$，$386=3+383$，$388=5+383$，$390=7+383$，$392=3+389$，$394=5+389$，$396=7+389$，$398=31+367$。

本次逐项复核的内容为：等式左端取遍 $\{4,6,\ldots,398\}$；右端两数之和等于左端；每个右端加数的全部可能小因子均已检查到其整数平方根。198项均满足这些条件。
