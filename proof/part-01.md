## 1. 相关定义

### 定义 1.1（素数集、冯·曼戈尔特函数与对数归一化四阶 Riesz 平滑算子）
设正整数集为 $\mathbb{Z}^+$，素数集定义为 $\mathbb{P} = \{p \in \mathbb{Z}^+ : p > 1 \text{ 且除 } 1, p \text{ 外无其他正因子}\}$。
定义冯·曼戈尔特函数（von Mangoldt function）$\Lambda: \mathbb{Z}^+ \to \mathbb{R}_{\ge 0}$ 为：
$$\Lambda(k) = \begin{cases} \log p, & \text{若 } k = p^m \ (p \in \mathbb{P},\ m \in \mathbb{Z}^+,\ m \ge 1), \\ 0, & \text{其他情形}. \end{cases}$$
定义切比雪夫函数 $\psi(x) = \sum_{k \le x} \Lambda(k)$，第一切比雪夫函数 $\vartheta(x) = \sum_{p \le x} \log p$，模 $q$ 算术级数切比雪夫函数 $\psi(x; q, r) = \sum_{\substack{k \le x \\ k \equiv r \pmod q}} \Lambda(k)$，以及模 $q$ 狄利克雷特征 $\chi$ 下的扭转切比雪夫函数 $\psi(x, \chi) = \sum_{k \le x} \Lambda(k)\chi(k)$。
对任意实数 $t \ge 400$ 及基准尺度 $n \ge 400$，定义带对数归一化因子的四阶 Cesàro–Riesz 平滑切比雪夫算子（其中 $4! = 24$）为：
$$\psi_4(t, \chi) = \frac{1}{4!\,(\log n)^5} \sum_{k \le t} \left(1 - \frac{k}{t}\right)^4 \Lambda(k)\chi(k) = \frac{1}{2\pi i (\log n)^5} \int_{\kappa - i\infty}^{\kappa + i\infty} \left(-\frac{L'}{L}(s, \chi)\right) \frac{t^s}{s(s+1)(s+2)(s+3)(s+4)} ds.$$

### 定义 1.2（加法特征、指数生成函数与次弧傅里叶相关核）
对任意实数 $\alpha \in \mathbb{R}$，记标准加法特征为 $e(\alpha) = e^{2\pi i \alpha}$，记 $\alpha$ 到最近整数的距离为 $\|\alpha\| = \min_{z \in \mathbb{Z}} |\alpha - z|$。
对于给定的偶数 $n > 2$，定义素数加权指数生成函数 $S(\alpha)$ 为：
$$S(\alpha) = \sum_{k=1}^{n} \Lambda(k) e(k\alpha).$$
对于任意可测子集 $\mathcal{A} \subset \mathbb{R}/\mathbb{Z}$、整数频率 $h \in \mathbb{Z}$ 及整数移位 $m \in \mathbb{Z}$，定义集合 $\mathcal{A}$ 的特征傅里叶核 $W_{\mathcal{A}}(h)$ 与次弧相关核 $a_{\mathcal{A}}(m)$ 分别为：
$$W_{\mathcal{A}}(h) = \int_{\mathcal{A}} e(h\alpha) d\alpha, \qquad a_{\mathcal{A}}(m) = \int_{\mathcal{A}} S(\alpha) e(-m\alpha) d\alpha = \sum_{k=1}^{n} \Lambda(k) W_{\mathcal{A}}(k-m).$$

### 定义 1.3（哥德巴赫加权表示数与纯素数对表示数）
对于大于 $2$ 的偶数 $n$，定义冯·曼戈尔特加权表示数 $R(n)$、纯素数对加权表示数 $R_{\text{prime}}(n)$ 以及高次素数幂余项 $R_{\text{pow}}(n)$ 如下：
$$R(n) = \sum_{\substack{a+b=n \\ a,b \ge 1}} \Lambda(a)\Lambda(b), \qquad R_{\text{prime}}(n) = \sum_{\substack{p+q=n \\ p,q \in \mathbb{P}}} (\log p)(\log q),$$
$$R_{\text{pow}}(n) = R(n) - R_{\text{prime}}(n) = \sum_{\substack{a+b=n \\ a=p^k (k \ge 2) \text{ 或 } b=q^m (m \ge 2)}} \Lambda(a)\Lambda(b).$$
由于对一切素数 $p, q \in \mathbb{P}$ 有 $(\log p)(\log q) > 0$，若对偶数 $n > 2$ 有 $R_{\text{prime}}(n) > 0$，则必存在素数对 $(p,q) \in \mathbb{P} \times \mathbb{P}$ 使得 $n = p + q$。

### 定义 1.4（基于阈值 $N_0 = 400$ 的主弧、双区段次弧完备覆盖剖分与窄扇区自洽性验证）
设解析分界阈值为 $N_0 = 400$（当 $n \ge 400$ 时，$\sqrt{n} \ge 20$，$n^{1/3} \ge 7.368063$，$n^{1/4} \ge \sqrt{20} \approx 4.472136$，$n^{1/6} \ge 2.714418$，$\log n \ge \ln(400) \approx 5.9914645$，$(\log n)^2 \ge 35.89765$，$(\log n)^3 \ge 215.0795$，$(\log n)^4 \ge 1288.641$）。取主弧模数上限参数 $P$、对偶尺度参数 $Q$、单弧半宽 $\eta_q$ 以及与 $q$ 无关的统一外层相位阈值 $\eta_*$ 分别为：
$$P = \frac{1}{2}\log n, \qquad Q = \frac{n}{(\log n)^2}, \qquad \eta_q = \frac{1}{qQ} = \frac{(\log n)^2}{qn}, \qquad \eta_* = \frac{\log n}{2n}.$$
对任意满足 $1 \le q \le P$ 且 $\gcd(a,q)=1$ 的整数对 $(q,a)$（其中 $1 \le a \le q$），定义单个主弧邻域为：
$$\mathfrak{M}(q,a) = \left\{ \alpha \in \left[\frac{1}{Q}, 1+\frac{1}{Q}\right] : \left|\alpha - \frac{a}{q}\right| \le \eta_q = \frac{(\log n)^2}{qn} \right\}.$$
定义主弧集合 $\mathfrak{M}$ 与次弧集合 $\mathfrak{m}$ 分别为：
$$\mathfrak{M} = \bigcup_{1 \le q \le P} \bigcup_{\substack{1 \le a \le q \\ \gcd(a,q)=1}} \mathfrak{M}(q,a), \qquad \mathfrak{m} = \left[\frac{1}{Q}, 1+\frac{1}{Q}\right] \setminus \mathfrak{M}.$$

**【次弧 $\mathfrak{m} = \mathfrak{m}_1 \sqcup \mathfrak{m}_2$ 的集合论无遗漏完备覆盖剖分】**：
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

**【窄扇区半宽 $\eta_{q,\nu}^{(2)}$ 的几何自洽性验证】**：
对全部三个参数 $H_\nu \ge H_{\text{I}} = n(\log n)^7$ 及一切 $n \ge 400, q > n^{1/2}$：
1. **严格包含关系 $\eta_{q,\nu}^{(2)} < \eta_q$**：因 $q > n^{1/2}$ 且 $\sqrt{H_\nu} \ge n^{1/2}(\log n)^{7/2}$，分母满足 $q^2 \sqrt{H_\nu} > q n (\log n)^{7/2} > q n$，故 $\eta_{q,\nu}^{(2)} = \frac{(\log n)^2}{q^2 \sqrt{H_\nu}} < \frac{(\log n)^2}{q n} = \eta_q$。
2. **各窄扇区 $\mathcal{I}_{q,a}^{(\nu)}$ 两两互不相交**：设 $\frac{a}{q} \ne \frac{a'}{q'}$ 为 $\mathfrak{m}_2$ 中任意两个不同既约有理点（$q, q' > n^{1/2}$）。在 $n \ge 400$ 上，因 $\frac{(\log n)^2}{\sqrt{H_\nu}} \le \frac{1}{n^{1/2}(\log n)^{3/2}} \le \frac{1}{20 \times (5.99146)^{3/2}} < \frac{1}{293}$，故：
   $$\eta_{q,\nu}^{(2)} + \eta_{q',\nu}^{(2)} = \frac{(\log n)^2}{\sqrt{H_\nu}}\left(\frac{1}{q^2} + \frac{1}{q'^2}\right) < \frac{1}{293}\left(\frac{1}{q^2} + \frac{1}{q'^2}\right) < \frac{1}{qq'} \le \left|\frac{a}{q} - \frac{a'}{q'}\right|.$$

---

## 2. 使用的引理

### 引理 1（加法特征正交性、$n \ge 400$ 主弧无交性严格验算与圆法积分分解）
对任意整数 $m \in \mathbb{Z}$ 及单位长度周期区间 $I = [\eta, \eta+1]$，有正交关系 $\int_I e(m\alpha)d\alpha = \delta_{m,0}$。当 $n \ge N_0 = 400$ 时，定义 1.4 中的各主弧邻域 $\mathfrak{M}(q,a)$ 两两严格互不相交，主弧总勒贝格测度满足 $\mu(\mathfrak{M}) \le \frac{(\log n)^3}{n}$，且 $R(n)$ 精确分解为：
$$R(n) = I_{\mathfrak{M}}(n) + I_{\mathfrak{m}}(n), \quad \text{其中 } I_{\mathfrak{M}}(n) = \int_{\mathfrak{M}} S(\alpha)^2 e(-n\alpha)d\alpha,\quad I_{\mathfrak{m}}(n) = \int_{\mathfrak{m}} S(\alpha)^2 e(-n\alpha)d\alpha.$$

**【证明】**
1. **正交性与全区间恒等式**：当 $m = 0$ 时 $\int_\eta^{\eta+1} 1\,d\alpha = 1$；当 $m \in \mathbb{Z}\setminus\{0\}$ 时，$\int_\eta^{\eta+1} e^{2\pi i m \alpha} d\alpha = \frac{e^{2\pi i m \eta}(e^{2\pi i m} - 1)}{2\pi i m} = 0$。将 $S(\alpha)^2 = \sum_{a=1}^n \sum_{b=1}^n \Lambda(a)\Lambda(b)e((a+b)\alpha)$ 代入区间 $I = [1/Q, 1+1/Q]$ 的积分，由正交关系得：
   $$\int_{1/Q}^{1+1/Q} S(\alpha)^2 e(-n\alpha) d\alpha = \sum_{a=1}^n \sum_{b=1}^n \Lambda(a)\Lambda(b) \delta_{a+b, n} = \sum_{a+b=n} \Lambda(a)\Lambda(b) = R(n).$$
2. **主弧互不相交性在 $n \ge 400$ 处的精确算术验算**：
   设 $\frac{a}{q} \ne \frac{a'}{q'}$ 为两个不同既约分数，且 $1 \le q, q' \le P = \frac{1}{2}\log n$。
   考察函数 $h(u) = e^u - u^3$，其中 $u = \log n \ge \ln(400) \approx 5.9914645$。
   在左端点 $u_0 = \ln(400)$ 处，$e^{u_0} = 400$，而 $u_0^3 \approx (5.9914645)^3 = 215.0795 < 400 = e^{u_0}$。
   对 $u \ge \ln(400)$ 求导得 $h'(u) = e^u - 3u^2 \ge 400 - 3 \times 35.8977 = 292.307 > 0$，故 $h(u)$ 在 $[\ln 400, \infty)$ 上严格单调递增，即对一切 $n \ge 400$ 恒有 $(\log n)^3 < n$，从而 $\frac{1}{2P} = \frac{1}{\log n} > \frac{(\log n)^2}{n} = \frac{1}{Q}$。
   因此，任意两个不同有理中心的距离满足：
   $$\left|\frac{a}{q} - \frac{a'}{q'}\right| \ge \frac{1}{qq'} = \frac{1}{2qq'} + \frac{1}{2qq'} \ge \frac{1}{2qP} + \frac{1}{2q'P} > \frac{1}{qQ} + \frac{1}{q'Q} = \eta_q + \eta_{q'}.$$
   故各闭区间 $\mathfrak{M}(q,a)$ 两两严格无交。其总测度为：
   $$\mu(\mathfrak{M}) = \sum_{q=1}^{\lfloor P \rfloor} \sum_{\substack{1 \le a \le q \\ \gcd(a,q)=1}} 2\eta_q = \sum_{q=1}^{\lfloor P \rfloor} \varphi(q) \frac{2(\log n)^2}{qn} \le 2 P \frac{(\log n)^2}{n} = \frac{(\log n)^3}{n}.$$
   由测度可加性得 $R(n) = I_{\mathfrak{M}}(n) + I_{\mathfrak{m}}(n)$。证毕。 $\blacksquare$

---

### 引理 2（奇异级数的欧拉乘积展开、周期零均值部分和界与精确算术上下界）
定义拉马努金和 $c_q(n) = \sum_{\substack{1 \le a \le q \\ \gcd(a,q)=1}} e\left(-\frac{an}{q}\right)$，以及哥德巴赫奇异级数 $\mathfrak{S}(n) = \sum_{q=1}^{\infty} \frac{\mu(q)^2}{\varphi(q)^2} c_q(n)$。则：
1. **周期零均值正交性与部分和上界**：对任意整数 $q \ge 2$ 及任意整数起点 $h_0 \in \mathbb{Z}$，$c_q(h)$ 在长度为 $q$ 的完整周期上的和恒等于 $0$，且任意长度 $L \ge 1$ 的部分和满足一致上界：
   $$\sum_{r=1}^{q} c_q(h_0 + r) = 0, \qquad \left|\sum_{r=1}^{L} c_q(h_0 + r)\right| \le \frac{1}{2} q \varphi(q) \quad (\forall L \ge 1).$$
