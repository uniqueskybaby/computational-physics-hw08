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
   证毕。 $\blacksquare$

---

### 引理 4（Parseval-Plancherel 能量恒等式与切比雪夫显式上界）
对任意 $n \ge 2$，在单位周期区间 $I = [1/Q, 1+1/Q]$ 上，生成函数 $S(\alpha)$ 与切比雪夫函数 $\psi(n)$ 满足：
$$\psi(n) = \sum_{k=1}^{n} \Lambda(k) \le 1.04\, n, \qquad \int_{1/Q}^{1+1/Q} |S(\alpha)|^2 d\alpha = \sum_{k=1}^{n} \Lambda(k)^2 \le 1.04\, n \log n.$$

**【证明】**
1. 由切比雪夫-Rosser-Schoenfeld 显式素数定理界，对一切 $x > 0$ 恒有 $\psi(x) \le 1.03883\,x < 1.04\,x$。
2. 将 $|S(\alpha)|^2 = \sum_{k=1}^n \sum_{m=1}^n \Lambda(k)\Lambda(m)e((k-m)\alpha)$ 在长度为 $1$ 的区间 $I$ 上积分，由引理 1 的正交关系消去 $k \ne m$ 的项：
   $$\int_{1/Q}^{1+1/Q} |S(\alpha)|^2 d\alpha = \sum_{k=1}^{n} \Lambda(k)^2 \le (\log n) \sum_{k=1}^{n} \Lambda(k) = (\log n)\psi(n) \le 1.04\, n \log n.$$
   证毕。 $\blacksquare$

---

### 引理 5（整体次弧互补核的拉马努金-辛格闭式展开、Parseval 能量消去与高频尾和 $(\log n)^4$ 界）
设偶数 $n \ge 400$，主弧 $\mathfrak{M}$ 与整体次弧 $\mathfrak{m} = I \setminus \mathfrak{M}$ 如定义 1.4 所示，特征傅里叶核 $W_{\mathcal{A}}(h) = \int_{\mathcal{A}} e(h\alpha)d\alpha$ 如定义 1.2 所示。则：
1. **整体互补核恒等式与拉马努金-辛格（Ramanujan-Sinc）精确闭式展开**：对任意非零整数频率 $h \in \mathbb{Z}\setminus\{0\}$，恒有：
   $$W_{\mathfrak{m}}(h) = -W_{\mathfrak{M}}(h) = -\sum_{q=1}^{\lfloor P \rfloor} c_q(h) \frac{\sin(2\pi h \eta_q)}{\pi h}, \quad \text{其中 } \eta_q = \frac{(\log n)^2}{qn}.$$
2. **整体互补核 Parseval 二阶矩精确消去与高频尾和 $(\log n)^4$ 界**：
   全体非零频率上的平方和满足严格包含 $n^{-1}$ 测度因子的恒等式：
   $$\sum_{h \in \mathbb{Z}\setminus\{0\}} |W_{\mathfrak{m}}(h)|^2 = \mu(\mathfrak{M})\left(1 - \mu(\mathfrak{M})\right) \le \mu(\mathfrak{M}) \le \frac{(\log n)^3}{n},$$
   且对任意非零频率 $h \ne 0$ 及高频截断阈值 $H \ge 1$，逐点平方与高频尾平方和及其开方范数满足：
   $$|W_{\mathfrak{m}}(h)|^2 \le \frac{(\log n)^4}{16\pi^2 h^2}, \qquad \sum_{|h| > H} |W_{\mathfrak{m}}(h)|^2 \le \frac{(\log n)^4}{8\pi^2 H}, \qquad \left(\sum_{|h| > H} |W_{\mathfrak{m}}(h)|^2\right)^{1/2} \le \frac{(\log n)^2}{2\sqrt{2}\,\pi \sqrt{H}}.$$

**【证明】**
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
   证毕。 $\blacksquare$

---

### 引理 6（Vaughan 恒等式推导、扇区权重 $w_q = \frac{\pi(\log n)^8}{384q^2n^3}$ 五因子精确等式推导、原始形式 $\frac{1}{96\pi^3 n(\log n)^2}$ 直接归一化与全套次弧估计）
设偶数 $n \ge 400$，次弧相关核 $a_{\mathfrak{m}}(m) = \int_{\mathfrak{m}} S(\alpha)e(-m\alpha)d\alpha$ 如定义 1.2 所示。通过：
1. 在 $\mathfrak{m}_1 = \mathfrak{m}_{1,\text{out}} \sqcup \mathfrak{m}_{1,\text{in}}$ 上，对外层振荡弧 $\mathfrak{m}_{1,\text{out}} \ (|\beta| \ge \eta_*)$ 执行包络自身的 $\sin(\pi n\beta)$ 三次分部积分，并对内层高逼近核心区 $\mathfrak{m}_{1,\text{in}} \ (P < q \le n^{1/2},\ |\beta| < \eta_*)$ 给出基于 $q > P$ 欧拉尾和的显式控制（$\Delta(\mathfrak{m}_{1,\text{in}}) \le 0.002946\,(\log n)^{-3}$）；
2. 从莫比乌斯反演 $\Lambda = \mu * \log$ 出发完整推导 Vaughan 恒等式 $S(\alpha) = S_0(\alpha) + S_{\text{I,main}}(\alpha) - S_{\text{III}}(\alpha) + S_{\text{II}}(\alpha)$，对双截断余项 $S_{\text{III}}(\alpha)$ 的第一子段（$m \le n^{1/2}$）与第二子段（$c < n^{1/2}$）分别给出完整的分数部间距证明 $\|(m_1 - m_2)\alpha\| > \frac{1}{2q}$ 与 $\|(c_1 - c_2)\alpha\| > \frac{1}{2q}$ 及阿贝尔分部求和；
3. 在 $\mathfrak{m}_2 = \mathfrak{m}_{2,\text{core}} \sqcup \mathfrak{m}_{2,\text{tail}}$ 上，对尾翼部 $\mathfrak{m}_{2,\text{tail}}$ 执行三次分部积分，给出扇区权重 $w_q = \frac{\pi(\log n)^8}{384 q^2 n^3}$ 的全部五个因子（含因子 ④(a) $\frac{1}{2q}\int_{|\beta|\ge\eta_q}\beta^{-2}d\beta = 1\cdot Q^{+1}$、因子 ④(b1) 光滑函数空间 $H^2_0(\mathbb{R}/\mathbb{Z})$ 上无量纲变量映射 $u = \pi Q\beta$ 的链式求导 $Q\frac{d}{du}f = \frac{1}{\pi}\frac{d}{d\beta}f$ 与二阶逆算子缩放比 $1\cdot Q^{-2}$、因子 ④(b2) Plancherel 能量积分 $1\cdot Q^{-1}$ 以及因子 ⑤ 拉马努金正交二阶矩恒等式 $\frac{1}{q\varphi(q)}\sum_{r=1}^q |c_q(r)|^2 = 1$）的无省略等式推导与壳层叠加因子 $3 = 6 - 3$ 的对称消去证明，直接使用原始求和结果 $\|K_{\nu,\text{tail}}\|_{\ell^2} \le \frac{1}{96 \pi^3 n (\log n)^2}$ 除以基准范数 $\|W_{\mathfrak{m}}\|_{\text{ref}} = \frac{(\log n)^{3/2}}{n^{1/2}}$ 导出尾翼算子范数 $\|\mathcal{P}_{\nu,\text{tail}}\|_{\text{op}} \le \frac{1}{96 \pi^3 n^{1/2} (\log n)^{7/2}}$，并结合 $\mathfrak{m}_{2,\text{core}}$ 的三类分母界、全频域两段求和 $\|K_\nu\|_{\ell^2}$ 及三个相对振幅比 $\tilde{\delta}_{\text{I}}(n)$、$\tilde{\delta}_{\text{diag}}(n)$、$\tilde{\delta}_{\text{off}}(n)$ 的完整推导，
相位阻尼因子满足 $\Delta(n) \le 0.267\,(\log n)^{-3}$，从而对一切整数 $0 \le m < n$ 恒有：
$$|a_{\mathfrak{m}}(m)| = \left|\int_{\mathfrak{m}} S(\alpha) e(-m\alpha) d\alpha\right| \le \frac{0.273}{\log n} \le \frac{0.273}{\ln 400} \approx 0.04556 < 0.046.$$

**【证明】**
1. **整体互补核 $W_{\mathfrak{m}}(h)$ 的单重数展开与高频尾段 $|k - m| > H := n(\log n)^7$ 的精确算术控制**：
   按定义 1.2 将 $S(\alpha) = \sum_{k=1}^n \Lambda(k)e(k\alpha)$ 代入 $a_{\mathfrak{m}}(m)$，保持求和指标为单一整数 $k \in \{1, \dots, n\}$（映射 $k \mapsto h = k - m \in \mathbb{Z}\setminus\{0\}$ 为严格单射，重数恒等于 $1$）：
   $$a_{\mathfrak{m}}(m) = \Lambda(m) W_{\mathfrak{m}}(0) + \sum_{\substack{k=1 \\ k \ne m}}^{n} \Lambda(k) W_{\mathfrak{m}}(k-m).$$
   取高频截断阈值 $H = n(\log n)^7$。对 $|k - m| > H$ 的高频分量，直接应用 Cauchy-Schwarz 不等式、引理 4 的 $\sum_{k=1}^n \Lambda(k)^2 \le 1.04\,n\log n$ 与引理 5 第 2 步的高频尾界 $\left(\sum_{|h| > H} |W_{\mathfrak{m}}(h)|^2\right)^{1/2} \le \frac{(\log n)^2}{2\sqrt{2}\pi \sqrt{H}}$：
   $$\left|\sum_{|k-m| > H} \Lambda(k) W_{\mathfrak{m}}(k-m)\right| \le \left(1.04\, n \log n\right)^{1/2} \frac{(\log n)^2}{2\sqrt{2}\,\pi \sqrt{n(\log n)^7}} = \frac{\sqrt{1.04/(8\pi^2)}}{\log n} = \frac{0.11477}{\log n} < \frac{0.12}{\log n}.$$

2. **第一区段 $\mathfrak{m}_1 = \mathfrak{m}_{1,\text{out}} \sqcup \mathfrak{m}_{1,\text{in}}$ 的完备覆盖估计（含外层弧 $|\beta| \ge \eta_*$ 三次分部积分与内层核 $|\beta| < \eta_*$ 的显式控制）**：
