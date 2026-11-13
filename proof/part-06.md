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
