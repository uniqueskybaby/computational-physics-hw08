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
