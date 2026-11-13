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
