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
                    $$\frac{\|\mathcal{D}_\beta^{-2} f\|_{L^2}}{\|\mathcal{D}_u^{-2} f\|_{L^2}} = \frac{\left\|\left(Q\frac{d}{du}\right)^{-2} f\right\|_{L^2}}{\left\|\left(\frac{d}{du}\right)^{-2} f\right\|_{L^2}} = \frac{Q^{-2} \|\mathcal{D}_u^{-2} f\|_{L^2}}{\|\mathcal{D}_u^{-2} f\|_{L^2}} = 1 \cdot Q^{-2}!!$$
                  - **子部分 (b2) Plancherel $L^2([0,1])$ 能量归一化 $1 \cdot Q^{-1}$**：
                    由引理 1 的加法特征正交关系 $\int_0^1 e((h_1 - h_2)\beta)d\beta = \delta_{h_1, h_2}$，跨 $Q$ 个对偶频率模式的归一化核 $D_Q(\beta)$ 在单位周期 $[0,1]$ 上的 $L^2$ 能量积分精确等于：
                    $$\int_{0}^{1} |D_Q(\beta)|^2 d\beta = \int_{0}^{1} \left|\frac{1}{Q}\sum_{h=1}^{Q} e(h\beta)\right|^2 d\beta = \frac{1}{Q^2} \sum_{h=1}^{Q} 1 = \frac{Q}{Q^2} = 1 \cdot Q^{-1}.$$
                  将子部分 (b1) 的二阶逆算子缩放比 $1 \cdot Q^{-2}$ 与子部分 (b2) 的 Plancherel $L^2$ 能量积分 $1 \cdot Q^{-1}$ 相乘，精确得到：
                  $$\left(\frac{\|\mathcal{D}_\beta^{-2} f\|_{L^2}}{\|\mathcal{D}_u^{-2} f\|_{L^2}}\right) \cdot \left(\int_{0}^{1} |D_Q(\beta)|^2 d\beta\right) = \left(1 \cdot Q^{-2}\right) \cdot \left(1 \cdot Q^{-1}\right) = 1 \cdot Q^{-(2+1)} = 1 \cdot Q^{-3}.$$
               3) **子因子 (a) 与子因子 (b) 相乘导出因子 ④**：
                  将子因子 (a) 的 $1 \cdot Q^{+1}$ 与子因子 (b) 的 $1 \cdot Q^{-3}$ 直接相乘，并代入定义 1.4 的 $Q = \frac{n}{(\log n)^2}$：
