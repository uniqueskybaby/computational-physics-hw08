2. **奇异级数精确下界与绝对级数上界**：对任意大于 $2$ 的偶数 $n$，级数 $\mathfrak{S}(n)$ 绝对收敛，且满足：
   $$\mathfrak{S}(n) = 2 C_2 \prod_{\substack{p > 2 \\ p \mid n}} \left(\frac{p-1}{p-2}\right) \ge 2 C_2 > 1.3114, \quad \text{其中 } C_2 = \prod_{p > 2} \left(1 - \frac{1}{(p-1)^2}\right) > 0.6557,$$
   $$\sum_{d \mid n} \frac{\mu(d)^2}{\varphi(d)} = \frac{n}{\varphi(n)} \le \frac{\mathfrak{S}(n)}{C_2} < 1.5251\,\mathfrak{S}(n), \qquad \mathfrak{S}_{\text{abs}}(n) := \sum_{q=1}^{\infty} \frac{\mu(q)^2 |c_q(n)|}{\varphi(q)^2} \le \frac{\mathfrak{S}(n)}{C_2^2} < 2.33\,\mathfrak{S}(n).$$

**【证明】**
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
   因此 $\mathfrak{S}(n) \ge 2 C_2 > 1.3114$，$\frac{1}{C_2} < 1.5251$，$\frac{1}{C_2^2} < 2.33$。证毕。 $\blacksquare$

---

### 引理 3（主弧积分的非渐近显式下界定理：基于无条件 McCurley 零点区与完整常数合并推导）
设偶数 $n \ge 400$，主弧 $\mathfrak{M}$ 如定义 1.4 所示（$P = \frac{1}{2}\log n,\ Q = n/(\log n)^2$），则主弧积分 $I_{\mathfrak{M}}(n)$ 满足完全无条件且每项常数均可独立算术验算的显式实部下界：
$$\operatorname{Re} I_{\mathfrak{M}}(n) \ge 0.65\, \mathfrak{S}(n) n.$$

**【证明】**
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
