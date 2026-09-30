+++
title = "高效的机器翻译 Reranking 方法"
date = 2026-09-30
draft = false
summary = "机器翻译 Reranking 加速方法综述。"
tags = ["Machine Translation", "LLM Decoding", "Reranking"]
+++

## LLM Decoding 和 Reranking

LLM Decoding 可以分成两种 [[1]](#ref-1)：

1. 设计 $f$，作用于每一步的 token 分布 $P(y_i \mid X, y_{<i})$，根据 $f\big(P(y_i \mid X, y_{<i})\big)$ 逐个选出 token。Greedy、Beam Search、Top-$k$/Top-$p$ 采样都属于这一类。
2. 设计 $g$，作用于完整输出 $y$ 上，例如 $g\big(P(y \mid X)\big)$，其中 $P(y \mid X) = \prod_i P(y_i \mid X, y_{<i})$。这个过程分成两步：先让 LLM 生成很多候选，然后用 $g$ 选择出最好的一个。所以这种解码也常被称为 Reranking。

## 机器翻译中的应用

机器翻译任务中，我们定义源句：$s$，模型生成翻译：$y$，参考翻译：$y^{\mathrm{ref}}$，翻译质量评估：$q(y, y^{\mathrm{ref}})$。

在机器翻译中做 Reranking 有一个很好的直觉：多生成一些候选，相比于只生成一个翻译，总能挑出更好的。
这里我们记候选集为 $\mathcal{C}$。
在机器翻译中，最常被应用的是 MBR Decoding 和 QE Reranking。

在 WMT General Machine Translation Shared Task 中，这两种 Reranking 方法多次被优胜方案采用 [[2]](#ref-2)。

### MBR Decoding

这个方法采用了贝叶斯决策：选择与其他候选“平均最一致”的那个候选。衡量一致程度的函数称为 **utility**，记为 $u(h, h_s)$，表示把 $h_s$ 当作参考时 $h$ 的好坏。选择规则如下：

$$
\hat{y}_{\mathrm{MBR}}
= \mathop{\mathrm{argmax}}\limits_{h\in\mathcal{C}}
\frac{1}{|\mathcal{C}|}
\sum_{h_s\in\mathcal{C}} u(h,h_s).
$$

这里真实的参考翻译 $y^{\mathrm{ref}}$ 在解码时是拿不到的，所以把其他候选 $h_s$ 当作伪参考（pseudo-reference）。因此 utility 通常直接使用 reference-based metric（例如 COMET、BLEURT、chrF++），即把前面的 $q(y, y^{\mathrm{ref}})$ 中的参考翻译换成伪参考：$u(h, h_s) = q(h, h_s)$。上式的求和平均就是 $h$ 的期望 utility，选出期望 utility 最高的候选，就能迎合该 metric 的偏好，选择出更好的输出。

下文中凡是涉及 MBR 的地方，“utility” 都指这个 reference-based metric，“utility 调用”指调用一次该 metric。

### QE Reranking

这个方法直接使用 reference-free 的模型（例如 CometKiwi）进行打分，然后选择：

$$
\hat{y}_{\mathrm{QE}}
= \mathop{\mathrm{argmax}}\limits_{h\in\mathcal{C}}
r_{\mathrm{QE}}(h;s).
$$

### 计算量分析

从输入到输出，以端到端的视角来看，对于每一个 source，需要的计算是：

- 候选生成：$|\mathcal{C}|$ 次 LLM 生成
- 如果使用 MBR Decoding：$|\mathcal{C}|^2$ 次 utility（reference-based metric）调用
- 如果采用 QE Reranking：$|\mathcal{C}|$ 次 QE 模型（reference-free metric）调用

所以加速可以从两个方向入手：减少候选生成的次数，或者减少 utility / QE 模型的调用次数。

## 加速方法

| 方法 | 适用 | 减少的计算 | 核心思路 |
| --- | --- | --- | --- |
| Quit [[3]](#ref-3) | MBR / QE | 候选生成 + 打分 | 增量生成，最优分数稳定后提前停止 |
| PMBR [[4]](#ref-4) | MBR | utility 调用 | 只算部分分数矩阵，用低秩矩阵补全恢复其余 |
| PruneMBR [[5]](#ref-5) | MBR | utility 调用 | 把“选最高分”变成“估计每个候选最终胜出的概率”，用 bootstrap 估计并剪枝 |
| CBMBR [[6]](#ref-6) | MBR | utility 调用 | 对伪参考聚类，用簇中心代替全部伪参考；等权的簇中心还能缓解多峰偏差 |
| BayesOpt [[7]](#ref-7) | QE | QE 模型调用 | 用高斯过程建模分数，只挑值得打分的候选 |

### 综合方法：Quit

现有的加速方法大多只针对 MBR 和 QE 的打分部分，而候选生成本身基本没有被处理。然而在当下基于 LLM 的翻译系统中，候选生成往往比 Reranking 更贵。Quit（Quantifying Uncertainty for Incremental Termination）[[3]](#ref-3) 同时减少了候选生成和 Reranking 的计算，把整个“生成—重排”流程看作一系列决策问题。

**增量生成。** 每一步生成 $k$ 个新候选（论文中 $k=8$），第 $i$ 步后候选集为 $\mathcal{C}_i$，$|\mathcal{C}_i| = ik$（对 MBR，伪参考集也随之增长）。记第 $i$ 步 reranker 选出的输出为 $\hat{y}_i$，其分数（即当前最优分数）为 $R_i(s)$；用满预算 $N_{\max}$ 时共 $N_{\max}/k$ 步。

**理想目标：提前停止的风险。** 在第 $i$ 步停下，损失的是它与满预算输出之间的真实质量差：

$$
\mathcal{L}_i(s) = \big|\, q(\hat{y}_i, y^{\mathrm{ref}}) - q(\hat{y}_{N_{\max}/k}, y^{\mathrm{ref}}) \,\big|.
$$

$\mathcal{L}_i(s)$ 越小，说明提前停止越不损失质量。但推理时既没有参考翻译，也没有满预算的输出，无法计算。

**近似一：用 reranker 分数代替真实质量。** 去掉对参考翻译的依赖：

$$
\mathcal{L}^{\mathrm{proxy}}_i(s) = \big|\, R_i(s) - R_{N_{\max}/k}(s) \,\big|.
$$

但 $R_{N_{\max}/k}(s)$ 要等到用满预算才知道，仍然无法计算。

**近似二：用最近窗口内的波动代替与终点的差距。** 把步骤划分成大小为 $w$ 的不重叠窗口 $\mathcal{B}_b$，只看窗口内最优分数的变化范围：

$$
\Delta_b R(s) = \max_{j\in\mathcal{B}_b} R_j(s) - \min_{j\in\mathcal{B}_b} R_j(s).
$$

$\Delta_b R(s)$ 大说明最优分数最近还在明显提升；小说明已经局部稳定，继续生成的收益有限。论文用 PRR 验证了它与真实风险 $\mathcal{L}_i(s)$ 相关。

**停止规则。** 当波动低于阈值 $\alpha$（论文中 $10^{-3}$）时停止：

$$
T = \min\{\, b\cdot w : \Delta_b R(s) \le \alpha \,\}.
$$

在 3 个 NMT 模型、19 个语言对上，Quit 对 MBR 带来 1.47–2.66 倍、对 QE Reranking 带来 3.43–4.12 倍的端到端加速，同时能够保持翻译质量。

### MBR Decoding：减少 utility 调用

MBR 的瓶颈在于要计算一个 $|\mathcal{C}| \times |\mathcal{C}|$ 的分数矩阵 $M$，其中 $M_{ij} = u(h_i, h_j)$。下面三种方法分别从“少算矩阵元素”“少评估候选”“少用伪参考”三个角度来减少计算。

#### PMBR

PMBR [[4]](#ref-4) 观察并实验验证分数矩阵 $M$ 具有低秩结构，因此可以把 MBR 看成一个矩阵补全问题：

1. 随机只计算 $M$ 中的一部分元素；
2. 用 Alternating Least Squares（ALS）把 $M$ 分解为 $M \approx UV^\top$，补全缺失的元素；
3. 在补全后的矩阵上按行求平均，选出 MBR 输出。

在 WMT22（en↔de、en↔ru）上，只需要原始 MBR $1/16$ 的 utility 计算，就能达到相同的 COMET22 分数，并优于其他近似基线。

#### PruneMBR

PruneMBR [[5]](#ref-5)（EMNLP 2023 Best Short Paper）的核心思想：**不必精确算出每个候选的期望 utility，只要知道它还有没有可能胜出。**

记候选集为 $H$，期望 utility 为 $U(y, \mathcal{Y}) = \mathbb{E}_{\hat{y}\sim\mathcal{Y}}[u(y,\hat{y})]$。MBR 的真正目标是 $\mathop{\mathrm{argmax}}_{y\in H} U\big(y, p_\theta(\cdot\mid x)\big)$，标准 MBR 用伪参考 $R$ 近似 $p_\theta$，需要 $|H||R|$ 次 utility 调用。PruneMBR 把“选最高分”改为估计每个候选成为真正赢家的概率：

$$
p\Big(\bigwedge_{\bar{y}\in H} U\big(y, p_\theta(\cdot\mid x)\big) \ge U\big(\bar{y}, p_\theta(\cdot\mid x)\big)\Big).
$$

**近似一：bootstrap。** $p_\theta$ 未知，就对当前伪参考 $R_t$ 有放回重采样，用“在 bootstrap 样本上胜出”来近似：

$$
\mathbb{E}_{\hat{R}_t\sim\mathrm{boot}(R_t)}\,
\mathbb{1}\Big(\bigwedge_{\bar{y}\in H_t} U(y,\hat{R}_t) \ge U(\bar{y},\hat{R}_t)\Big).
$$

**近似二：只和当前赢家比。** $H_t$ 很大时上式方差很大，于是只和当前赢家 $\bar{\bar{y}} = \mathop{\mathrm{argmax}}_{\bar{y}\in H_t} U(\bar{y}, R_t)$ 比较：

$$
\mathbb{E}_{\hat{R}_t\sim\mathrm{boot}(R_t)}\,
\mathbb{1}\big(U(y,\hat{R}_t) \ge U(\bar{\bar{y}},\hat{R}_t)\big).
$$

“胜过所有人”必然“胜过 $\bar{\bar{y}}$”，所以这是上一式的上界，且与 $|H_t|$ 无关。若它小于 $1-\alpha$，就剪掉 $y$。bootstrap 只是重新组合已算好的 utility，不需要额外调用。

**算法。** 伪参考数量按 schedule 逐步翻倍（COMET 从 8 开始，到 256 为止）；每一步补充伪参考，只为剩下的候选计算新的 utility，再做 bootstrap 剪枝；只剩一个候选或到达最后一步时停止。初始伪参考太少会让 bootstrap 偏差大、误剪赢家，因此 $r_1$ 决定了速度与精度的权衡。

结果：与标准 MBR 在统计上无显著差异，utility 调用在 chrF++ 下减少至少 7 倍，在 COMET 下至少 15 倍。

#### CBMBR

CBMBR [[6]](#ref-6) 的思想很巧妙：**伪参考只是特征空间里的向量，相近的伪参考可以合并成一个代表点。**

COMET 先用编码器 $f$ 分别编码源句、候选和参考，再由输出层 $\phi$ 打分，即 $u(h, h_s) = \phi\big(f(s), f(h), f(h_s)\big)$。编码只要 $|\mathcal{C}|$ 次，瓶颈是 $|\mathcal{C}|^2$ 次 $\phi$ 的计算。它依赖两个假设：

1. **语义相近的句子，向量也相近。** 在 STS-B 上，COMET 句向量与人工相似度的 Pearson 相关达 73.6，甚至高于 LaBSE（72.7）。
2. **$\phi$ 近似线性。** 于是一个簇 $\mathcal{K}_i$ 内“先打分再平均”可以换成“先平均再打分”：

$$
\frac{1}{|\mathcal{K}_i|}\sum_{h_s\in\mathcal{K}_i}\phi\big(f(s), f(h), f(h_s)\big)
\approx
\phi\big(f(s), f(h), c_i\big),\quad c_i = \frac{1}{|\mathcal{K}_i|}\sum_{h_s\in\mathcal{K}_i} f(h_s).
$$

具体做法是用 kmeans++ 与 k-means 把伪参考向量聚成 $k$ 类（论文中 $k=64$），然后

$$
\hat{y}_{\mathrm{CBMBR}}
= \mathop{\mathrm{argmax}}\limits_{h\in\mathcal{C}}
\frac{1}{k}\sum_{i=1}^{k} \phi\big(f(s), f(h), c_i\big),
$$

复杂度从 $O(|\mathcal{C}|^2)$ 降到 $O(|\mathcal{C}|\,k)$。

**为什么反而能超过标准 MBR？** 上式对簇中心**等权**平均。标准 MBR 平等对待每个伪参考，候选呈多峰分布时（例如来自多个系统），样本多的译法会占主导。CBMBR 让每个簇只算一票，对这种偏差更鲁棒。对照实验也支持这一点：按簇大小加权、更接近标准 MBR 的 CBMBR$_{\mathrm{cnt}}$，得分和标准 MBR 几乎一样（86.6 vs 86.7），不如 CBMBR（87.0）。

**结果。** 在多样化候选上，与标准 MBR 的差距在 0.1 COMET 以内，期望 utility 计算加速 5.7 倍（端到端 1.8 倍）。在多系统候选上，最多比标准 MBR 高 0.5 COMET。局限是只适用于能单独编码每个句子的 metric。

### QE Reranking：减少 QE 模型调用

#### BayesOpt

QE Reranking 的计算只有 $|\mathcal{C}|$ 次调用，但 QE 模型（如 CometKiwi）本身也很贵。Cheng 等人 [[7]](#ref-7) 把 Reranking 看成一个贝叶斯优化（Bayesian Optimization）问题：目标是找到 $\mathop{\mathrm{argmax}}_{h\in\mathcal{C}} r_{\mathrm{QE}}(h;s)$，但不必给每个候选都打分。

1. 用高斯过程（Gaussian Process）对未打分候选的分数建模，核函数是候选句向量（生成时顺带得到的 mean-pooled 表示，几乎零额外开销）上的 RBF 核；
2. 每一步用 acquisition function（Expected Improvement）在 exploration（和已打分候选差异大的）与 exploitation（和高分候选相似的）之间权衡，选出下一个要打分的候选；
3. 打分预算用完后，返回已打分候选中分数最高的那个。

在 200 个候选的设置下，平均只需约 70 次 CometKiwi 调用，就能达到随机选 180 个候选打分的效果。论文还提出了 multi-fidelity 版本：先用更便宜的小 QE 模型打分作为代理，进一步提升效率。

## 参考

1. <span id="ref-1"></span>*CMU Advanced NLP Spring 2025 (7): Decoding Algorithms*. [YouTube](https://www.youtube.com/watch?v=cN8yX_ZZWJw&list=PLqC25OT8ZpD3WxQ0FwWMGPS_BcWdcKyZy&index=8)
2. <span id="ref-2"></span>Tom Kocmi, Ekaterina Artemova, Eleftherios Avramidis, et al. *Findings of the WMT25 General Machine Translation Shared Task: Time to Stop Evaluating on Easy Test Sets*. WMT 2025. [ACL Anthology](https://aclanthology.org/2025.wmt-1.22/)
3. <span id="ref-3"></span>Guangyu Chen, Boxuan Lyu, Hidetaka Kamigaito, Kotaro Funakoshi, Manabu Okumura. *Quit While You're Ahead: Quit for Efficient Candidate Generation in Machine Translation Reranking*. arXiv:2609.00588, 2026. [arXiv](https://arxiv.org/abs/2609.00588)
4. <span id="ref-4"></span>Firas Trabelsi, David Vilar, Mara Finkelstein, Markus Freitag. *Efficient Minimum Bayes Risk Decoding using Low-Rank Matrix Completion Algorithms*. arXiv:2406.02832, 2024. [arXiv](https://arxiv.org/abs/2406.02832)
5. <span id="ref-5"></span>Julius Cheng, Andreas Vlachos. *Faster Minimum Bayes Risk Decoding with Confidence-based Pruning*. EMNLP 2023. [ACL Anthology](https://aclanthology.org/2023.emnlp-main.767/)
6. <span id="ref-6"></span>Hiroyuki Deguchi, Yusuke Sakai, Hidetaka Kamigaito, Taro Watanabe, Hideki Tanaka, Masao Utiyama. *Centroid-Based Efficient Minimum Bayes Risk Decoding*. Findings of ACL 2024. [ACL Anthology](https://aclanthology.org/2024.findings-acl.654/)
7. <span id="ref-7"></span>Julius Cheng, Maike Züfle, Vilém Zouhar, Andreas Vlachos. *A Bayesian Optimization Approach to Machine Translation Reranking*. NAACL 2025. [ACL Anthology](https://aclanthology.org/2025.naacl-long.145/)
