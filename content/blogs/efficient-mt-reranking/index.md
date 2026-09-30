+++
title = "Efficient Machine Translation Reranking"
date = 2026-09-30
draft = false
summary = "A summarization of efficient machine translation reranking methods."
tags = ["Machine Translation", "LLM Decoding", "Reranking"]
+++

## LLM Decoding and Reranking

LLM decoding methods fall into two families [[1]](#ref-1):

1. Design a function $f$ that acts on the per-step token distribution $P(y_i \mid X, y_{<i})$, and pick tokens one at a time according to $f\big(P(y_i \mid X, y_{<i})\big)$. Greedy decoding, beam search, and top-$k$/top-$p$ sampling all belong to this family.
2. Design a function $g$ that acts on the complete output $y$, e.g. $g\big(P(y \mid X)\big)$, where $P(y \mid X) = \prod_i P(y_i \mid X, y_{<i})$. This is a two-step process: first let the LLM generate many candidates, then use $g$ to select the best one. For this reason, this kind of decoding is often called reranking.

## Application to Machine Translation

In machine translation, we denote the source sentence by $s$, the model translation by $y$, the reference translation by $y^{\mathrm{ref}}$, and the translation quality metric by $q(y, y^{\mathrm{ref}})$.

Reranking in machine translation rests on a simple intuition: if we generate several candidates instead of just one, we can always pick a better one.
We denote the candidate set by $\mathcal{C}$.
The two most widely used reranking methods in machine translation are MBR decoding and QE reranking.

In the WMT General Machine Translation Shared Task, both reranking methods have repeatedly been adopted by winning systems [[2]](#ref-2).

### MBR Decoding

This method is based on Bayesian decision theory: it selects the candidate that "agrees most, on average" with the other candidates. The function that measures this agreement is called the **utility**, written $u(h, h_s)$, which scores how good $h$ is when $h_s$ is treated as the reference. The decision rule is:

$$
\hat{y}_{\mathrm{MBR}}
= \mathop{\mathrm{argmax}}\limits_{h\in\mathcal{C}}
\frac{1}{|\mathcal{C}|}
\sum_{h_s\in\mathcal{C}} u(h,h_s).
$$

The true reference $y^{\mathrm{ref}}$ is not available at decoding time, so the other candidates $h_s$ are used as pseudo-references. The utility is therefore usually a reference-based metric (e.g. COMET, BLEURT, chrF++): we simply replace the reference in $q(y, y^{\mathrm{ref}})$ with a pseudo-reference, $u(h, h_s) = q(h, h_s)$. The average in the formula is the expected utility of $h$; choosing the candidate with the highest expected utility caters to the preferences of the metric and yields a better output.

In the rest of this post, whenever MBR is involved, "utility" refers to this reference-based metric, and a "utility call" means one call to that metric.

### QE Reranking

This method directly scores each candidate with a reference-free model (e.g. CometKiwi) and selects:

$$
\hat{y}_{\mathrm{QE}}
= \mathop{\mathrm{argmax}}\limits_{h\in\mathcal{C}}
r_{\mathrm{QE}}(h;s).
$$

### Computational Cost

From an end-to-end perspective, the computation required for each source is:

- Candidate generation: $|\mathcal{C}|$ LLM generations
- With MBR decoding: $|\mathcal{C}|^2$ utility (reference-based metric) calls
- With QE reranking: $|\mathcal{C}|$ QE model (reference-free metric) calls

### Why Does Acceleration Matter?

Take $|\mathcal{C}|=256$: translating a single sentence means generating 256 candidates, and MBR then needs $256^2 = 65{,}536$ utility calls. On top of that, today's utilities are usually neural models, so calling them this many times costs a lot of compute. The quality gains of reranking are well established, but at this cost it fits competitions and offline evaluation far better than real systems with latency and compute budgets.

I care less about squeezing out another fraction of a point on a leaderboard than about whether a method can actually be used. Reranking already works well; what holds it back is mostly cost, so this post collects several ways to bring that cost down.

There are two directions for acceleration: reduce the number of candidate generations, or reduce the number of utility / QE model calls.

## Acceleration Methods

| Method | Applies to | Computation reduced | Key idea |
| --- | --- | --- | --- |
| Quit [[3]](#ref-3) | MBR / QE | Candidate generation + scoring | Generate incrementally; stop early once the best score stabilizes |
| PMBR [[4]](#ref-4) | MBR | Utility calls | Compute only part of the score matrix; recover the rest via low-rank matrix completion |
| PruneMBR [[5]](#ref-5) | MBR | Utility calls | Turn "pick the highest score" into "estimate each candidate's probability of winning"; estimate it with bootstrap and prune |
| CBMBR [[6]](#ref-6) | MBR | Utility calls | Cluster the pseudo-references and use cluster centroids instead of all of them; equally weighted centroids also reduce multimodal bias |
| BayesOpt [[7]](#ref-7) | QE | QE model calls | Model scores with a Gaussian process; only score candidates worth scoring |

### A Unified Method: Quit

Most existing acceleration methods target only the scoring stage of MBR and QE, leaving candidate generation itself largely untouched. Yet in today's LLM-based translation systems, candidate generation is often more expensive than reranking. Quit (Quantifying Uncertainty for Incremental Termination) [[3]](#ref-3) reduces the cost of both candidate generation and reranking by treating the whole "generate-then-rerank" pipeline as a sequence of decisions.

**Incremental generation.** At each step, generate $k$ new candidates ($k=8$ in the paper); after step $i$ the candidate set is $\mathcal{C}_i$ with $|\mathcal{C}_i| = ik$ (for MBR, the pseudo-reference set grows along with it). Let $\hat{y}_i$ be the output selected by the reranker at step $i$ and $R_i(s)$ its score, i.e. the current best score; the full budget $N_{\max}$ corresponds to $N_{\max}/k$ steps.

**Ideal objective: the early-stopping risk.** Stopping at step $i$ costs the true quality gap between its output and the full-budget output:

$$
\mathcal{L}_i(s) = \big|\, q(\hat{y}_i, y^{\mathrm{ref}}) - q(\hat{y}_{N_{\max}/k}, y^{\mathrm{ref}}) \,\big|.
$$

The smaller $\mathcal{L}_i(s)$, the less quality is lost by stopping early. But at inference time neither the reference nor the full-budget output is available, so it cannot be computed.

**Approximation 1: replace true quality with the reranker score.** This removes the dependence on the reference:

$$
\mathcal{L}^{\mathrm{proxy}}_i(s) = \big|\, R_i(s) - R_{N_{\max}/k}(s) \,\big|.
$$

But $R_{N_{\max}/k}(s)$ is only known once the full budget is used, so this is still not computable.

**Approximation 2: replace the gap to the end with recent fluctuation.** Partition the steps into non-overlapping windows $\mathcal{B}_b$ of size $w$ and look only at the range of the best score within a window:

$$
\Delta_b R(s) = \max_{j\in\mathcal{B}_b} R_j(s) - \min_{j\in\mathcal{B}_b} R_j(s).
$$

A large $\Delta_b R(s)$ means the best score is still improving noticeably; a small one means it has locally stabilized and further generation has limited benefit. The paper verifies with PRR that it correlates with the true risk $\mathcal{L}_i(s)$.

**Stopping rule.** Stop once the fluctuation falls below a threshold $\alpha$ ($10^{-3}$ in the paper):

$$
T = \min\{\, b\cdot w : \Delta_b R(s) \le \alpha \,\}.
$$

Across 3 NMT models and 19 language pairs, Quit achieves end-to-end speedups of 1.47–2.66× for MBR and 3.43–4.12× for QE reranking while preserving translation quality.

### MBR Decoding: Reducing Utility Calls

The bottleneck of MBR is computing a $|\mathcal{C}| \times |\mathcal{C}|$ score matrix $M$ with $M_{ij} = u(h_i, h_j)$. The three methods below reduce this cost from three angles: computing fewer matrix entries, evaluating fewer candidates, and using fewer pseudo-references.

#### PMBR

PMBR [[4]](#ref-4) observes, and verifies experimentally, that the score matrix $M$ has a low-rank structure, so MBR can be cast as a matrix completion problem:

1. Compute only a random subset of the entries of $M$;
2. Use Alternating Least Squares (ALS) to factorize $M \approx UV^\top$ and fill in the missing entries;
3. Average each row of the completed matrix and select the MBR output.

On WMT22 (en↔de, en↔ru), it matches the COMET22 score of vanilla MBR with only $1/16$ of the utility computations, and outperforms other approximation baselines.

#### PruneMBR

The core idea of PruneMBR [[5]](#ref-5) (EMNLP 2023 Best Short Paper): **we don't need each candidate's exact expected utility, only whether it can still win.**

Write the candidate set as $H$ and the expected utility as $U(y, \mathcal{Y}) = \mathbb{E}_{\hat{y}\sim\mathcal{Y}}[u(y,\hat{y})]$. The true MBR objective is $\mathop{\mathrm{argmax}}_{y\in H} U\big(y, p_\theta(\cdot\mid x)\big)$; standard MBR approximates $p_\theta$ with pseudo-references $R$, costing $|H||R|$ utility calls. PruneMBR replaces "pick the highest score" with estimating each candidate's probability of being the true winner:

$$
p\Big(\bigwedge_{\bar{y}\in H} U\big(y, p_\theta(\cdot\mid x)\big) \ge U\big(\bar{y}, p_\theta(\cdot\mid x)\big)\Big).
$$

**Approximation 1: bootstrap.** Since $p_\theta$ is unknown, resample the current pseudo-references $R_t$ with replacement and approximate it by "winning on a bootstrap sample":

$$
\mathbb{E}_{\hat{R}_t\sim\mathrm{boot}(R_t)}\,
\mathbb{1}\Big(\bigwedge_{\bar{y}\in H_t} U(y,\hat{R}_t) \ge U(\bar{y},\hat{R}_t)\Big).
$$

**Approximation 2: compare only with the current winner.** When $H_t$ is large this has high variance, so compare only with the current winner $\bar{\bar{y}} = \mathop{\mathrm{argmax}}_{\bar{y}\in H_t} U(\bar{y}, R_t)$:

$$
\mathbb{E}_{\hat{R}_t\sim\mathrm{boot}(R_t)}\,
\mathbb{1}\big(U(y,\hat{R}_t) \ge U(\bar{\bar{y}},\hat{R}_t)\big).
$$

Beating everyone implies beating $\bar{\bar{y}}$, so this is an upper bound of the previous quantity and independent of $|H_t|$. If it falls below $1-\alpha$, $y$ is pruned. Bootstrap only recombines utilities already computed, so it needs no extra calls.

**Algorithm.** The number of pseudo-references doubles according to a schedule (from 8 for COMET, up to 256). Each step adds pseudo-references, computes new utilities only for the surviving candidates, and prunes via bootstrap; it stops when one candidate remains or the last step is reached. Too few initial pseudo-references bias the bootstrap and risk pruning the winner, so $r_1$ sets the speed–accuracy trade-off.

Results: statistically indistinguishable from standard MBR, with at least 7× fewer utility calls for chrF++ and 15× for COMET.

#### CBMBR

CBMBR [[6]](#ref-6) rests on an elegant idea: **pseudo-references are just vectors in feature space, so nearby ones can be merged into a single representative point.**

COMET encodes the source, candidate, and reference separately with an encoder $f$, then scores them with an output layer $\phi$: $u(h, h_s) = \phi\big(f(s), f(h), f(h_s)\big)$. Encoding takes only $|\mathcal{C}|$ passes; the bottleneck is the $|\mathcal{C}|^2$ evaluations of $\phi$. The method relies on two assumptions:

1. **Semantically similar sentences have nearby vectors.** On STS-B, COMET's sentence vectors reach a Pearson correlation of 73.6 with human similarity, even above LaBSE (72.7).
2. **$\phi$ is approximately linear.** Within a cluster $\mathcal{K}_i$, "score, then average" can then be swapped for "average, then score":

$$
\frac{1}{|\mathcal{K}_i|}\sum_{h_s\in\mathcal{K}_i}\phi\big(f(s), f(h), f(h_s)\big)
\approx
\phi\big(f(s), f(h), c_i\big),\quad c_i = \frac{1}{|\mathcal{K}_i|}\sum_{h_s\in\mathcal{K}_i} f(h_s).
$$

Concretely, the pseudo-reference vectors are clustered into $k$ clusters with kmeans++ and k-means ($k=64$ in the paper), and

$$
\hat{y}_{\mathrm{CBMBR}}
= \mathop{\mathrm{argmax}}\limits_{h\in\mathcal{C}}
\frac{1}{k}\sum_{i=1}^{k} \phi\big(f(s), f(h), c_i\big),
$$

reducing the complexity from $O(|\mathcal{C}|^2)$ to $O(|\mathcal{C}|\,k)$.

**Why can it beat standard MBR?** The formula averages the centroids with **equal weights**. Standard MBR treats every pseudo-reference equally, so when candidates are multimodal (e.g. from several systems), the translation with the most samples dominates. CBMBR gives each cluster one vote and is more robust to this bias. An ablation supports this: CBMBR$_{\mathrm{cnt}}$, which weights centroids by cluster size and thus approximates standard MBR more closely, scores almost the same as standard MBR (86.6 vs 86.7) and below CBMBR (87.0).

**Results.** With diverse candidates, it stays within 0.1 COMET of standard MBR while speeding up the expected-utility computation 5.7× (1.8× end to end). With multi-system candidates, it beats standard MBR by up to 0.5 COMET. Its limitation is that it only applies to metrics that encode each sentence independently.

### QE Reranking: Reducing QE Model Calls

#### BayesOpt

QE reranking needs only $|\mathcal{C}|$ calls, but QE models (e.g. CometKiwi) are themselves expensive. Cheng et al. [[7]](#ref-7) cast reranking as a Bayesian optimization problem: the goal is to find $\mathop{\mathrm{argmax}}_{h\in\mathcal{C}} r_{\mathrm{QE}}(h;s)$ without scoring every candidate.

1. Model the scores of unscored candidates with a Gaussian process, using an RBF kernel over candidate sentence vectors (mean-pooled representations obtained as a by-product of generation, at almost no extra cost);
2. At each step, use an acquisition function (Expected Improvement) to balance exploration (candidates that differ from those already scored) and exploitation (candidates similar to high-scoring ones), and pick the next candidate to score;
3. When the scoring budget is exhausted, return the highest-scoring candidate among those scored.

With 200 candidates, it needs only about 70 CometKiwi calls on average to match the result of scoring a random subset of 180 candidates. The paper also proposes a multi-fidelity variant that uses a cheaper, smaller QE model as a proxy to further improve efficiency.

## References

1. <span id="ref-1"></span>*CMU Advanced NLP Spring 2025 (7): Decoding Algorithms*. [YouTube](https://www.youtube.com/watch?v=cN8yX_ZZWJw&list=PLqC25OT8ZpD3WxQ0FwWMGPS_BcWdcKyZy&index=8)
2. <span id="ref-2"></span>Tom Kocmi, Ekaterina Artemova, Eleftherios Avramidis, et al. *Findings of the WMT25 General Machine Translation Shared Task: Time to Stop Evaluating on Easy Test Sets*. WMT 2025. [ACL Anthology](https://aclanthology.org/2025.wmt-1.22/)
3. <span id="ref-3"></span>Guangyu Chen, Boxuan Lyu, Hidetaka Kamigaito, Kotaro Funakoshi, Manabu Okumura. *Quit While You're Ahead: Quit for Efficient Candidate Generation in Machine Translation Reranking*. arXiv:2609.00588, 2026. [arXiv](https://arxiv.org/abs/2609.00588)
4. <span id="ref-4"></span>Firas Trabelsi, David Vilar, Mara Finkelstein, Markus Freitag. *Efficient Minimum Bayes Risk Decoding using Low-Rank Matrix Completion Algorithms*. arXiv:2406.02832, 2024. [arXiv](https://arxiv.org/abs/2406.02832)
5. <span id="ref-5"></span>Julius Cheng, Andreas Vlachos. *Faster Minimum Bayes Risk Decoding with Confidence-based Pruning*. EMNLP 2023. [ACL Anthology](https://aclanthology.org/2023.emnlp-main.767/)
6. <span id="ref-6"></span>Hiroyuki Deguchi, Yusuke Sakai, Hidetaka Kamigaito, Taro Watanabe, Hideki Tanaka, Masao Utiyama. *Centroid-Based Efficient Minimum Bayes Risk Decoding*. Findings of ACL 2024. [ACL Anthology](https://aclanthology.org/2024.findings-acl.654/)
7. <span id="ref-7"></span>Julius Cheng, Maike Züfle, Vilém Zouhar, Andreas Vlachos. *A Bayesian Optimization Approach to Machine Translation Reranking*. NAACL 2025. [ACL Anthology](https://aclanthology.org/2025.naacl-long.145/)
