+++
title = "機械翻訳リランキングの高速化"
date = 2026-09-30
draft = false
summary = "機械翻訳におけるリランキング高速化手法のまとめ。"
tags = ["Machine Translation", "LLM Decoding", "Reranking"]
+++

## LLM のデコーディングとリランキング

LLM のデコーディングは大きく 2 種類に分けられます [[1]](#ref-1)。

1. 各ステップのトークン分布 $P(y_i \mid X, y_{<i})$ に作用する関数 $f$ を設計し、$f\big(P(y_i \mid X, y_{<i})\big)$ に基づいてトークンを 1 つずつ選ぶ方法。Greedy、ビームサーチ、Top-$k$/Top-$p$ サンプリングはいずれもこちらに属します。
2. 出力全体 $y$ に作用する関数 $g$ を設計する方法。例えば $g\big(P(y \mid X)\big)$ で、$P(y \mid X) = \prod_i P(y_i \mid X, y_{<i})$ です。この過程は 2 段階からなります。まず LLM に多数の候補を生成させ、次に $g$ で最良の 1 つを選びます。そのため、この種のデコーディングはリランキング（Reranking）とも呼ばれます。

## 機械翻訳への応用

機械翻訳タスクにおいて、原文を $s$、モデルの翻訳を $y$、参照訳を $y^{\mathrm{ref}}$、翻訳品質の評価指標を $q(y, y^{\mathrm{ref}})$ と定義します。

機械翻訳でリランキングを行う背景には、シンプルな直感があります。翻訳を 1 つだけ生成するよりも、候補を複数生成しておけば、より良いものを選べるはずだ、というものです。
ここでは候補集合を $\mathcal{C}$ と書きます。
機械翻訳で最もよく使われるのは MBR Decoding と QE Reranking です。

WMT General Machine Translation Shared Task では、この 2 つのリランキング手法が上位システムに繰り返し採用されています [[2]](#ref-2)。

### MBR Decoding

この手法はベイズ決定理論に基づき、他の候補と「平均的に最もよく一致する」候補を選びます。この一致度を測る関数を **utility** と呼び、$u(h, h_s)$ と書きます。これは $h_s$ を参照とみなしたときの $h$ の良さを表します。選択規則は次のとおりです。

$$
\hat{y}_{\mathrm{MBR}}
= \mathop{\mathrm{argmax}}\limits_{h\in\mathcal{C}}
\frac{1}{|\mathcal{C}|}
\sum_{h_s\in\mathcal{C}} u(h,h_s).
$$

デコード時には本当の参照訳 $y^{\mathrm{ref}}$ は手に入らないため、他の候補 $h_s$ を擬似参照（pseudo-reference）として用います。そのため utility には通常、reference-based metric（例：COMET、BLEURT、chrF++）をそのまま使います。つまり、前述の $q(y, y^{\mathrm{ref}})$ の参照訳を擬似参照に置き換えた $u(h, h_s) = q(h, h_s)$ です。上式の平均は $h$ の期待 utility であり、期待 utility が最も高い候補を選ぶことで、その metric の好みに合った、より良い出力が得られます。

以下、MBR に関する箇所での「utility」はこの reference-based metric を指し、「utility 呼び出し」はその metric を 1 回呼び出すことを指します。

### QE Reranking

この手法は reference-free なモデル（例：CometKiwi）で各候補を直接スコアリングし、次のように選択します。

$$
\hat{y}_{\mathrm{QE}}
= \mathop{\mathrm{argmax}}\limits_{h\in\mathcal{C}}
r_{\mathrm{QE}}(h;s).
$$

### 計算量の分析

入力から出力までをエンドツーエンドで見ると、1 つの原文あたりに必要な計算は次のとおりです。

- 候補生成：$|\mathcal{C}|$ 回の LLM 生成
- MBR Decoding の場合：$|\mathcal{C}|^2$ 回の utility（reference-based metric）呼び出し
- QE Reranking の場合：$|\mathcal{C}|$ 回の QE モデル（reference-free metric）呼び出し

### なぜ高速化が重要なのか？

$|\mathcal{C}|=256$ を例にとると、1 文を翻訳するだけで 256 個の候補を生成し、MBR ではさらに $256^2 = 65536$ 回の utility 呼び出しが必要になります。しかも現在の utility は一般にニューラルネットワークであり、大量に呼び出すには高価な計算資源が必要です。リランキングによる品質向上は繰り返し示されていますが、このコストでは、コンペティションやオフライン評価には向いていても、レイテンシや計算資源に制約のある実システムにそのまま組み込むのは困難です。

リーダーボードでスコアをもう少し伸ばすことよりも、私はある手法が実際に使えるかどうかに関心があります。リランキングの効果はすでに十分で、足かせになっているのは主にコストです。そこでこの記事では、そのコストを下げるいくつかの手法をまとめました。

高速化には 2 つの方向があります。候補生成の回数を減らすか、utility / QE モデルの呼び出し回数を減らすかです。

## 高速化手法

| 手法 | 対象 | 削減する計算 | 基本アイデア |
| --- | --- | --- | --- |
| Quit [[3]](#ref-3) | MBR / QE | 候補生成 + スコアリング | 候補を逐次生成し、最良スコアが安定したら早期停止 |
| PMBR [[4]](#ref-4) | MBR | utility 呼び出し | スコア行列の一部だけを計算し、残りを低ランク行列補完で復元 |
| PruneMBR [[5]](#ref-5) | MBR | utility 呼び出し | 「最高スコアを選ぶ」を「各候補が最終的に勝つ確率を推定する」に置き換え、bootstrap で推定して枝刈り |
| CBMBR [[6]](#ref-6) | MBR | utility 呼び出し | 擬似参照をクラスタリングし、全擬似参照の代わりにクラスタ中心を使う。等しい重みの中心は多峰性によるバイアスも抑える |
| BayesOpt [[7]](#ref-7) | QE | QE モデル呼び出し | ガウス過程でスコアをモデル化し、スコアリングする価値のある候補だけを評価 |

### 統合的な手法：Quit

既存の高速化手法の多くは MBR と QE のスコアリング部分だけを対象としており、候補生成そのものはほとんど扱われていません。しかし現在の LLM ベースの翻訳システムでは、候補生成のほうがリランキングよりもコストが高いことが多いです。Quit（Quantifying Uncertainty for Incremental Termination）[[3]](#ref-3) は候補生成とリランキングの両方の計算を削減し、「生成→リランキング」のパイプライン全体を一連の意思決定問題として扱います。

**逐次生成。** 各ステップで $k$ 個の新しい候補を生成します（論文では $k=8$）。ステップ $i$ の後の候補集合を $\mathcal{C}_i$、$|\mathcal{C}_i| = ik$ とします（MBR の場合は擬似参照集合も一緒に増えます）。ステップ $i$ でリランカーが選んだ出力を $\hat{y}_i$、そのスコア（その時点の最良スコア）を $R_i(s)$ とし、予算 $N_{\max}$ を使い切るまでのステップ数は $N_{\max}/k$ です。

**理想的な目的：早期停止のリスク。** ステップ $i$ で止めたときの損失は、その出力と予算を使い切ったときの出力との真の品質差です。

$$
\mathcal{L}_i(s) = \big|\, q(\hat{y}_i, y^{\mathrm{ref}}) - q(\hat{y}_{N_{\max}/k}, y^{\mathrm{ref}}) \,\big|.
$$

$\mathcal{L}_i(s)$ が小さいほど、早期停止で品質が失われません。しかし推論時には参照訳も予算を使い切った出力もないため、計算できません。

**近似 1：真の品質をリランカーのスコアで置き換える。** これで参照訳への依存がなくなります。

$$
\mathcal{L}^{\mathrm{proxy}}_i(s) = \big|\, R_i(s) - R_{N_{\max}/k}(s) \,\big|.
$$

ただし $R_{N_{\max}/k}(s)$ は予算を使い切るまで分からないため、まだ計算できません。

**近似 2：終点との差を直近の変動で置き換える。** ステップを大きさ $w$ の重ならないウィンドウ $\mathcal{B}_b$ に分け、ウィンドウ内での最良スコアの変動幅だけを見ます。

$$
\Delta_b R(s) = \max_{j\in\mathcal{B}_b} R_j(s) - \min_{j\in\mathcal{B}_b} R_j(s).
$$

$\Delta_b R(s)$ が大きければ最良スコアはまだ明らかに向上しており、小さければ局所的に安定していて、生成を続ける効果は小さいと考えられます。論文では PRR により、これが真のリスク $\mathcal{L}_i(s)$ と相関することを確かめています。

**停止規則。** 変動が閾値 $\alpha$（論文では $10^{-3}$）を下回ったら停止します。

$$
T = \min\{\, b\cdot w : \Delta_b R(s) \le \alpha \,\}.
$$

3 つの NMT モデル・19 言語対において、Quit は翻訳品質を維持しつつ、MBR で 1.47–2.66 倍、QE Reranking で 3.43–4.12 倍のエンドツーエンド高速化を達成しています。

### MBR Decoding：utility 呼び出しの削減

MBR のボトルネックは、$|\mathcal{C}| \times |\mathcal{C}|$ のスコア行列 $M$（$M_{ij} = u(h_i, h_j)$）を計算する必要がある点です。以下の 3 つの手法は、それぞれ「計算する行列要素を減らす」「評価する候補を減らす」「使う擬似参照を減らす」という観点から計算量を削減します。

#### PMBR

PMBR [[4]](#ref-4) は、スコア行列 $M$ が低ランク構造を持つことを観察し、実験的に確かめたうえで、MBR を行列補完問題として捉えます。

1. $M$ の要素のうち、ランダムに選んだ一部だけを計算します。
2. Alternating Least Squares（ALS）で $M \approx UV^\top$ と分解し、欠けている要素を補完します。
3. 補完後の行列で行ごとに平均を取り、MBR の出力を選びます。

WMT22（en↔de、en↔ru）では、通常の MBR の $1/16$ の utility 計算だけで同等の COMET22 スコアを達成し、他の近似手法も上回っています。

#### PruneMBR

PruneMBR [[5]](#ref-5)（EMNLP 2023 Best Short Paper）の基本アイデアは、**各候補の期待 utility を正確に求める必要はなく、まだ勝つ可能性があるかどうかさえ分かればよい**、というものです。

候補集合を $H$、期待 utility を $U(y, \mathcal{Y}) = \mathbb{E}_{\hat{y}\sim\mathcal{Y}}[u(y,\hat{y})]$ とします。MBR の本来の目的は $\mathop{\mathrm{argmax}}_{y\in H} U\big(y, p_\theta(\cdot\mid x)\big)$ で、標準的な MBR は擬似参照 $R$ で $p_\theta$ を近似するため $|H||R|$ 回の utility 呼び出しが必要です。PruneMBR は「最高スコアを選ぶ」を、各候補が真の勝者である確率の推定に置き換えます。

$$
p\Big(\bigwedge_{\bar{y}\in H} U\big(y, p_\theta(\cdot\mid x)\big) \ge U\big(\bar{y}, p_\theta(\cdot\mid x)\big)\Big).
$$

**近似 1：bootstrap。** $p_\theta$ は未知なので、現在の擬似参照 $R_t$ から復元抽出し、「bootstrap サンプル上で勝つ」で近似します。

$$
\mathbb{E}_{\hat{R}_t\sim\mathrm{boot}(R_t)}\,
\mathbb{1}\Big(\bigwedge_{\bar{y}\in H_t} U(y,\hat{R}_t) \ge U(\bar{y},\hat{R}_t)\Big).
$$

**近似 2：現在の勝者とだけ比べる。** $H_t$ が大きいと上式は分散が大きいため、現在の勝者 $\bar{\bar{y}} = \mathop{\mathrm{argmax}}_{\bar{y}\in H_t} U(\bar{y}, R_t)$ とだけ比較します。

$$
\mathbb{E}_{\hat{R}_t\sim\mathrm{boot}(R_t)}\,
\mathbb{1}\big(U(y,\hat{R}_t) \ge U(\bar{\bar{y}},\hat{R}_t)\big).
$$

全員に勝つなら必ず $\bar{\bar{y}}$ にも勝つので、これは前の式の上界であり、$|H_t|$ に依存しません。これが $1-\alpha$ を下回れば $y$ を枝刈りします。bootstrap は計算済みの utility を組み直すだけなので、追加の呼び出しは不要です。

**アルゴリズム。** 擬似参照の数はスケジュールに従って倍々に増やします（COMET では 8 から始めて 256 まで）。各ステップで擬似参照を追加し、残った候補についてだけ新しい utility を計算して、bootstrap で枝刈りします。候補が 1 つになるか最終ステップに達したら停止します。初期の擬似参照が少なすぎると bootstrap のバイアスが大きく、勝者を誤って枝刈りしやすいため、$r_1$ が速度と精度のトレードオフを決めます。

結果：標準的な MBR と統計的な有意差はなく、utility 呼び出しを chrF++ で 7 分の 1 以下、COMET で 15 分の 1 以下に削減しました。

#### CBMBR

CBMBR [[6]](#ref-6) のアイデアは巧妙です。**擬似参照は特徴空間上のベクトルにすぎないので、近いもの同士は 1 つの代表点にまとめられる**、というものです。

COMET はエンコーダ $f$ で原文・候補・参照を個別にエンコードし、出力層 $\phi$ でスコアを出します。つまり $u(h, h_s) = \phi\big(f(s), f(h), f(h_s)\big)$ です。エンコードは $|\mathcal{C}|$ 回で済み、ボトルネックは $|\mathcal{C}|^2$ 回の $\phi$ の計算です。この手法は 2 つの仮定に基づきます。

1. **意味の近い文はベクトルも近い。** STS-B で、COMET の文ベクトルは人手の類似度と Pearson 相関 73.6 を示し、LaBSE（72.7）さえ上回ります。
2. **$\phi$ はほぼ線形である。** するとクラスタ $\mathcal{K}_i$ 内では「スコアを計算してから平均」を「平均してからスコアを計算」に置き換えられます。

$$
\frac{1}{|\mathcal{K}_i|}\sum_{h_s\in\mathcal{K}_i}\phi\big(f(s), f(h), f(h_s)\big)
\approx
\phi\big(f(s), f(h), c_i\big),\quad c_i = \frac{1}{|\mathcal{K}_i|}\sum_{h_s\in\mathcal{K}_i} f(h_s).
$$

具体的には、kmeans++ と k-means で擬似参照ベクトルを $k$ 個のクラスタにまとめ（論文では $k=64$）、次のように選択します。

$$
\hat{y}_{\mathrm{CBMBR}}
= \mathop{\mathrm{argmax}}\limits_{h\in\mathcal{C}}
\frac{1}{k}\sum_{i=1}^{k} \phi\big(f(s), f(h), c_i\big).
$$

計算量は $O(|\mathcal{C}|^2)$ から $O(|\mathcal{C}|\,k)$ に下がります。

**なぜ標準的な MBR を上回れるのか？** 上式はクラスタ中心を**等しい重み**で平均します。標準的な MBR は各擬似参照を平等に扱うため、候補が多峰的な場合（例えば複数システム由来）にはサンプル数の多い訳が支配的になります。CBMBR は各クラスタに 1 票ずつしか与えないので、このバイアスに頑健です。比較実験もこれを裏付けています。クラスタの大きさで重み付けし、標準的な MBR により近い CBMBR$_{\mathrm{cnt}}$ は、標準的な MBR とほぼ同じスコア（86.6 vs 86.7）で、CBMBR（87.0）には及びませんでした。

**結果。** 多様な候補では標準的な MBR との差が 0.1 COMET 以内で、期待 utility の計算が 5.7 倍（エンドツーエンドで 1.8 倍）高速化されました。複数システムの候補では、標準的な MBR を最大 0.5 COMET 上回りました。制約は、各文を個別にエンコードできる metric にしか使えない点です。

### QE Reranking：QE モデル呼び出しの削減

#### BayesOpt

QE Reranking の計算は $|\mathcal{C}|$ 回の呼び出しだけですが、QE モデル（例：CometKiwi）自体のコストが高いです。Cheng ら [[7]](#ref-7) はリランキングをベイズ最適化（Bayesian Optimization）問題として定式化しました。目標は $\mathop{\mathrm{argmax}}_{h\in\mathcal{C}} r_{\mathrm{QE}}(h;s)$ を見つけることですが、すべての候補をスコアリングする必要はありません。

1. ガウス過程（Gaussian Process）で未スコアの候補のスコアをモデル化します。カーネルは候補の文ベクトル（生成時についでに得られる mean-pooled 表現で、追加コストはほぼゼロ）上の RBF カーネルです。
2. 各ステップで獲得関数（Expected Improvement）を使い、exploration（スコア済みの候補と大きく異なるもの）と exploitation（高スコアの候補に似ているもの）のバランスを取りながら、次にスコアリングする候補を選びます。
3. スコアリングの予算を使い切ったら、スコア済みの候補のうち最もスコアの高いものを返します。

候補数 200 の設定では、平均約 70 回の CometKiwi 呼び出しだけで、ランダムに選んだ 180 候補をスコアリングした場合と同等の結果が得られます。論文ではさらに multi-fidelity 版も提案されており、より安価な小型 QE モデルのスコアを代理として取り入れることで、効率をさらに高めています。

## 参考文献

1. <span id="ref-1"></span>*CMU Advanced NLP Spring 2025 (7): Decoding Algorithms*. [YouTube](https://www.youtube.com/watch?v=cN8yX_ZZWJw&list=PLqC25OT8ZpD3WxQ0FwWMGPS_BcWdcKyZy&index=8)
2. <span id="ref-2"></span>Tom Kocmi, Ekaterina Artemova, Eleftherios Avramidis, et al. *Findings of the WMT25 General Machine Translation Shared Task: Time to Stop Evaluating on Easy Test Sets*. WMT 2025. [ACL Anthology](https://aclanthology.org/2025.wmt-1.22/)
3. <span id="ref-3"></span>Guangyu Chen, Boxuan Lyu, Hidetaka Kamigaito, Kotaro Funakoshi, Manabu Okumura. *Quit While You're Ahead: Quit for Efficient Candidate Generation in Machine Translation Reranking*. arXiv:2609.00588, 2026. [arXiv](https://arxiv.org/abs/2609.00588)
4. <span id="ref-4"></span>Firas Trabelsi, David Vilar, Mara Finkelstein, Markus Freitag. *Efficient Minimum Bayes Risk Decoding using Low-Rank Matrix Completion Algorithms*. arXiv:2406.02832, 2024. [arXiv](https://arxiv.org/abs/2406.02832)
5. <span id="ref-5"></span>Julius Cheng, Andreas Vlachos. *Faster Minimum Bayes Risk Decoding with Confidence-based Pruning*. EMNLP 2023. [ACL Anthology](https://aclanthology.org/2023.emnlp-main.767/)
6. <span id="ref-6"></span>Hiroyuki Deguchi, Yusuke Sakai, Hidetaka Kamigaito, Taro Watanabe, Hideki Tanaka, Masao Utiyama. *Centroid-Based Efficient Minimum Bayes Risk Decoding*. Findings of ACL 2024. [ACL Anthology](https://aclanthology.org/2024.findings-acl.654/)
7. <span id="ref-7"></span>Julius Cheng, Maike Züfle, Vilém Zouhar, Andreas Vlachos. *A Bayesian Optimization Approach to Machine Translation Reranking*. NAACL 2025. [ACL Anthology](https://aclanthology.org/2025.naacl-long.145/)
