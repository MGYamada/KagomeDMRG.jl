# Literature-motivated VBC preparation motifs

作成: 2026-09-21。実装は [vbc_motif_tools.jl](../../examples/vbc_motif_tools.jl)。
これは有限円筒に加える一時的なexchange誘導の定義であり、文献の波動関数・相・energyの再現ではない。
既存の `period9` / `period27` 試行状態とは異なる構成で、名前の置換は行っていない。

## 目視確認した資料と適用範囲

- [Fang et al., arXiv:2306.09563v1](https://arxiv.org/pdf/2306.09563v1),
  Fig.1(b), PDF p.2: cyanのhourglass領域にあるbowtieの6辺を選ぶ。図のbond色による強弱分類は再現しない。
  共有頂点を持つ2三角形の5site clusterであり、並進胞は9site。
  図は4種類のbondを区別するが、本実装はhourglassの6辺だけを選択する。
  他の強いred bondを含む**全4種類の強弱模様は再現していない**。
- [Cheng and Li, arXiv:2507.20308v1](https://arxiv.org/pdf/2507.20308v1),
  Fig.5, PDF p.9: blackの中心hexagon、redの6翼、greenのdumbbell、grayの残りを対応付けた。
  取得PDFの先頭は2025-07-27のv1。出版情報は
  [PRB 113, 085136 (2026)](https://journals.aps.org/prb/abstract/10.1103/5tvd-253q)
  と照合したが、出版版の図がv1と完全に同一かは未確認。

両PDFの該当ページをPopplerで画像化し、拡大した図も目視した。
hourglassは4class全体の読み取りを完了したという主張を避け、明確なbowtie領域に限定した。
windmillの図中の小数は丸められたcorrelationで、以下では誘導strengthを選ぶためだけに使用する。
実測値のtarget・収束判定・energy誤差棒には使用しない。

## 座標と並進原点

`a1=(1,0), a2=(1/2,sqrt(3)/2)`、`A=0, B=a1/2, C=a2/2`。
以下の `A(x,y)` 等は無限平面上のcell座標で、負の座標も許す。
原点 `origin=(ox,oy)` は表の全siteを `ox*a1+oy*a2` だけ動かす。
その後に軸方向を `0 <= x < Lx` で切り、`Ly*a2` で円周を同一視する。
端の追加bondは作らない。既存のbond順序・向き・`wy`を保持する。

`Ly` は3の倍数を要求する。N54比較は `(Lx,Ly)=(6,3), Q=6`。
N9 `(1,3)` はHamiltonian照合用の切り出しとして許すが、完全なVBC clusterを含むとは限らない。
以下の並進対称性は無限模様の性質であり、軸OBC全系の対称性ではない。

## Hourglass: 9site胞中の共有C頂点bowtie

原点の中心を `C(0,0)` とする。次の各三角形の3辺、計6辺を選ぶ。

| 三角形 | 3頂点 |
| --- | --- |
| lower | `C(0,0), A(0,0), B(0,0)` |
| upper | `C(0,0), A(0,1), B(-1,1)` |

中心を `x-y = 0 (mod 3)` のC siteへ並進する。
基本並進は `(1,1), (-1,2)`、determinantは3、したがって9site胞となる。
胞中18 NN bondのうち6辺がbowtie、残り12辺がbackgroundである。
選択した辺のraw scoreを1、他を0とし、周期胞平均を引いて最大正値を1へ規格化する。
したがって **bowtie weight=1、background weight=-1/2**。
これは図の一部を探索可能にする選択的誘導であり、文献のbond energy順位全体をencodeしていない。

## Windmill: 27site胞と18site cluster

中心hexagonの中心は `a1*(1/2)+a2*(1/2)`、Cartesianでは `(3/4,sqrt(3)/4)`。
並進は `(3,0),(0,3)`。hexagonの6siteと、6翼の12siteが18site clusterを作る。
残る9siteは3本のdumbbellである。以下を同じ27site胞の周期像として扱う。

中心hexagonは次のcycleの隣接6辺である。

`B(0,0) - A(1,0) - C(1,0) - B(0,1) - A(0,1) - C(0,0) - B(0,0)`

| 赤いwing | endpoint 1 | endpoint 2 |
| --- | --- | --- |
| 1 | `B(1,0)` | `C(2,-1)` |
| 2 | `C(1,-1)` | `A(1,-1)` |
| 3 | `A(0,0)` | `B(-1,0)` |
| 4 | `B(-1,1)` | `C(-1,1)` |
| 5 | `C(0,1)` | `A(0,2)` |
| 6 | `A(1,1)` | `B(1,1)` |

| 緑のdumbbell | 3site chain; 隣接する2辺を選ぶ |
| --- | --- |
| 1 | `B(-2,2) - A(-1,2) - B(-1,2)` |
| 2 | `A(-1,0) - C(-1,0) - A(-1,1)` |
| 3 | `C(0,2) - B(0,2) - C(1,1)` |

図の表示窓には隣接する周期像も含まれるため、図に見えるdumbbell数と胞内の数を混同しない。

| class | 胞内bond数 | raw score `s` | 最終weight |
| --- | ---: | ---: | ---: |
| black ring | 6 | 0.3722 | 0.398760668 |
| red wing | 6 | 0.6180 | 1 |
| green dumbbell | 6 | 0.4466 | 0.580746861 |
| gray background | 36 | 0.0743 | -0.329917921 |

`sbar=0.2091777777777778`, `w=(s-sbar)/(0.6180-sbar)`。
これは無限周期模様での平均0・最大正値1の規格化で、OBCで欠けた辺に合わせて再中心化はしない。
両templateの同じlambdaは最大の強化量を合わせるが、同じ準備biasの大きさや最適化費用を保証しない。

## APIと除去

```julia
include("examples/vbc_motif_tools.jl")
lattice = kagome_cylinder(6, 3)
template = vbc_bond_template(lattice; kind=:windmill, origin=(0,0))
prepared = vbc_prepared_lattice(lattice, template, 0.3)
unperturbed = vbc_prepared_lattice(lattice, template, 0.0)
```

`Jxy_b=Jz_b=1+lambda*w_b` とし、全NN bondを残す。
base格子は毎回canonicalなisotropic `J=1` を渡す。前段のprepared格子へ重ね掛けしない。
lambdaは有限かつ非負、全couplingは正を要求する。
テンプレートはgeometry・bond順序・windingと照合する。
lambda=0では元のNN模型へ厳密に戻る。固定Qの状態準備・MPS charge・複素性は呼び出し側が担当する。
最終比較のHamiltonianに誘導を残さず、warm start時にHamiltonian-dependent environmentを再構築する。

## 実際に行った幾何確認

2026-09-21にPythonの整数座標列挙で3x3周期平面のNN54辺と対応付けた。
hourglassは54辺中18辺を選択し、mod3並進stabilizerは `(0,0),(1,1),(2,2)`。
windmillはring/wing/dumbbellが各6辺、background36辺、stabilizerは `(0,0)` のみだった。
別agentの独立Cartesian確認でも、選択辺の距離二乗1/4、classの無重複、
windmill18siteとdumbbell9siteの非重複な27site被覆が一致した。
この文書の幾何確認でDMRGは実行していない。Hamiltonianの独立spin-basis照合と
誘導除去・保存の検査は [validate_vbc_preparation.jl](../../examples/validate_vbc_preparation.jl)
および比較campaignの報告に分ける。

残る解釈上の制限: hourglassの全4bond classは未再現、文献波動関数は未再現、
一つの原点と点群domainの準備だけでは初期条件依存を網羅しない。
誘導除去後に模様が残っても、自発秩序・基底状態・相をそれだけで確定しない。
