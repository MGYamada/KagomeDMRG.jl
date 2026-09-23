# 最近接1/9の静的研究: 規約と判定

更新日: 2026-09-23。現在の優先順位は[ROADMAP](../../ROADMAP.md)、
競合文献と方針の根拠は[研究方向の再検討](research_direction_20260920.md)を参照する。
本書は数値規約・精度基準・記録への入口に限定する。
当面はコードベース整備を優先し、追加の静的研究とCSL計算は保留する。以下の規約は再開時にも維持する。

## 模型と幾何

`H(h)=H0−hΣSz`、`H0=Σ_NN Si·Sj`、J=1、θ=0。
整数chargeは `Q=2ΣSz`、1/9では `Q0=N/9`、`N%9=0`。
一様磁場は固定Q内で定数なので、同じQをhの格子上で再最適化しない。
文献の `h/J≈0.35–0.42` は[特定手法の報告](https://arxiv.org/abs/2306.09563)であり、合格区間ではない。

軸はOBC、wrapは `Ly*a2`、順序はx、y、A/B/C。文献のYC名をLyから推定しない。
9-site周期の一方向を `A1=a1+a2, A2=−a1+2a2` と取ると `3a2=A1+A2`。
27-site周期 `3a1,3a2` もLyが3の倍数で整合するが、周期整合性は十分なbulkの保証ではない。
N27は27-site模様の一周期、N54も狭幅系である。端除外幅・長さ・幅への依存を分ける。

[既存のperiod9/27初期状態](order_seed_design.md)は文献のVBC波動関数ではない。
結合変調・局所磁場による準備は別Hamiltonianとして記録し、除去後の同じNN模型で状態を比較する。
固定Qで `〈S+〉=0` でも磁気秩序は否定できない。Sz・bond分布と縦横相関を併用する。

## 固定Qから磁場区間を求める

全交換energyを `E(Q)` とし、`FQ(h)=E(Q)−hQ/2` の下側包絡を取る。
energy/siteをそのまま次の式へ代入しない。隣接sectorだけなら

```text
h− = E(Q0) − E(Q0−2)
h+ = E(Q0+2) − E(Q0)
Δh = E(Q0+2) + E(Q0−2) − 2E(Q0)
```

比較したsector集合全体では

```text
h_lower = max_{Q<Q0} 2[E(Q0)−E(Q)]/(Q0−Q)
h_upper = min_{Q>Q0} 2[E(Q)−E(Q0)]/(Q−Q0)
```

`h_lower>=h_upper` なら、そのenergy集合には有限幅の安定区間がない。
隣接Qだけでは遠方sectorへの飛び越しを除外できず、未探索範囲を記録する。
θ=0の等方模型はspin反転により `E(Q)=E(−Q)`。既知の鏡像sectorを重複計算しない。
Qを変える場合は別の初期MPSを使い、保存状態のchargeを変更して再開しない。

変分上界同士の差から厳密な境界の上下限は得られない。
追加磁化 `Sz_i(Q0±2)−Sz_i(Q0)` の位置を調べ、端の励起と内部への磁化注入を区別する。
有限系の磁化区間、bulk plateau、中性gap、秩序、Hall応答は別々の結論である。

## 精度と解釈

既存の静的比較の目安を維持する。値は判定対象の実測量であり、設定cutoffで代用しない。

| 量 | 上限 |
| --- | ---: |
| 親との差・solverとの差・取得済みの末尾sweep差の最大絶対値/N | 1e−6 J |
| 最大Sz変化 | 1e−4 |
| 最大bond energy変化 | 1e−4 J |
| 分散/N | 1e−5 J² |
| 最終sweepの最大実測切断誤差 | 1e−6 |

固定χの定常性とχ変更の応答を区別する。χと追加sweepが同時に変わる比較は純粋なχ依存ではない。
1-sweep比較のenergy値は `max(abs(Enew−Eparent), abs(E_MPO−Esolver))/N` とする。
1-sweep呼出しに内部の二sweep差を補わず、最終solver energyとMPO期待値の差も別に記録する。
5条件通過は基底状態・相同定の十分条件ではなく、同じχ・sweep数も精度一致を意味しない。
候補差とsweep・χ・準備条件による変動を比較し、必要な状態に精度改善を集中する。
分散はHamiltonianの分散であり、基底energyの誤差棒ではない。

未達結果は探索と次ケース選択に使えるが、収束した優劣・境界・相とは主張しない。
保存・電荷・Hermiticity等の必須整合性は維持し、違反した状態を物理的証拠として扱わない。
高価な全観測は代表状態へ集中させ、欠測と未達を保存する。基準を結果に合わせて緩めない。

## 実行と輸送への接続

ケース・総wall上限・χ上限・測定頻度を物理比較ごとに事前設定する。
Julia/BLAS各1 thread、数値process一つを出発点に、途中保存と外側の停止監視を使う。
旧checkpointは元のsource・依存環境で厳密に読み込み、新しい由来へ付け替えない。
時間切れではworker全体の完了と保存済みbatchを区別する。

研究再開後、既知CSLは別模型 `J2=J3=0.5,Q=0` の測定対照として並行する。
静的NN研究や限定した追跡診断を待たせないが、1/9の量子化ポンプ解釈には、既知非零・零応答、
枝連続性、複数cut、実空間/Schmidt一致、χ・刻み・サイズ依存が必要である。
詳細は[flux設計](../flux_insertion_design.md)と[CSL対照](p3_csl_control_design.md)を参照する。

## 研究記録

過去の設定・数値・時間切れ・検証範囲は各記録を正本とし、新しい結果で上書きしない。

| 研究 | 記録 |
| --- | --- |
| N18の静的基準 | [Q0/2/4](p4_static18_validation.md)、[Q6追加](p4_static18_q6_validation.md)、[中性二重項](p4_static18_neutral_validation.md) |
| N27の準備と再開 | [初回停止](p4_static27_pilot.md)、[保存状態](p4_static27_progress.md)、[累積4](p4_static27_refine.md) |
| 隣接Q・周期・長さ・CSL準備 | [段階的キャンペーン](p4_p3_staged_campaign.md)、[隣接Q追加](p4_static27_sector_refinement.md) |
| Q5の精度依存とCSL診断 | [固定χ](p4_q5_csl_fixed_chi_followup.md)、[同一親χ比較](p4_q5_matched_chi_comparison.md) |
| 最新のN27 χ512比較 | [Q1追加・Q5固定χ](p4_static512_sector_followup.md) |
| N54・Q6の結合準備比較 | [3枝の除去後構造と精度](p4_vbc54_preparation_comparison.md)、[文献motif対応](vbc_motif_mapping.md) |
| N54・Q6の同一親χ比較 | [random／windmillのχ128/256・追加2 sweeps](p4_vbc54_matched_parent_comparison.md) |
| 既知CSLの資源判断 | [保存結果の再監査とN72準備案](p3_csl_resource_decision_20260920.md)、[N72実行可能性監査](p3_csl72_feasibility_20260922.md) |
