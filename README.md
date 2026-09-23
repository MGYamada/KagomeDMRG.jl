# KagomeDMRG.jl

スピン1/2の最近接kagome反強磁性Heisenberg模型を、U(1)円筒DMRGで調べる
研究用Juliaパッケージです。主対象は磁化 `M/Msat = 1/9` です。

**当面は研究計算の追加より、コードベースの充実を優先します。**
実行・保存後診断の共通化、APIとデータ契約、テスト・CI、実測に基づく性能改善の順に整えます。
保存後診断の共通APIと小系driverを実装・検証しました。次に設定・実行記録の契約を整えます。
N54の固定χ256追加計算と既知CSLのN72準備は保留します。
開発の完了条件と研究再開時の候補は [ROADMAP](ROADMAP.md) にまとめています。

## 実装と現在の限界

- ITensors.jl + ITensorMPS.jlによる複素U(1)二サイトDMRG、固定磁化sector、格子・flux gauge。
- 局所磁化、結合energy、相関、chirality、Schmidt診断と、実測切断誤差・分散。
- 完了点のcheckpoint保存・再開と、診断に応じたflux刻みの調整・状態復元。
- [保存後の静的診断](docs/static_diagnostics.md)：直接測定と保存状態に共通のAPI、小系のsolve／diagnose分離例。

小系の独立ED照合と保存・測定の検証を研究基盤としています。
[N54・Q6の結合準備比較](docs/research/p4_vbc54_preparation_comparison.md)を実装・実行しました。
続く[random／windmillの同一親χ比較](docs/research/p4_vbc54_matched_parent_comparison.md)も4ケースを完了しました。
全4子で5精度条件は未達で、縮小する中央差から異なる安定した構造を主張できません。
bulk plateau・中性gap・微視的な相は未確立で、既知CSLの非零ポンプ校正も未完了です。
研究再開後も、1/9の量子化ポンプの解釈には既知CSLの校正が必要です。

## 最小実行例

Juliaupを利用する場合、リポジトリのルートでJulia 1.13と研究用環境を準備します。

```sh
juliaup add 1.13
julia +1.13 --project=research --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia +1.13 --project=research --startup-file=no --threads=1
```

Juliaのプロンプトで9サイトを計算します。第2引数はflux角 `theta`（ラジアン）です。

```julia
using KagomeDMRG, LinearAlgebra, ITensors
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()

lattice = kagome_cylinder(1, 3)
result = run_dmrg(lattice, 0.37; seed=11)
@show result.energy
@show result.sz sum(result.sz)
@show result.max_truncation_errors
```

この例は `Q=2Sz_total=1`、物理的な `Sz_total=1/2` です。
小系の動作確認用で、bulk輸送を測る例ではありません。
研究計算には局所版NDTensorsを固定した `research/` の環境を使います。
過去checkpointの再開には、保存時のsource・依存環境との一致が必要です。
[環境と再開の詳細](docs/getting_started.md)を参照してください。

## 資料

- [開発・研究ロードマップ](ROADMAP.md)：コード整備の優先順位、完了条件、保留中の研究。
- [Getting started / API (English)](docs/getting_started.md)：環境、使用方法、物理規約。
- [物理・数値設計](docs/flux_insertion_design.md)：電荷・gauge・ポンプの定義と解釈条件。
- [研究方針の根拠](docs/research/research_direction_20260920.md)：競合する文献、既存の証拠と限界。
- [テストガイド](test/README.md)：目的に応じた検証の実行方法。

研究上の議論は日本語、公開API・docstring・保存データのキーは英語を基本とします。
