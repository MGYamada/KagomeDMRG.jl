# KagomeDMRG.jl

スピン1/2の最近接kagome反強磁性Heisenberg模型を、U(1)円筒DMRGで調べる
研究用Juliaパッケージです。主対象は磁化 `M/Msat = 1/9` です。

**研究の中心は、競合する秩序・液体候補を、初期状態・端・有限幅の影響から切り分けることです。**
次は長い同一円筒で、中央の結合模様・磁化・相関の準備依存と境界感度を比較します。
数値精度や計算資源は、その物理的な違いを判別するために配分します。
具体的な次の研究と判断基準は [ROADMAP](ROADMAP.md) にまとめています。

## 実装と現在の限界

- ITensors.jl + ITensorMPS.jlによる複素U(1)二サイトDMRG、固定磁化sector、格子・flux gauge。
- 局所磁化、結合energy、相関、chirality、Schmidt診断と、実測切断誤差・分散。
- 完了点のcheckpoint保存・再開と、診断に応じたflux刻みの調整・状態復元。

小系の独立ED照合と保存・測定の検証を研究基盤としています。
長い円筒での競合状態比較はこれから整備する計画です。
bulk plateau・中性gap・微視的な相は未確立で、既知CSLの非零ポンプ校正も未完了です。
静的研究を先行し、1/9の量子化ポンプの解釈は校正後に進めます。

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

- [研究ロードマップ](ROADMAP.md)：優先する問い、次の比較、進め方の判断基準。
- [Getting started / API (English)](docs/getting_started.md)：環境、使用方法、物理規約。
- [物理・数値設計](docs/flux_insertion_design.md)：電荷・gauge・ポンプの定義と解釈条件。
- [研究方針の根拠](docs/research/research_direction_20260920.md)：競合する文献、既存の証拠と限界。
- [テストガイド](test/README.md)：目的に応じた検証の実行方法。

研究上の議論は日本語、公開API・docstring・保存データのキーは英語を基本とします。
