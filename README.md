# KagomeDMRG.jl

Kagome Heisenberg 模型の U(1) cylinder DMRG と flux insertion による
スピンポンプ測定を開発するための研究用 Julia パッケージです。
主対象はスピン 1/2 の最近接反強磁性模型の `M/Msat = 1/9` です。

- [実装ロードマップとバックエンドの比較](ROADMAP.md)
- [flux insertion の物理・数値設計](docs/flux_insertion_design.md)
- [Getting started (English)](docs/getting_started.md)

## 現在の実装

格子・有向 bond と winding、固定磁化 sector、seam/uniform gauge、
ITensors.jl + ITensorMPS.jl による複素 U(1) 二サイト DMRG を実装しています。
局所 Sz、相関、実空間の移送量、各 sweep の実測切断誤差を取得できます。
独立したスピン基底 ED との小系照合を含む 2,205 件のテストが、Julia 1.12.7、
ITensors 0.9.31、ITensorMPS 0.4.1、KrylovKit 0.10.4 で通過しました。
これは現在の依存関係での検証であり、互換範囲の全バージョンを検証したものではありません。
[検証記録](docs/research/p0_reference_validation.md)に数値と限界を記載します。
縮退した QN 切断境界に upstream の制限があり、空状態を検出して停止します。
一般的な縮退境界の切断誤差は引き続き校正が必要です。

checkpoint、枝の連続性を診断する flux continuation、既知 CSL のポンプ検証は
未実装です。単一 flux 点での最適化や移送量の読み出しだけでは、量子化や相の同定は
主張しません。SU(3)₁、Hall 応答を持つ／持たない D(Z₃) は候補として扱います。

## 最小実行例

リポジトリのルートで依存関係を準備し、テストを実行します。

```sh
julia --project=. --startup-file=no -e 'using Pkg; Pkg.instantiate()'
julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
julia --project=. --startup-file=no --threads=1
```

Julia のプロンプトで、9 サイトの小系を計算します。`theta` はラジアンです。

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

この系は `Q=2Sz_total=1`、`Sz_total=1/2` です。
`max_truncation_errors` は各 sweep で測定した最大値であり、設定値 `cutoff` とは
区別します。9 サイト例は実装検証用で、バルクのポンプ測定用ではありません。
独立 ED と照合し、選択したメタデータと数値を保存する再現用スクリプトもあります。

```sh
julia --project=. --startup-file=no examples/validate_small_system.jl
```

## 言語とリリース方針

現時点の議論・研究計画は日本語で進めます。公開 API、docstring、機械可読な結果の
キーは英語とし、英語の導入文書を維持します。英語リリースまでに研究文書の英訳と
リリース用の再現性・検証確認を進めます。
