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
Julia 1.13.0 でも同じ主要パッケージ版で全 2,205 件が通過しました。
checkpoint と Schmidt 診断の追加時は全 3,782 件が通過しました。
適応的 continuation と切断 guard の追加後は、Julia 1.13.0 で全 4,283 件が通過しています。
独立 sparse ED の縮退・収束検査を追加後は、全 4,348 件が通過しました。
これらは記録した依存関係での検証であり、互換範囲の全バージョンを検証したものではありません。
[検証記録](docs/research/p0_reference_validation.md)に数値と限界を記載します。
[保存・診断の検証](docs/research/p1_restart_schmidt_validation.md)も別途記録しています。
縮退・近接した QN 切断境界に upstream の制限があり、空状態と
報告スペクトル／保持次元の不整合を検出して停止します。
この検出は上流カーネルの修正ではなく、一般的な切断誤差の校正は引き続き必要です。

完了した DMRG 点の checkpoint 保存・再開と、Schmidt 電荷・エンタングルメント診断も
実装しています。[保存・診断 API](docs/getting_started.md#checkpoint-restart-and-schmidt-diagnostics)を参照してください。
診断値による flux 点の受理・棄却、刻み半減、checkpoint 復元を行う
`continue_flux` を追加しました。零応答対照と独立 ED による初期検証の範囲は
[P2 研究記録](docs/research/p2_continuation_validation.md)に記載します。
18 サイトの最近接相互作用系では、切断なしの χ=512 で独立 ED と照合し、
刻み・初期状態・正負 flux を変えた `0→±0.37→0` の四経路を完走しました。
χ=128 の初期点棄却と χ=256 の切断 guard 停止も
[相互作用系の検証記録](docs/research/p2_interacting18_validation.md)に残しています。
既知 CSL のポンプ検証は未実装です。
単一 flux 点での最適化や移送量の読み出しだけでは、量子化や相の同定は
主張しません。SU(3)₁、Hall 応答を持つ／持たない D(Z₃) は候補として扱います。

## 最小実行例

Juliaup を使い、リポジトリのルートで Julia 1.13 の環境を準備します。
`override` はこのディレクトリ以下で起動する Julia に適用されます。

```sh
juliaup add 1.13
juliaup override set 1.13
julia +1.13 --project=. --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia +1.13 --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
julia +1.13 --project=. --startup-file=no --threads=1
```

Julia 1.13 用の依存関係は `Manifest-v1.13.toml` に保存しています。
既存の `Manifest.toml` は Julia 1.12.7 の検証環境として保持します。
`Project.toml` の互換範囲は Julia 1.13 も含みます。

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
julia +1.13 --project=. --startup-file=no examples/validate_small_system.jl
```

別 Julia process での再開と独立 ED・Schmidt 診断の照合は次で実行できます。
受理済み／試行 snapshot を分離し、元の状態を上書きしません。
同一 Julia・依存関係・ソースの信頼済みローカル計算用で、途中 sweep の再開ではありません。

```sh
julia +1.13 --project=. --startup-file=no --threads=1 examples/validate_restart.jl
```

適応的追跡の零応答・往復・再開、および縮退した初期状態による追跡停止は次で再現できます。
採否の閾値は実験ごとに明示し、`completed` はその診断を通過したという意味です。

```sh
julia +1.13 --project=. --startup-file=no --threads=1 examples/validate_continuation.jl
```

18 サイトの相互作用系で、bond dimension・刻み・初期状態・正負 flux と独立 ED を
比較する研究用スクリプトもあります。磁化は 1/9 の固定 sector であり、
plateau の成立や量子化の検証ではありません。

```sh
julia +1.13 --project=. --startup-file=no --threads=1 examples/validate_interacting18.jl
```

既知 CSL の正応答対照については、[原論文・補足の監査と実装条件](docs/research/p3_csl_control_design.md)を参照してください。

## 言語とリリース方針

現時点の議論・研究計画は日本語で進めます。公開 API、docstring、機械可読な結果の
キーは英語とし、英語の導入文書を維持します。英語リリースまでに研究文書の英訳と
リリース用の再現性・検証確認を進めます。
