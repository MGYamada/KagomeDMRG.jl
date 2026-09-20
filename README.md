# KagomeDMRG.jl

Kagome Heisenberg 模型の U(1) cylinder DMRG と flux insertion による
スピンポンプ測定を開発するための研究用 Julia パッケージです。
主対象はスピン 1/2 の最近接反強磁性模型の `M/Msat = 1/9` です。

- [実装ロードマップとバックエンドの比較](ROADMAP.md)
- [1/9 plateau の研究方針・既知領域・次の限定計算](docs/research/p4_static_plateau_strategy.md)
- [flux insertion の物理・数値設計](docs/flux_insertion_design.md)
- [Getting started (English)](docs/getting_started.md)

次の優先課題は、最近接等方模型 `J1=1,J2=J3=0` の静的な磁化境界と競合秩序です。
文献の `h/J ≈ 0.35–0.42` を出発点に、固定Q間の交換エネルギーを比較します。
既知CSLの校正は並行して進め、1/9ポンプを解釈する前の条件とします。
この方針は研究計画であり、当リポジトリでbulk plateauを確認したという意味ではありません。
N18のQ=0,2,4,6を比較した範囲では、1/9 sectorの区間は
`0.2585201582<h/J<0.4869821958`です。Q6の追加で区間は変わりませんでした。
[追加比較の精度・出典・未探索範囲](docs/research/p4_static18_q6_validation.md)を記録しています。
保存済みN18Q2では、約`0.00186467503 J`上の中性二重項が円周運動量`±2π/3`を持ち、
密度・NN bondの両方に非零の遷移強度を示しました。
[有限系の構造診断と独立検証](docs/research/p4_static18_neutral_validation.md)として保存し、
bulk gapやVBCの同定とは区別しています。
N27ではQ1/Q3のχ256・累積8を固定し、Q5だけ累積10 sweepsへ進めると、有限trialの区間は
`0.2586414244<h/J<0.4837321046`となりました。Q5追加2 sweepsで幅は約0.50%狭まり、
分布変化は小さくなったものの、精度5条件は引き続き未達です。
[固定χの追加比較と検証範囲](docs/research/p4_q5_csl_fixed_chi_followup.md)を参照してください。

## 現在の実装

格子・有向 bond と winding、固定磁化 sector、seam/uniform gauge、
ITensors.jl + ITensorMPS.jl による複素 U(1) 二サイト DMRG を実装しています。
整数電荷 `Q=2Sz_total` を明示でき、省略時は従来どおり 1/9 磁化です。
例えば `run_dmrg(kagome_cylinder(1, 4), 0.37; Q=0)` は12サイトの零磁化を選びます。
checkpoint と flux 継続は保存状態の Q を引き継ぎ、指定値との不一致を拒否します。
NN模型のN12Q0・θ=0,0.37とN9Q=−1,1,3の5点を独立EDと照合しました。
精度・時間・未検証範囲は[明示Qの研究記録](docs/research/p3_explicit_charge_validation.md)を参照してください。
別名の `kagome_j1j2j3_cylinder` で六角形内のJ2/J3を追加し、`bond_energies` と
結合種類の保存にも対応しました。`J1=1,J2=J3=.5` のN12Q0・零/非零fluxが独立EDと一致しました。
[拡張模型の検証記録](docs/research/p3_extended_model_validation.md)に数値と範囲を記載しています。
N18Q0のθ=0,.37でも独立EDとの残差≤`9.60e-12`で一致し、bond energy・chirality・
実空間/Schmidt移送を照合しました。[N18の測定検証](docs/research/p3_extended18_validation.md)。
この小系の移送は数値精度内で零であり、非零CSLポンプの検証は未完了です。
局所 Sz、相関、実空間の移送量、各 sweep の実測切断誤差を取得できます。
独立したスピン基底 ED と小系を照合しています。切断修正・保存履歴の修正時に
行った独立校正と、実行範囲・中断を含む詳細は
[検証記録](docs/research/p1_qn_truncation_calibration.md#テスト実行の範囲)を参照してください。
Julia 1.13.0 の記録した依存関係での結果であり、互換範囲の全版を検証したものではありません。
過去の [P0 小系照合](docs/research/p0_reference_validation.md)と
[保存・診断の検証](docs/research/p1_restart_schmidt_validation.md)は、当時の環境の記録として保持します。
縮退・近接した QN 切断境界の不整合に対し、局所版 NDTensors `0.4.31+1` を固定し、
実際の保持状態と報告 spectrum を同じ選択規則に揃えました。
空状態と報告 rank の不一致を検出する guard も維持しています。
[切断の独立校正と有限 χ 比較](docs/research/p1_qn_truncation_calibration.md)に、
修正範囲、残る精度制限、保存する実行履歴の修正を記録します。

完了した DMRG 点の checkpoint 保存・再開と、Schmidt 電荷・エンタングルメント診断も
実装しています。[保存・診断 API](docs/getting_started.md#checkpoint-restart-and-schmidt-diagnostics)を参照してください。
診断値による flux 点の受理・棄却、刻み半減、checkpoint 復元を行う
`continue_flux` を追加しました。零応答対照と独立 ED による初期検証の範囲は
[P2 研究記録](docs/research/p2_continuation_validation.md)に記載します。
18 サイトの最近接相互作用系では、切断なしの χ=512 で独立 ED と照合し、
刻み・初期状態・正負 flux を変えた `0→±0.37→0` の四経路を完走しました。
χ=128 の初期点棄却と χ=256 の切断 guard 停止も
[相互作用系の検証記録](docs/research/p2_interacting18_validation.md)に残しています。
修正版では停止を解消し、18 サイトの χ=128,256,512、θ=0,0.37 を計 10 点比較しました。
χ=512 は ED と一致し、χ=128,256 は sweep 数を増やしても精度不足が残りました。
既知 CSL のポンプ再現は未検証です。
単一 flux 点での最適化や移送量の読み出しだけでは、量子化や相の同定は
主張しません。SU(3)₁、Hall 応答を持つ／持たない D(Z₃) は候補として扱います。

## 最小実行例

Juliaup を使い、リポジトリのルートで Julia 1.13 の環境を準備します。
`override` はこのディレクトリ以下で起動する Julia に適用されます。

```sh
juliaup add 1.13
juliaup override set 1.13
julia +1.13 --project=research --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia +1.13 --project=research --startup-file=no -e 'using Pkg; Pkg.test("KagomeDMRG")'
julia +1.13 --project=research --startup-file=no --threads=1
```

研究用環境は `research/Project.toml` です。`research/Manifest-v1.13.toml` は Julia 1.13、
`research/Manifest.toml` は Julia 1.12.7 用で、両者は本checkoutと局所版 NDTensors を
相対pathで指定します。研究スクリプトは `--project=research` で実行してください。
過去の検証を再現する場合は、記録に対応する過去の checkout と依存環境を使います。
`Project.toml` の互換範囲は Julia 1.13 も含みます。

研究用環境では `Pkg.test("KagomeDMRG")` で全回帰テストを実行します。分野指定には
`Pkg.test("KagomeDMRG"; test_args=["truncation", "checkpoint"])` を使います。
ライブラリ開発ではrootの `Project.toml` も使えます。
`julia --project=. --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.test()'`
は互換範囲からローカル環境を解決します。rootに生成する `Manifest*.toml` はGit管理外で、
固定した研究用依存関係は `research/` のmanifestを参照します。
[テストの構成と実行方法](test/README.md)を参照してください。

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
julia +1.13 --project=research --startup-file=no examples/validate_small_system.jl
```

別 Julia process での再開と独立 ED・Schmidt 診断の照合は次で実行できます。
受理済み／試行 snapshot を分離し、元の状態を上書きしません。
同一 Julia・依存関係・ソースの信頼済みローカル計算用で、途中 sweep の再開ではありません。

```sh
julia +1.13 --project=research --startup-file=no --threads=1 examples/validate_restart.jl
```

適応的追跡の零応答・往復・再開、および縮退した初期状態による追跡停止は次で再現できます。
採否の閾値は実験ごとに明示し、`completed` はその診断を通過したという意味です。

```sh
julia +1.13 --project=research --startup-file=no --threads=1 examples/validate_continuation.jl
```

18 サイトの相互作用系で、bond dimension・刻み・初期状態・正負 flux と独立 ED を
比較する研究用スクリプトもあります。磁化は 1/9 の固定 sector であり、
plateau の成立や量子化の検証ではありません。

```sh
julia +1.13 --project=research --startup-file=no --threads=1 examples/validate_interacting18.jl
```

既知 CSL の正応答対照については、[原論文・補足の監査と実装条件](docs/research/p3_csl_control_design.md)を参照してください。

N27の隣接磁化sector、9/27サイト周期の初期試行状態、CSL対照の状態準備には、
configでケースを固定して一つずつ実行する `examples/run_research.py` を用意しました。
各走行は起動・保存・診断込み最大600秒、Julia/BLAS各1 threadです。
保存したtrialの数値整合性と基底状態への収束を区別します。
[ケースと事前判定基準](docs/research/p4_p3_staged_campaign.md)、
[周期初期状態の定義・校正](docs/research/order_seed_design.md)を参照してください。
以前のcheckpointを再開する際は、そのsource・依存関係・solver設定との厳密な一致が必要です。
時間切れの記録では外側の `execution.toml` と完了したbatchを併読します。

## 言語とリリース方針

現時点の議論・研究計画は日本語で進めます。公開 API、docstring、機械可読な結果の
キーは英語とし、英語の導入文書を維持します。英語リリースまでに研究文書の英訳と
リリース用の再現性・検証確認を進めます。
