# テスト

固定した研究用環境から、全回帰テスト、分野指定、一覧表示を次のコマンドで実行する。

```sh
julia --project=research --startup-file=no -e 'using Pkg; Pkg.test("KagomeDMRG")'
julia --project=research --startup-file=no -e 'using Pkg; Pkg.test("KagomeDMRG"; test_args=["truncation", "checkpoint"])'
julia --project=research --startup-file=no test/runtests.jl --help
```

`research/Manifest.toml` はJulia 1.12用、`research/Manifest-v1.13.toml` は1.13用で、
`research/Project.toml` が本checkoutとvendored NDTensorsを参照する。
ライブラリ開発のroot環境では従来どおり
`julia --project=. --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.test()'`
を使える。そこで生成するroot manifestはGit管理外であり、研究用の固定環境とは区別する。
`test/Project.toml`はテスト依存とtestからの相対sourceを明記し、どちらの起動環境でも
vendored NDTensorsを同じ場所から読む。
以下の過去の実行時間・assertion数は当時の環境の記録で、今回の配置変更の検証結果ではない。

[AGENTS.md](../AGENTS.md#practical-test-time-budget) の時間予算に従い、日常の
重点確認は約60秒、全体は約3分を目安とする。5分を超える場合は高コストの原因を
調べる。統合変更の全体確認は1回とし、その後の修正では影響した分野だけを再実行する。
起動・precompile が支配的な場合は数値計算の時間と分けて記録する。

| 分野 | 主な検証 |
| --- | --- |
| `lattice` | 幾何・電荷・winding・局所／円周 flux、J1/J2/J3族 |
| `ed` | 独立 spin-basis Hamiltonian、平面六角形の拡張結合、ゲージ同値、疎固有値と縮退 |
| `dmrg` | MPO/ED、残差・密度・相関、解析的切断対照 |
| `observables` | 累積移送、bond energy、Schmidt電荷、三spin chiralityと周期画像・ゲージ |
| `checkpoint` | 原子的保存、Hamiltonian 更新、破損・不整合の拒否 |
| `truncation` | QN 保持 rank、実損失、cutoff・noise・停止 guard |
| `continuation` | 逐次更新・逆走、棄却・復元、各停止分岐 |
| `provenance` | source・active環境の変更拒否、root Manifestなし、cached moduleでの環境再捕捉 |
| `diagnostics` | 静的診断の直接／別process一致、独立複素状態参照、欠測・失敗・精度未達、入力保存不変性 |

テストの拡張モードは廃止し、重複した組合せ自体を削除した。
QN 校正は17例、低水準の切断選択は6例とし、格子・Schmidt・flux追跡も
異なる境界条件や不具合を検出する代表例へ絞った。数値許容差は変更していない。

独立 ED と dense 参照は production の MPO・切断選択を呼ばない。
保存ファイルと固定場対照は `checkpoint_helpers.jl`、解析的少数spin対照は
`truncation_helpers.jl` にまとめる。continuation の主経路は実 DMRG を使い、
制御分岐の対照は既知の固定場固有状態を使う。

18サイトの χ 比較などの研究計算は `examples/` から実行する。
過去の計算・テストの記録は `docs/research/` に当時の範囲として保持する。

重複削減時には1,211件が4分13.0秒で成功した。その後の明示Q対応では、既存の
DMRG・固有値計算を非既定Qへ差し替え、整数電荷・保存契約の安価な検証を追加した。
数値最適化の本数は増やしていない。

明示Q対応の統合実行は1,263件が4分17.3秒で成功した。拡張交換模型では、独立幾何・
不等なJ2/J3の小行列MPO照合・結合エネルギーとfamily保存の検査を追加し、
通常suiteのDMRG・固有値計算本数は維持した。

2026-09-19、Julia 1.13.0、Julia/BLAS 各1 threadで拡張模型の統合 `Pkg.test()` を
1回実行し、**1,322 / 1,322件が成功**した。Test計測は4分17.7秒、
コマンド全体は266.78秒、package precompileは約3秒だった。
3分の目安は未達、5分の見直し基準内である。重点確認を日常の基本とし、
研究driverの点比較はこのsuiteへ追加しない。

同日のDMRG進行callback追加後は `dmrg checkpoint truncation` の重点検証を実行し、
**532 / 532件が成功**した。Test計測は1分39.3秒、package precompileは約3秒。
callback有無での同一初期状態の結果一致、phase/bond/sweep通知、variance省略、例外の伝播、
callbackを保存しないことを追加確認した。60秒の目安は超えたため、ここから通常suite全体を
重複実行せず、影響した分野の検証を終えて研究driverへ進んだ。
上記1,322件は変更前の全体結果として保持する。

同日のchirality追加後は `lattice observables checkpoint` を1回実行し、
**736 / 736件が成功**した。Test計測94.0秒、外側Julia起動を除く
`Pkg.test` 呼出103.058秒、package precompile約3秒（内数）。
独立Pauli三spin行列、周期三角形の距離参照、N18の合成chiral状態で符号・
ゲージ・正規化・入力保持を検証した。[数値とsource hash](../docs/research/p3_chirality_validation.md)。
初回MPO/MPS・Schmidt等のJITを含む3群の実行で60秒目安を超えたため、
全体suiteは重複実行していない。chirality単独の測定費用は未分離である。

## Manifest分離後の確認

研究環境から `Pkg.test("KagomeDMRG"; allow_reresolve=false)` を実行し、
Julia 1.13.0・Julia/BLAS各1 threadで **1,470 / 1,470件が成功**した。
Test計測は4分17.4秒、`Pkg.test`呼出全体は277.019秒、事前precompileは約14秒（内数）。
root Manifestがないpackage copy、active Project/Manifest変更、同名ファイルのhash区別、
既存cacheを必須とする別processでの環境再捕捉、旧環境checkpointの拒否を含む。

Julia 1.12.7では依存再解決、import、active環境hashと実行例の構文を確認した。
1.12のprecompileは約19秒で、数値suite全体を同版で再実行したとは扱わない。
両研究manifestのpackage versionとtree hashは配置変更前と同一で、local pathのみ移動した。
4 Python launcherはASTとmock起動で研究環境・cwd・出力先・1 thread・時間上限を確認した。
研究scanと過去の数値結果の再計算は行っていない。

## 共通静的診断追加後の確認（2026-09-23）

研究環境から `Pkg.test("KagomeDMRG"; allow_reresolve=false)` を1回実行し、
Julia 1.13.0・Julia/BLAS各1 threadで **1,855 / 1,855件が成功**した。
Test計測は5分0.3秒、`Pkg.test`呼出は305.561秒（外側Julia起動を除く）。
別processの保存後測定は、実DMRGのN9、非零flux・uniform gauge・Q=-1・磁場付き複素状態、
幾何Schmidt cutを持つLx2解析状態を一つのchild processにまとめている。
独立spin-basis照合、参照reportの由来とhash、未実施・失敗・精度未達、元ファイル不変性を含む。

3分目標は未達で、5分の見直し基準にも達した。別process検証には新しいJuliaでの
測定コードのコンパイルが必要であり、別途実行した小系driverではsolve/diagnose呼出の
98.67%／99.39%をコンパイルが占めた（suite全体の内訳を同率と推定するものではない）。
初版の重点197件は96.4秒だったが、参照由来・3ケース保存後照合を含む最終版の重点群だけの
時間は未分離。全体再実行は重複せず、今後は影響した群を選ぶ。起動・JIT費用の削減は次の
検証環境整備に残し、独立参照・strict再開・許容差は維持する。
Julia 1.12で今回の追加機能を検証したとは扱わない。
[実行例・数値・保存契約](../docs/static_diagnostics.md)を参照する。
