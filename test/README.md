# テスト

全回帰テスト、分野指定、一覧表示は次のコマンドで実行する。

```sh
julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
julia --project=. --startup-file=no -e 'using Pkg; Pkg.test(; test_args=["truncation", "checkpoint"])'
julia --project=. --startup-file=no test/runtests.jl --help
```

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
| `provenance` | 実行コードとの不一致拒否、precompile 更新 |

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
