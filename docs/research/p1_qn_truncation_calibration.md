# P1: QN 切断の修正と有限 χ 校正

実施日: 2026-09-19。対象は CPU の複素 U(1) 二サイト DMRG。
最近接 spin-1/2 kagome Heisenberg 模型の `N=18,Q=2,M/Msat=1/9` を
小系参照とする。量子化ポンプ、plateau、相同定はこの検証の対象ではない。

## 原因と限定した修正

従来の NDTensors 0.4.31 では、全 sector の重みから決めた保持 rank を、
各 block で scalar threshold に変換して再判定していた。
閾値付近の同値・近接した重みで二つの判定が食い違い、実際の保持次元が
報告 spectrum より小さくなる。閾値の微小補正を外すだけでは完全縮退を処理できない。

`vendor/NDTensors` に元の registry tree
`311282faab99a123706a738fcceb9caba81225ae` をライセンスとともに保持し、
局所版 `0.4.31+1` として固定した。修正は block-sparse SVD と Hermitian eigen の
保持 rank 選択であり、dense 分解、DMRG sweep、Hamiltonian と電荷規約は変更しない。
これはこのリポジトリの研究用修正であり、上流で公開された release ではない。
実装の参照は [上流 block 分解](https://github.com/ITensor/ITensors.jl/blob/main/NDTensors/src/blocksparse/linearalgebra.jl)、
公開契約は [ITensor 分解 API](https://docs.itensor.org/ITensors/stable/ITensorType.html)。

全 block の `(重み, block, block 内 index)` をまとめ、一度だけ選択する。
既定の `min_blockdim=0` は上流の cutoff・mindim・maxdim と誤差の正規化規則を使う。
完全に等しい重みは block 座標、block 内 index の順に決める。
`maxdim` は上限であり、縮退空間の全状態を保持するとは限らない。
正の `min_blockdim` では sector ごとの最低保持数も同じ上限に算入し、
最低数だけで上限を超える設定は例外とする。
返す spectrum と誤差は実際に選ばれた状態に対応する。
空状態・非有限 norm・報告 rank と保持次元の不一致を検出する guard は維持する。

## 18 サイトの失敗入力

元のコード revision `59c7466` を隔離した出力ディレクトリへ復元し、
元の manifest と NDTensors 0.4.31 で `χ=256, seed=11` を再実行した。
`examples/capture_qn_boundary.jl` は専用プロセスで `replacebond!` の入力・出力を
観測するだけであり、元の分解は変更しない。

| 零 flux、sweep 4、左からの bond 9 更新 | 元の値 |
| --- | ---: |
| 報告保持 rank | 256 |
| 実際の保持 rank | 255 |
| 報告した捨てた確率 | 4.311879671383606e-8 |
| 保持状態への射影から測った捨てた確率 | 4.485943627674516e-8 |
| 独立 dense SVD による最良 rank 256 の損失 | 4.311879671383731e-8 |
| 境界の第 256 重み | 1.7406398874382037e-9 |
| 境界の第 257 重み | 1.7406398867654784e-9 |

これは四サイトの近接重み例と同じ原因であり、元の報告値は実損失を約 3.9% 過小評価した。
保存した入力の SHA-256 は
`f108c66616189a1ca441b66d2725186edfe5df2f987e4731442cf27cf62cd01d`。
元環境の source/manifest hash、設定、入力 checksum は
[捕捉記録](data/p1_qn_upstream_capture.toml)に残す。
この失敗を新しい成功結果で上書きしない。

同じ入力を局所版の SVD/eigen、それぞれ左右の factorize に渡すと、
四通りすべてで実保持 rank と報告 rank が 256 になった。

| 分解・方向 | 報告した損失 | 実保持状態の射影損失 |
| --- | ---: | ---: |
| SVD・left | 4.311879671383606e-8 | 4.311879671383600e-8 |
| SVD・right | 4.311879671383606e-8 | 4.311879671383602e-8 |
| eigen・left | 4.311879668873603e-8 | 4.311879671401754e-8 |
| eigen・right | 4.311879670622746e-8 | 4.311879671406266e-8 |

dense 最良 rank 256 との損失差と報告／実測差は、事前の絶対許容差 `1e-12` 内。
loss は差分 norm から計算しており、二つのほぼ等しい norm の差し引きを避ける。

## 保存する実行履歴の修正

レビューで、計算後に source を編集してから保存すると、実行したコードではなく
保存時の disk 上のコードの hash が記録される問題を確認した。
KagomeDMRG と vendored backend のロード時に immutable な source identity を確定し、
計算結果と baseline に持たせる。実行開始・終了、保存、読み込み、継続で
disk とロード済みコードを照合し、不一致なら Julia の再起動を要求する。
backend を先にロードしてからその source が変わる場合も照合する。
source・manifest の precompile 依存関係を登録し、新プロセスで古い hash が残ることを防ぐ。

チェック対象は明示した本体 source、vendored backend の Project/src/ext、
Julia minor 版に対応する manifest、および記録対象の runtime/package version である。
任意の外部環境・全メソッドの動的書き換えを証明する機構ではない。
旧 source の checkpoint を新 source として読み込む移行は行わず、schema 1 の
source/runtime 一致条件を保つ。今回の失敗入力 replay は、checksum を確認した
限定的な分解入力の検証であり、旧 checkpoint の再開ではない。

実行履歴を大きな NamedTuple のまま各結果に含めると、checkpoint の統合検証で
LLVM の aggregate 展開に長時間を要した。最終実装は不変な非 parametric struct に
hash tuple を保持し、結果の型へ全 hash の個数を展開しない。hash 内容・照合条件・
保存 schema は同じで、実行履歴の値そのものを可変にしたわけではない。
変更後の atomic checkpoint 検証は **109/109 assertions** が成功し、
単独 testset は 98.2 秒で完了した。変更前の統合実行中の同 testset は 13 分 51.7 秒を
要したが、別の実行条件・並行負荷も含むため厳密な速度比の benchmark とはしない。

## 独立な校正と数値比較

`test/test_qn_calibration.jl` は二・四 spin の係数行列を独立に構成し、
LinearAlgebra の dense SVD/eigen と照合する。完全縮退・近接重み・複素 block 内回転、
左右分割、cutoff の相対／絶対規則、入力 norm の変更を含む。
noise ありでは既知の正半定値摂動を密度行列に加え、その trace 損失と
元の波動関数の射影損失を別々に測る。noise を加えた spectrum の誤差を
元の波動関数の捨てた確率とは呼ばない。
この独立校正は **3,783 assertions** が成功した。
vendor 側の保持数・制約・実 SVD/eigen 回帰は **346 assertions** が成功した。

18 サイト比較の設定は結果を見る前に固定した。
幾何は `(Lx,Ly)=(2,3)`、wrap は `3a2`、軸 OBC、円周 PBC、
既定の cell/sublattice 順序、`Jxy=Jz=1,hz=0,Q=2`、seam gauge。
各 case の θ は `0→0.37` の逐次計算とし、その case 自身の零 flux 状態を warm start と
密度 baseline に使う。これは受理済み continuation 経路ではなく、すべて trial 保存とする。

固定 case は `(χ,sweeps)=(128,6),(128,12),(256,6),(256,12),(512,6)`。
共通条件は `seed=11,cutoff=0,noise=0,eigsolve_tol=1e-11,
eigsolve_krylovdim=40,eigsolve_maxiter=20`、Julia/BLAS 各 1 thread、variance 計測あり。
中央 bond 9 の左右すべての更新で、dense 参照と実保持状態の射影損失を再測定する。
全点を独立 spin-basis ED と比較し、energy/site `1e-8`、residual `2e-6`、
局所 Sz `1e-6`、基底空間からの漏れ norm `2e-5` を事前の比較閾値とする。
この閾値を超えた結果も保存し、切断処理の修正と物理状態の収束を分ける。

### 二つの実行記録

最初の実行では `(128,6),(128,12),(256,6)` の 6 点が完了した。
統合テストのコンパイル遅延を修正するため次の case の途中で中断し、
[中断時の記録](data/p1_finite_chi_interrupted.toml)と source 一式を保存した。
中断は数値判定の不合格ではなく、未完了 case は `interrupted` として残している。

ユーザーの検証圧縮方針に従い、完了した 6 点は再計算しない。
残る `(256,12),(512,6)` は新規の零 flux 状態から独立に計算する。
二つの実行で production source の差は `src/provenance.jl` の履歴表現だけであり、
全数値 source、backend、manifest、独立 ED helper の hash は一致することを照合した。
driver には固定 protocol から case を選ぶ引数だけを追加し、各記録に選択 case を保存する。
旧 checkpoint の再開条件を緩める処理は行っていない。
残る 4 点も完了し、[完了側の記録](data/p1_finite_chi_completion.toml)に保存した。
合計 10 点の trial checkpoint、payload checksum、各 source identity を再照合した。

二つの実行の結果は、それぞれの source identity とともに読む。
最初の結果を最終コードの hash で付け替えない。
バイナリ入力・MPS・中断時の source snapshot は `outputs/` に保持し、
数値表と hash を含む TOML はこの研究記録とともに version 管理する。

### 18 サイトの結果

`ΔE/N` は独立 ED からのエネルギー差の絶対値を site 数で割った値、leakage は
ED 基底空間への射影からの距離 norm。局所 Sz はその射影を正規化した状態と比べる。
最後の列は最終 sweep の全 bond・両方向にわたる最大実測切断誤差である。

| χ | sweeps | θ | ΔE/N | 全系 residual | max ΔSz | leakage | 最終切断誤差 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 128 | 6 | 0 | 6.017e-6 | 2.469e-2 | 3.346e-3 | 1.314e-2 | 7.838e-6 |
| 128 | 6 | 0.37 | 6.931e-6 | 2.653e-2 | 3.226e-3 | 1.509e-2 | 8.433e-6 |
| 128 | 12 | 0 | 5.980e-6 | 2.467e-2 | 8.691e-4 | 5.866e-3 | 7.591e-6 |
| 128 | 12 | 0.37 | 6.933e-6 | 2.654e-2 | 4.229e-3 | 1.941e-2 | 8.436e-6 |
| 256 | 6 | 0 | 1.771e-8 | 1.566e-3 | 3.986e-7 | 2.118e-4 | 4.254e-8 |
| 256 | 6 | 0.37 | 2.348e-8 | 1.808e-3 | 6.877e-7 | 2.415e-4 | 5.591e-8 |
| 256 | 12 | 0 | 1.771e-8 | 1.566e-3 | 3.986e-7 | 2.118e-4 | 4.254e-8 |
| 256 | 12 | 0.37 | 2.348e-8 | 1.808e-3 | 6.877e-7 | 2.415e-4 | 5.591e-8 |
| 512 | 6 | 0 | 4.441e-16 | 7.174e-12 | 1.304e-12 | 6.349e-12 | 0 |
| 512 | 6 | 0.37 | 2.467e-16 | 8.612e-13 | 8.670e-12 | 4.605e-11 | 0 |

χ=512 の二点は四つの事前 ED 比較閾値を満たす。最大 link rank は 386、
すべての sweep の切断誤差は 0 であり、この固定電荷小系の参照では切断していない。
χ=128,256 の八点は精度基準を満たさない。χ=256 は局所 Sz の基準には入るが、
energy・residual・基底空間からの漏れが残る。
6→12 sweep で χ=256 の残差はほぼ変わらず、χ=128 の密度には変化が残った。
この比較では sweep 数だけを増やしても誤差は解消せず、χ を増やす効果が大きい。
単一 seed の結果であり、有限 χ の大域的な変分最小値を保証するものではない。

中央 bond の **168 更新**で、実保持 rank と報告 rank が一致した。
報告損失と実測射影損失の最大差は `5.30e-21`、dense 最良 rank の損失との最大差は
`4.49e-20`。この照合は中央 bond に対して行い、他の bond は既存 guard と
sweep 最大誤差で監視する。noise ありの校正は小さい解析対照の範囲である。

全 10 点の全電荷誤差は `2.67e-15` 以下、左右移送の不一致は `3.22e-15` 以下、
実空間と Schmidt 移送の差は `1.88e-15` 以下だった。
variance と独立 residual² の差は最大 `1.40e-13`。
χ=512 の variance はこの丸め誤差の床にあり、小さな符号付き値をそのまま保存した。

θ=`0→0.37` の右側移送は χ=512 で `−2.656502326040776e-4`、
χ=256 では約 `−2.656821835e-4`。χ=128 は 6 sweep で `−2.621283033e-4`、
12 sweep で `−2.015338406e-4` だった。これは短い有限区間の密度変化であり、
2π のポンプや量子化の測定ではない。

### テスト実行の範囲

コードレビュー時の変更前 `Pkg.test()` は 4,348 assertions が成功した。
今回の変更後の全体実行では、格子・ED/MPO・DMRG・観測量・Schmidt・checkpoint・
切断校正まで **8,059 assertions** が成功した後、continuation のコンパイル中に
明示的に中断した。全体実行の最終 status は成功ではなく中断である。

ユーザーの方針により、以後は通過済みの組合せを繰り返さず、履歴表現変更の影響を受ける
保存・再開・source 変更検出を重点確認する。切断の数値 source はその後変更していない。
全体 suite の新しい成功件数として 8,059 を扱わない。
最終表現では checkpoint **109 件**に加え、隔離した N=9 の零結合・固定場対照で
source 変更拒否、正常再開・再保存、短い `continue_flux` 再開の **20 件**が成功した。
後者は変更した模型による保存契約の検証であり、最近接模型の研究結果には含めない。
Julia 1.13.0 で実行した。Julia 1.12.7 は manifest の再解決までで、
今回の追加検証の成功を同版に拡張しない。

## 再現と研究上の範囲

通常の回帰はリポジトリの project 環境で実行する。

```sh
julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
julia --project=. --startup-file=no --threads=1 examples/validate_finite_chi.jl outputs/finite-chi-new
```

固定 protocol のうち未完了だった二つだけを実行するコマンドは次の通り。

```sh
julia --project=. --startup-file=no --threads=1 examples/validate_finite_chi.jl outputs/finite-chi-completion-new outputs/qn-boundary-upstream-oys7v4ta/capture 256:12,512:6
```

最後の引数を省略すると五つの case すべてを新規計算する。

同じ捕捉入力を比較する場合は、第二引数に信頼できる capture ディレクトリを渡す。
元の失敗を新規捕捉するには、旧 revision の source/manifest を別ディレクトリへ復元し、
その project を有効にして現行 `examples/capture_qn_boundary.jl` を実行する。
いずれも新規出力先を要求し、失敗・未収束結果を保持する。

この系には軸方向 cut が一つだけで内部 bulk column がない。
有限 χ の residual、密度誤差、移送の依存性を調べるための小系であり、
Hall 応答・非零 quantized pump・topological branch の検証にはならない。
既知 CSL、複数 cut の相互作用系、長さ・幅・初期状態の検証は次の段階とする。
GPU、generic 非 Hermitian 分解、大規模 DMRG の性能は検証していない。
