# P2 初期検証: 適応的 flux 追跡と未解決な枝の保存

実行日: 2026-09-19。Julia 1.13.0、ITensors 0.9.31、ITensorMPS 0.4.1、
NDTensors 0.4.31。Julia と BLAS はそれぞれ 1 thread とした。

今回の goal は、受理・棄却・刻み半減・checkpoint 復元を行う driver を実装し、
独立 ED と幾何学的 cut を持つ零応答対照で検証することだった。
`FluxPolicy` と `continue_flux` を追加した。ここで示す結果は測定・追跡基盤の
初期検証であり、最近接模型の量子化ポンプや相同定の結果ではない。

## 実装と採否

- θ は unwrapped に保持し、一つの trajectory の点は逐次計算する。
  各試行の直前に最後の受理済み snapshot を読み直し、同じ site indices と
  solver 設定で新しい Hamiltonian を構築する。棄却状態を次の初期値に使わない。
- 各試行で overlap 振幅、指定したサイトの Sz、Schmidt entropy・左 Sz 平均の変化、
  全系 variance、最後の sweep の実測切断誤差、最終二 sweep の energy 差と
  最終期待値との差を確認する。初期状態にも品質判定を適用する。
- 累積移送は最初の零 flux MPS の実測密度を基準にする。
  左右の和、幾何学 cut `c` に対応する bond `3Ly*c` の Schmidt 電荷変化、
  指定した cut 間の右移送の差を記録し、明示された有限の閾値と照合する。
- 閾値は実験ごとに指定する。期待する 2/3 や 0 は採否条件に含めない。
  刻みは棄却時に半減し、下限では停止する。受理後に自動では増やさない。
  逆方向の目標列と、非零 θ の受理 snapshot からの再開を扱う。
- 有効な棄却状態は `trial/` に保存し、受理点は別の `accepted/` に公開する。
  最小刻み・試行予算・不正な診断は `unresolved` として journal に残す。
  `completed` は指定した有限の診断を通過した意味で、断熱性の証明ではない。

現 driver は seam gauge、noise=0、全点での variance、二 sweep 以上に限定する。
θ に依存する uniform-gauge 変換を含んだままの raw overlap は使わない。
variance の負値は `100eps(Float64)*max(1,E²)` の丸め範囲と明示された閾値の
両方を満たす場合だけ認める。全 sweep の energy と切断誤差も保存する。
公開 API と実行例は [Getting started](../getting_started.md#adaptive-flux-continuation) を参照。

## 保存した数値結果

再現コマンド:

```sh
julia --project=. --startup-file=no --threads=1 examples/validate_continuation.jl
```

今回の出力は `outputs/p2-continuation-yxp5sJ/validation.toml` と、同じディレクトリ内の
各 trajectory journal／snapshot。集約記録を
[p2_continuation_validation.toml](data/p2_continuation_validation.toml) に同一内容で保存した。
モデルの全 bond、幾何、電荷、solver、policy、観測領域、全試行、依存バージョン、
git revision と選択したソース／manifest の SHA-256 を含む。
実行前後のソース hash は一致し、集約の検証状態は `passed`。
局所のネットワーク設定や環境変数の dump は保存していない。

### 27 サイト: 解析的零応答と往復・再開

`Lx=3,Ly=3,N=27,Q=3`、wrap `3a2`、軸 OBC・円周 PBC、順序は `x,y,A/B/C`。
対照模型は `Jxy=Jz=0` とし、最初の 15 サイトに `hz=+1`、残りの 12 サイトに
`hz=-1` を加えた `H=-Σ hz_i Sz_i`。一意な積状態 `Up^15 Dn^12` が
エネルギー `-13.5` の基底状態で、Hamiltonian は flux に依存しない。
これは最近接 Heisenberg 模型とは異なる模型である。

最初に 0→2π を計算し、保存点から別 journal で 4π→6π→0 を追跡した。
全受理列は `0,2π,4π,6π,4π,2π,0`。step は 2π、下限 π/32、
二 sweep、maxdim=2、cutoff=0。二つの幾何学 cut は bond 9 と 18、
監視する中央サイトは 10:18 とした。overlap 下限は `1-1e-10`。

| 指標 | 実測 |
| --- | ---: |
| エネルギーの解析値からの最大誤差 | 0 |
| 左右の実空間移送の最大絶対値 | 0 |
| Schmidt 電荷による累積移送の最大絶対値 | 0 |
| 二つの cut の右移送の差 | 0 |
| 隣接受理点の overlap 振幅 | 1 |

絶対左 Sz 平均はそれぞれ 4.5 と 6、entropy は 0。
大きい θ 刻みが許されるのはこの対照模型が flux に依存しないためであり、
相互作用系に同じ刻みを推奨するものではない。

### 9 サイト: 独立 ED と追跡 driver の照合

最近接等方的 `Jxy=Jz=1,hz=0`、`Lx=1,Ly=3,N=9,Q=1`、seed=11。
0→0.17→0.37 を計算し、距離探索から独立に構成した spin-basis Hamiltonian と比較した。
零 flux の二重縮退は全基底空間への射影で扱う。

| 指標 | 最大値 |
| --- | ---: |
| エネルギー誤差 / site | `1.49e-16` 未満 |
| ED residual norm | `1.06e-13` 未満 |
| 基底空間外の振幅 norm | `3.62e-13` 未満 |
| 局所 Sz 誤差 | `1.83e-13` 未満 |

この走行は driver と ED の整合を調べるため、overlap 下限を明示的に 0 とした。
最初の overlap は約 0.71376、次は約 0.99668 であり、連続な枝の認定には用いない。
`Lx=1` には軸方向の内部 cut がなく、bond 4 は補助的な MPS cut である。
この結果から axial pump は読み出せない。

### 縮退した始点: 刻み半減で解消しない棄却

同じ 9 サイト模型で、零 flux の二次元基底空間内の実数 ED 固有ベクトルを
初期 MPS として零 flux DMRG に与えた。overlap 振幅の下限を 0.95、
最初の刻み 0.01、下限 0.0025 とすると次の結果になった。

| 試行 θ | 零 flux 状態との overlap 振幅 | 判定 |
| --- | ---: | --- |
| 0.01 | 0.7070306006 | `low_overlap` |
| 0.005 | 0.7070698490 | `low_overlap` |
| 0.0025 | 0.7070886048 | `low_overlap` |

受理列は `[0]` のまま `unresolved / min_step` となり、三つの棄却 MPS を保存した。
これらも ED の基底エネルギーと一致し、記録した点の最大 residual は `2.17e-14` 未満。
低エネルギー・小さい残差でも所与の始点からの連続性を満たすとは限らない。
独立レビューでも θ を `1e-2` から `1e-5` に減らした際、選んだ零 flux ベクトルとの
**二乗** overlap は約 0.5 に、零 flux 全基底空間への重みは 1 に近づいた。
[独立 ED 記録](data/p2_degeneracy_ed.toml)を保存した。
これは有限系の縮退に伴う枝選択の問題であり、トポロジカル縮退の証拠ではない。

## 縮退 QN 切断の独立監査

NDTensors 0.4.31 の installed source と独立した二・四サイト状態を調べた。
`src/truncate.jl` は切断境界の重みが相対 `1e-3` 未満で近接すると block 選択閾値を
上げるが、`src/blocksparse/linearalgebra.jl` が返す spectrum と切断誤差は元の
rank に対応したままになる。SVD と eigen の両方で次を再現した。

| Schmidt 確率 | 要求 maxdim | 実保持 rank | 報告損失 | 実際の損失 |
| --- | ---: | ---: | ---: | ---: |
| 0.6, 0.2, 0.2 | 2 | 1 | 0.2 | 0.4 |
| 0.6, 0.2001, 0.1999 | 2 | 1 | 0.1999 | 0.4 |
| 0.6, 0.21, 0.19 | 2 | 2 | 0.19 | 0.19 |
| 0.6, 0.2, 0.2 | 3 | 3 | 0 | 0 |

正規化後も状態は非零なので norm 検査だけでは検出できない。
`length(eigs(spec)) == dim(linkind(psi,bond))` を毎更新で照合する guard を追加した。
実際の `replacebond!` と独立な state overlap に基づく損失を比較する回帰テストを含む。
正常な二・四サイト DMRG、noise=0 と `1e-5` の計 32 更新では次元が一致した。
[監査数値記録](data/p2_qn_truncation_audit.toml)にバージョンと上流ソース hash を保存した。
分解 API の参照先は [ITensor 公式文書](https://docs.itensor.org/ITensors/stable/ITensorType.html#Decompositions)。
上の機構説明は installed source と実験による確認である。

この guard は既知の不整合を検出するもので、上流の切断処理を修正してはいない。
次元が一致することは一般的な誤差保証ではない。失敗が続く場合は縮退群を保持する
maxdim の増加や根本修正が必要で、θ 刻み半減だけでは解決しない。

## 検証・残課題

独立担当による focused tests は continuation 417 assertions、切断 guard 84 assertions が成功。
最終的に次の標準 package check も完了し、全 **4,283 assertions** が成功した
（既存 3,782 + continuation 417 + 切断 guard 84）。今回の拡張を含む全 suite は
Julia 1.12.7 では再実行していない。

```sh
julia --project=. --startup-file=no --threads=1 -e 'using Pkg; Pkg.test()'
```

状態の破損を伴う数値例外、意図的な品質棄却、最小刻み、予算切れ、初期品質、
負の variance、非有限 variance、再開時の基準保持を含む。
棄却後に以前の snapshot が不変であり、新しい試行が正規化済みの受理状態から始まることを検査した。

独立レビューで、NaN variance による保存拒否が journal を `running` のままにする
不具合を再現し修正した。修正後は `unresolved / trial_validation_error`、試行は
`invalid_trial` となり、例外の型だけを保存する。再検証は 9 assertions 成功。
再現スクリプトの初回には集約 TOML の辞書型の不備があり、数値 journal を保持したまま
書出し処理を修正して新しい出力先で再実行した。上の記録は修正後の成功走行である。

P2 の全研究検証や P3 は未完了である。今後は相互作用する長い cylinder で、
刻み半減・bond dimension・長さ・幅・初期状態依存を比較する。
既知 CSL の正応答対照は原論文の結合と幾何を照合して実装する必要がある。
局所 Krylov の収束フラグ、一般的な切断誤差保証、途中 sweep 再開、
自動的な solver 精度増加、uniform gauge の共通規約による overlap、
大規模 I/O の性能評価は今回の実装範囲に含めない。
SU(3)₁、Hall-active / Hall-inactive D(Z₃)、VBC、磁気秩序、gapless な状態は
引き続き候補であり、量子化や相の識別を主張しない。
