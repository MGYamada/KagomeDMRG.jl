# P3a: 拡張模型N18Q0の独立ED・測定照合

日付: 2026-09-19。状態: **θ=0,0.37の2点で全照合基準を通過**。
対象は `J1=1,J2=J3=0.5,Q=0` の測定校正であり、最近接1/9模型とは別に扱う。
[N12の拡張交換検証](p3_extended_model_validation.md)と
[chirality演算子の校正](p3_chirality_validation.md)に続く研究単位である。
Hamiltonian・DMRG・切断・continuation・checkpointのproduction実装は変更しない。

## 問いと事前に固定した条件

独立なspin-basis EDとDMRGの状態が一致する範囲で、拡張交換、局所Sz、bond energy、
方向付きchirality、実空間移送とSchmidt電荷を同時に照合する。
`(Lx,Ly)=(2,3),N=18`、軸OBC、wrap `3a2`、既存の全単位胞を残す端と
`x,y,A/B/C` 順を使う。Q0の基底数は48,620、J1/J2/J3のbond数は30/24/12、
CCW基本三角形は9個。軸cutはsite9と10の間の一つで、内部bulk列はない。

模型は[Gong–Zhu–Shengの原著](https://www.nature.com/articles/srep06317)で用いられた
六角形内の拡張交換を参考にする。今回の18サイトは原著の288サイトのポンプ走行と
異なる有限円筒であり、原著の端・入力を厳密に再現したとはしない。
原著との対応と不足情報は[CSL対照設計](p3_csl_control_design.md)を参照する。

| 項目 | 固定値 |
| --- | --- |
| flux・電荷 | seam gauge、θ=0,0.37 radを順番に計算、Q=0、hz=0 |
| DMRG | seed11、初期linkdim4、χ512、6 sweeps、cutoff=noise=0 |
| 局所Krylov | tol=`1e-11`、krylovdim40、maxiter20 |
| ED | 独立疎行列、BlockLanczos、nev=blocksize=4、seed5678、krylovdim80、maxiter200、tol=`1e-11` |
| ground cluster | 最低値から`1e-8 J`以内。clusterより高い要求準位も存在すること |
| 再試行 | 比較基準未達のθに限り、同じ状態から追加6 sweepsを最大1回 |
| 資源 | 数値process一つ、Julia/BLAS各1 thread、block-sparse threadingなし |
| 時間 | 起動・compile・診断・保存を含め600秒。既存の外側監視が単一solve中も停止 |

N18のχ512は中央二分割の全Hilbert空間を含める上限であり、大系の有限χ収束設定ではない。
θ=.37は比較基準を満たしたθ=0のMPSと同じsite indicesから始める。
全保存状態はtrialで、`continue_flux`の受理軌道としては登録しない。
期待する移送値やchiralityの符号を採否に使わない。

事前基準はenergy/site差≤`1e-8 J`、独立残差・ground cluster外norm≤`2e-6`、
射影して対応付けた同じED状態とのSz・bond energy・chirality差≤`1e-6`。
EDは全行列残差≤`1e-11 J`と直交性を別計算する。
有限個の低準位からground spaceの完全性やbulk gapは証明しない。

同じDMRG状態に対する独立chiralityとゲージ間の差、Schmidt/密度・左右移送の
整合性には`1e-10`、EDとの移送差には`1e-6`、Q/norm誤差には`1e-12`を要求する。
varianceは独立残差の二乗との差≤`1e-9 J²`を確認し、
`100eps(Float64)*max(1,E²)`を差し引きの丸めscaleとして保存する。
設定cutoffや局所solver toleranceを誤差推定値とは呼ばない。

## 独立参照と保存

[研究driver](../../examples/validate_extended18.jl)は、productionのbond/MPOを使わずに
平面の距離・chordless六角形から構築する既存の
`reference_extended_hamiltonian(...;sparse_matrix=true)`と
`_reference_low_eigensystem`を組み合わせる。NN専用の`reference_eigensystem`は呼ばない。

新しい[chirality参照](../../examples/extended18_chirality_reference.jl)は
Cartesian Pauli行列とLevi-Civitaから8×8の三spin行列を作り、固定Q基底へ埋め込む。
productionのladder展開・MPO構築を参照しない。
seamでは`alpha=-image_y*theta`による`U C U†`を使う。
DMRG状態そのものと射影ED状態の双方で比較し、uniform gaugeではMPSと演算子を共に回す。
三角形の幾何・CCW方向は前回の独立校正を用いる。

解析対照`±sqrt(3)/4`、向き反転、一般角のゲージ回転、入力正規化、
N18Q0へ埋め込んだ継ぎ目上の三spin状態を確認する。
実模型のchiralityが零でも、独立演算子検証はこれら非零の対照で符号を固定する。

```sh
python3 examples/run_extended18.py outputs/p3-extended18-NEW --wall-seconds 600
```

既存の時間監視・資源計測を再利用し、研究driver/参照helperのhashを明示的に追加保存する。
実行途中でsourceが変わった場合は成功にしない。ED固有ベクトルを原子的に保存し、
trial checkpointの状態と零flux baselineは読込み後にも照合する。
過去のN27/N18記録・checkpointのsource identityは変更しない。

## 実測結果

両点とも初回6 sweepsで通過し、追加sweepは不要だった。最終保持次元は512、
全12 sweepsの実測切断誤差は0。設定・全profile・180個のproduction/manifest hash、
9個の研究driver/参照hashは[数値記録](data/p3_extended18_validation.toml)、
外側の終了・資源量は[実行記録](data/p3_extended18_execution.toml)を正本とする。
ローカルの固有ベクトルとtrial checkpointは`outputs/p3-extended18-20260919/`にある。

| θ | ED最低E/J | energy/site差 | 独立残差/J | ground cluster外norm | 最大Sz差 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | −8.041823140097150 | 7.89e−16 | 9.60e−12 | 5.89e−12 | 5.25e−13 |
| 0.37 | −8.038044018455132 | 9.87e−17 | 4.96e−13 | 3.01e−13 | 1.72e−14 |

要求した4準位では両点のground clusterは1状態で、次準位との差はそれぞれ
`0.1540075466 J`と`0.1531875824 J`だった。全4準位のED残差は最大`2.74e−12 J`。
これは固定Qの有限クラスター準位間隔であり、bulkの中性gapではない。
別のblock seedによる低準位の完全性の追加監査は行っていない。

| θ | 最大bond energy差/J | 同じMPSの独立chiralityとの差 | 射影EDのchiralityとの差 | seam/uniformのchirality差 |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 1.35e−12 | 1.10e−17 | 3.23e−14 | 0 |
| 0.37 | 3.16e−14 | 3.57e−17 | 1.20e−14 | 2.38e−17 |

9三角形のchiralityの最大絶対値はθ=0で`2.97e−14`、θ=.37で`2.09e−15`。
非零の合成対照では`±sqrt(3)/4`が再現され、N18Q0の継ぎ目の対照誤差は
最大`1.11e−16`だった。実模型で零だった値を、非零chiralityの観測としては記述しない。

θ=.37の零flux baselineからの右移送は`+2.87777e−13`、左移送は`−2.88674e−13`。
左右和の絶対値は`8.97e−16`、実空間/Schmidt移送差は`8.42e−16`、
射影EDとの移送差は`2.92e−13`以下だった。基準profileの再読込み差は両点とも0。
これらは有限区間の零応答であり、2πのポンプや受理済みcontinuation軌道ではない。

raw varianceは順に`+9.95e−14 J²`、`−1.42e−14 J²`だった。
独立残差の二乗は約`9.21e−23 J²`、`2.46e−25 J²`で、
差し引きの丸めscale約`1.44e−12 J²`より十分小さい。
負のraw varianceは丸め範囲内にあり、絶対値を取り直した誤差推定には使わない。

全工程は**149.397秒**、worker CPUは147.414秒、peak RSSは3,379,904,512 bytes
（約3.148 GiB）。Julia 1.13.0、ITensors 0.9.31、ITensorMPS 0.4.1、
NDTensors 0.4.31+1、KrylovKit 0.10.4で実行した。
peak RSSはOSのchild high-waterであり、全process treeの同時peakとは異なる。

| 工程 | θ=0 / 秒 | θ=.37 / 秒 |
| --- | ---: | ---: |
| 独立ED | 4.655 | 3.572 |
| DMRG（energy・variance・Szを含む） | 66.431 | 23.485 |
| 診断全体 | 33.331 | 1.174 |
| うち9三角形のchirality | 1.262 | 0.115 |

最初のED開始まで約13.6秒を要した。初回の各時間には読込・JITが含まれ、
特に最初の診断は呼出側のcompileも含むため、個別計測の和とは一致しない。
起動・compileを外側149.397秒から除外しておらず、JIT時間だけの厳密な分離は未実施。
次点のchirality測定費用をN36や原著規模へそのまま外挿しない。

研究driver・物理規約は別担当が独立レビューし、完了後の数値・閾値・保存hashも
別担当と照合した。完了後の監査は保存値の再集計・hash・コードの確認であり、
追加のMPS収縮やED対角化を再実行したものではない。
新しい参照の解析selfcheckと実模型2点の検証を実行した。
productionと既存testを変更していないため、既通過のpackage suiteは再実行していない。
既存の736件等を今回の新しい全体テスト結果として数え直さない。

## 解釈と次の研究単位

零磁場でchirality seedを含まないspin-twist Hamiltonianは、物理的時間反転Tで不変である。
`T i T⁻¹=-i` と `T S± T⁻¹=-S∓` により、twistを含む二つの横交換項が入れ替わる。
局所z回転UはTと可換なので、dressed chiralityもTで符号を反転する。
したがって、Q0の非縮退な有限系基底状態では局所Szとdressed chiralityが零となる。
これは本模型の式からの推論であり、小系の零移送をCSLの否定や正応答校正の完了とは扱わない。

CSL側の後続候補はN36 `(3,4),Q=0` の零flux状態準備と資源測定である。
複素枝・seed・χを判断に必要な点だけ比較し、既知非零ポンプは別単位とする。
実CSLの軌道には、現行continuationの密度・Schmidt等に加え、
chiralityとbond profileの記録を組み込む必要がある。

最近接1/9側はN27追加収束の保留を維持できる。
次の小さい静的候補はN18・θ=0・Q=6の1点であり、
既存のQ=0,2,4からQ=6への飛び越しが上側磁場境界を縮めるかを調べる。
`h26=(E6-E2)/2`を比較し、Q≥8の未探索と内部bulkがない制限は残す。
このQ6計算とN36状態準備は本記録では未実行であり、校正と一括して完了扱いにしない。
