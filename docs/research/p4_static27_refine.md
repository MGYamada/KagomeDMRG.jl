# P4段階B: 保存したN27状態から2 sweepsを追加する

2026-09-19。[前回の2-sweep状態](p4_static27_progress.md)を起点に、
同じχ上限で最適化を進めたときのenergy・variance・密度・bond profileと費用を比較する。
以下は実行前に固定した範囲である。

**結果: 累積4 sweepsまで189.979秒で完了。実保持次元256の新trialと全診断を保存した。**
2-sweep状態からEは0.474536 J下がり、varianceは0.315463→0.00849796 J²へ低下した。
列磁化と結合パターンも大きく変わり、前回の初期状態の偏りは最適化に対して安定していない。
3→4 sweepにも0.040673 Jのenergy変化が残るため、基底状態収束は未確立である。

## 問いと固定条件

最近接・等方Heisenberg模型の `(Lx,Ly)=(3,3),N=27,Q=3,θ=0`。
wrap `3a2`、軸OBC・円周PBC、順序 `x,y,A/B/C`、完全な単位胞を残す端、48 bonds。
`Jxy=Jz=1,J2=J3=hz=0`、seam gauge、総Sz=1.5。
元の試行は複素ランダムMPS・seed11から2 sweepsを行い、実保持最大次元64、
`E/J=−10.825816813361072`、variance `0.31546287414180085 J²` だった。

今回は[元outcome](data/p4_static27_progress_outcome.toml)で指定されたtrialを検証して読み、
同じsites/psi0/Qで **2 sweepsを追加、累積2→4** とする。
maxdim256、cutoff=noise=0、Krylov tolerance `1e-11`、dimension40、最大反復20を維持する。
seed11を設定に残すが、初期状態を新たに乱数生成する走行ではない。
最適化呼出しではvarianceを省略し、完了状態の保存後に別計時する。

```sh
python3 examples/run_static27_refine.py outputs/p4-static27-refine-new --wall-seconds 600
```

[launcher](../../examples/run_static27_refine.py)は既存の上限停止とCPU/RSS計測を再利用する。
起動・親checkpoint検証・JIT・診断・保存・最大2秒の停止猶予を含めwall600秒、
Julia/BLAS各1thread、数値process一つ。追加のχ・seed・Qや自動追加sweepは行わない。

## 再開と比較の契約

[worker](../../examples/validate_static27_refine.jl)は親outcome・validationのSHAと、
checkpointのmetadata/payload/checksumsを照合する。
`load_checkpoint(...;status=:trial)` によりsource/runtime・模型・Q・site identity・
θ・solver設定・energy・Szを検証してから、`run_dmrg` へsites/psi0/Qを明示する。
productionは変更せず、元trialの受理ラベルやvariance未測定設定を書き換えない。

各bond/sweep通知に今回のsweep番号と累積番号を残す。
完了した新状態を、自身を零flux baselineとする別のtrialとして先に保存・再読込する。
親状態は比較の由来として保存し、最適化による密度変化をflux pumpと解釈しない。
前回と同じSchmidt・相関・bond energy・varianceの式と読出し許容差を用い、
各phaseで得られた値を保存してから整合性を判定する。

比較量は全energy、variance、局所Szと列別Sz、bond energy、Schmidt量、
親状態との正規化overlap、実保持次元と切断誤差、工程時間・CPU・peak RSS。
varianceの低下や保持次元256への到達は期待値で合格判定せず、実測事項とする。
4 sweepsや小さいenergy変化だけで基底収束とは判定しない。
Q3のみでは磁場区間を求めず、三列からbulk秩序・ギャップ・相を同定しない。

## 検証と結果

実行前に親原記録・production179ファイル・旧追加4ファイルのSHA、checkpointの
checksum・設定・模型・runtimeを独立担当と主担当が確認した。
新driverはJuliaで読込み2.18秒、親の実MPS再読込・検証24.99秒で通過した。
親の188ファイルの不変性、E/J=−10.825816813361072、最大保持次元64も一致した。
この事前確認では新たなDMRG最適化は行っていない。

launcherはJuliaを起動しないmockでworker指定、単一呼出し、絶対出力先、600秒・各thread1、
不正wall値と既存出力への上書き拒否を確認した（約0.034秒）。
独立担当が再開・保存・比較のコードと本プロトコルを読み、実行前の必須修正はなかった。
production・通常テスト・依存は前回から変更していないため、直前の重点532件の結果を維持し、
package suiteは再実行しない。研究走行の観測整合性と保存は別に検証する。

## 累積2→4 sweepsの比較

実行先は `outputs/p4-static27-refine-20260919`。正本の
[validation](data/p4_static27_refine_validation.toml)、
[execution](data/p4_static27_refine_execution.toml)、
[progress](data/p4_static27_refine_progress.toml)をbyte単位でコピーし、
主要値・SHA・次の判断を[outcome](data/p4_static27_refine_outcome.toml)に保存した。
workerは `completed_refinement_pilot_accuracy_unestablished`、launcherは終了確認付きexit 0。
元snapshotを保持し、新しい `trial/checkpoint-iYdlHH` を保存・再読込した。

| 観測量 | 累積2 sweeps | 累積4 sweeps |
| --- | ---: | ---: |
| 全交換energy E/J | −10.825816813361072 | −11.300352744875930 |
| E/(NJ) | −0.4009561782726323 | −0.4185315831435530 |
| variance/J² | 0.31546287414180085 | 0.008497959029512003 |
| 最大保持次元 | 64 | 256 |
| 最終sweepの最大実測切断誤差 | 0 | 1.359629631238424e−5 |
| cut9のSchmidt entropy | 0.3693990511 | 1.4848921292 |
| cut18のSchmidt entropy | 1.5473701032 | 1.6301858196 |

varianceは前回の0.02694倍（約37.1分の1）になったが、基底energyの誤差が37.1分の1に
なったという意味ではない。varianceの丸め床は `2.84e-12 J²`、エネルギー標準偏差は
0.0921844 Jであり、固有状態精度や基底空間への近さを独立EDで保証した値でもない。

| 今回のsweep | 累積sweep | 末尾の局所solver energy/J | 最大保持次元 | 最大実測切断誤差 |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 3 | −11.259680103462461 | 256 | 0 |
| 2 | 4 | −11.300352744875937 | 256 | 1.359629631238424e−5 |

ここで初めて設定上限256に達し、累積4 sweep目には上限での切断が実測された。
局所solverの末尾energyと、上表の再計算した全MPS期待値は区別して保存している。
3→4の末尾energy差は−0.0406726414 J。追加sweepによる変動と有限χの影響はまだ分離できない。
切断誤差を局所Szや磁場境界の誤差限界には用いない。

| 列別Sz | 累積2 sweeps | 累積4 sweeps | 変化 |
| --- | ---: | ---: | ---: |
| 第1列 | 1.4268392033 | 0.4732156572 | −0.9536235460 |
| 第2列 | −0.0360801767 | 0.3985618761 | +0.4346420528 |
| 第3列 | 0.1092409734 | 0.6282224667 | +0.5189814933 |

総Sz=1.5を保ちながら、前回の片側への強い偏りが大きく変わった。
最大の局所Sz変化は0.387862、最大bond energy変化は0.387083 J。
親状態との正規化overlapは0.066867で、2-sweep状態から波動関数も大きく変化している。
これはθ=0の最適化の比較であり、fluxによる移送・枝の断絶・相転移の証拠ではない。
累積4 sweepsにも不均一性は残るが、端局在・VBC・磁気秩序を同定しない。

## 費用と保存・観測整合性

| 項目 | 実測 |
| --- | ---: |
| 起動込みwall time | 189.979006秒 |
| 子のCPU時間 | 187.993604秒 |
| 子のpeak RSS | 2,486,943,744 bytes = 2.316147 GiB |
| 親trialの検証・再読込 | 23.5054秒 |
| 追加2 sweepsとenergy/Sz読出し | 142.1317秒 |
| 新checkpoint保存・再読込 | 1.1026秒 |
| Schmidt量 | 0.6068秒 |
| 相関とbond profile | 2.4533秒 |
| 別工程variance・状態照合 | 10.7050秒 |
| うちH†H収縮 | 10.6149秒 |
| callback処理合計 | 0.1341秒 |
| 新checkpoint payload | 5,810,386 bytes |

時刻はJIT・検証・callback費用を含む。親再読込と新checkpoint工程は、純粋なI/Oではなく
energy・Sz等の再照合も含む。新trialの最大保持次元256に対してvarianceは約10.6秒で測れた。
今回の費用では追加の最適化が大きいが、記録のない初回8-sweep試行の停止工程は特定しない。

総磁化誤差は `5.78e-15` 以下、checkpointのoverlap誤差は `3.56e-15` 以下。
bond energy和の差は `1.07e-14 J` 以下、固定Qのzz恒等式・Schmidt分散等も許容差を通過した。
別variance計測の前後でnorm=1、保存状態への正規化overlapは `1.0000000000000036` で、
1をわずかに超える分は丸め誤差の範囲である。Qとsite indicesも一致した。
これらは保存・読出しの整合性検査であり、基底状態収束を保証しない。
元trialのvariance未測定設定は維持し、後段の値は研究診断として別記録に置く。

production179ファイル・今回の追加6ファイルは走行時のSHAと一致し、
元outcome・validation・checkpointと旧sourceの不変性も終了時に確認した。
独立担当も新旧checkpoint・同一site identity・模型・runtime・settingsを照合し、
116通知／104 bond更新、Schmidt確率と電荷、相関の恒等式、48 bondの和、variance、
前回との全比較差を再集計して不一致がなかった。
overlapはworkerが計算した保存値の整合確認であり、独立担当はbinary MPSを再計算していない。
研究サイズのDMRG・EDも追加実行していない。

## 当時提案した次の限定研究単位

次はこの4-sweep trialを同一source/runtime・模型・Q・site indicesで検証してから、
同じmaxdim256でさらに2 sweeps（累積6）を追加する。起動込み600秒、各thread1の別予算とする。
3→4にenergy変化が残り、密度も大きく変わったため、まず固定χでのsweep依存を調べる。
その変動が小さくなった段階でχ512比較を選び、有限χの影響を分離する。
この次走行は未実施であり、今回の予算を延長した追加計算ではない。
Q1/Q5の比較、磁場区間、別seedや競合秩序は、代表点の精度と費用を踏まえて別途選ぶ。

**同日追記:** この追加収束は保留した。保存済みtrialを維持し、
ロードマップの独立した[P3a chirality実装・校正](p3_chirality_validation.md)を完了した。
次の作業順は[現行ロードマップ](../../ROADMAP.md)を参照する。
source変更前の179ファイルと6個のdriverはhashを照合して退避し、旧trialを新コードの
結果として読み替えていない。この追記は上記の実測値・当時の計算条件を変更しない。
後に分離したコミット `48d15e2b566b4d4b50a6c5295b84348cc5f30f8d` のGit treeも、
この179個のproduction/manifestと6個のdriverについて全SHAが一致することを確認した。
