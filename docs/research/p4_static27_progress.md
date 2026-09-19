# P4段階B: N27の進行記録と診断費用の分離

2026-09-19。[初回8-sweep試行](p4_static27_pilot.md)は600秒上限内で停止し、
`run_dmrg` 内のどこまで進んだか不明だった。本記録では工程と局所更新の進行を残し、
完了状態の保存をvariance計測より先に行う。以下は実行前に固定した研究単位である。

**結果: 57.812秒で2 sweepsと全診断を完了し、初めてN27のtrial状態を保存した。**
設定上限はχ256だが、実際の最大bond次元は16→64。
エネルギー変化とvarianceが大きく、基底状態収束は未確立である。
段階Bの診断用完了点は得られたが、隣接Qを比較するための精度評価は継続する。

## 問いと実行範囲

同じ最近接・等方Heisenberg模型の `(Lx,Ly)=(3,3),N=27,Q=3,θ=0`、
wrap `3a2`、軸OBC・円周PBC、順序 `x,y,A/B/C`、完全単位胞の端、48 bonds。
`Jxy=Jz=1,J2=J3=hz=0`、seam gauge、総Sz=1.5。
χ256、seed11、初期link dimension 4、2 sweeps、cutoff=noise=0、
Krylov tolerance `1e-11`、dimension 40、最大反復20に固定する。
2 sweepsを収束条件とせず、新しい複素ランダムMPSから代表1点だけを計算する。

```sh
python3 examples/run_static27_progress.py outputs/p4-static27-progress-new --wall-seconds 600
```

[launcher](../../examples/run_static27_progress.py)は既存の停止監視・資源測定を再利用する。
起動・JIT・全計測・保存・最大2秒の停止猶予を含むwall上限600秒、
Julia/BLAS各1thread、数値process一つ。χ変更・追加sweep・別seed・隣接Qは実行しない。
今回のN27は独立EDを行わない。旧sourceのcheckpointからの再開もしない。

## 進行記録と状態保存

`run_dmrg(...; progress_callback=...)` に任意の通知を追加する。
callbackはscalar-onlyのNamedTupleを受け取り、MPS/MPOや可変のsolver設定を受け取らない。
既定値は `nothing`。solver設定・切断guard・Hamiltonian・保存schemaは変更しない。
callback例外はそのまま停止として伝播させる。

- 初期状態、MPO、sweep、全energy、variance（指定時のみ）、Szの開始・終了。
- 各bond更新のsweep番号・half-sweep・bond・局所energy・実測切断誤差・保持次元。
- 各sweepの完了と局所energy・当該sweepの最大切断誤差・最大bond次元。

[worker](../../examples/validate_static27_progress.jl)は各通知を `progress.toml` に原子的に保存する。
局所energyはupstream solverの値であり、全MPSの最終期待値とは区別する。
callbackの累積時間には通知の保存・標準出力が含まれ、全経過時間にも計上される。
工程の時刻には初回JITが含まれる。初回8-sweep試行の主要原因を今回の時間差だけで断定しない。

実行順は次のとおり。

1. `measure_variance=false` で2 sweepsを実行し、全energy・Sz・実測切断誤差を保存する。
2. norm・総磁化・診断の有限性を確認し、完了状態をtrial checkpointへ保存・再読込する。
   保存・読込はMPO再構築とenergy/Sz照合を含むため、純粋なI/O時間とは呼ばない。
3. cut 9/18のSchmidt量を測り、平均左Szを密度和と照合する。
4. zz/+-相関を一度だけ測り、θ=0のbond energyを導出する。
   bond和を全energy、Schmidt分散を左領域zz相関と照合する。
5. `real(inner(H,psi,H,psi))−E²` を別計時する。
   符号付きvariance、丸め床 `100eps(Float64)*max(1,E²)`、H†Hの虚部を保存し、
   前後のnorm・保存状態への正規化overlap・Q・site indicesを確認する。

norm/全磁化/Schmidt確率和の許容差は `1e-12`、密度・相関の恒等式は `1e-10`、
bond energy和・Schmidt分散は `1e-9`。これらは観測読出しの整合性であり、基底精度の基準ではない。
Schmidt平均はcut間の一致を要求しない。checkpointのvariance未測定ラベルは書き換えず、
後段のvarianceは別の研究診断として保存する。

## 判定範囲

今回の完了条件は、限定した試行を実施し、完了または停止の工程と得られた量を
独立監査付きで保存して次の研究単位を決めること。上限で終了したworkerの
`running` は最後のsnapshotにすぎず、実行状態はlauncherの記録と照合する。
2 sweepsが返っても追加sweep・χ・seed依存は未評価であり、基底状態の収束とは扱わない。
状態が得られなければ未取得量を明記し、段階Bの数値的完了とも扱わない。
Q3のみでは磁場区間を求められず、三列はbulkの保証にならない。
phase同定・Hall応答・競合秩序の優劣は今回の判定対象外である。

## 検証と結果

実行前に `Pkg.test(; test_args=["dmrg","checkpoint","truncation"])` を実行し、
532/532件が成功した。Test計測99.3秒、package precompile約3秒。
同じ初期MPSのcallback有無によるenergy・Sz・overlap・sweep診断の一致、phase順、
variance省略、例外時に偽の完了を出さないこと、callback非保存を確認した。
全package suiteは重複実行していない。最初のsandbox内起動はJulia launcherのlock作成で
失敗し、数値計算開始前だった。承認済み実行で上記重点検証を完了した。

研究workerの読込み・固定設定・列磁化・複素相関配列化を別に軽量確認し、
読込みは2.19秒だった。launcherはJuliaを起動しないmockで、worker指定、単一呼出し、
1thread設定、600秒上限、不正上限と既存出力への上書き拒否を確認した（約0.97秒）。
共通supervisorの停止・CPU/RSS換算は前回の確認を再利用した。

独立担当がcallback・保存・観測式・判定範囲をレビューした。
Schmidt診断の不合格値も保存されるよう、保存を判定より先へ移した。
その後、研究に使うsourceを固定して走行する。初回N27のdriver・launcherは変更せず、
前回の保存SHAと一致する。productionのcallback追加は新sourceとして記録し、
過去の状態やsource監査を現在の版として読み替えない。

## 実測結果

出力は `outputs/p4-static27-progress-20260919`。
[validation](data/p4_static27_progress_validation.toml)、
[execution](data/p4_static27_progress_execution.toml)、
[progress](data/p4_static27_progress_progress.toml)を元ファイルとbyte単位で一致する形で保存した。
SHAと主要結果・次の判断は[別のoutcome記録](data/p4_static27_progress_outcome.toml)にまとめる。
workerは `completed_diagnostic_pilot_accuracy_unestablished`、launcherは終了確認付きexit 0である。

| 量 | 実測 |
| --- | ---: |
| 起動込みwall time | 57.812239秒 |
| 子のCPU時間 | 56.919411秒 |
| 子のpeak RSS | 1,559,904,256 bytes = 1.452774 GiB |
| 全交換エネルギー E/J | −10.825816813361072 |
| E/(NJ) | −0.4009561782726323 |
| 総Sz | 1.5000000000000002 |
| variance/J² | 0.31546287414180085 |
| variance/(NJ²) | 0.011683810153400032 |
| エネルギー標準偏差/J | 0.5616608177021082 |
| checkpoint payload | 514,274 bytes |
| callback合計 | 0.120882秒 |

| sweep | 局所solverの末尾energy/J | 実際の最大bond次元 | 最大実測切断誤差 |
| --- | ---: | ---: | ---: |
| 1 | −10.001215900841197 | 16 | 0 |
| 2 | −10.825816813361083 | 64 | 0 |

2回のsweep末energyは約0.824601 J変化している。
**切断誤差0は収束を意味しない。** 今回は設定上限256へ達する前の初期最適化であり、
2 sweepsの間に到達した変分空間が十分かどうかを示すものではない。
varianceは丸め床 `2.60e-12 J²` を大幅に上回る。
独立ED比較・追加sweep・χ/seed比較は行っておらず、基底energyの誤差保証はない。

| 工程 | wall秒 |
| --- | ---: |
| 2 sweeps、準備、最終energy/Sz | 43.9005 |
| checkpoint保存・再読込・整合性照合 | 0.7983 |
| Schmidt cut 9/18 | 0.3104 |
| zz/+-相関とbond profile | 2.2894 |
| 別工程のvarianceと状態照合 | 2.5309 |

variance工程内のH†H収縮自体は2.4830秒だった。
104回のbond更新と2回のsweep完了を保存した。
`run_dmrg` 内の進行時計では初期化19.320秒、MPO 2.879秒、sweep工程13.553秒、
全energy 1.415秒、Sz 2.225秒である。初回JITやcallback処理を含む値で、純粋なkernel時間ではない。
外側の工程時間との残差は呼出し前のJIT・dispatch等も含み、資源の二重計上はしない。
今回の低い保持次元でvarianceが短かったことから、**旧8-sweep試行の停止原因は断定できない**。
最大保持次元256へ達した状態での費用は未測定である。

列別Szは、左から `1.4268392033, −0.0360801767, 0.1092409734`。
中央列はsite 10:18だが、この強い非一様性は未収束の初期状態の観測であり、
端励起・VBC・磁気秩序の確定した性質として扱わない。
cut 9/18のSchmidt entropyはそれぞれ0.36939905、1.54737010（自然対数）。
左右で同じ値になる条件は課していない。

## 保存・整合性と次の研究単位

trial checkpoint `trial/checkpoint-fG0rfi` の保存・再読込が成功し、
元の状態とのoverlap誤差は `1.45e-15` 以下。
checkpointはvariance未測定のまま保持し、別診断の前後もnorm=1、
保存状態への正規化overlapは `0.9999999999999986` だった。
Q・site indices、総磁化、相関の恒等式とSchmidt量も整合した。
bond energy和の誤差は `1.07e-14 J` 以下、Schmidt分散の照合誤差は `1.78e-15` 以下。
これらは読出し・保存の検証であって、基底状態収束の判定ではない。

保存後にproduction 179ファイルと追加4ファイルのSHA、rawコピー、
progressと最終sweep診断、checkpointのpayload・metadataのSHAを確認した。
独立担当もsource・checkpoint・設定を照合し、116通知の順序、104 bond更新、
sweepごとの最大切断誤差を再構成した。保存相関から48 bondのenergy、固定Q恒等式、
Schmidt確率から平均・分散・entropy、H†Hからvarianceを再集計し、不一致はなかった。
この独立監査は保存数値の再集計であり、binary MPSの再読込・DMRG再計算は追加していない。

次はこのtrialを**同じsource・模型・Q・site indicesで検証してから**再開し、
同じmaxdim256で2 sweepsを追加する別の限定研究単位を選ぶ。
合計4 sweepsへ進め、保持次元の成長、切断誤差、全energy、variance、列磁化の変動と費用を見る。
既定でacceptedを読む `resume_dmrg` は用いず、`load_checkpoint(...;status=:trial)` 後に
同じsites/psi0/Qを明示した `run_dmrg` を使う。元snapshotのSHAと累積sweep数2→4を保存し、
元trialをacceptedへ付け替えない。保持次元256への到達も次回に観測する事項である。
これは今回実行したものではなく、次回も起動前に600秒上限を設定する計画である。
固定したχの上限に達した段階の観測を得るまで、χ512・別seed・隣接Qへは拡大しない。
未収束が続けば追加sweepとχの影響を分けて判断し、期待する秩序や磁場幅に合わせない。

後続の[累積2→4 sweeps比較](p4_static27_refine.md)を実施し、実保持次元256の新trialを
取得した。本記録の2-sweep状態と診断は変更せず、比較の起点として保持する。
