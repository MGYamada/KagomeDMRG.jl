# N27 Q5とCSL零flux準備の固定χ追加sweep比較

実施日: 2026-09-20。状態: **NNは全診断完了・精度未達。CSLは2 sweepsと状態保存後、分散測定中に上限停止**。

後続の[同一親χ比較・CSL分散補足](p4_q5_matched_chi_comparison.md)で、このCSL保存trialの
分散0.24237056 J²を追加DMRGなしで取得した。本記録の未完了状態や欠測欄は当時のまま保持する。

Q5の固定χ256・累積8→10 sweepsでは、Szとbondの変化は前回の6→8より小さくなったが、
5精度条件はいずれも未達だった。Q1/Q3を固定した有限trial区間は
`0.258641424389<h/J<0.483732104557`へ変わった。
Q5は10 sweeps、Q1/Q3は8 sweepsなので、同じsweep数や精度での比較ではない。
CSL対照も交換energyは低下したが、切断とsweep energyの基準を超えている。
新しい分散は未取得であり、非零flux・非零ポンプは未実行のまま保持する。

## 1. NN Q5の固定χ結果

全交換energyはJ、分散はJ²単位。全状態は未収束trialのまま保持する。

| 量 | 累積8 | 累積10 | 10−8 |
| --- | ---: | ---: | ---: |
| E/J | −10.819403413718224 | −10.820525943702517 | −0.001122529984293 |
| 分散/J² | 0.008554313414 | 0.008202898286 | −0.000351415128 |
| 最終sweep切断誤差 | 5.4125737e−5 | 5.8985646e−5 | +4.8599091e−6 |

最終9→10のsweep energy差は0.000138190803 J。
分散はわずかに改善した一方、切断誤差は増えた。
密度・bondの最大変化は8→10でそれぞれ0.02968740、0.06857957 Jで、
前回6→8の0.12851927、0.24319599 Jより小さい。
9 sweep時点のprofileは測定していないため、この空間変化を9→10だけには帰属しない。

| 精度条件 | 上限 | 新しい測定値 | 判定 |
| --- | ---: | ---: | --- |
| energy変動/N（J） | 1e−6 | 4.1575185e−5 | 未達 |
| 最大Sz変化 | 1e−4 | 0.02968740 | 未達 |
| 最大bond変化（J） | 1e−4 | 0.06857957 | 未達 |
| 分散/N（J²） | 1e−5 | 3.0381105e−4 | 未達 |
| 最終sweep実測切断誤差 | 1e−6 | 5.8985646e−5 | 未達 |

SzのRMS差は0.00819834（ℏ=1）、bondは0.01842424 J。
円周並進y=0,1,2を照合し、両観測量とも移動なしが最小だった。
この差を剛体的な円周位置ずれだけで除くことはできない。x方向の並進は適用していない。

Q1のE/J=−11.562899472648773、Q3のE/J=−11.304258048259314を固定した。

| 比較 | h−/J | h+/J | 幅Δh/J |
| --- | ---: | ---: | ---: |
| Q5累積8 | 0.258641424389 | 0.484854634541 | 0.226213210152 |
| Q5累積10 | 0.258641424389 | 0.483732104557 | 0.225090680167 |

上限と幅の変化はともに−0.001122529984 J。幅は約0.496%狭まった。
これは追加最適化への感度であり、真の磁場境界の誤差棒ではない。

| Q3→Q5の列別追加Sz | x=0 | x=1 | x=2 |
| --- | ---: | ---: | ---: |
| Q5累積8 | 0.14385723 | 0.24073260 | 0.61541017 |
| Q5累積10 | 0.13721303 | 0.23302850 | 0.62975846 |

各行の和は1で、最大誤差8.22e−15。
同Q5の8→10差は`(−0.00664419,−0.00770410,+0.01434829)`で和は零となる。
右端列への偏りは残るが、三列・未収束の比較から端励起やbulk plateauを識別しない。

## 2. CSL対照

J1=1、J2=J3=0.5、N36・Q0・θ=0の別模型で、χ256の追加2 sweepsを完了した。
以前の準備4＋追加2＋flux開始用零flux2の履歴に今回の2を足すと累積10となる。
今回のworkerが直接記録するのは`2 sweeps in this invocation`であり、
累積8→10という説明はこの履歴に基づく。新しくfluxを印加した結果ではない。

| 量 | 親の零flux trial | 追加2 sweeps後 | 開始基準との比較 |
| --- | ---: | ---: | --- |
| E/J | −16.04147694678745 | −16.042073560344676 | 差−0.000596613557 J |
| 分散/J² | 0.240731514746 | **未取得** | 1e−3以下かは不明 |
| 最終sweep切断誤差 | 6.4077604e−4 | 7.2742117e−4 | 1e−5を超え未達 |
| 最終sweep energy差/J | 0.005864713767 | 0.000124208093 | 1e−4を超え未達 |
| 三角形chiralityの平均 | 0.001283092360 | 0.000405446232 | 受理条件には使わない |

外側600秒上限に従い、598.299秒でworkerの終了を確認した。
worker経過約507.25秒で2 sweepsが完了し、checkpoint保存と厳密reloadは約511.28秒までに終わった。
Schmidt・bond・chiralityを測定し、worker経過516.54秒に入ったH†Hによる分散測定が
上限までに終わらなかった。進行記録のworker時計と外側時計は独立なので、
各時刻の差を厳密な個別処理時間とは扱わない。

保存trialは`checkpoint-N8JX0z`。solver設定は`measure_variance=false`であり、
後処理の分散も未取得のままである。元validationの`status="running"`と
`active_phase="variance"`、外側executionの`status="timed_out"`をそのまま保存する。
途中のrecordにworker終了時のsource/parent確認フラグ、全診断のintegrity判定、
親MPSとのoverlapはない。これらを完了・通過・零値へ置き換えない。
保存したMPSのreload overlap誤差は1.78e−15で、親とのoverlapとは異なる。

取得済みの切断・energy差の二条件だけでも非零flux開始には不適格であり、
分散の値を推測して補う必要はない。chiralityの減少も相の同定や否定には使わない。
NNの磁場区間や1/9ポンプへこの別模型の結果を混ぜない。

親との最大profile差はSz=1.49320e−4、bond=0.00827773 J、chirality=0.00117431。
RMS差はそれぞれ5.70896e−5、0.00234008 J、0.00092971だった。
二つのSchmidt cutでentropy変化は+0.00723121、−0.00307312（自然対数）。
全てθ=0同士の再最適化であり、密度やSchmidt電荷の差をfluxによるポンプには数えない。

## 問いと事前の計算予算

NNのN27・Q5では、前回の累積6→8 sweepsでSz・bond分布が大きく変わった。
同じχ256で累積8→10へ進め、その後の固定χ変動を調べる。
Q1/Q3の既存χ256・累積8 trialは固定し、Q5だけを更新したときの
有限trial磁場区間と列別追加磁化の変化を測る。
前回の結果は[隣接sector比較](p4_static27_sector_refinement.md)に残す。

別模型のCSL対照はJ1=1、J2=J3=0.5、N36=(3,4)、Q=0、θ=0とし、
検証済みχ256 trialから2 sweepsを1 batch追加する。
variance・実測切断誤差・sweep energy差・密度・bond・chirality・Schmidt量を比較する。
今回は零flux準備の限定比較とし、非零fluxのtrajectoryは実行しない。

数値workerは同時に一つ、Julia/BLAS各1 thread、各走行は起動・保存・診断込み600秒上限。
両caseとも2 sweeps×1 batch、cutoff=noise=0、eigsolve_tol=1e-11、
Krylov次元40・最大反復20を維持する。上限で未完了ならその状態と未取得量を保存する。
同χ反復、χ拡大、seed・サイズ・sector走査へ自動拡大しない。

NNはwrap=3a2、CSLはwrap=4a2、軸OBC・円周PBC、胞内A/B/C順を用いる。
両端は範囲外xへの結合だけを除く。整数Q=2ΣSz、複素MPS、seam gaugeを維持する。
過去のcheckpointは対応するsource/runtimeの保存backendで厳密に再開し、
現行analysis driverのhashを別に残す。sourceやcheckpointの由来を新しい版へ付け替えない。

## 判定

NNの固定χ基準はenergy変動/N≤1e−6 J、最大Sz/bond変動≤1e−4、
variance/N≤1e−5 J²、最終sweep切断誤差≤1e−6。
energy変動は親からの変化と最終batch内のsweep・期待値差の最大を用いる。
各条件を個別に判定し、保存・観測量の整合性と精度を区別する。

Q1/Q3固定でh−=E3−E1、h+=E5−E3とする。
変分上界同士の差は厳密な境界の誤差区間ではない。
三列の密度差から端励起・domain wall・bulk応答を識別しない。

CSLの零flux事前基準はvariance≤1e−3 J²、最終切断誤差≤1e−5、
最後のsweep・期待値energy差≤1e−4 Jを維持する。
零fluxの診断だけで受理済みcontinuation、CSL成立、非零ポンプ校正を主張しない。
いずれも未収束・悪化・零chiralityを含めて保存する。

## 保存・費用・検証

NNは起動・保存・診断込み306.188秒、DMRG呼出し258.518秒、保存後の観測診断13.447秒。
worker CPUは304.062秒、peak RSSは2.663 GiBだった。
初回load/JITと親照合は包含時間に含み、独立費用へ分離していない。
trial `checkpoint-DRzoIB` を保存し、厳密reloadとSz・複素相関・bond和・Schmidt電荷/分散の
整合性を通過した。H†Hの虚部は2.84e−15 J²、丸め尺度2.60e−12 J²以下だった。
この整合性の通過と、上記5精度条件の未達は区別する。

CSLは停止までのworker CPU593.529秒、peak RSS3.520 GiB。
完了したDMRG呼出しは461.577秒、checkpoint保存/再読込みは約1.23/0.94秒、
Schmidt/bond/chiralityは約0.55/2.49/2.18秒だった。
分散処理は未完了なので所要時間の確定値を記入しない。

[NN測定](data/p4_q5_fixed_chi_8to10_validation.toml)と
[実行記録](data/p4_q5_fixed_chi_8to10_execution.toml)は元出力と同一byteで保存した。
[入力hash・全比較](data/p4_q5_fixed_chi_8to10_comparison.json)は
[専用解析](../../examples/analyze_q5_fixed_chi.py)で作成した。
完成状態・親checkpoint・config・保存source・runtime・指定累積sweepを検査してから比較し、
未完了recordと既存出力の上書きを拒否する。

CSLは保存時の180 source hashに一致するbackendを別directoryへ復元した。
177ファイルは現checkoutの一致byte、残る3ファイルは保存revision
`ba0391f248e606fd810218df8a93e66c552f2364`のgit blobから復元し、全hashを照合した。
[復元の由来と全source hash](data/p3_csl36_fixed_chi_backend_snapshot.json)を保存する。
これは現在の保存契約への移行ではなく、過去の実装・依存環境の再現である。
checkpoint本体とsource snapshotはローカル成果物で、文書だけから厳密再開はできない。

[CSLの未完了worker記録](data/p3_csl36_fixed_chi_8to10_validation.toml)と
[外側の時間切れ記録](data/p3_csl36_fixed_chi_8to10_execution.toml)は元出力と同一byteで保存した。
[CSL比較・欠測・入力hash](data/p3_csl36_fixed_chi_8to10_comparison.json)は
[CSL解析](../../examples/analyze_csl_preparation.py)の`--timed-out-variance`で作成した。
通常モードはこのrecordを拒否する。時間切れ専用モードは、終了確認済みの分散処理中停止と
2 sweeps・保存・既取得診断を検査し、分散、全診断のintegrity判定、worker終了フラグ、
親とのoverlapを`null`で保持する。保存状態のreload overlapは別の量として記録する。

[独立監査](data/p4_q5_csl_fixed_chi_provenance_audit.json)は、
実ファイルhash・metadata・source/runtime/config・親への結合に加え、
NNの区間・全5条件・列差・円周並進と、CSLの保存確率からのSchmidt量・列Sz・bond和を
集計scriptのhelperを呼ばずに再計算した。CSLの累積履歴4+2+2+2も親系列から確認した。
終了後のhash一致は独立した外部監査であり、未実行のworker終了時guardとは区別する。
この監査はPythonでのbyte・metadata・算術検査で、新しいMPS収縮・EDではない。
前回の監査では固定Q3のconfig原本archiveを発見できず、現configのhashと保存driverを照合した。
後続の上記比較で、旧`analysis-sources/examples/configs/`配置にある原本を発見しhashを確認した。
独立再計算と比較JSONの差はNNで0、CSLで最大2.22e−15で、保存値との不一致はなかった。

数値kernelは変更せず、package suiteと小系EDは再実行していない。
新解析の構文、通常モードによる未完了入力の拒否、既存出力の保持を確認し、実結果に適用した。
CSL専用モードは不正な時間切れ入力10件を拒否し、実結果でも欠測と未達の区別を確認した。

## 再現コマンド

次の2本を逐次実行した。再実行時は新しい出力directoryを選ぶ。

```sh
python3 examples/run_research.py static examples/configs/nn27_q5_chi256_resume8_to10.toml outputs/p4-q5-fixed-chi-8to10-20260920 --backend outputs/source-snapshots/n27-refine-before-chirality-20260919 --wall-seconds 600
python3 examples/run_research.py csl examples/configs/csl36_chi256_resume8_to10.toml outputs/p3-csl36-fixed-chi-8to10-20260920 --backend outputs/source-snapshots/csl36-pre-env-split-20260920 --wall-seconds 600
```

比較は以下で再作成できる。既存JSONの上書きは拒否するため、新しい出力名を選ぶ。
CSLの`--new`と`--execution`には、保存checkpointへの相対参照を検証するため元の出力を指定する。

```sh
python3 examples/analyze_q5_fixed_chi.py outputs/q5-fixed-chi-new.json --fixed-q1 docs/research/data/p4_sector_refine_nn27_q1_chi256_validation.toml --fixed-q3 docs/research/data/p4_campaign_nn27_q3_chi256_validation.toml --old-q5 docs/research/data/p4_sector_refine_nn27_q5_chi256_validation.toml --new-q5 docs/research/data/p4_q5_fixed_chi_8to10_validation.toml --new-execution docs/research/data/p4_q5_fixed_chi_8to10_execution.toml
python3 examples/analyze_csl_preparation.py outputs/csl-fixed-chi-new.json --old docs/research/data/p3_campaign_csl36_flux_readiness_chi256_validation.toml --new outputs/p3-csl36-fixed-chi-8to10-20260920/validation.toml --execution outputs/p3-csl36-fixed-chi-8to10-20260920/execution.toml --config examples/configs/csl36_chi256_resume8_to10.toml --timed-out-variance
```

旧backendはいずれもroot環境を使う。修正済みlauncherはこれを保持し、
現在のcheckoutで新規計算する場合は`research/Project.toml`を選ぶ。
今回のrunに使ったanalysis driverは各出力の`analysis-sources/`へ別hash付きで保存した。

## 次の限定比較

NNでは固定χの反復だけで精度を満たせなかった。次の候補は、同じQ5累積10の親から
χ256とχ512を同じ追加sweep数だけ進め、sweep数とχの効果を分離する比較である。
各走行の費用を別に予算化し、Q1/Q3の未収束や未探索Q≥7を残したままの比較とする。
周期9/27の競合状態、N54の長さ依存も未解決であり、このQ5の改善だけで優劣を決めない。

CSLは2 sweepsを保存できた一方、同じ600秒枠に分散診断が収まらなかった。
次は保存trialの分散を別途計時する計画と、χ変更を含む状態準備の予算を検討する。
これは未実行の次候補であり、今回の時間切れを後から完走扱いにしない。
数値上の開始基準、十分な長さ・幅、非零ポンプ校正は引き続き未達である。
