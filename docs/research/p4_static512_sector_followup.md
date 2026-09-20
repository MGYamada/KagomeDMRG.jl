# N27のχ512追加最適化と隣接sector比較

実施日: 2026-09-20。状態: **Q1・Q5とも計算・保存・別診断を完了、全5精度条件は未達**。

Q1のχ拡大を伴う追加最適化で分散は約8.454分の1となった。
Q5の固定χ512・追加1 sweepでは分散が6.286%減った一方、切断誤差は38.675%増えた。
既存Q3の完了trialと組み合わせた有限trial区間は
`0.258104370198<h/J<0.483570840219`。未収束の交換energyの差であり、bulk plateauの証拠とはしない。

## Q1の結果

Q1のχ256・累積8からχ512・累積9へ進め、全energyは0.000951728335 J低下した。
分散は約8.454分の1、実測切断誤差は約13.17分の1になったが、全5精度条件を超えている。
χ拡大と追加sweepを合わせた効果であり、純粋なχ依存の分離ではない。

| Q1の量 | 親χ256・累積8 | 新χ512・累積9 |
| --- | ---: | ---: |
| E/J | −11.562899472648773 | −11.563851200983553 |
| 分散/J² | 0.004353531611 | 0.000514962258 |
| 最終sweep実測切断誤差 | 3.4764758e−5 | 2.6395149e−6 |

| 個別の精度条件 | 上限 | 新しい値 | 判定 |
| --- | ---: | ---: | --- |
| energy変動/N（J） | 1e−6 | 3.5249198e−5 | 未達 |
| 最大Sz変化 | 1e−4 | 0.00876031 | 未達 |
| 最大bond変化（J） | 1e−4 | 0.03516845 | 未達 |
| 分散/N（J²） | 1e−5 | 1.9072676e−5 | 未達 |
| 最終sweep切断誤差 | 1e−6 | 2.6395149e−6 | 未達 |

親とのoverlap絶対値は0.99700284。χ変更を含むため固定χstationarityは評価対象外とした。
保存・reload、再測定energy/Sz/norm/Q、相関・Schmidt・bond和の整合性は通過した。
solve・状態保存は571.108秒、独立診断は87.186秒で、それぞれ900秒／180秒上限内だった。
DMRG呼出しは534.376秒、診断中のH†H収縮は47.377秒。

Q1/Q3/Q5をχ512・累積9で選ぶと、比較したtrialの区間は
`0.258104370198<h/J<0.483607503261`、幅0.225503133063 Jとなった。
同じχ・累積sweep数であり、精度一致や基底収束を意味しない。
Q5を累積10へ進めると上側だけがさらに0.000036663042 J下がる。

[Q1 solve測定](data/p4_static512_q1_solve_validation.toml)・
[実行](data/p4_static512_q1_solve_execution.toml)、
[別診断](data/p4_static512_q1_diagnostics_validation.toml)・
[診断実行](data/p4_static512_q1_diagnostics_execution.toml)は元出力と同一byteで保存した。

## Q5の固定χ追加1 sweep

χ512の累積9から10へ進め、energyは0.000036663042 J低下した。
分散の改善と切断誤差の悪化が同時に起きており、固定χでの収束は確認できない。
親とのoverlap絶対値0.99943508も、局所分布の定常性を保証しない。

| Q5の量 | 親χ512・累積9 | 新χ512・累積10 |
| --- | ---: | ---: |
| E/J | −10.822139327524651 | −10.822175990566798 |
| 分散/J² | 0.000827696370 | 0.000775670951 |
| 最終sweep実測切断誤差 | 4.3042434e−6 | 5.9688929e−6 |

| 個別の精度条件 | 上限 | 新しい値 | 判定 |
| --- | ---: | ---: | --- |
| energy変動/N（J） | 1e−6 | 1.3578904e−6 | 未達 |
| 最大Sz変化 | 1e−4 | 0.00368541 | 未達 |
| 最大bond変化（J） | 1e−4 | 0.00804630 | 未達 |
| 分散/N（J²） | 1e−5 | 2.8728554e−5 | 未達 |
| 最終sweep切断誤差 | 1e−6 | 5.9688929e−6 | 未達 |

保存と観測integrityは通過、固定χstationarityは未達と別に記録した。
[Q5 solve測定](data/p4_static512_q5_solve_validation.toml)・
[実行](data/p4_static512_q5_solve_execution.toml)、
[別診断](data/p4_static512_q5_diagnostics_validation.toml)・
[診断実行](data/p4_static512_q5_diagnostics_execution.toml)を元出力と同一byteで保存した。

## 隣接sectorによる有限trial区間

全交換energyから `h−=E3−E1, h+=E5−E3` を計算する。
整数chargeが2ずつ異なるため、物理磁化差は1である。
既存Q3・χ512・累積9の値はE/J=−11.305746830785633、分散/J²=0.000740651514、
切断誤差3.8014601e−6。これも精度未達である。

| 選択したtrial | h−/J | h+/J | 幅/J |
| --- | ---: | ---: | ---: |
| Q1/Q3/Q5ともχ256・累積8 | 0.258641424389 | 0.484854634541 | 0.226213210152 |
| 前回: Q5だけχ512・累積9 | 0.258641424389 | 0.482118720735 | 0.223477296345 |
| 前回の選択からQ3もχ512・累積9へ | 0.257152641863 | 0.483607503261 | 0.226454861398 |
| Q1/Q3/Q5ともχ512・累積9 | 0.258104370198 | 0.483607503261 | 0.225503133063 |
| 今回: Q5だけさらに累積10 | 0.258104370198 | 0.483570840219 | 0.225466470021 |

前回の選択からQ3を差し替えると、下側は0.001488782526 J下がり上側は同量上がる。
次にQ1を差し替えると下側が0.000951728335 J上がる。Q5の追加sweepは上側を
0.000036663042 J下げる。境界の変化は複数sectorの誤差を含み、変分上界同士の差から
厳密な境界の上下限や誤差棒は得られない。Q≥7との競合も未確認である。

## 密度・bond分布

| 親→新trial | Sz差の最大値 / RMS | bond差の最大値 / RMS（J） |
| --- | ---: | ---: |
| Q1: χ256累積8→χ512累積9 | 0.00876031 / 0.00438949 | 0.03516845 / 0.01392610 |
| Q5: χ512累積9→10 | 0.00368541 / 0.00104934 | 0.00804630 / 0.00244515 |

円周方向の並進y=0,1,2を全27site・48bondで比較すると、両ケースともSzとbondの
RMS最小はy=0だった。Q5のSz RMSは順に0.00104934、0.00119127、0.00122501、
bond RMSは0.00244515、0.00275800、0.00284532である。
有限trialの円周方向の位置比較であり、9/27-site VBCとの完全なtemplate照合や相の選択ではない。

最新のsector差をx列ごとに合計すると、以下になる。各列は9site。
値は符号付きSz差で、絶対値化・負値の切捨て・再正規化をしていない。

| 高Q−低Q | x=0 | x=1 | x=2 | 全体 |
| --- | ---: | ---: | ---: | ---: |
| Q3−Q1 | 0.35814517 | 0.19084541 | 0.45100942 | 1 |
| Q5−Q3 | 0.14235784 | 0.23412076 | 0.62352140 | 1 |

丸め前の総和誤差は4e−15未満。追加磁化は一様には分布しないが、3列の有限trialから
端局在とbulk応答を分離できない。幅・長さ依存と収束を確かめる前に秩序・plateauと解釈しない。

## 実測費用

| ケース・段階 | wall秒 | worker CPU秒 | peak RSS GiB | wall上限秒 |
| --- | ---: | ---: | ---: | ---: |
| Q1 solve・保存 | 571.108 | 569.196 | 2.747 | 900 |
| Q1 保存状態の診断 | 87.186 | 85.944 | 3.853 | 180 |
| Q5 solve・保存 | 668.185 | 666.009 | 2.770 | 900 |
| Q5 保存状態の診断 | 87.743 | 86.752 | 4.028 | 180 |

Q5のDMRG呼出しは631.513秒、別processでの全観測は51.745秒、うちH†Hは46.237秒。
4段階のwall合計は1414.223秒。起動・読込み・再測定・保存等の費用を含め、各上限内で終了した。
RSSは各worker終了後のOS high-water値であり、launcherや同時process木全体のpeakではない。
分離によって完了状態と診断の成否を区別できたが、起動費用も増すため高速化を主張しない。

## 問いとケース

前回の[同一親χ比較](p4_q5_matched_chi_comparison.md)では、Q5のχ256→512で
分散と切断誤差が改善したが、精度基準は未達だった。今回は次の2ケースに限定する。

| ケース | 親の状態 | 新しい計算 | 解く問い |
| --- | --- | --- | --- |
| Q5 | χ512・累積9 | χ512で1 sweep、累積10 | χ拡大後の固定χ変動と精度 |
| Q1 | χ256・累積8 | χ512で1 sweep、累積9 | 下側境界を決める隣接sectorの追加最適化への感度 |

Q1の変化にはχ変更と追加1 sweepが両方含まれる。純粋なχ依存を分離する比較ではない。
既存のQ3・χ512・累積9 trialと今回の完了trialを選んでQ1/Q3/Q5の有限区間を比較する。
同じχ・sweep数を同じ精度や収束の証拠にはしない。
旧Q3走行は次のbatchで時間切れになったため、完了済み最初のbatchだけを由来・hash付きで
利用し、旧worker全体の完走とは区別する。未完了batchは比較に入れない。

模型は最近接・等方Heisenberg、J=1、θ=h=0、N27=(3,3)、wrap=3a2。
軸OBC・円周PBC、胞内A/B/C順、範囲外xに出る結合のみ除いた端を維持する。
整数Q=2ΣSz、複素MPS、seam gauge、cutoff=noise=0、eigsolve_tol=1e−11、
Krylov次元40・最大反復20、seed=11を保持する。

## 事前の計算予算と保存

前回のQ5・χ512はsolve呼出し495.529秒、全診断50.300秒、起動等を含め579.474秒だった。
今回、各ケースの **solve・状態保存は900秒以内、別processの観測診断は180秒以内** とする。
どちらも起動・読込み・保存を含み、外側supervisorで単一solve中も上限停止する。
数値processは同時に一つ、Julia/BLAS各1 thread、2ケース×1 sweepで打ち切る。
この予算は研究用の限定走行で、通常package testの時間上限を変更しない。

solveでは元の保存backend/runtimeから厳密再開し、完了trialを原子的に保存・reloadする。
観測診断は独立した出力へ記録し、solve時の未測定分散や未実行integrityを後から書き換えない。
途中停止の場合は保存の有無と完了段階を区別し、未測定値は補間しない。
Hamiltonian・DMRG kernel・切断法・checkpoint契約は変更しない。

親の正本はQ1が[p4_sector_refine_nn27_q1_chi256_validation.toml](data/p4_sector_refine_nn27_q1_chi256_validation.toml)、
Q5が[p4_q5_matched_chi512_8to9_validation.toml](data/p4_q5_matched_chi512_8to9_validation.toml)。
新configでrecordとmetadataのSHA-256を固定する。

## 判定

energy変動/N≤1e−6 J、最大Sz/bond変化≤1e−4、variance/N≤1e−5 J²、
最終sweep実測切断誤差≤1e−6を維持する。Q5では固定χの変化を調べる。
Q1の親比較はχ変更なのでstationarityは評価対象外とし、各閾値への比較のみ行う。
一回の新しい呼出しは1 sweepのため、その呼出し内の最後の二sweep差は未測定。
親→新状態のenergy差と最終solver energy対MPO期待値の差は別に保存する。

比較は全交換energyで `h−=E3−E1, h+=E5−E3` を求め、Q1/Q3/Q5をどれに差し替えたかを
明示して区間変化を分解する。変分上界同士の差は厳密な境界の上下限や誤差棒ではない。
列別の追加Sz、Sz/bond分布の最大差・RMS、円周並進y=0,1,2を比較する。
Q≥7、長さ・幅・seed依存、端とbulkの分離、neutral gapと相同定は未解決のまま残す。

CSLは別模型の[前回の分散補足](p4_q5_matched_chi_comparison.md)で開始精度未達と判定済み。
今回はNNの二つの問いに計算資源を割り当て、CSLの新たなDMRGや非零fluxは含めない。

## 実行と検証の方法

[分離worker](../../examples/research_static_split.jl)は既存のconfig読込み・厳密な親読込み・
観測診断を再利用する。[launcher](../../examples/run_static_split.py)は既存の外側supervisorを
変更せず、各段階に明示的な上限を渡す。solveの完了は`completed_solve_diagnostics_deferred`、
診断の完了は`completed_saved_state_diagnostics_accuracy_separate`として別出力へ保存する。
診断はsolveのrecord/execution/checkpoint/config/解析sourceのhashと正常終了を照合し、
新規MPOによるenergy、Sz、norm、Qの再測定後に相関・分散・Schmidt量を測る。
Q5の固定χ判定とQ1のχ変更を区別し、全観測integrityと数値精度も別に判定する。

起動前に別担当が式・保存由来・設定・時間監視を確認し、Julia/Pythonの構文とCLIを確認した。
旧Q3の保存済みχ512・累積9はH†Hの実部・虚部も記録しており、新たな分散収縮は要求しない。
今回の新しい走行は以下を逐次実行した。再現時は新しい出力directoryを指定する。

```sh
python3 examples/run_static_split.py solve examples/configs/nn27_q1_chi512_split8_to9.toml outputs/p4-static512-q1-8to9-20260920-solve --wall-seconds 900
python3 examples/run_static_split.py diagnose examples/configs/nn27_q1_chi512_split8_to9.toml outputs/p4-static512-q1-8to9-20260920-diagnostics --solve-output outputs/p4-static512-q1-8to9-20260920-solve --wall-seconds 180
python3 examples/run_static_split.py solve examples/configs/nn27_q5_chi512_split9_to10.toml outputs/p4-static512-q5-9to10-20260920-solve --wall-seconds 900
python3 examples/run_static_split.py diagnose examples/configs/nn27_q5_chi512_split9_to10.toml outputs/p4-static512-q5-9to10-20260920-diagnostics --solve-output outputs/p4-static512-q5-9to10-20260920-solve --wall-seconds 180
```

backendは親と同じ`outputs/source-snapshots/n27-refine-before-chirality-20260919`。
origin revisionは`48d15e2b566b4d4b50a6c5295b84348cc5f30f8d`、実使用sourceは179ファイルの
SHA-256表で固定した。Julia 1.13.0、ITensors 0.9.31、ITensorMPS 0.4.1、NDTensors 0.4.31+1。
実行configと研究driver一式のコピーを各出力の`analysis-sources/`へ保存した。
元のsource revisionと新しい研究driverを混同せず、各recordのhashで識別する。

[集計script](../../examples/analyze_static512_followup.py)と
[機械可読な全比較](data/p4_static512_sector_followup_comparison.json)には、7状態の選択理由、
実保持次元、全5区間、分布差、符号付き列差と円周並進を保存した。
docsへ保存した診断recordからも再解析し、数値・由来hash・設定・主張の範囲が一致した。
再解析は以下をrepository rootで実行する（出力は存在しないpathにする）。

```sh
python3 -B examples/analyze_static512_followup.py /tmp/static512-reanalysis.json \
  --q1-diagnostics docs/research/data/p4_static512_q1_diagnostics_validation.toml \
  --q5-diagnostics docs/research/data/p4_static512_q5_diagnostics_validation.toml
```

[独立監査](data/p4_static512_sector_followup_independent_audit.json)は集計helperをimportせず、
保存byte・入力hash・runtime・設定・親子関係と式を別実装で照合した。
571項目を通過、不一致0、最大算術残差2.16e−14。8recordの同一byte保存、4段階の正常終了と
時間上限、相関からのbond再計算、固定Q共分散、Schmidt charge、precision flags、区間・分布差を含む。
監査本体は0.277秒で、埋込みPython sourceをリード担当も再実行して同一判定・入力hashを確認した。
この監査は保存結果の算術・由来の検証であり、新たなMPS収縮や独立ED照合ではない。
旧Q3 χ256のraw complex H†H欠測、旧Q3 χ512の第2batch中断を保持している。
旧χ512の完了第1batchについては既存のpost-run hash auditと今回のbyte照合を利用した。

研究driverと集計だけの変更で、Hamiltonian・DMRG kernel・checkpoint契約は変更していない。
Python/Julia構文・CLI確認、実際の4段階走行、厳密reload・保存状態の再測定、再解析と
独立監査を実施した。package全suiteと既存ED計算は再実行していない。

## 次の判断

**追記（同日）:** 以下は計算直後の候補だったが、ユーザーの指摘を受けた
[研究方針の再検討](research_direction_20260920.md)でQ3追加の優先度を下げた。
現在は長い同一円筒の競合する中央構造・境界依存の比較を先に置く。以下は判断履歴として残す。

今回の限定goalは完了したが、N27の収束やbulk plateauの研究目標は未達である。
次の候補は、両磁場境界に現れるQ3・χ512・累積9からの固定χ追加1 sweepと別診断とする。
旧走行は第2batch中断のため、この固定χ変動はまだ測れていない。
同じ保存先行方式と900秒／180秒上限、各1 threadを候補予算とし、起動前に最新費用を確認する。
現在の分離workerはQ1/Q5のみを許可する。Q3の実行には、case許可と正確な親hash、
旧走行の完了第1batch・部分完了監査の照合を追加し、再レビューする必要がある。

Q5はenergy変化が小さくなっても分布・分散・切断の基準を満たしていない。
固定χ反復のみで改善すると仮定せず、さらにχを増す前には時間・メモリを見積もり、
問いを一つに絞った比較を選ぶ。今回の2ケースを機械的に全sector・全seedへ拡大しない。
Q≥7のsector skipping、別初期状態、長さ・幅依存、neutral gap、端とbulkの分離は引き続き未解決。
既知CSLの非零ポンプ校正も別模型で未達のまま残す。
