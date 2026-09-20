# N27 Q5の同一親・同一sweep数によるχ比較

実施日: 2026-09-20。状態: **NNの同一親・各1 sweep比較とCSL欠測分散回収を完了。両模型とも精度未達**。

後続の[Q1追加χ512・Q5固定χ512比較](p4_static512_sector_followup.md)を別記録に保存した。
本記録の同一親比較とCSL分散補足は、実施時点の結果として保持する。

同じ累積8の親から1 sweep進めると、χ512はχ256より全energyが0.001751574625 J低く、
分散は約9.815分の1、実測切断誤差は約13.44分の1となった。
χ256を9→10へさらに1 sweep進めた既存結果では、energy低下が0.000138190803 Jにとどまり、
分散と切断誤差は増えた。今回の限定比較はχ拡大を調べる有効性を示すが、
χ512も全5精度条件を超えており、基底状態・磁場境界・秩序の収束は未確立である。

## NNの結果

全状態はθ=0の未収束trial。energyは全交換energyである。

| 量 | 共通の親 χ256・8 sweeps | χ256・9 sweeps | χ512・9 sweeps | 既存χ256・10 sweeps |
| --- | ---: | ---: | ---: | ---: |
| E/J | −10.819403413718224 | −10.820387752899613 | −10.822139327524651 | −10.820525943702517 |
| 分散/J² | 0.008554313414 | 0.008123807838 | 0.000827696370 | 0.008202898286 |
| 最終sweep実測切断誤差 | 5.41257e−5 | 5.78302e−5 | 4.30424e−6 | 5.89856e−5 |

`δχE=−0.001751574625 J`、`δsE=−0.000984339181 J`。
新しいχ256・累積9のMPO期待値は、既存8→10走行の最初のsolver報告energyと数値上一致した。
期待値とsolver報告は区別して保存し、このenergy一致からMPS全体の同一性は主張しない。
χ512・9とχ256・10の比較はχとsweep数が異なる補助比較に限る。

| 精度条件 | 上限 | χ256・9（親→新） | χ512・9（親→新） |
| --- | ---: | ---: | ---: |
| energy変動/N（J） | 1e−6 | 3.64570e−5 | 1.01330e−4 |
| 最大Sz変化 | 1e−4 | 0.02263730 | 0.02750969 |
| 最大bond変化（J） | 1e−4 | 0.05289354 | 0.06562477 |
| 分散/N（J²） | 1e−5 | 3.00882e−4 | 3.06554e−5 |
| 最終sweep実測切断誤差 | 1e−6 | 5.78302e−5 | 4.30424e−6 |

両列とも全条件を超える。χ512への親比較はχ変更を含むため、固定χstationarityは
評価対象外とし、個別閾値への比較だけを示す。energy低下が大きいこと自体は失敗ではなく、
精度と観測量の安定をまだ確認できていないという意味である。
両走行の最後のsolver energyとMPO期待値の差は0、最後の二sweep差は未測定。

χ256・9→χ512・9の最大Sz差は0.00761733、RMSは0.00261190。
bondの最大差は0.01490427 J、RMSは0.00502549 Jだった。
円周並進y=0,1,2では両観測量とも移動なしが最小で、単なる円周位置ずれでは差を除けない。
x方向の並進やVBCの波動関数テンプレートへの照合は行っていない。
χ256・9→10の最大Sz/bond変化は0.00705010／0.01568603 Jであり、
今回のχ変更による空間差と同程度に残っている。

### Q1/Q3を固定した有限trial磁場区間

Q1/Q3はχ256・累積8のE/J=−11.562899472648773／−11.304258048259314を維持した。

| Q5の選択 | h−/J | h+/J | 幅Δh/J |
| --- | ---: | ---: | ---: |
| 共通の親χ256・8 | 0.258641424389 | 0.484854634541 | 0.226213210152 |
| χ256・9 | 0.258641424389 | 0.483870295360 | 0.225228870970 |
| χ512・9 | 0.258641424389 | 0.482118720735 | 0.223477296345 |
| 既存χ256・10 | 0.258641424389 | 0.483732104557 | 0.225090680167 |

χ256・9→χ512・9で上限と幅は0.001751574625 J低下した。
Q1/Q3は未収束のままで、χ512列はQ5だけχを増やした比較である。
同じ精度でのsector比較、真の境界の誤差棒、bulk plateauとは扱わない。

| Q3→Q5の列別追加Sz | x=0 | x=1 | x=2 |
| --- | ---: | ---: | ---: |
| χ256・9 | 0.13685313 | 0.23331666 | 0.62983021 |
| χ512・9 | 0.14075457 | 0.23370040 | 0.62554503 |

各行の和は1、最大和則誤差4.45e−15以下。
χ変更の同Q差は `(0.00390144,0.00038374,−0.00428518)` で和は零。
右端列への偏りは残るが、三列・未収束の結果から端励起とbulk応答は分離できない。

### 保存・費用

| 量 | χ256 | χ512 |
| --- | ---: | ---: |
| 起動・保存・全診断込み秒 | 186.069 | 579.474 |
| DMRG呼出し秒 | 136.606 | 495.529 |
| 全観測診断秒 | 13.315 | 50.300 |
| うちH†H分散秒 | 10.420 | 45.984 |
| worker CPU秒 | 181.036 | 577.715 |
| worker peak RSS GiB | 2.808 | 3.539 |

初回load/JIT等は包含時間に入り、独立時間には分離していない。
両trialは保存・厳密reload、Sz・複素相関・bond和・Schmidt平均/分散の整合性を通過した。
χ512は600秒上限に近く、今後の費用を同じ時間内で保証できない。

[χ256測定](data/p4_q5_matched_chi256_8to9_validation.toml)・
[実行](data/p4_q5_matched_chi256_8to9_execution.toml)、
[χ512測定](data/p4_q5_matched_chi512_8to9_validation.toml)・
[実行](data/p4_q5_matched_chi512_8to9_execution.toml)は元出力と同一byteで保存した。
[比較JSON](data/p4_q5_matched_chi_8to9_comparison.json)は
[解析script](../../examples/analyze_q5_matched_chi.py)で生成し、両分岐の同一親・設定・
runtime/source、checkpoint checksum、実測solver設定、保存driver/configを照合した。
元のcheckpoint本体とbackend snapshotはローカル成果物である。

## 別模型CSLの欠測分散を回収

前回のN36・Q0・J1=1,J2=J3=0.5・θ=0・χ256の保存trial
`checkpoint-N8JX0z`を元のsource/runtimeで厳密reloadし、MPO期待値と複素H†Hを再測定した。
追加DMRG、ED、flux印加は行っていない。

| 量 | 親の零flux状態 | 前回保存trialの今回後処理 |
| --- | ---: | ---: |
| E/J | −16.04147694678745 | −16.042073560344676 |
| 分散/J² | 0.240731514746 | 0.242370558567 |
| 最終sweep切断誤差（前回測定） | 6.4077604e−4 | 7.2742117e−4 |
| 最終sweep energy差/J（前回測定） | 0.005864713767 | 0.000124208093 |

分散は0.001639043821 J²、約0.681%増加した。今回取得した分散は開始基準1e−3 J²を超え、
前回測定の切断1e−5・sweep差1e−4 Jの基準も未達のままである。
固定χ256でenergyが下がっても分散が減るとは限らず、同じχの反復だけで準備が解決したとはいえない。
非零CSLポンプの校正は引き続き未達で、NNの区間とは別の結果である。

新しいH†Hの虚部は−2.026e−14 J²、丸め尺度5.714e−12 J²以内。
再測定energyと保存energyの差は複素数の絶対値で6.731e−16 J、Sz profile差は0。
normは前後とも1、全Sz誤差1.91e−16、測定前後の全tensorの最大変化は0だった。
10個の診断整合性条件と入力・解析source・backend不変guardを通過した。
測定前後overlap誤差1.78e−15は、異なる累積sweep状態間のoverlapではない。

300秒上限に対し起動・保存込み123.524秒、H†H収縮92.134秒、strict load23.878秒。
worker CPU122.280秒、peak RSS3.639 GiBだった。
前回は保存後の残り時間で分散処理を終えられなかったが、今回の独立後処理では完了した。
元の`running` validationと`timed_out` executionはそのまま保持する。
この新しい診断は、旧workerの未実行の終了guard・全観測integrity・親とのoverlapを
遡って実行済みにするものではない。

[分散補足の測定](data/p3_csl36_saved_moment_validation.toml)・
[実行記録](data/p3_csl36_saved_moment_execution.toml)を元出力と同一byteで保存した。
[後処理worker](../../examples/audit_csl_saved_moment.jl)と
[300秒launcher](../../examples/run_csl_saved_moment.py)が、元の停止位置、入力6ファイルの
固定hash、保存driver、model/source/runtime/settingsを検証してから測定する。

## 検証範囲

数値kernelは変更していない。全package suiteや小系EDは重複実行していない。
今回の実走行はNN2件とCSL後処理1件で、全て保存・診断まで上限内に完了した。
新解析の構文、未完了record・誤χ・変更した親hash・既存出力の拒否、
旧レイアウトのconfig/driver検証を行った。CSL workerのJulia構文とPython launcherの
構文・CLIも確認し、実データに対する10診断を通過した。

固定Q3のconfig原本は古い保存配置
`outputs/p4-campaign-nn27-q3-chi256/analysis-sources/examples/configs/nn27_q3_chi256_resume.toml`
にあり、今回hashを確認した。過去の監査で未発見だった原本の所在を補足する。
Q3の旧recordに`backend_unchanged`欄がない事実は変更していない。

[独立監査](data/p4_q5_matched_chi_csl_provenance_audit.json)では解析helperを呼ばず、
実ファイルhash・source/runtime・設定・親のつながり、保存相関からのbond/電荷和則、
Schmidt確率からの平均・分散・entropy、磁場区間とprofile/列差・円周並進を再計算した。
NNの最大算術残差は1.42e−14、CSLの分散再計算差は0、不一致は0件だった。
leadも独立監査と集計JSONのscalar/profile/円周並進/精度判定58項目を照合し、差は0だった。
これは保存値の独立算術・byte監査であり、新しいEDやMPS収縮ではない。
再現用Pythonソースは監査JSON内に保存した。

## 次の判断

Q5のχ256追加反復をそのまま継続せず、χ512での固定χ変動を次の候補にする。
一方、磁場区間の解決にはQ1/Q3側の精度も必要であり、Q5だけの改善を収束した区間へ
読み替えない。既存Q3・χ512 trialと並ぶQ1の限定χ比較を候補として残す。
χ512の実測が上限に近いため、追加計算前に1 sweepの保存と後処理の予算を分けて設定する。
今回の単位から広いχ/seed/sector走査へ自動拡大しない。
CSL側も固定χ256での分散・切断の改善が見られないため、次の準備はχ依存の限定比較を
別予算で検討する。今回確保した後処理の分離を利用し、非零fluxへはまだ進めない。

## 事前に固定した問いと計算予算

最近接・等方Heisenberg模型のN27=(3,3)、Q=5、θ=0で、同じ累積8-sweep・χ256の
trialからχ256とχ512をそれぞれ1 sweepだけ進める。親は
`outputs/p4-sector-refine-nn27-q5-chi256-20260920/trial/checkpoint-fdv5cd`。
既存のχ256・累積8→10走行には9 sweep時点のdensity/bond/varianceがないため、
今回のχ256分岐で累積9の観測量を揃える。

比較する量は `δχE=E(512,9)−E(256,9)` と `δsE=E(256,9)−E(256,8)`。
同じ最適化1 sweepに対するχ上限変更の効果を比べるもので、十分なsweep数で収束した
状態間のχ誤差を分離する実験ではない。既存のχ256・累積10は補助比較として保持する。
Q1/Q3はχ256・累積8の既存trialに固定し、`h−=E3−E1, h+=E5−E3` と
列別の追加Szを再計算する。未探索Q≥7、端励起、bulk plateau、相同定は未判定のままとする。

各NN走行は起動・保存・診断込み600秒、1 sweep×1 batch。数値processは同時に一つ、
Julia/BLAS各1 thread。cutoff=noise=0、eigsolve_tol=1e−11、Krylov次元40・最大反復20、
複素MPS・seam gauge・整数Q=2ΣSzを維持する。wrap=3a2、軸OBC・円周PBC、
胞内A/B/C順、範囲外xに出る結合のみ除いた端を使う。
旧source/runtimeの保存backendを厳密に再使用し、親の由来は付け替えない。

精度基準はenergy変動/N≤1e−6 J、最大Sz/bond変化≤1e−4、
variance/N≤1e−5 J²、実測切断誤差≤1e−6を維持する。
1 sweepのため「最後の二sweep差」は欠測とし、親→新状態のenergy変動と
最終solver energy対MPO期待値の差を区別して保存する。
χ512への親比較はχ変更なので固定χのstationarityとは呼ばない。
円周並進y=0,1,2も調べ、profile差が単なる円周位置ずれで説明できるか確認する。

別模型CSL N36=(3,4)、Q0、J1=1・J2=J3=0.5、θ=0では、前回保存した
`checkpoint-N8JX0z`を再読込み、追加DMRGなしで欠測H†H分散を測定する。
この後処理は独立した300秒上限とし、入力hash・source/runtime・期待値・norm・電荷・
保存profileとの一致を確認する。元の時間切れrecordは変更しない。
新しい分散が取得できても、既取得の切断・sweep差はCSL開始基準未達であり、
非零flux受理や非零ポンプ校正の達成とはしない。

## 再現コマンド

それぞれ新しい出力directoryを用い、以下を逐次実行する。

```sh
python3 examples/run_research.py static examples/configs/nn27_q5_chi256_matched8_to9.toml outputs/p4-q5-matched-chi256-8to9-20260920-r1 --backend outputs/source-snapshots/n27-refine-before-chirality-20260919 --wall-seconds 600
python3 examples/run_research.py static examples/configs/nn27_q5_chi512_matched8_to9.toml outputs/p4-q5-matched-chi512-8to9-20260920 --backend outputs/source-snapshots/n27-refine-before-chirality-20260919 --wall-seconds 600
python3 examples/run_csl_saved_moment.py outputs/p3-csl36-saved-moment-20260920 --wall-seconds 300
```

初回χ256起動はsandbox内でJulia launcherのlockfile作成に失敗し、数値workerに到達しなかった。
そのexecutionを保持し、許可された再実行を別の`-r1`出力へ分離した。
[起動失敗の記録](data/p4_q5_matched_chi256_launch_failure.toml)は研究結果に含めない。

比較は次のコマンドで再作成できる。既存出力の上書きを拒否するため新しいJSON名を選ぶ。

```sh
python3 examples/analyze_q5_matched_chi.py outputs/q5-matched-chi-new.json --fixed-q1 docs/research/data/p4_sector_refine_nn27_q1_chi256_validation.toml --fixed-q3 docs/research/data/p4_campaign_nn27_q3_chi256_validation.toml --parent-q5 docs/research/data/p4_sector_refine_nn27_q5_chi256_validation.toml --chi256-q5 docs/research/data/p4_q5_matched_chi256_8to9_validation.toml --chi512-q5 docs/research/data/p4_q5_matched_chi512_8to9_validation.toml --previous25610 docs/research/data/p4_q5_fixed_chi_8to10_validation.toml
```
