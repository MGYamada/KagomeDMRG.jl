# N27隣接sectorの固定χ追加sweep比較

実施日: 2026-09-20。状態: **両sectorの追加2 sweeps・保存・診断を完了。全精度基準は未達**。

Q1/Q5をχ256のまま累積6→8へ進めると、Q3固定での有限trial磁場区間は
`0.258641424389<h/J<0.484854634541`となった。前回より幅は約3.99%狭まった。
Q5は7→8でもenergyが大きく変わり、6→8の空間分布も再配列した。
7 sweep時点の密度は未測定なので、分布変化を最終sweepだけには帰属しない。
両sectorで分散と切断誤差は増えた。
同χ・同累積sweep数の比較は得られたが、数値精度の一致や基底収束は得られていない。

## 追加計算の結果

全交換energyをJ単位、分散をJ²単位で示す。全て未収束のtrialである。
Q3の累積8の既存状態は再最適化せず、同じ比較基準として使った。

| Q | E（累積6） | E（累積8） | ΔE（8−6） | 分散（6→8） | 最終sweep切断誤差（8） |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | −11.562622938931 | −11.562899472649 | −0.000276533718 | 0.003935862→0.004353532 | 3.47648e−5 |
| 5 | −10.810283686912 | −10.819403413718 | −0.009119726806 | 0.004854255→0.008554313 | 5.41257e−5 |

固定したQ3はE=−11.304258048259314、分散0.007326738807591937、
最終切断誤差4.982766955179003e−5。Q3も基底収束の基準には達していない。

| 比較集合 | h−/J | h+/J | 幅Δh/J |
| --- | ---: | ---: | ---: |
| 旧: Q1/5累積6、Q3累積8 | 0.258364890672 | 0.493974361348 | 0.235609470676 |
| 新: 全Qでχ256・累積8 | 0.258641424389 | 0.484854634541 | 0.226213210152 |
| 新−旧 | +0.000276533718 | −0.009119726806 | −0.009396260524 |

Q5の7→8のenergy変化は0.007817118209 Jであり、6→8の改善の大部分を占める。
Q1は7→8でも4.7876117e−5 J動いた。前回の小さなsweep差だけから
その後も境界が安定すると仮定できない。これらは追加最適化への感度で、厳密誤差棒ではない。

## 精度と空間構造の変化

事前基準を変えず、新しい測定から各判定を再計算した。**Q1/Q5とも下表の全5条件が未達**。
energy欄は`max(|E8−E6|, 最終batch内sweep差と最終期待値差の最大)/N`。

| 量 | 上限 | Q1 | Q5 |
| --- | ---: | ---: | ---: |
| energy変動/N（J） | 1e−6 | 1.02420e−5 | 3.37768e−4 |
| 最大Sz変化 | 1e−4 | 0.01571197 | 0.12851927 |
| 最大bond変化（J） | 1e−4 | 0.07425098 | 0.24319599 |
| 分散/N（J²） | 1e−5 | 1.61242e−4 | 3.16826e−4 |
| 最終sweep実測切断誤差 | 1e−6 | 3.47648e−5 | 5.41257e−5 |

旧→新のSz/bond profileのRMS差はQ1で0.00850953/0.02691915、
Q5で0.05819208/0.09605324だった。円周並進y=0,1,2を比較すると、
両sector・両観測量とも移動なしが最小で、今回の変化は剛体的な円周位置合わせで除けない。
x方向の移動は行っていない。これは秩序の同定ではなく、有限trialの再配列の診断である。

追加磁化の列別差は符号を保った`Sz(higher Q)−Sz(lower Q)`として求めた。

| 比較 | 状態 | x=0 | x=1 | x=2 |
| --- | --- | ---: | ---: | ---: |
| Q1→Q3 | 旧 | 0.35791624 | 0.18556656 | 0.45651720 |
| Q1→Q3 | 新 | 0.35933568 | 0.18981086 | 0.45085345 |
| Q3→Q5 | 旧 | 0.30636158 | 0.39890723 | 0.29473119 |
| Q3→Q5 | 新 | 0.14385723 | 0.24073260 | 0.61541017 |

各行の和は1で、最大誤差8.22e−15。Q5では右端列への配分が大きく変わったが、
三列・未収束の比較から端励起、domain wall、bulk plateauを識別できない。

## 保存・費用・検証

| ケース | 起動・保存・診断込み（秒） | DMRG呼出し（秒） | 保存後の観測診断（秒） | worker CPU（秒） | peak RSS（GiB） |
| --- | ---: | ---: | ---: | ---: | ---: |
| Q1再試行 | 322.767 | 273.069 | 13.510 | 318.751 | 2.887 |
| Q5 | 332.652 | 286.089 | 13.236 | 331.618 | 2.713 |

初回load/JIT、親照合、checkpoint保存・reload等は包含時間に含むが、独立の費用へ分離していない。
DMRG呼出しもその内部のcompileを含み得る。両走行は逐次で、各600秒上限内に正常終了した。
checkpointはQ1が`checkpoint-q5FuXN`、Q5が`checkpoint-fdv5cd`で、いずれもtrialのまま保持する。

- 両親状態と新状態の実MPSを旧backendで厳密loadし、site identity・電荷・norm・energyを照合した。
  新状態ではSz、複素相関、bond energy和、Schmidt電荷・分散の整合性と保存再読込が通過した。
- 今回測定したH†Hの虚部はQ1で1.49e−14、Q5で−5.25e−15 J²で、各丸め尺度以下だった。
  分散の増加を丸め誤差で置き換えない。Q3の旧欠測欄は変更せず、
  [以前の保存状態監査](data/p4_campaign_moment_audit_validation.toml)と区別する。
- 実行workerはsource/backend/親の不変性を確認した。別担当もconfig・親・保存状態・
  実行source・179-file backendのhashとruntimeを照合し、磁場差・全精度判定・
  円周並進と列別磁化差を独立に再計算した。
  [独立監査](data/p4_static27_sector_refinement_audit.json)の範囲はPythonでのhash・算術・metadata検査である。
  別のDMRG、MPS再収縮、独立EDを実行したとは扱わない。
  Q3旧走行はconfigのarchiveだけが未保存だった。現行configと記録hash、
  保存済みの4 driver、backendとcheckpointを照合し、欠けたarchiveは未保存のまま明記した。
- 既存集計selfcheckは0.00555秒で通過した。新しい解析scriptの構文、
  未完了のQ5 recordが拒否されること、実結果への適用を確認した。
  production sourceは変更せず、通常package suiteは再実行していない。

[Q1測定](data/p4_sector_refine_nn27_q1_chi256_validation.toml)・
[Q1実行記録](data/p4_sector_refine_nn27_q1_chi256_execution.toml)、
[Q5測定](data/p4_sector_refine_nn27_q5_chi256_validation.toml)・
[Q5実行記録](data/p4_sector_refine_nn27_q5_chi256_execution.toml)は元recordと同一byteで保存した。
[全比較と入力・解析hash](data/p4_static27_sector_refinement_comparison.json)は、
[専用解析](../../examples/analyze_sector_refinement.py)で再作成できる。
比較に使うbatchは指定した累積sweep数で固定し、energyの最小値では選んでいない。

## 問いと事前の計算予算

最近接等方Heisenberg模型、θ=0、`(Lx,Ly)=(3,3),N=27`で、
Q1/Q5のχ256・累積6 sweepsのtrialから各2 sweepsを追加する。
Q3の既存χ256・累積8 sweepsを固定して、隣接sectorの最適化だけによる
有限trial磁場区間とSz・bond profileの変化を調べる。
各sectorは1 batch、各600秒以内、数値workerは逐次一つ、Julia/BLAS各1 thread。
同じχ・累積sweep数でも数値精度が一致するとは仮定しない。
この単位では追加χ、seed、Q、CSL計算へ自動拡大しない。

wrapは`3a2`、x方向OBC・y方向PBC、胞内A/B/C順、端は範囲外xへの結合だけを除く。
Q=2ΣSzであり、1/9の対象はQ3、隣接Q1/5の磁化差は各1。
旧backend snapshotと同一runtimeから再開し、元のtrial・証拠は変更しない。

事前基準は[前回の走行](p4_p3_staged_campaign.md)から変更しない。
固定χでのenergy変動/N≤1e−6 J、最大Sz/bond変動≤1e−4、
variance/N≤1e−5 J²、最終sweep切断誤差≤1e−6を個別に確認する。
観測整合性と保存状態の厳密reloadの成功は、この精度基準の成功と分ける。
上限で未完了ならその状態を記録し、全診断が完了したbatchだけを比較する。

## 再現コマンド

```sh
python3 examples/run_research.py static examples/configs/nn27_q1_chi256_resume6_to8.toml outputs/p4-sector-refine-nn27-q1-chi256-20260920-retry --backend outputs/source-snapshots/n27-refine-before-chirality-20260919 --wall-seconds 600
python3 examples/run_research.py static examples/configs/nn27_q5_chi256_resume6_to8.toml outputs/p4-sector-refine-nn27-q5-chi256-20260920 --backend outputs/source-snapshots/n27-refine-before-chirality-20260919 --wall-seconds 600
```

二つを逐次実行する。再実行時は新しい出力directoryを使う。
snapshotの元commitは`48d15e2b566b4d4b50a6c5295b84348cc5f30f8d`。
実行時Git照会のHEADと保存した179 source hashを区別する。
現行analysis driverは別hashとして保存し、親recordのhashもconfigに固定した。
checkpoint本体とbackend snapshotはローカル成果物で、文書だけで別環境の厳密再開はできない。
この旧backendは現行のactive環境hash結合を追加する前の版であり、
その新しい保存契約を今回の再開へ遡って適用したとは扱わない。
保存backendを作業directoryにして`--project=.`で起動し、当時のsource・manifest・runtimeと
模型・site identity・実MPSを照合する。

Q1の最初の起動は、Julia launcherのlockfile作成がsandboxに拒否され、0.187秒で停止した。
Julia workerのvalidationは未作成で、DMRG結果には数えない。
元の`outputs/p4-sector-refine-nn27-q1-chi256-20260920/execution.toml`を保持し、
同じconfigと600秒上限で、権限を拡張した上記の別出力へ再試行した。
[起動失敗の実行記録](data/p4_sector_refine_nn27_q1_chi256_launch_failure.toml)も保存した。

新しい比較出力は次で再作成する。既存出力の上書きは拒否する。

```sh
python3 examples/analyze_sector_refinement.py outputs/p4-sector-refinement-new.json --old-q1 docs/research/data/p4_campaign_nn27_q1_chi256_validation.toml --old-q5 docs/research/data/p4_campaign_nn27_q5_chi256_validation.toml --fixed-q3 docs/research/data/p4_campaign_nn27_q3_chi256_validation.toml --new-q1 docs/research/data/p4_sector_refine_nn27_q1_chi256_validation.toml --new-q5 docs/research/data/p4_sector_refine_nn27_q5_chi256_validation.toml
```

## 比較式と解釈の範囲

全交換energyを使い、`F_Q=E_Q−hQ/2`、`h−=E3−E1`、`h+=E5−E3`とする。
Q3固定なので`δh−=−δE1`、`δh+=δE5`。
前回のQ1/5累積6・Q3累積8の区間は`0.258364890672<h/J<0.493974361348`だった。
新旧の差は最適化条件に対する感度であり、変分上界同士の差から得る厳密誤差区間ではない。
Q7以上、χ・seed・幅・長さの収束、bulk plateau・VBC・相同定は未判定とする。
三列のsector間磁化差は端局在とbulk応答を確定する証拠にはしない。
別模型のCSL非零ポンプ校正も未達のまま扱う。

## 次の限定比較

今回の最も大きな未解決量はQ5の最終sweepのenergy変動と、6→8の密度再配列である。
次の候補は、この保存状態から同じχ256でさらに2 sweepsを1 batchだけ追加し、
再配列後の固定χ変動を測ることである（未実行、別出力・各1 thread・上限600秒）。
分散と切断の停滞が残る場合は、同χの反復を全面走査にせず、χ拡大を別比較として選ぶ。
Q1/Q3のχ依存、Q7以上、N54の精度・端/中央の分離も未解決として保持する。
この比較単位の完了は、研究全体の収束や相同定の完了ではない。

上記候補のQ5累積8→10比較を、[次の限定研究](p4_q5_csl_fixed_chi_followup.md)で実施した。
変動は小さくなったが全5精度条件は未達で、同χ反復の次にはχ依存を分離する比較を候補とした。
