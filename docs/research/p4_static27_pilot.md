# P4段階B: 27サイト・1/9磁化の資源・精度pilot

実行日: 2026-09-19。[静的研究方針](p4_static_plateau_strategy.md)の段階B。
最近接・等方Heisenberg模型で、磁化を固定した代表1点の費用と診断を測る。
この段階では磁化境界、量子化ポンプ、bulk plateau、競合秩序の優劣は判定しない。

**結果: 600秒の予算内で時間上限停止。段階Bの数値的完了は未達。**
起動込み598.239秒、CPU 595.894秒、peak RSS 2.08034 GiBを記録したが、
`run_dmrg` は戻らず、最終状態・エネルギー・観測量は取得できなかった。
完了sweep数も不明であり、この試行からχ256の精度を判定しない。

## 実行前に固定した範囲

`(Lx,Ly)=(3,3),N=27,Q=3,ΣSz=1.5,Nup=15,Ndown=12`。
wrap `3a2`、軸OBC・円周PBC、順序 `x,y,A/B/C`、完全な単位胞を残す端で、
最近接bondは48本。`Jxy=Jz=1,J2=J3=hz=θ=0`、seam gauge。
中央領域はsite 10:18の一列、幾何cutはMPS bond9と18。
9-siteと27-siteの競合周期を収容できるが、十分なbulkや秩序収束は仮定しない。

複素ランダムMPSのseed=11、χ=256、8 sweeps、cutoff=noise=0。
局所Krylov tolerance `1e-11`、dimension 40、最大反復20、varianceを計測する。
追加sweep、χ変更、別seed、磁化セクターの追加はこの走行では行わない。
N27の独立EDや全Hilbert空間へのMPS展開も行わない。

```sh
python3 examples/run_static27.py outputs/p4-static27-new --wall-seconds 600
```

外側から起動・import・compileを含む600秒を監視し、停止猶予の最大2秒も上限に含める。
Julia/BLAS各1thread、数値process一つ。新規/空の出力ディレクトリだけを受け付ける。
[N27 launcher](../../examples/run_static27.py)は既存N18の停止監視を変更せず再利用する。
workerと独立した `execution.toml` へ終了状態、CPU時間とpeak RSSを保存する。

メモリは二つの方法で記録する。

- workerの `Sys.maxrss()` は、Julia processの高水位をbytesで各phase前後に保存する。
  起動・JIT・DMRG・診断・保存を含み、phase固有の割当量や増分ではない。
- launcherの `resource.getrusage(RUSAGE_CHILDREN)` は終了後の子のCPU使用時間と
  RSS高水位を保存する。Darwinではbytes、LinuxではKiBをbytesへ変換する。
  過去の子の実績がある場合はRSSを差し引かず、取得不能を明示する。
  Python親や同時process tree全体のピークを測るものではない。

## 診断と完了条件

[worker](../../examples/validate_static27.jl)はDMRGが返したSz profileとvarianceを保存し、
trial checkpointを先に保存・再読込する。その後、全zz/+-相関、bond energy、
二つのcutのSchmidt量を追加計測する。
pm相関は複素のままreal/imagを分けて保存する。相関の測定は一度とし、θ=0のbond energyを
`Jz*real(zz[i,j])+Jxy*real(pm[i,j])` から算出して、和をMPO期待値と照合する。

観測整合性は以下で調べる。これらは状態の読出しの検証であり、基底状態精度の合格基準ではない。

- normと全磁化の誤差 `≤1e-12`。
- `diag(zz)=1/4`、`diag(pm)=1/2+Sz`、zzの実対称性、pmのHermitian性、
  固定磁化の恒等式 `Σj zz[i,j]=1.5 Sz[i]` の誤差 `≤1e-10`。
- bond energy和とMPO全energyの差 `≤1e-9`。
- cutごとのSchmidt平均と左領域Szの差 `≤1e-10`、Schmidt分散と
  `Σleft zz−(Σleft Sz)²` の差 `≤1e-9`、確率和の誤差 `≤1e-12`。
  異なるcutの平均やentropyが一致する条件は課さない。
- raw varianceは符号付きで保存し、負値は `100eps(Float64)*max(1,E²)` の丸め幅を確認する。
  各sweepのenergyと実測切断誤差が有限で、切断誤差が0から1の範囲内にあることを確認する。

全phaseが完了して整合性が通っても、状態は **単一χ・単一seedで精度未確立** と表示する。
varianceから計算した `sqrt(max(variance,0))` はMPOのエネルギー分散に由来し、
独立ED残差や基底エネルギーの誤差保証ではない。丸め床以下の値をゼロ誤差とは扱わない。
局所Krylovの収束フラグはbackendから取得していない。

最後のsweep間energy変化、variance、保持rankと捨てた重み、実時間・メモリを合わせて、
次にχ512比較、追加sweep、別seedのどれが必要かを判断する。期待する秩序や磁場幅に
近いことを採否条件にしない。途中停止・未達の状態も元の診断とともに保持する。

## 実行前の確認

独立担当がQ・幾何・観測式・保存・停止監視と、精度未確立の扱いを読み取り監査した。
Juliaでdriverの読込み、列磁化の集計、複素行列のreal/imag配列化、固定したサイズ・χ・
sweep数を確認した。最初の解析確認で `sum(fill(1/18,...))` の厳密等値を要求したため
丸め誤差でassertionが失敗した。確認式を絶対許容差 `1e-15` の比較に直して成功した。
研究の数値許容差やsolver設定は変更していない。

launcherは軽量な子processでCPU/RSS取得、Darwin/Linuxの単位換算、先行する子がいる
場合の取得不能表示、終了未確認、出力上書き拒否、worker記録の維持を確認した。
CLIのhelp・600秒超拒否・1thread設定も確認した。N18のdriver・launcher・参照helperは
過去の保存SHAと一致し、今回の追加によって変更されていない。

production・依存・通常テストは変更していないため、全 `Pkg.test()` は再実行していない。

## 実測と停止結果

出力先は `outputs/p4-static27-20260919`。終了したlauncherの
[execution原記録](data/p4_static27_execution.toml)、workerが停止前に原子的に保存した
[partial原記録](data/p4_static27_partial.toml)をbyte単位でコピーした。
両方のSHAと、未取得量を含む結果の解釈は[別のoutcome記録](data/p4_static27_outcome.toml)に保存する。

| 項目 | 実測・状態 |
| --- | --- |
| 外側の時間上限 | 600秒。停止猶予を含み、延長なし |
| SIGTERM送信 | 起動後598.000395秒 |
| 子の終了確認までのwall time | 598.239454秒、exit code `−15` |
| 子のCPU時間 | user 590.683291秒 + system 5.210483秒 = 595.893774秒 |
| 子のpeak RSS | 2,233,745,408 bytes = 2.080337524 GiB |
| 最後に保存されたphase | `DMRG_8_sweeps` の開始時 |
| 完了sweep数・完了状態 | 不明・未取得 |
| checkpoint・最終energy・Sz/bond profile・Schmidt量 | 未取得 |
| variance・実測切断誤差・観測整合性検査 | 未取得・未実行 |

launcherは子の終了後にprocess groupの残存processを掃除し、停止結果を保存した。
workerの `status="running"` は停止前の最後のsnapshotであり、現在も実行中という意味ではない。
実行終了の正本はlauncherの `status="timed_out"` と終了確認である。
出力先には原記録2ファイルだけがあり、取得できなかった量を補間・推定していない。

phase名は `DMRG_8_sweeps` だが、一つの `run_dmrg` 呼出しにMPO準備、JIT、
最適化sweep、energy・variance・Szの読出しが含まれる。
したがって「8 sweeps自体が終わらなかった」「varianceが支配的だった」のいずれも、
今回の記録からは特定できない。χ256の精度不足やメモリ不足とも判定しない。
Q3の完了エネルギーがなく、隣接Qも未実施なので、新しい磁場区間は得られていない。

## 保存後の監査と次の判断

独立担当が時間・CPU・RSS換算、終了状態、rawコピーのbyte一致、outcomeのSHAを確認した。
開始時記録のproduction 179ファイルと追加3ファイルは、停止後の現在のsource SHAと一致した。
workerの `finally` に置いた終了時sourceチェックは未実行であり、この外部照合と区別する。
予定した観測整合性検査は完了状態がないため未実行である。
今回確認した範囲は、実行前の軽量検証、実際の上限停止と資源記録、保存後の整合性である。

今回の10分の予備実験と結果評価を閉じ、**段階Bは未完了のまま**とする。
次の限定研究単位は、同じN27・Q3・χ256・seed11でまず2 sweepsを設定し、
sweepごとの進行・energy・実測切断誤差を記録する。
`measure_variance=false` で完了状態を先にcheckpointへ保存し、varianceを別に計時する。
これは次回の計画であり、今回のdriverを変更して実行したものではない。
次の時間上限も起動前に固定し、今回の予算は延長しない。
2 sweepsを収束条件とはせず、取得した状態と費用を見て追加sweep・χ・seed比較を選ぶ。
χ512、隣接Q、長さ拡張の走行は、進行と診断の費用を分けて測るまで保留する。

後続の[2-sweep進行・診断研究](p4_static27_progress.md)では、57.812秒でtrial保存と
全診断を完了した。実保持次元64の未収束状態であり、本記録の8-sweep試行の
停止結果・未取得量を変更するものではない。
