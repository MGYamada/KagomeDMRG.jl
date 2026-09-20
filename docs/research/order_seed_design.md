# 固定Qの9／27サイト周期の初期試行状態

2026-09-19。対象はNN等方模型の `Q=N/9`、零fluxでの初期状態比較。
実装は [`examples/order_seed_tools.jl`](../../examples/order_seed_tools.jl)。
これは明示的な **density/bond-biased trial** の構築であり、文献のVBCの
結合模様・波動関数の再現や、最適化後のVBCを主張するものではない。
乱数seedの変更を秩序候補の構築とは扱わない。

## 構成と電荷

各原始単位胞 `A(x,y),B(x,y),C(x,y)` のうち二サイトを選び、
サイト順 `i<j` で次の状態を作る。

```
|pair(phi)> = (|Up_i Dn_j> - exp(i phi)|Dn_i Up_j>)/sqrt(2)
phi = +0.17 (odd seed), -0.17 (even seed)
```

残る一サイトは `2Sz=+1` または `-1` の積状態とする。
pairの合計Szは厳密に0だが、非零phiではsingletとは呼ばない。
`<S+_i S-_j>=-exp(i phi)/2` は非零の虚部を持つので、複素性は
全波動関数のglobal phaseだけではない。これは初期試行状態へのバイアスであり、
Hamiltonianにchirality項・pinning場・追加結合を加えていない。
最適化後に残るchiralityや密度模様は別に測り、この初期値を証拠として流用しない。

`seed` は乱数生成には使わず、原点
`(ox,oy)=(seed mod 3, floor(seed/3) mod 3)` とphiの符号だけを指定する。
以下で `u=x-ox,v=y-oy,c−=(u−v) mod 3,c+=(u+v) mod 3` とする。

| 規則 | period9 | period27 |
| --- | --- | --- |
| spectatorの整数電荷 | `c−=0,1,2` に `+1,+1,-1` | 同じ |
| pairの選択 | `c−=0,1,2` に AB,BC,AC | `c+=0,1,2` に AB,BC,AC |
| 原始並進（a1,a2成分） | `(1,1),(-1,2)` | `(3,0),(0,3)` |
| supercell内のサイト数 | 9 | 27 |

3個の原始単位胞あたりspectatorがUp二つ・Down一つなのでQ=1。
period27の9個の原始単位胞ではUp六つ・Down三つでQ=3となる。
`Ly%3=0` とし、period27では完全な模様を含む `Lx%3=0` も要求する。
したがって `(Lx,Ly)=(3,3),N27,Q3` と `(6,3),N54,Q6` の両方で
目標Qを厳密に満たし、同じwrap `3a2` と同じ切り方の開放端を使える。

period9の並進は `dx−dy=0 mod 3` で、面積3の√3×√3周期となる。
period27は同時に `dx+dy=0 mod 3` が必要である。
3を法として2が可逆なので `dx=dy=0 mod 3` となり、面積9の3×3周期を持つ。
有限円筒のx方向にはOBCがある。ここでいう周期は背景となる無限模様の周期であり、
有限Hamiltonianのx並進対称性を仮定するものではない。

## APIと保存

```julia
include("examples/order_seed_tools.jl")
lattice = kagome_cylinder(3, 3)
sites = spin_sites(lattice)
psi, metadata = make_order_seed(lattice, sites; kind=:period9, Q=3, seed=11)
# 同じsitesをrun_dmrg(...; sites, psi0=psi, Q=3, ...)へ渡す。
```

canonicalな `x,y,A,B,C` 順序と整数QN `Sz=2*physical Sz` を検査する。
pairの準備にはU(1)を保存する二spin unitaryを使い、非隣接AC pairの処理後も
サイト順を戻す。各pairは一原始単位胞に収まり、異なる単位胞は積になっている。
そのため厳密なSchmidt rankは最大2で、構成時の `cutoff=0,maxdim=2` は
このansatzの近似切断ではない。DMRGのmaxdimやcharge sectorを2へ固定してはいけない。

metadataは許可リスト方式で、geometry/Q/seed、並進・位相・pair/spectatorの
具体的規則とサイト列、解析的なSzと無重みの零flux `Si·Sj`、測定norm・保持次元を記録する。
モデル・solver・source/dependency hash・親checkpointの由来は研究driver側が保存する。
環境変数やローカルネットワーク情報は取得・記録しない。

## 測定・比較指標

最適化前の解析値は、spectatorで `Sz=±1/2`、pair上では `Sz=0`。
選択pairで `<Si·Sj>=−1/4−cos(phi)/2`、それ以外のbondは
`<Si·Sj>=<Sz_i><Sz_j>` となる。これは初期構成の校正値であり、
DMRGの収束先を選ぶ受理基準には使わない。

最適化後には同じQ・模型で全energy、variance、実測切断誤差、末尾sweep変化とともに、
次の量を比較する。

- 全サイトSzと全NNの無重み `Si·Sj`。J=1の本研究では後者はbond energyに等しい。
- 列ごとの平均・RMS変調。端と中央を分け、N27には十分なbulk領域がないと明記する。
- 各sublatticeの平均を除いた密度と、同じ幾何学的bond種別の平均を除いた結合模様。
  両trialに共通する「up triangle内のpairを強くした」バイアスを、周期の違いと混同しない。
- 解析テンプレートとの規格化内積およびRMS差。許される原点の平行移動も比較し、
  ドメインの位置差と異なる秩序を分ける。テンプレートへの近さをenergy収束条件にしない。
- cell座標に対する密度・bondの複素Fourier振幅。period9の `±(2π/3,−2π/3)` と、
  period27で必要になる独立な3×3波数格子を区別する。OBCのx方向のFourier解析は
  模様の記述であり、保存された運動量量子数ではない。
- 縦・横spin相関、Schmidt量、零flux chirality。密度一様性やchirality単独でCSLとはしない。

N27→N54は同じterminationと原点規則を保ち、全列のprofileを保存する。
N27の中央一列とN54の中央列を比較しても、長さ依存が十分小さい・端から十分遠い
ことは別途確立する必要がある。これらのtrialはいずれも低い初期entanglementと
強い局所バイアスを持ち、競合VBCを網羅していない。

今回の集計実装は全Sz/bond profile、列と明示した領域の統計、
副格子平均を除いたSzの3×3波数格子、保存Czzからの連結相関の領域統計を含む。
Fourier振幅は各副格子のサイト数で割り、原始単位胞の座標で位相を定義する。
N54での3×3波数格子は3-cell周期の成分を取り出す部分的な記述であり、
6列に沿う全波数・全profile変動を尽くすものではない。
相関が未記録の場合は欠損として残す。原点を走査したテンプレート照合、
bond種別ごとのFourier解析、NN保存状態のchiralityは今回の集計では未実施である。

## 検証範囲

`order_seed_calibration()` は4×4のpair unitary、3を法とする並進stabilizer、
N27/N54それぞれの両trialについて、厳密Q、norm、χ≤2、全Sz・NN bondと
複素pair coherenceの解析値との一致を確認する。DMRG・ED対角化は行わない。

2026-09-19、Julia 1.13.0、Julia/BLAS各1 threadで実行して成功した。
外側に60秒のprocess timeoutを設け、起動・import込み30.729561秒、
関数本体27.460337秒だった。最初のsandbox内launcher試行はlockfile作成権限により
Julia起動前に停止したため、許可後に同じ時間制限付きコマンドを再実行した。

| ケース | 全Szの最大誤差 | NN bondの最大誤差 | 複素pair coherenceの最大誤差 |
| --- | ---: | ---: | ---: |
| N27Q3・period9 | 1.11e−16 | 2.00e−15 | 1.35e−15 |
| N27Q3・period27 | 1.11e−16 | 2.00e−15 | 1.35e−15 |
| N54Q6・period9 | 1.11e−16 | 3.89e−15 | 2.65e−15 |
| N54Q6・period27 | 1.11e−16 | 3.89e−15 | 2.65e−15 |

全四ケースでnorm=1、保持次元2、整数QNは指定Qに一致した。
4×4 gateのunitarity誤差は1.59e−16、指定pair状態との差は1.39e−17。
3を法とした並進stabilizerはperiod9で `(0,0),(1,1),(2,2)`、
period27で `(0,0)` のみとなり、想定した最小周期を確認した。
両状態間の最大差は両サイズともSzで0.5、NN bondで0.7427923834547798であり、
乱数seed名だけが異なる同じ状態ではない。

この検証でDMRG・ED対角化・全package suiteは実行していない。
χや初期バイアスを変えた最適化後の安定性、文献VBCとの模様の照合、bulk秩序は未判定。
構成コードはこの校正後に変更していない。研究workerが利用するときは、
このhelperのsource hashも当該走行の由来として保存する。
