# N72既知CSL対照の零flux準備: 実行可能性監査

2026-09-22。状態: **ソース・保存記録の読取と計画の具体化のみ。N72の実装・DMRG・MPS収縮は未実行**。
対象は `J1=1, J2=J3=0.5, Q=0` の測定対照であり、最近接1/9のN54静的比較とは別である。
本記録は[9月20日の資源判断](p3_csl_resource_decision_20260920.md)を実装可能な単位へ分解する。
非零ポンプ校正は未達のままで、静的比較の数値workerと同時に走らせない。

## 結論

模型・DMRG・checkpointの共通APIはN72を扱える構造になっている。
一方、現行CSL研究workerはN36専用であり、config変更だけではN72を走らせられない。
必要なのは新しいHamiltonianやbackendではなく、**N72の幾何選択、段階別χ、診断分離、
比較全体のwall予算を持つ零flux専用driver**である。
まず新規複素seed11の最初の2 sweepsを8時間準備単位の第1段階として保存し、
実χ・時間・メモリ・保存の費用を測る。同じ初期化を再実行する独立pilotにはしない。

N36の費用外挿だけではN72・χ512の実行可能性も収束も確定できない。
今回は使用可能メモリを確認していない。
peak RSSの過去記録はあるが、現行supervisorにメモリ上限の強制機能はない。

## 固定する模型と幾何

`(Lx,Ly)=(6,4), N=72, Q=2ΣSz=0, Nup=Ndown=36, θ=0`。
軸OBC、wrap=`4a2`、順序はx、y、A/B/C、整数indexは `3(4x+y)+s`。
範囲外xへ出るbondだけを除き、端の結合変更、Zeeman項、chirality Hamiltonianは加えない。
J3は六角形対向頂点のみで、同距離の直線鎖上結合は含めない。

| 対象 | N72での選択 |
| --- | --- |
| 列Sz | x=0,…,5の各12 site |
| 中央領域 | x=2,3、site25:48 |
| 幾何cut | c=2,3,4、右領域はx≥c |
| 対応Schmidt bond | 24,36,48 |
| 全bond数 | J1:136、J2:128、J3:64、計328 |
| 中央の両端が含まれるbond | J1:40、J2:32、J3:16、計88 |
| CCW elementary triangle | 全44、中央に全頂点が含まれるもの12 |

cut2/3/4は両側に少なくとも2列を持つが、それだけでbulkとは呼ばない。
中央選択はsite番号の定数ではなく、`lattice.sites` のx座標から生成する。
bondは両端、triangleは全頂点が中央に属するものを選ぶ。
Schmidt bondは `3Ly*c` から生成し、幾何cutと同じ領域を測ることを検査する。
ここに示す数は文書化された座標式・bond数式の算術確認であり、N72の新規数値校正ではない。

## 再利用できるものと必要な変更

| 現行の実装 | 再利用・制約 |
| --- | --- |
| `src/extended_model.jl` | `kagome_j1j2j3_cylinder(6,4;J1=1,J2=.5,J3=.5)`、family分類は幾何に依存する汎用実装 |
| `src/dmrg.jl` | 明示Q、複素固定電荷random MPS、warm start、任意maxdim、実測切断、`measure_variance=false`、進行callbackを利用可能 |
| `src/checkpoint.jl` | trialの原子的保存、model/site/Q/source/runtime/environment照合、再loadを利用可能。旧N36を移行する必要なし |
| `src/observables.jl`, `src/schmidt.jl`, `src/chirality.jl` | Sz・bond・Schmidt・triangle診断はN36専用でない |
| `examples/run_static18.py`, `run_static27.py` | 汎用supervisor関数は明示wallを受け、process group終了・実行記録・終了後RSSを記録する。600秒制限は旧CLI側にある |
| `examples/research_csl.jl` | `(3,4)`、N36、cuts1/2、bond12/24、center13:24、3列集計、600秒上限が固定。各batchでH†Hを実行する |
| `examples/run_research.py` | 現行research環境を選べるが、CLIのwallは600秒以下。8時間単位の実行は不可 |
| `examples/research_static_split.jl` | solveと診断の記録を分離し、完了したsolve・入力hash・snapshotを厳密に結ぶ設計を参照できる。ただし旧NN backend・1-sweep契約のまま流用は不可 |
| `examples/audit_csl_saved_moment.jl` | 保存trialを変更せずH†Hを補足する設計を参照できる。固定された旧N36入力hashとbackendは流用不可 |

既存N36の履歴と専用driverを維持し、新しい零flux準備driverで必要な契約だけを実装する案が小さい。
新しいdriver・config・supervisorソースをallowlistで保存し、現行 `research` 環境とvendored
NDTensorsのidentityを全段階で一致させる。実行途中にpackage sourceや依存環境を変更しない。
段階ごとの出力は新規directoryにし、前段の完了したsolve記録とtrialのhashを結ぶ。
診断が未了のsolveに `integrity_passed` やflux-readyを付けない。

## 一つの準備単位の実施案

全段階を逐次実行し、Julia/BLAS各1 thread、数値process一つ。
初期化seed11、initial linkdim4、Krylov tolerance1e−11、dimension40、maxiter20、
noise0、cutoff0をN36対照から引き継ぐ。configured cutoffを誤差評価に使わない。

1. **χ128、2 sweeps。** freshなN72・Q0の複素MPSから開始。
   最初の段階のinclusive wall上限を1800秒とし、起動・compile・保存・再loadも含める。
   energy、各sweepの実測切断とenergy、実maxlinkdim、全Sz、norm/Q、保存overlap、
   時間・peak RSSを記録する。H†Hは行わない。
2. **χ256、2 sweeps。** 第1段階の完了・identity一致・trial再loadを確認して継続。
   同じ最小診断を保存する。第1段階がtimeoutしたとき、途中checkpointだけから自動再開しない。
3. **χ512、2 sweeps、最大でもさらに2 sweeps。** 同じ親系列を継続し、
   両χ512段階で中央bond/chirality・cut2/3/4のSchmidt・全Szを取得する。
   追加2 sweepsはこの比較単位に含む上限であり、終了後の同条件自動反復はしない。
4. **最高χで最後に完了したtrialの独立診断。** fresh MPOのenergy照合後に複素H†Hを測り、
   raw variance・虚部・roundoff scaleを保存する。全bond和とfamily和も照合する。
   元checkpointの `measure_variance=false` を書き換えず、別記録にhash付きで結ぶ。

総inclusive wallは28800秒、そのうち少なくとも3600秒を最終診断へ予約する。
したがってsolve・保存・途中のprofile測定を合わせて25200秒を越えない。
第1段階の1800秒もこの総量に含め、後段の開始時に残時間を引き渡す。
各段階のPython launcherを別processにすれば、既存のRUSAGE_CHILDREN high-waterの
「先行childがあると段階別peak RSSを得られない」という制約を避けて記録できる。
外側の集計は複数workerのRSS最大値を加算せず、各実行のhigh-waterとして保存する。
これは実装案であり、8時間を守る新規orchestratorはまだ存在しない。

第1段階で数値integrity・保存・資源に問題があれば、その記録を成果として止める。
正常なら、実測時間とメモリ余裕からχ512までの進行可否を判断する。
χ512に入る前に数値が悪くても、それだけで研究上の失敗ではないが、
残予算に入りそうにないsweepを始めて診断枠を使い切らない。

## 資源見積りの限界と終了判断

保存記録から再確認したN36・χ256の追加2 sweepsは461.576714秒。
別processの保存状態H†H補足のpeak RSSは3.639252 GiBだった。
`t ∝ Nχ^p, p=2,3` という既存の感度scenarioでは、N72・χ512の4 sweepsは
2.05145–4.10290時間となる。これは信頼区間・上限ではなく、solver反復・charge block分布・
起動・保存・観測費用を予測しない。メモリを同じ式で外挿して実行可と判断しない。

N36のfreshな最初の2 sweepsはconfigured χ128に対し実maxlinkdim64だった。
N72でも最初のstage時間をχ128飽和時の性能とみなさず、各段階の実χを必ず併記する。
χ128→256→512は前状態を成長させる準備であり、同一成熟度の性能benchmarkではない。

零flux開始の数値条件は既存のvariance≤1e−3 J²、最終実測切断≤1e−5、
sweep energy差≤1e−4 Jを維持し、norm/Q/Schmidt-density/bond和のintegrityも確認する。
これらに達しても単一seed・一幅・一長さでCSLや量子化を主張しない。
未達ならχ512追加反復やfluxへ自動拡張せず、実測費用から次の資源判断を行う。

既存 `csl_flux!` はvarianceなしtrialを受け取ると、新しいvariance-enabled零flux
最適化を行ってbaselineを作る。今回の分離した診断だけでtrialを既存flux APIの
accepted状態へ昇格できる設計ではない。最初の短い非零continuationも別予算・別実装判断となる。
ポンプの期待値やchiralityの符号・大きさは、準備の受理や枝選択の条件にしない。

## 実際に行った確認

- AGENTS、ROADMAP、flux設計、CSL設計・資源判断、上表の実装を読取。
- PythonでN72のsite/cut/Schmidt indexと文書化されたbond・triangle数を算術確認。
- 保存TOMLからN36の時間・RSS・最初の実χを再読し、上の費用scenarioを再計算。
- メモリ容量の読取は実行環境に拒否され、未取得。
- DMRG、ED、MPS load/収縮、Julia/package試験は実行していない。コード変更のない監査であり、
  既存の小系校正を繰り返していない。

未解決は、新規driverの実装とfocused検査、N72の実時間・メモリ、零flux精度、
複素seed依存、非零校正、長さ・幅依存である。
