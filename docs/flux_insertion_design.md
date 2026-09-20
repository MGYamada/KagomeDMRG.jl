# スピン flux 挿入の物理・数値設計

作成日: 2026-09-19。[ロードマップ](../ROADMAP.md)の仕様詳細。
以下の設定値は検証を始めるための提案であり、計算結果や収束保証ではない。
P0 と P1 の小系基準計算は実装済みで、実測結果は
[検証記録](research/p0_reference_validation.md)に分けて記載する。
完了点の checkpoint と Schmidt charge 診断を追加した。
continuation の受理・棄却・刻み半減・復元を追加し、対照系で初期検証した。
既知 CSL のポンプ・研究サイズの枝追跡・相同定は引き続き研究計画である。

## 1. 測定する U(1) 応答

`ħ=1` とし、局所演算子 `S+` が運ぶ物理的な Sz を 1 とする。
Hall 係数はここでは無次元化して

\[
q_\ell=t^{\mathsf T}K^{-1}\ell\pmod 1,\qquad
\sigma_{xy}^{S^z}=t^{\mathsf T}K^{-1}t,
\qquad \Delta S_R^z(2\pi)=\sigma_{xy}^{S^z}
\]

と定義する。最後の等式は、向きの規約を揃え、バルクが適切に gapped で、
目的の枝を追跡できる極限での予測である。
`j_x^z = κxy E_y^z` の係数を用いる規約なら `κxy=σxy^Sz/(2π)`。
分数スピンは局所励起を加える自由度のため modulo 1 だが、数値で測る移送量は
最初から modulo 1 に丸めず、連続な実数として保存する。
K 行列による応答の枠組みは
[Wen–Zee](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.46.2290)を参照。

ユーザーが指定した Hall 応答を持つ D(Z₃) 候補は

\[
K=\begin{pmatrix}0&3\\3&0\end{pmatrix},\quad
K^{-1}=\begin{pmatrix}0&1/3\\1/3&0\end{pmatrix},\quad t=(1,1)^{\mathsf T}.
\]

したがって `q_l=(l1+l2)/3 mod 1`, `σxy^Sz=2/3`,
`c₋=signature(K)=0`, `|det K|=9`。
熱 Hall 応答を決める c₋ がゼロでも、保存 Sz の Hall 応答は非零になり得る。
これは候補の有効理論上の計算であり、最近接 kagome 模型での実現を意味しない。
追加の可逆相の積み重ねは今回の代表例には含めない。

2π flux が挿入する anyon の型は、向きの符号を除いて `ℓflux=t mod K Z²`。
今回の表の 3 例では、最小の正整数 n で `n K⁻¹t` が整数ベクトルとなる値が 3。
よって `0→2π→4π→6π` を追跡し、非零 Hall の 2 例では累積移送 `0,2/3,4/3,2` を検証する。
ゼロ Hall の例 `t=(1,0)` でも flux は非自明な anyon を挿入し得る。

**6π で戻るのはバルクの topological sector。有限 cylinder の端には移送された
整数スピンが残り得るので、全波動関数や端密度の復帰を要求しない。**
D(Z₃) の 9 セクターのうち、一種類の flux が辿る軌道は 3 セクターにとどまる。
3 周期が見えても縮退数 3 と結論しない。

## 2. モデルと磁化規約

\[
H(\theta)=\sum_{b=(i,j)}\left[
J_b^z S_i^z S_j^z+
\frac{J_b^{xy}}2\left(e^{i A_b(\theta)}S_i^+S_j^-
+e^{-i A_b(\theta)}S_i^-S_j^+\right)\right]
-h\sum_i S_i^z.
\]

本研究の主対象では最近接ボンドに `Jz=Jxy=J>0`、それ以外はゼロ。
XXZ、遠距離結合、chirality 項は既知相の対照実験や追加研究として区別する。
一つの物理ボンドは一度だけ列挙し、その中で共役な 2 項を作る。

この規約で [Fang et al.のEq.(1), Fig.1](https://arxiv.org/pdf/2306.09563)は
`0.35<h/J<0.42` に1/9 plateauを報告している。`h/J=0.38` を比較用の代表点にするが、
有限円筒の磁場境界をこの値に合わせない。まずθ=0の静的NN研究を優先し、
既知CSLのポンプ校正を並行する。文献の競合候補・予算・段階別の判定は
[静的研究方針](research/p4_static_plateau_strategy.md)に記載する。

\[
M_{\rm sat}=N/2,\quad M=N/18,\quad Q=2M=N/9,
\quad N_\uparrow=5N/9,\quad N_\downarrow=4N/9.
\]

- 既定の 1/9 では `N % 9 == 0` を要求する。N が偶数なら M は整数、奇数なら半整数でよい。
- ITensor の spin QN `Sz` の整数ラベルは物理的な `2Sz`。
  目標値は `QN("Sz",N÷9)` であり、物理 Sz=N/9 と取り違えない。
  これはゼロ磁化の交互 Up/Dn 初期状態からは到達しない。
  [公式 QN DMRG チュートリアル](https://docs.itensor.org/ITensorMPS/stable/tutorials/QN_DMRG.html)。
- 固定 M では `-hM` は定数なので、flux 走行中の h 調整は不要。
  plateau 自体の確認は別の磁化 sector も用いる。
- 隣接 sector の全交換エネルギー `E0(M)` から得る
  `h−=E0(M)−E0(M−1)`, `h+=E0(M+1)−E0(M)` は最初の目安。
  より遠い sector を含めた凸包、端の磁化、サイズ依存を調べて plateau を確認する。
  磁化 plateau があっても、同じ M 内の中性励起まで gapped とは限らない。
  有限 cylinder の最小励起が端に局在する可能性もあるため、励起の空間分布を確認し、
  端の gap をそのまま bulk gap としない。

Q表記では差分をQ0±2で取る。遠方Qも比較した場合は
`max_{Q<Q0} 2[E0(Q0)−E0(Q)]/(Q0−Q)` と
`min_{Q>Q0} 2[E0(Q)−E0(Q0)]/(Q−Q0)` が、比較したsector内での下限・上限となる。
変分エネルギーの差は厳密な境界の誤差限界ではなく、χ・初期状態等による変動を併記する。
未探索sectorによる飛び越しは未判定とし、追加磁化の端への局在とbulkの応答を分ける。

この静的比較の最初の実例として、NNのN18・θ=0・Q=0,2,4を独立ED/DMRGで照合した。
比較したsector内の区間は `0.2585201582<h/J<0.4869821958`。
二列で内部bulkはなく、遠方Qは未探索なのでbulk plateauの検証ではない。
数値・列磁化・精度・10分の外側時間監視は[段階Aの研究記録](research/p4_static18_validation.md)に分ける。
続くQ6の1点追加では、`(E6-E2)/2=0.6027813016 J`となり、
比較集合`{0,2,4,6}`の区間も同じだった。新規Q6は残差2.05e−12以下でEDと一致し、
起動込み101.142秒で完了した。[出典と精度の確認](research/p4_static18_q6_validation.md)。
Q≥8は未探索として残す。保存済みQ2の近接中性二重項については、密度・bondの
両方に非零の遷移を測定した。基底との差0.00186467503 J、円周運動量±2π/3、
局所分散の和への比率4.01%／8.40%は有限系の構造診断に限る。
[由来・残差・並進・和則と解釈](research/p4_static18_neutral_validation.md)を別記した。
続くN27・Q3・χ256の代表試行は時間上限で停止し、完了状態・観測量は未取得。
資源量のみを[N27の停止記録](research/p4_static27_pilot.md)に保存し、新しい磁場区間や
有限χの精度結果としては扱わない。
その後の2-sweep試行は進行・保存・別計時診断を完了したが、設定上限256に対し実保持次元64、
variance 0.31546 J²の初期状態である。[新しい研究記録](research/p4_static27_progress.md)と分け、
この完了状態も磁化境界・bulk収束の証拠とはしない。
同じtrialから累積4 sweepsへ進めると、実保持次元256、variance 0.008498 J²となり、
列磁化も大きく変わった。[固定χのsweep比較](research/p4_static27_refine.md)として保存し、
まだ収束した秩序や磁場区間とは解釈しない。

共通 API は `target_sector(N; Q=...)`、`initial_mps(...; Q=...)`、
`run_dmrg(...; Q=...)` で別の固定 sector を明示できる。
整数 `Q` の範囲 `−N≤Q≤N` と N との偶奇一致を要求し、`Nup=(N+Q)/2` とする。
省略時は warm start があっても 1/9 のままで、入力 MPS の Q 不一致を拒否する。
保存・再開・continuation は検証済みの保存状態または始点の Q を継承し、
任意に渡した Q は一致を要求する。途中で Q を変える操作ではない。
既存の schema 1 の model/state 電荷欄を使い、baseline・MPS・電荷欄を照合する。
旧ソースの snapshot は従来どおり source identity 不一致で拒否し、自動移行しない。

## 3. Cylinder とゲージ

初期実装では `x=0,…,Lx−1` を開境界、`y=0,…,Ly−1` を周期境界とする。
`a1=(1,0)`, `a2=(1/2,√3/2)`、副格子座標 `A=0, B=a1/2, C=a2/2`、
`N=3LxLy` を採用する。Ly は原始単位胞数であり、文献の YC 表記とは自動的に同一視しない。
wrap vector は `Ly*a2`。MPS の順序は x、次に y、最後に A,B,C で、
`index=3(x Ly+y)+s`（s=1,2,3）とする。端は範囲外の x へ出る結合だけを除く。
幾何学的 cut c は完全な単位胞列の間に置き、右領域を `x>=c` とする。
これは Cartesian の水平座標だけで分類する切断ではない。

競合する√3×√3と3×3 VBCを許す比較には `Ly%3=0` を用いる。
前者の周期基底を `a1+a2,−a1+2a2` と取ると和は `3a2`、後者は `3a1,3a2` である。
Lxは開境界なので3の倍数を必須としないが、同じ端の `Lx=3→6` を長さ比較の候補にする。
N27は3×3の一周期であり、周期整合性だけで十分なbulkがあるとはみなさない。

実装した最近接ボンドのテンプレートは各 `(x,y)` に対し
`A-B`, `A-C`, `B-C` の胞内 3 本と、
`A(x,y)-B(x−1,y)`, `A(x,y)-C(x,y−1)`,
`B(x,y)-C(x+1,y−1)` の 3 本。
x が範囲外のボンドを除き、y を wrap する前の winding を保存する。
この切り方で Ly≥3 のとき bond 数は `(6Lx−2)Ly`、バルク配位数は 4。
Ly=1,2 の特殊な周期同定は初期の研究対象から外す。

別名の `kagome_j1j2j3_cylinder` は六角形内のJ2と対向頂点間のJ3を加える。
第三近接距離だけでJ3を選ばず、直線鎖上の同距離結合を除く。
端は無限格子で定義した各bondの両端だけで切り、全六角形が残ることは要求しない。
全familyに同じ位相規約を適用し、零結合も配列には残す。
独立幾何・EDの構成と照合は[拡張模型の検証](research/p3_extended_model_validation.md)に分ける。

ボンドは少なくとも `(i,j,Jxy,Jz,wy)` を保持する。
`wy` は i から j への向きで周期境界を跨ぐ符号付き巻き数で、向きを反転すると符号も反転する。
MPS index の大小は winding の判定に使わない。

seam gauge では `A_b(θ)=wy*θ`。これを初期実装の既定とする。

- `SzSz` の係数は twist しない。境界を跨ぐ斜めボンドも漏らさない。
- 三角形・六角形など可縮ループの有向位相和は 0、円周を一周するループは ±θ。
  各 plaquette に磁束を入れる問題と混同しない。
- `H(θ+2π)=H(θ)` が同じ基底で成り立つ。
- 局所回転 `U=exp(i Σi χi Szi)` に対し、
  `U H(A) U† = H(Aij+χi−χj)`。
  uniform gauge との行列・スペクトル照合を小系で行う。
- uniform gauge は `ηi=yi+δsublattice,C/2`、`χi=−θηi/Ly` とし、
  `Auniform=wy θ+χi−χj` を実装して seam gauge との行列照合に使う。
  θ=0 の 18 サイト ED では円周並進と低準位の運動量を照合した
  （[監査記録](research/data/p2_ed18_translation.toml)）。
  一般の θ でのゲージ補正を含む並進と MPS の運動量診断は未実装である。
  2π 後の比較には大きなゲージ変換が必要で、未変換の MPS overlap を比較しない。
- 長距離結合や多体項を追加する場合も、同じゲージ規則から twist を導出する。
  chirality seed を残す場合にはその境界項も含める。

## 4. ITensor による最小構成

現在の API では `using ITensors, ITensorMPS` を用いる。
`siteinds("S=1/2",N; conserve_qns=true)` を一度作り、全 θ で再利用する。
対象 sector の Up/Dn 個数を満たす初期状態を複数用意し、複素 MPS を扱う。
以下は項の規約を示す設計用の抜粋。
公開 API は `twisted_exchange_mpo(sites,lattice,theta; gauge=:seam,hz=nothing)` であり、
実行可能な最小例は [`examples/validate_small_system.jl`](../examples/validate_small_system.jl)。

```julia
using ITensors, ITensorMPS

function twisted_exchange_mpo(sites, bonds, theta)
    terms = OpSum()
    for b in bonds
        a = (b.Jxy / 2) * cis(b.wy * theta)
        terms += b.Jz, "Sz", b.i, "Sz", b.j
        terms += a, "S+", b.i, "S-", b.j
        terms += conj(a), "S-", b.i, "S+", b.j
    end
    return MPO(ComplexF64, terms, sites)
end
```

固定 M のため定数となる Zeeman 項はこの抜粋から省略した。
一般の θ では必ず複素 Hermitian Hamiltonian として扱う。
API は [OpSum/MPO](https://docs.itensor.org/ITensorMPS/stable/OpSum.html)、
[MPS](https://docs.itensor.org/ITensorMPS/stable/MPSandMPO.html)に基づく。
P0 で Julia 1.12.7、ITensors 0.9.31、ITensorMPS 0.4.1 を使って照合し、
解決した依存関係を当時のroot `Manifest.toml` に保存した。
現在の研究用環境は `research/Project.toml` に分離し、Julia 1.12用を
`research/Manifest.toml`、1.13用を `research/Manifest-v1.13.toml` に保持する。
研究例は `--project=research` で起動し、実際に選択された環境のhashも保存する。
rootの `Project.toml` はライブラリ開発用にも使えるが、そのローカルmanifestはGit管理外とする。
過去の数値・checkpoint・環境hashは当時の記録として保持し、新配置へ書き換えない。

初期は各 θ の MPO を構築し直し、`dmrg(H, psi_previous; ...)` で最適化する。
実際の切断誤差、energy/variance、各観測量の収束を測り、`cutoff` の設定値だけを
誤差として報告しない。Krylov solver の収束も確認する。
現実装の `run_dmrg` は単一点を最適化し、両 half-sweep の実測切断誤差の最大値、
再計算した最終エネルギー、`⟨H†H⟩−⟨H⟩²` を返す。
noise が非零なら切断誤差は摂動した密度行列のものであり、波動関数の捨てた確率と同一視しない。
upstream API が局所 Krylov の convergence info を公開しないため、それを確認済みとはしない。
小系では独立した全系 residual も検証する。warm start の入力はコピーし、site indices と Q を照合する。
任意の `progress_callback(event)` は工程の開始・終了、guard通過後のbond更新、sweep完了を
scalar-onlyのNamedTupleで通知する。callbackへMPS/MPOは渡さず、設定やcheckpointにも保存しない。
局所・sweep末のenergyはsolver報告値であり、最終MPO期待値と区別する。
通知の経過時間はcallback処理も含む。varianceを省略した場合、その工程の通知は出さない。
静的な資源測定では `measure_variance=false` の完了trialを先に保存し、varianceを別診断として
計測できる。これは全点varianceを要求するcontinuationの受理条件を緩める変更ではない。
実装の小系照合と研究走行は[N27の進行記録](research/p4_static27_progress.md)に分ける。
従来の upstream NDTensors 0.4.31 では、等しい Schmidt 重みを異なる QN sector で切る境界に
空状態と誤った切断誤差を返す事例がある。近接した重みでも、非零状態を残しつつ
切断誤差を過小報告する条件を四サイトの解析状態で確認した。observer は正準中心の
norm がゼロ・非有限の場合と、報告 spectrum の長さが保持 bond 次元と異なる場合に停止する。
現環境は局所版 NDTensors `0.4.31+1` を固定し、全 sector の保持状態を一度選んで
block rank と spectrum に共用する。`maxdim` は上限であり、同値の重みは block 座標・
block 内 index の順に選ぶため縮退空間全体を保持するとは限らない。
正の `min_blockdim` も全体の上限に算入し、両立しない制約は例外とする。
SVD と Hermitian density の経路を修正し、一般非 Hermitian eigen の経路は従来通りとする。
再現条件は [P2 の監査記録](research/p2_continuation_validation.md)、
修正・独立校正は [P1 切断検証](research/p1_qn_truncation_calibration.md)を参照。
必要なら固定 3 成分の MPO 和を使う最適化を比較する。
ITensorMPS は MPO の配列による和を受け取れるが、速度向上は実測で判断する。
[DMRG API](https://docs.itensor.org/ITensorMPS/stable/DMRG.html)。

## 5. Adiabatic continuation と枝の診断

小さい flux ごとの静的 DMRG は実時間発展ではなく、変分状態の continuation。
各 θ で最小エネルギー状態を探し直すことも、前段の MPS を渡すことだけも、
目的の topological branch を保証しない。
有限系の一意な基底状態を常に完全収束させると、2π で元の状態に戻って
移送が消える場合がある。バルクに断熱的な応答と、微小な端・セクター混成の
avoided crossing をどう通過したかを区別する。

1. θ=0 の候補を複数の固定 Q 初期状態から求める。エネルギーだけでなく、
   バルク密度、bond energy、chirality、EE を比較し、異なる候補を保存する。
2. θ は unwrapped の値で保持し、まず `δθ=π/12` 程度から試す。
   受理済み状態を保護した上で、同じ site indices と前段状態から次の θ を計算する。
3. 正規化 overlap、バルクの局所量、EE/Schmidt charge、最大切断誤差、
   エネルギー・残差を調べる。overlap は系の長さに依存するため、
   全サイズ共通の固定閾値だけでは採否を決めない。
4. 急変・未収束なら checkpoint に戻り刻みを半減し、必要に応じ bond dimension、
   sweeps、Krylov 精度を増す。解消しなければ枝の追跡失敗として記録する。
5. `0→2π→4π→6π`、負の θ、同じ枝を戻る走行を比較する。
   逆走の不一致は履歴依存・枝飛びの診断であり、反対 chirality の枝への切替とは区別する。
6. 適応刻み以外に `δθ/2` の全走行も比較する。
   期待する 2/3 に近づくことを収束条件や枝の選択基準に使わない。

現在の `continue_flux` は上記の一部を実装する。seam gauge、noise=0、全点の variance、
少なくとも二 sweep を必須とし、`FluxPolicy` に有限の閾値を明示する。
初期状態を含め、overlap の絶対値（振幅）、指定サイトの密度変化、Schmidt entropy・
左 Sz 平均の変化、最後の sweep の実測切断誤差、variance、最後の二 sweep の energy 差と
最終期待値との差を判定する。負の variance は `100eps(Float64)*max(1,E²)` の
丸め範囲と明示した variance 閾値の両方を満たす必要がある。
累積移送の左右和、Schmidt 電荷との一致、指定 cut 間の差も記録・照合する。
全サイズに共通の overlap 閾値や期待 pump 値は内蔵していない。

各試行は最後の受理 checkpoint を読み直して始め、棄却時は刻みを半減する。
最小刻み・試行数制限まで解消しなければ `unresolved` として記録する。
数値的に有効な棄却状態は trial snapshot に残す。不正な診断・状態は受理せず、
型と分類した理由だけを保存する。`completed` は指定した有限の診断を通過した意味であり、
物理的な断熱性・相の同定を意味しない。零 flux の縮退では、正しい基底状態でも
刻みを減らして overlap が 1 に近づくとは限らない。小系の実例と制限は
[検証記録](research/p2_continuation_validation.md)に示す。
18 サイト相互作用系では、半分の刻みをあらかじめ設定した別走行との比較を行った。
bond dimension や sweeps の自動増加、uniform gauge で共通規約に変換した overlap、
研究サイズでの刻み・サイズ収束の比較は今後の課題である。

本模型は h≠0 でも θ=0 の Sz 基底で Hamiltonian が実数なので、複素共役に関する
対称性が残る。実数演算だけでは複素の枝を探索できない。
複素初期状態や微小な chirality seed を診断に使い、seed の符号・強さ・除去後の安定性を調べる。
seed を残した結果は別 Hamiltonian の結果として表示し、最近接模型の結果と混ぜない。

枝が安定に追えない場合は、複数の低エネルギー候補、端の制御、iDMRG を検討する。
実時間 TDVP による緩やかな flux ramp は独立の追加検証になり得るが、
bulk gap と小さな端の splitting に対する時間尺度の検討が必要であり、初期必須機能にはしない。

## 6. 観測量

中央付近の幾何学的な切断 c に対して

\[
P_c(\theta)=\sum_{i\in R_c}
\left[\langle S_i^z\rangle_\theta-\langle S_i^z\rangle_0\right]
\]

を右向きの移送量とする。固定 Q のため全系の変化はゼロだが、左右それぞれの期待値は
分数だけ変化できる。0-flux の実測 profile を基準にし、端の初期磁化を 1/18 と仮定しない。

- 左の変化と右の変化の和がゼロであることを確認する。
- 中央の複数 cut で `Pc` が一致し、密度変化が主に端に蓄積するかを確認する。
  バルク全体の再配列や domain wall の移動を量子化輸送と取り違えない。
- 列ごとの密度変化、中央の bond energy、scalar chirality を保存する。
- 対応する MPS cut で `Σα λα² SzL,α` を計算し、左領域の直接和と照合する。
  QN label の向き・符号・定数 offset を既知の積状態で校正し、物理 Sz に換算する。
  sector ごとの切断確率と entanglement spectrum を追跡する。
- `⟨H²⟩−⟨H⟩²` は可能なサイズ・代表点で評価する。
  全 θ で高価なら頻度を制御するが、energy の安定だけを正しさの根拠にしない。

`bond_energies(psi,lattice,theta;gauge)` は各bondの交換エネルギーを同じ位相規約で返す。
和に `−dot(hz,sz_profile(psi))` を加えると全エネルギーとなる。
この初版は全相関行列を再利用するためO(N²)の保存量を使い、大系での費用は未計測。
`oriented_triangles(lattice)` は周期画像を持つCCWの基本三角形を返し、
`triangle_chiralities(psi,lattice,theta;gauge)` はその順に正規化期待値を返す。
`theta=0` では `Si·(Sj×Sk)`、非零fluxでは局所移送で定義したdressed演算子とする。
各頂点の画像をmとしてseamの角は `alpha=−m*theta`、uniformでは既存の
`gauge_angles`を加える。三角形辺の位相と一致し、状態と観測演算子を共に変換する。
非零fluxで継ぎ目を跨ぐ裸の三スピン積とは異なる。
向き・符号・独立三spin対照は[chiralityの仕様・校正記録](research/p3_chirality_validation.md)にまとめる。
初版は三角形ごとにMPOを収縮し、大系での測定費用は未計測。
静的なSz・bond patternの比較は先に進められる。固定Qでは横磁化の一体期待値がゼロなので、
磁気秩序の検討には縦相関と `⟨S+i S−j⟩` も使い、端で誘起された変調の長さ依存を調べる。

`schmidt_diagnostics(psi,b)` はこのうち MPS prefix `1:b` の確率・絶対電荷・
entropy・左物理 Sz の平均と分散を実装した。幾何学 cut `c` は `b=3Ly*c` に対応する。
MPS のコピーを正準化して切断なし SVD を行い、左テンソル群の flux の和から
境界 link の向き付き QN を引いて整数 `q_left` を得る。最後に 2 で割って物理 Sz にする。
実空間密度を使って offset を合わせることはしない。
積状態、非零全電荷の独立スピン基底、link QN の定数 shift・arrow 反転、
27 サイトの合成状態での複数 cut の左右移送により校正した。
これは読み出しの校正であり、Hamiltonian の flux 応答の検証ではない。

静的 DMRG の各点で得た平衡電流を、架空の時間刻みで積分して pump としない。
特に `∂E/∂θ` は円周方向の twist に共役な応答であり、軸方向の移送量そのものではない。
時間積分による検証をするなら、実時間発展と切断を横切る電流演算子を別途実装する。

## 7. 検証とサイズの拡大

| 検証 | 対象と判定 |
| --- | --- |
| 代数・格子 | 二サイト交換の固有値、配位数・bond 数、ボンド反転と共役、U(1) 保存、Wilson loop |
| 小系 ED | 例: `Lx=1,Ly=3,N=9` は sector 内 dense、`Lx=2,Ly=3,N=18` は sparse ED。非零 θ で MPO と照合 |
| ゲージ | seam/uniform のユニタリ同値、2π 周期、π 以外も含む一般の θ での Hermiticity |
| 負の対照 | Jxy=0 と局所縦磁場で一意に選択した固定 Q 積状態で、flux 非依存かつゼロ pump。これは読み出し検証で相同定ではない |
| 正の対照 | 既知の kagome CSL における 2π 当たり ±1/2 の移送と sector flow |
| 収束 | bond dimension、刻み、sweeps、長さ、円周、端の処理、初期状態を変えた誤差評価 |
| 保存・再開 | 同一 checkpoint からの再開と連続走行で、受理した θ 列と物理量が一致 |

実装した N=9、Q=1 の sector は 126 次元で、θ=0 の最近接模型では基底状態が二重縮退する。
DMRG 状態の ED 基底空間への射影と全系 residual を確認し、射影後の同じ状態の密度・相関を比較する。
θ=0.37 では一意な基底状態との比較も行う。この有限系の縮退をトポロジカル縮退とは解釈しない。
N=18、Q=2 は `binomial(18,10)=43758` 次元である。
独立な sparse BlockLanczos 参照を追加し、複素初期 block、収束情報、全行列 residual、
直交性を確認する。9 サイトの dense 対角化で縮退空間も照合した。
18 サイトの θ=0,0.185,0.37 で低い三準位を計算した
（[数値記録](research/data/p2_ed18_reference.toml)）。有限個の Ritz 対の収束だけでは
全基底空間の完全性や bulk gap を保証しない。
9 サイト系には内部軸方向 cut がなく、その照合から軸方向 pump の検証はできない。

18 サイトでは χ=512 の `0→±0.37→0` を、事前に決めた二種類の刻みと二つの seed で
比較した。四経路の延べ 24 保存点は独立 ED と一致し、実空間と Schmidt 移送も整合した。
χ=128 の初期点は棄却され、χ=256 は零 flux の sweep 中に切断 guard が停止した。
その後の局所版による [有限 χ 校正](research/p1_qn_truncation_calibration.md)では、
χ=256 の停止は解消した。θ=0,0.37 の 10 trial 点を比較し、χ=512 は ED 基準内、
χ=128,256 は 12 sweep でも精度不足だった。中央 bond の 168 更新で実保持 rank と
報告損失を独立に照合した。この比較は新しい受理済み continuation 経路ではない。
この幾何には cut が一つしかなく内部 bulk column はない。
[相互作用系の検証記録](research/p2_interacting18_validation.md)に、有限区間の応答、
失敗と検証の限界を保存した。2π pump や magnetization plateau の成立は未検証である。

正の対照には Gong–Zhu–Sheng の拡張 kagome Heisenberg 模型を用いる。
論文は `J′=0.5`、`3×24×4` cylinder の U(1) DMRG で 2π 当たり 1/2 の移送を報告している。
最近接模型の 1/9 plateau とは異なる検証用モデルであり、元論文の結合図と
端・円周の定義を確認して再現する。距離だけで J2/J3 を一括定義しない。
[原論文と flux insertion 手法](https://www.nature.com/articles/srep06317)。
本文・補足の模型、幾何、flux 規約を照合した
[CSL 対照の実装前監査](research/p3_csl_control_design.md)では、J3 を六角形の対向頂点に
限定したテンプレート案と Q=0 対応の変更点を整理した。
原著の全 site 列・整数 QN 入力・pinning・内部刻み履歴は未取得であり、
提案した Q=0 と端の切り方を原著入力そのものとは扱わない。

小系では一意な基底状態が 2π で戻ることが正常であり、ED の最小固有値追跡に
2/3 のポンプを要求しない。ED はまず Hamiltonian と観測量の独立した検証に使う。

今後の実行順序は [ロードマップ](../ROADMAP.md)を正本とする。
QN 切断の修正・独立校正と 18 サイトの有限 χ 比較を終え、
明示Q、拡張結合、bond energyとN12Q0の小系照合も完了した。
NNの `N18,θ=0,Q=0,2,4` の静的比較を完了した。
`Lx=3,Ly=3,N=27,Q=3` のχ256・8-sweep試行は時間上限で完了状態を得られず、
続く2-sweep試行で進行記録とcheckpoint先行保存、varianceの別計時を完了した。
同χ上限のtrialから累積4 sweepsへ進め、実保持次元256の状態と全診断を取得した。
energy・密度の変動が残るこの時点では、追加収束を一度保留し、P3aの独立した観測量・小系校正へ進んだ。
その後の[段階的研究](research/p4_p3_staged_campaign.md)ではN27を再開し、
Q3のχ256累積8・χ512累積9、隣接Q1/5、周期9/27の初期状態とN54を比較した。
[隣接sectorの固定χ追加sweep比較](research/p4_static27_sector_refinement.md)に続き、
Q5だけをχ256・累積10へ進めた[再配列後の比較](research/p4_q5_csl_fixed_chi_followup.md)も別記録に残す。
分布変動は小さくなったが、固定χの精度5条件は未達だった。
各trialの数値整合性と精度・基底収束は区別する。
初期状態の精度を満たした場合に短区間追跡へ進み、実測後に
`Lx=6,Ly=3,N=54` 等で端の分離と複数cutを調べる。
既知CSL側ではchiralityの独立三spin・合成状態による校正を完了した。
拡張模型N18Q0のθ=0,.37も独立ED/DMRGと照合し、実空間/Schmidt移送・chirality・
bond energyの基準を通過した。[149.397秒の実測・零応答の解釈](research/p3_extended18_validation.md)。
これは小系校正であり、既知非零ポンプは未検証である。その後のCSL側では
`Lx=3,Ly=4,N=36` の零flux状態をχ128→256で準備・診断したが、
初期精度gateを通過せず、非零fluxは未着手である。追加の零flux精度改善が必要となる。
さらにχ256の零flux2 sweepsを追加してtrialを保存・reloadしたが、保存後の分散測定で
600秒上限に達した。新分散の欠測と切断・sweep energy基準の未達を
[固定χ追加研究](research/p4_q5_csl_fixed_chi_followup.md)に残し、非零fluxへは進んでいない。
NN静的研究では、N27を一時保留していた間にN18・Q6の追加で既存Q2境界が不変であることを確認した。
保存済みQ2の近接中性二重項も、縮退基底に不変な密度・bond遷移強度を測定した。
この小系の追加再加工を自動継続せず、NNの静的sector・秩序比較を優先し、
別模型のCSL零flux精度改善も進める。各計算前に別の予算を固定し、全組合せ走査を条件にしない。

`Lx=12,Ly=3,N=108` や `Lx=18 or 24,Ly=6,N=324 or 432` はその先の候補であり、
現段階で着手サイズや所要時間を確約しない。文献の YC 幅と比較するときは格子図を照合する。
9-siteと27-siteの競合する秩序単位胞と両立するcylinderを含める。
非整合な円周だけで VBC を排除しない。長さは相関長・端の侵入長より十分大きく取る。

代表点の費用を測ってから、判断を変え得るχ・seedの比較を選ぶ。
全χ×全seedを必須にせず、未収束ならsweep数とχの影響を分ける。
その結果と時間・メモリを見て1024,2048,…が必要か判断する。
小系の初期照合目標は `|ED−DMRG|/N < 10⁻⁸ J`、局所 Sz の差 `<10⁻⁶`。
研究用 pump の暫定精度目標は `<10⁻²` とし、少なくとも刻み・bond dimension・長さの
変更による差を個別に提示する。設定 cutoff からこの精度を推定しない。

## 8. 独自 backend の設計上の注意

MPS 方式では、整数 charge とその縮退次元を持つ link、charge を保存する局所テンソル、
複素 Hermitian な有効 Hamiltonian、sector ごとの SVD を実装する。
伝統的 block 方式では、各 block の charge と縮約された `Sz,S+,S−` 接続演算子を管理し、
`QL+qsite1+qsite2+QR=Qtarget` を満たす部分空間で局所問題を解く。
両方式で同じ物理的な sector と Schmidt 重みが得られることを検証する。

`../SUNDMRG.jl` から参考にできるのは sweep の段階分け、作業領域再利用、
Krylov、checkpoint、MPI ownership/cleanup の設計。
SU(N) の irrep、Wigner/Racah 係数、multiplet 次元の重みは U(1) 用に置き換える。
`step_density.jl` の実数 `syrk/syev` と `Float64` の格納は、単に型名を変えるだけでは不十分。
複素共役を含む `Ψ Ψ†`、Hermitian 固有値分解、または SVD と整合するカーネルが必要。
転置と随伴を明確に区別する。

flux 非依存の構造は再利用してよいが、θ が変わった後に旧 Hamiltonian の数値環境を
そのまま使わない。3 成分の環境を別々に保持する場合も、MPS 更新に合わせて再縮約する。
切断の選択は全 sector の特異値を合わせた重みに基づけ、物理的な移送に必要な sector を
固定の個数制限で排除しない。

CPU は block-sparse threading と BLAS threading の競合を避けて比較する。
GPU は ComplexF64 の実問題で速度・精度を測り、転送と小 block の overhead を含めて判断する。
[ITensor の threading 文書](https://docs.itensor.org/ITensors/stable/Multithreading.html)。

## 9. ポンプから相同定へ

2/3 が安定に得られた場合も、SU(3)₁ 型か D(Z₃) 型かは未確定である。
Abelian な最小エンタングルメント状態に対し、幅依存の
`S=α Ly−γ` から `γ=½ log 3` と `log 3` を比較する。
ただし少数の細い cylinder の外挿だけで決めず、有限 bond dimension と秩序の影響を確認する。

さらに複数セクター、momentum-resolved entanglement、momentum polarization または
modular データを検討する。これらには並進対称性・ゲージ・セクター基準の追加設計が必要。
単一の U(1)-resolved entanglement spectrum や chirality の期待値だけから c₋ を決めない。
[Zaletel–Mong–Pollmann の原論文](https://arxiv.org/abs/1211.3733)。

ゼロ応答の D(Z₃) 候補では、非自明な sector flow や分数 charge など、pump 以外の証拠が必要。
ゼロ pump は、自明相、gapless 状態、枝追跡失敗でも起こり得る。
逆に、2/3 に近い一つの端密度差だけでは量子化の証拠としない。

1/9 plateau には異なる提案が存在する。
[2023 年の iPEPS/PESS 研究](https://arxiv.org/abs/2306.09563)は VBC と gapless な性質を報告し、
[2024 年の VMC 研究](https://arxiv.org/abs/2407.20629)は chiral Z₃ 候補を提案している。
さらに[Cheng–Li, PRB 113, 085136 (2026)](https://journals.aps.org/prb/abstract/10.1103/5tvd-253q)は
3×3 windmill VBCを提案している。[2025年の著者稿](https://arxiv.org/abs/2512.11670v1)には
gaplessなchiral spin density waveの提案もあるが、今回査読誌出版と定量的な規約対応は未確認。
各手法の主張と本研究の観測を区別し、密度の一様性やchiralityだけで候補を絞らない。
2024年VMCのspinon Chern numberと、ここで指定したK,tの物理的なSz応答は自動的には同じではない。
spinon に ±Θ の境界位相を与える場合、`S+=f†up fdown` の物理的 twist の大きさは 2Θ になるため、
flux の規格化も照合する。今回比較する応答はユーザー指定の K,t に基づき、
数値収束の判定は期待する応答値とは独立に行う。

## 10. 保存する成果物

各走行について格子定義・順序・couplings・Q・ゲージ・初期状態の seed・
採用した θ 列・solver 設定・パッケージバージョンを保存する。
各受理点にはエネルギー、可能なら variance、実測切断誤差、各 cut の pump、
Sz profile、chirality、EE、charge spectrum、前段 overlap、枝の診断結果を対応させる。
checkpoint は最後の受理点を原子的に保存し、再開時に格子・site indices・Q・設定を照合する。

現実装の `save_checkpoint` は完了した DMRG 点と零 flux の基準 MPS を保存する。
受理済み／試行を別ディレクトリに置き、新しい snapshot を一時ディレクトリから
同じ親への rename で公開する。既存 snapshot の上書きは行わない。
TOML metadata と payload のサイズ・SHA-256 を確認し、同一 Julia・依存バージョン・
ソース・実際のactive Projectと選択manifestの一致を要求してから Serialization payload を読む。
sourceはmodule評価時、環境・runtimeは`__init__`で捕捉し、precompile cacheを共有しても
他の実行環境へ由来を付け替えない。ロード後のactive環境変更も拒否する。
復元後は電荷・site identity・norm・密度・基準状態を照合し、保存した模型の
Hamiltonian を再構築して energy を確認する。`resume_dmrg` は受理済み状態の
solver 設定で新しい DMRG batch を開始し、Hamiltonian 環境は再構築する。
保存された path と受理ラベルは呼出側が指定するもので、自動の枝判定ではない。
`continue_flux` から保存する場合は明示的な診断 policy に基づく受理ラベルとなり、
試行・理由・全 unwrapped path は同じ走行の `trajectory.toml` に原子的に更新される。
受理点から再開する際も初期零 flux の基準を引き継ぐ。新しい policy・目標列・刻みを
指定した別 journal を作る方式であり、driver 内部の途中命令を復元するものではない。
同一環境の信頼済みローカル snapshot 用で、長期交換形式・途中 sweep の再開・
電源断に対する耐久性は保証していない。
[実測検証と再現手順](research/p1_restart_schmidt_validation.md)を参照。

最初の研究レポートは、`ΔSz_right(θ)` と誤差、列ごとの密度変化、
flux に対する entanglement/charge flow、収束比較、候補の判定範囲から構成する。
CPU/GPU 型番・計算時間・最大メモリは性能比較に使える。
環境変数の一括 dump、ローカルのネットワーク設定、接続先、認証情報は保存しない。
