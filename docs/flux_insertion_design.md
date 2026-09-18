# スピン flux 挿入の物理・数値設計

作成日: 2026-09-19。[ロードマップ](../ROADMAP.md)の仕様詳細。
以下の設定値は検証を始めるための提案であり、計算結果や収束保証ではない。

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

\[
M_{\rm sat}=N/2,\quad M=N/18,\quad Q=2M=N/9,
\quad N_\uparrow=5N/9,\quad N_\downarrow=4N/9.
\]

- `N % 9 == 0` を要求する。N が偶数なら M は整数、奇数なら半整数でよい。
- ITensor の spin QN `Sz` の整数ラベルは物理的な `2Sz`。
  目標値は `QN("Sz",N÷9)` であり、物理 Sz=N/9 と取り違えない。
  これはゼロ磁化の交互 Up/Dn 初期状態からは到達しない。
  [公式 QN DMRG チュートリアル](https://docs.itensor.org/ITensorMPS/stable/tutorials/QN_DMRG.html)。
- 固定 M では `-hM` は定数なので、flux 走行中の h 調整は不要。
  plateau 自体の確認は別の磁化 sector も用いる。
- 隣接 sector の交換エネルギー `E0(M)` から得る
  `h−=E0(M)−E0(M−1)`, `h+=E0(M+1)−E0(M)` は最初の目安。
  より遠い sector を含めた凸包、端の磁化、サイズ依存を調べて plateau を確認する。
  磁化 plateau があっても、同じ M 内の中性励起まで gapped とは限らない。
  有限 cylinder の最小励起が端に局在する可能性もあるため、励起の空間分布を確認し、
  端の gap をそのまま bulk gap としない。

## 3. Cylinder とゲージ

初期実装では `x=0,…,Lx−1` を開境界、`y=0,…,Ly−1` を周期境界とする。
`a1=(1,0)`, `a2=(1/2,√3/2)`、副格子座標 `A=0, B=a1/2, C=a2/2`、
`N=3LxLy` を採用する。Ly は原始単位胞数であり、文献の YC 表記とは自動的に同一視しない。
wrap vector、端の切り方、MPS の site ordering を明示して比較する。

最近接ボンドの候補テンプレートは各 `(x,y)` に対し
`A-B`, `A-C`, `B-C` の胞内 3 本と、
`A(x,y)-B(x−1,y)`, `A(x,y)-C(x,y−1)`,
`B(x,y)-C(x+1,y−1)` の 3 本。
x が範囲外のボンドを除き、y を wrap する前の winding を保存する。
この切り方で Ly≥3 のとき bond 数は `(6Lx−2)Ly`、バルク配位数は 4。
Ly=1,2 の特殊な周期同定は初期の研究対象から外す。

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
- uniform gauge は円周方向の並進を扱う次段階で使用する。
  2π 後の比較には大きなゲージ変換が必要で、未変換の MPS overlap を比較しない。
- 長距離結合や多体項を追加する場合も、同じゲージ規則から twist を導出する。
  chirality seed を残す場合にはその境界項も含める。

## 4. ITensor による最小構成

現在の API では `using ITensors, ITensorMPS` を用いる。
`siteinds("S=1/2",N; conserve_qns=true)` を一度作り、全 θ で再利用する。
対象 sector の Up/Dn 個数を満たす初期状態を複数用意し、複素 MPS を扱う。
以下は設計用の抜粋であり、完成した実行例ではない。

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
使用する Julia・依存パッケージの具体的な互換バージョンは P0 で解決して固定する。

初期は各 θ の MPO を構築し直し、`dmrg(H, psi_previous; ...)` で最適化する。
実際の切断誤差、energy/variance、各観測量の収束を測り、`cutoff` の設定値だけを
誤差として報告しない。Krylov solver の収束も確認する。
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

正の対照には Gong–Zhu–Sheng の拡張 kagome Heisenberg 模型を用いる。
論文は `J′=0.5`、`3×24×4` cylinder の U(1) DMRG で 2π 当たり 1/2 の移送を報告している。
最近接模型の 1/9 plateau とは異なる検証用モデルであり、元論文の結合図と
端・円周の定義を確認して再現する。距離だけで J2/J3 を一括定義しない。
[原論文と flux insertion 手法](https://www.nature.com/articles/srep06317)。

小系では一意な基底状態が 2π で戻ることが正常であり、ED の最小固有値追跡に
2/3 のポンプを要求しない。ED はまず Hamiltonian と観測量の独立した検証に使う。

試行サイズは `Lx=12,Ly=3,N=108`、次に `Lx=18 or 24,Ly=6,N=324 or 432` を候補とする。
メモリと収束に応じて調整し、文献の YC 幅と比較するときは格子図を照合する。
9-site 以上の秩序単位胞と両立する cylinder を含める。
非整合な円周だけで VBC を排除しない。長さは相関長・端の侵入長より十分大きく取る。

初期の bond dimension は例えば 256→512→1024 と段階的に増やす。
これで量子化が決まるとは想定せず、実測で 2048,4096,… が必要か判断する。
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
後者の spinon Chern number と、ここで指定した K,t の物理的な Sz 応答は自動的には同じではない。
spinon に ±Θ の境界位相を与える場合、`S+=f†up fdown` の物理的 twist の大きさは 2Θ になるため、
flux の規格化も照合する。今回比較する応答はユーザー指定の K,t に基づき、
数値収束の判定は期待する応答値とは独立に行う。

## 10. 保存する成果物

各走行について格子定義・順序・couplings・Q・ゲージ・初期状態の seed・
採用した θ 列・solver 設定・パッケージバージョンを保存する。
各受理点にはエネルギー、可能なら variance、実測切断誤差、各 cut の pump、
Sz profile、chirality、EE、charge spectrum、前段 overlap、枝の診断結果を対応させる。
checkpoint は最後の受理点を原子的に保存し、再開時に格子・site indices・Q・設定を照合する。

最初の研究レポートは、`ΔSz_right(θ)` と誤差、列ごとの密度変化、
flux に対する entanglement/charge flow、収束比較、候補の判定範囲から構成する。
CPU/GPU 型番・計算時間・最大メモリは性能比較に使える。
環境変数の一括 dump、ローカルのネットワーク設定、接続先、認証情報は保存しない。
