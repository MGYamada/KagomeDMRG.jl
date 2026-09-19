# P2: 18 サイト相互作用系の精度と flux 往復

実行日: 2026-09-19。対象は最近接等方的 kagome Heisenberg 模型、
`Jxy=Jz=1`、縦磁場・chirality seed なし、固定 `M/Msat=1/9` sector。
今回の goal は、軸方向 cut を持つ小系で独立 sparse ED と DMRG を照合し、
bond dimension・flux 刻み・初期状態への依存を記録することだった。
同時に [既知 CSL 対照の原著監査](p3_csl_control_design.md)を行った。

## 幾何と数値設定

`Lx=2,Ly=3,N=18,Q=2,Nup=10,Ndown=8`、物理的全 `Sz=1`。
wrap は `3a2`、軸 OBC・円周 PBC、完全な単位胞を残す端、site 順序は `x,y,A/B/C`。
有向最近接 bond は 30 本。seam gauge で `S+_i S-_j` の係数を
`(1/2)cis(wy*theta)` とし、`SzSz` は twist しない。

幾何学 cut は bond 9 の一つだけで、右領域は site 10:18。
二列しかないため内部の bulk column はなく、密度の監視対象は全 18 サイトとする。
cut 間の一致を検証できる形状ではない。零 flux の実測 profile を基準に、
左右の累積移送と絶対左 Schmidt 電荷の変化を保存する。

全 DMRG 点で六 sweep を設定し、noise=0、cutoff=0、局所固有値 tolerance `1e-11`、
Krylov dimension 40、最大 iteration 20、variance 測定あり。
中央 cut の厳密な固定 Q rank の上限は

```text
sum(min(binomial(9,n), binomial(9,10-n)) for n=1:9) = 386
```

であり、maxdim=512 はこの小系を切断なしで表現できる。
512 を研究サイズでも十分な bond dimension とする意味ではない。

追跡前に固定した policy は overlap 振幅 `≥0.90`、密度変化 `≤0.05`、
entropy 変化 `≤0.10`、左 Schmidt Sz 変化 `≤0.05`、variance 絶対値 `≤1e-8`、
最後の sweep の実測切断誤差 `≤1e-10`、energy の収束差 `≤1e-9`。
整合性 tolerance は `1e-9`。一つしかない cut の spread=0 は非自明な検査にならない。
期待 pump 値や ED 最低エネルギーは continuation の採否条件に入れていない。
独立 ED との比較は採否とは別に実行する。

## 独立 sparse ED

sector の次元は `binomial(18,10)=43758`。
[reference_eigensolve.jl](../../test/reference_eigensolve.jl)は、production の bond や MPO を
使わず距離探索から構成した疎行列を BlockLanczos で対角化する。
単一ベクトル Lanczos による縮退の取り逃しを避けるため複素初期 block を使い、
収束個数、全行列で再計算した残差、直交性を照合する。
[KrylovKit の一次文書](https://jutho.github.io/KrylovKit.jl/stable/man/eig/)に従った。
9 サイトの dense 固有空間と、既知の三重縮退行列で helper を検証した。

独立担当は seed=1234、blocksize=3、nev=3、tolerance `1e-11`、
Krylov dimension 80、最大 iteration 200 を用いた。
各点の三対は収束し、再計算 residual の最大は `7.38e-12` 未満だった。
[ED 参照記録](data/p2_ed18_reference.toml)に全準位・設定・ソース hash を保存した。

| θ | 計算した最低エネルギー | 最低二準位の間隔 |
| --- | ---: | ---: |
| 0 | −7.442339744470409 | 0.001864675030 |
| 0.185 | −7.442134635417066 | 0.001042166319 |
| 0.37 | −7.441523009592125 | 0.000175737922 |

主走行では別の seed=5678、blocksize=4 で各保存点を照合した。
部分スペクトルの residual 収束だけで完全性を数学的に保証するものではない。
この準位間隔を neutral bulk gap と呼ばず、励起の端への局在も未判定とする。
特に θ=0.37 では間隔が小さいため、variance 閾値の通過だけに頼らず
実際の residual と基底空間射影を記録する。

独立に円周方向の site 置換 `T:(x,y,s)→(x,y+1 mod 3,s)` を spin basis 上に構築した。
θ=0 で `T³=I` と `T†HT=H` は厳密に成立した。
計算した基底ベクトルの T 固有値は 1、第一励起の二重項は `exp(±2πi/3)`。
部分空間の translation leakage は `4.70e-12` 未満だった。
[対称性監査記録](data/p2_ed18_translation.toml)を保存した。
これは有限円筒の円周対称性であり、トポロジカル縮退の証拠ではない。

## bond dimension と未解決の結果

零 flux の seed=11 を同じ六 sweep の設定で比較した。
χ は設定した maxdim であり、実際に保持した次元とは区別する。

| χ | 最終 energy の ED 誤差 / site | 独立全系 residual | 最大 Sz 誤差 | 結果 |
| ---: | ---: | ---: | ---: | --- |
| 128 | `6.02e-6` | `2.47e-2` | `3.35e-3` | 初期点を棄却、flux 走行なし |
| 256 | 取得できず | 取得できず | 取得できず | sweep 4、往路、bond 9 で切断 guard が停止 |
| 512 | `9.38e-16` | `7.18e-12` | `5.23e-13` | 初期点を受理、実際の最大 link 次元 386 |

χ=128 の raw variance は `6.10e-4`、最後の sweep の実測切断誤差は
`7.84e-6`。variance・切断誤差・sweep energy 差が policy を満たさず、
`unresolved` とした。初期点の移送 0 は基準の定義によるもので、零応答の観測ではない。

χ=256 では `DMRG truncation spectrum does not match the retained bond dimension`
が発生した。既知の QN 切断整合性 guard が相互作用する最近接模型でも作動することを
確認したが、この走行では内部の不一致を捕捉する以上の根本原因の切り分けはしていない。
完了状態や checkpoint がないため、停止前の局所 sweep energy を収束値として使わない。
[失敗の原記録](data/p2_chi256_failure/validation.toml)と
[停止位置・設定・観測の解釈](data/p2_chi256_failure/failure_observation.toml)を保存した。
原記録の `caller_assigned_status` は完了時に予定していたラベルであり、実際の
`status="error"` と `truncation_guard_fired=true` がこの結果を表す。
零 flux の初期状態での失敗なので、continuation の刻み半減では解決しない。

χ=512 は全保存点で実測切断誤差 0、保持次元は最大 386 だった。
これは切断のない小系参照の成立であり、有限 χ での一般的な切断校正の完了ではない。

## 相互作用系の往復と精度比較

χ=512 の四経路は、設定した点をすべて受理して完走した。
初期点の再利用を含む延べ 24 保存点のすべてが、continuation とは独立した
ED 照合の閾値を満たした。試行の棄却・刻み半減はこの四経路では発生しなかった。
[全数値記録](data/p2_interacting18_validation.toml)に各訪問点、診断、採否を保存した。

| seed | 最大刻み | 経路 | 保存点 | 折返し点の ΔSz_right | θ=0 へ戻した ΔSz_right |
| ---: | ---: | --- | ---: | ---: | ---: |
| 11 | 0.185 | `0→0.37→0` | 5 | `−2.65650232642054e-4` | `1.62e-13` |
| 11 | 0.0925 | `0→0.37→0` | 9 | `−2.65650232658548e-4` | `8.31e-14` |
| 29 | 0.185 | `0→0.37→0` | 5 | `−2.65650232751803e-4` | `5.78e-14` |
| 11 | 0.185 | `0→−0.37→0` | 5 | `−2.65650232646266e-4` | `1.65e-13` |

全 χ=512 保存点で、energy/site の ED 誤差は `2.03e-15` 以下、独立全系 residual は
`7.83e-12` 以下、基底空間からの leakage は `6.09e-11` 以下、最大 Sz 誤差は
`1.26e-11` 以下だった。MPO と独立行列で測った同一状態の energy 差は `2.63e-13` 以下。
raw variance には差し引きの丸め誤差で負になる点もあり、独立 residual の二乗との差は
`2.63e-13` 以下だった。この床より小さい誤差を raw variance だけから評価しない。

隣接する受理状態の overlap 振幅の最小は粗い刻みで `0.99414`、細かい刻みで
`0.99853`。左右移送の打消し誤差は `3.81e-15` 以下、右移送と左 Schmidt Sz 変化の
符号反転との差は `4.34e-15` 以下だった。
往路・復路の同じ θ を別の訪問として対応させると、刻み半減による右移送の差は
最大 `4.21e-13`、seed 変更による差は最大 `1.15e-13`。
零 flux 帰還時の site profile の初期値との差は全経路で `1.18e-12` 以下だった。

θ=±0.37 の energy 差は `8.89e-15`、最大密度差は `2.51e-12`、
`|⟨ψ(−θ)|conj(ψ(θ))⟩|=0.9999999999999976`。
現在の模型では Sz 基底で `H(−θ)=conj(H(θ))` であり、得られた状態はこの関係と整合する。
正負で同符号の小さな密度再分配は、この有限円筒の観測として記録する。
固定された chiral 枝の Hall 応答や 2π pump を測った結果とはしない。

独立レビューでは χ=128 の初期点と、χ=512 の粗い刻みの初期点・θ=0.37・帰還点を
別の spin-basis 展開と 512×512 の全 Schmidt SVD で照合した。
最大密度差 `2.39e-15`、同一状態 energy 差 `8.71e-14`、Schmidt 確率差 `8.89e-16`、
右移送差 `5.14e-15` 以下だった。
[独立照合記録](data/p2_interacting18_independent_audit.toml)に結果と source hash を残した。
細かい刻み・seed=29・負 flux の経路全体を別担当が独立再走行したわけではない。

## 再現方法と保存

```sh
julia --project=. --startup-file=no --threads=1 examples/validate_interacting18.jl
```

Julia 1.13.0、ITensors 0.9.31、ITensorMPS 0.4.1、NDTensors 0.4.31、
KrylovKit 0.10.4 を使用し、各 process の Julia/BLAS は 1 thread。
主走行は `outputs/p2-interacting18-study/` に置いた。
既に同じ source・模型・solver で試走した seed=11 の初期状態を再利用する実行では、
第二引数として信頼済み preflight record を渡す。
読み込み時に checkpoint の source/runtime/geometry/charge と solver 設定を照合する。

```sh
julia --project=. --startup-file=no --threads=1 examples/validate_interacting18.jl \
  outputs/p2-interacting18-study outputs/p2-interacting-preflight/preflight.toml
```

どちらも出力先が既存なら上書きしない。
全試行の journal と MPS、独立 ED の数値、実際の採否理由を保持する。
機械可読記録の `completed_with_recorded_outcomes` は実験手順を終えた意味で、
すべての精度・枝が合格した意味ではない。
主走行前後の source・manifest hash は一致し、保存記録の `code.status="verified"` と
現在の対象ファイルの hash も照合した。上記の TOML は原出力から内容を変えずに複製した。
大きな MPS・ED vector と journal は `outputs/` に保持し、git の対象外とする。
χ=256 の追加走行の設定は失敗観測記録に列挙し、同じ `run_dmrg` の入力を再構成できる。
1/9 はここでは固定した磁化 sector を表す。隣接 Q のエネルギーや有限の磁場区間を
調べていないため、magnetization plateau の成立は主張しない。

## 検証範囲と次段階

標準 `Pkg.test()` は Julia 1.13.0 で全 **4,348 assertions** が成功した。
既存 4,283 に、独立 sparse 固有値・縮退・入力・非収束の 65 検査を追加した。
18 サイトの長めの研究走行は標準テストから分離した明示的な実行例であり、
Julia 1.12.7 で今回の拡張 suite を再実行した結果ではない。

P2 全体や P3 の完了ではない。今回の短い flux 区間から 2π の Hall 応答は評価できず、
長さ・幅への収束、複数 bulk cut、隣接磁化 sector、励起の局在は未検証である。
相互作用系でも作動した QN 切断 guard の根本修正・独立校正を、研究サイズへの拡大前に進める。
既知 CSL 対照には Q=0 の明示指定、六角形内 J2/J3、独立した拡張模型 ED と
chirality 診断が必要であり、実装条件を [P3 設計](p3_csl_control_design.md)にまとめた。
SU(3)₁、Hall-active/Hall-inactive D(Z₃)、秩序相・gapless 状態を識別する結論は出していない。
