# J1–J2–J3交換模型と結合エネルギーの小系検証

2026-09-19。明示Qの共通APIに続き、既知CSL対照に必要な拡張交換模型を実装する。
今回の単位は幾何、Hamiltonian、結合エネルギー、保存契約の小系照合である。
chirality観測量、N18Q0の移送、既知CSLのポンプ検証は後続とする。

## 模型と独立性

[原著 Eq.(1)](https://www.nature.com/articles/srep06317) のJ2/J3は六角形内の交換であり、
J3には同じ距離の直線鎖上の第三近接を含めない。
境界を跨ぐ全交換へのtwistは[著者稿のMethods](https://arxiv.org/pdf/1312.4519v2)とも照合した。
このcylinderの具体的なsite列と端は[先行する設計](p3_csl_control_design.md)の規約であり、
原著の有限cylinder入力を完全に復元したものではない。

`kagome_j1j2j3_cylinder(Lx,Ly; J1=1.0,J2=0.5,J3=0.5)` を別名で追加した。
既存の `kagome_cylinder` と1/9磁化の既定動作は維持する。
各familyで `Jxy=Jz=Ja` とし、零結合も明示的に保持する。
全結合の横交換は `(Ja/2)cis(wy*theta)`、縦交換はtwistしない。

productionは無限格子の代表bondを生成し、両端のxが保持範囲にあるものを残す。
独立EDのJ1/J2は距離判定、J3は整数座標の平面NNグラフからchordlessな6頂点cycleを
探索し、対向頂点を選んだ後に端を切り出す。productionのテンプレートは参照しない。
2cellの余白を含めるため、系内に六角形全体が残らなくても有効な端のbondを保持する。
periodicなsite番号と符号付き巻き数への変換は最後に行う。

familyごとの結合数は `((6Lx−2)Ly,(6Lx−4)Ly,(3Lx−2)Ly)`。
幾何参照は `1≤Lx≤3, 3≤Ly≤4` に限定し、行列参照はさらにN≤18、
dense次元≤4096を要求する。研究用の無制限な六角形探索器にはしない。
通常テストは幾何の代表2サイズと、`J1=.8,J2=.3,J3=.6` のN9複素MPO/EDを使う。
不等なJ2/J3はfamilyの取り違えを検出するためであり、追加のDMRG・固有値計算は行わない。

## 観測量と保存

`bond_families(lattice)` は実際の `(i,j,wy)` を幾何と照合し、bond順に
`:J1,:J2,:J3,:other` を返す。反転・結合値・零結合に依存しない。
モデル名を記憶するmutableなタグを追加せず、bond配列の編集も毎回反映する。

`bond_energies(psi,lattice,theta;gauge=:seam)` はbond順の規格化された期待値
`Jz*<Sz_i Sz_j> + Jxy*real(cis(A_b)*<S+_i S-_j>)` を返す。
場は含めず、全エネルギーは `sum(bond_energies(...))-dot(hz,sz_profile(psi))` となる。
初版は既存の2つの全相関行列を再利用するため、保存量はO(N²)。
大系での費用は未測定で、必要になればbondだけに絞る実装を計測して検討する。

checkpointは全bondの結合・向き・巻き数に加え、familyとその定義を保存する。
schema 1を維持し、読み込み側で実際の要求格子からfamilyを再計算して一致を確認する。
旧ソースのsnapshotは従来どおりsource identity不一致で拒否し、自動移行しない。
新しいproductionファイルもprovenanceの対象とする。

## 実行前に固定した研究計算

[validate_extended_model.jl](../../examples/validate_extended_model.jl) で
`(Lx,Ly)=(1,4),N=12,Q=0,J1=1,J2=J3=.5` のθ=0,0.37だけを実行する。
結合数は16/8/4、ED次元は924。全固有対を求め、最低準位から `1e-10` 以内の
全状態を含む基底空間にDMRG状態を射影する。

Julia/BLAS各1 thread、seed11、maxdim64、8 sweep、cutoff/noise=0、局所eigsolve
tolerance `1e-12`、Krylov次元30、最大反復10とする。
θ=0.37は零fluxの状態と同じsite indicesを使ったwarm startで逐次進める。
θ間で300秒の上限を確認して超過時は後続点を始めない。
12サイトの最大Schmidt rankは64なので、この設定は有限χの収束走査ではない。

基準はenergy/site誤差≤`1e-8`、独立残差と基底空間からの漏れ≤`2e-6`、
局所Sz・各bondエネルギー誤差≤`1e-6`、電荷/norm誤差≤`1e-12`。
独立seed33の複素MPSに対するMPO作用、2π周期、J2=J3=0のNN還元、
seam/uniform間の行列・状態・bondエネルギーの差は≤`1e-11`を要求する。
全Hilbert空間の4096² MPO行列は作らず、MPO作用を直接収縮する。
trialを保存・再読込し、Qとfamily、状態の一致とNN格子への誤った読み込みの拒否も確認する。

量子化pumpの期待値に近づくことを最適化・受理条件に使わない。
Lx=1には軸cutがなく、この2点は受理済みのflux trajectoryでもない。
小系の一致はCSL・量子化・磁化plateauを示さない。

## 実測結果

2点とも事前基準を全て通過した。数値・設定・依存版・source hashの正本は
[実行記録TOML](data/p3_extended_model_validation.toml)。trial snapshotは
`outputs/p3-extended-model-20260919/trial/` に保存した。

| θ | ED最低エネルギー | energy/site誤差 | 独立残差 | 基底空間からの漏れ | 最大bond energy誤差 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | −5.039405220964036 | 5.18e−16 | 3.62e−14 | 1.88e−14 | 1.09e−14 |
| 0.37 | −5.037683613046604 | 2.96e−16 | 4.19e−13 | 2.02e−13 | 5.06e−14 |

両点のQ0最低準位は非縮退だった。記録した次準位との差は有限クラスター・固定Qの
準位差であり、bulkの中性gapやCSLの証拠とは解釈しない。
局所Szの最大誤差は `4.96e−15` 未満だった。

familyごとの交換エネルギーの和は次の通り。正負も含めそのまま保存した。

| θ | J1の和 | J2の和 | J3の和 |
| ---: | ---: | ---: | ---: |
| 0 | −4.4478248412613866 | 0.08394619233591522 | −0.6755265720385608 |
| 0.37 | −4.447245711084801 | 0.0808783696611267 | −0.6713162716229298 |

全bondの和とHamiltonian期待値の差は≤`3.56e−15`。
両ゲージのMPO作用と独立EDの差は≤`1.48e−15`、θ=.37のゲージ行列差は `3.69e−15`、
各bond energyのゲージ間差は≤`5.56e−16`だった。
J2=J3=0のNN還元はMPO作用で差0、2π周期の差は≤`1.37e−15`だった。
全trialの再読込でQ・family・状態が一致し、NN格子を指定した誤読込を拒否した。

統合 `Pkg.test()` は1回で **1,322 / 1,322件成功**。
Test計測4分17.7秒、起動込み266.78秒で、通常suiteのDMRG・固有値計算本数は
増やしていない。3分目標は未達、5分の見直し基準内である。
研究driver全体は83.19秒、import後の内部計測は80.79秒。
最初のDMRGはJITを含め40.64秒、次点は0.49秒、密EDは各0.47–0.51秒だった。
内部時間には診断・保存とその初回compileも含む。
geometryと位相、観測量、保存契約は別担当が原著・式・コードを独立監査し、
その後に主担当が統合した実コードで以上を実測した。

今回完了したのは **J1–J2–J3の幾何・交換Hamiltonian・bond energy・保存とN12照合**。
次は方向付き三角形のchiralityを三spin行列で校正し、N18Q0の独立残差・
射影・実空間/Schmidt移送を検証する。P3a全体とCSLポンプ検証は未完了のままとする。
