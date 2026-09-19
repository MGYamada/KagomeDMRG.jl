# P3a: 方向付き scalar chirality の実装・校正

日付: 2026-09-19。状態: **観測量の実装と独立校正を完了**。
N27の追加sweepを保留し、[ロードマップ](../../ROADMAP.md)の観測量を整備する。
対象は観測演算子であり、chirality seed項をHamiltonianに追加する変更ではない。
既知CSLの非零ポンプや、最近接1/9の相同定はこの校正に含めない。

## 三角形と演算子の定義

`oriented_triangles(lattice)` は、現在の端・wrap `Ly*a2` における
最近接の基本三角形を返す。`KagomeTriangle.sites` は反時計回りの3頂点、
`image_y` は各頂点の周期画像、`kind` は `:up` または `:down`。
J2/J3結合を追加した格子でも、同じ基本三角形を使う。

| 種類 | 周期画像を巻く前の頂点順 | 範囲 |
| --- | --- | --- |
| up | `A(x,y),B(x,y),C(x,y)` | `0≤x<Lx,0≤y<Ly` |
| down | `A(x,y),B(x−1,y),C(x,y−1)` | `1≤x<Lx,0≤y<Ly` |

両者の符号付き面積の2倍は `sqrt(3)/8>0`、総数は `(2Lx−1)Ly`。
`y=0` のdown三角形は画像 `(0,0,−1)` を持ち、その辺のwindingは
`(0,−1,+1)`。和は0であり、局所的なfluxを囲まない。

物理スピン `Sz=±1/2`、`hbar=1` で

\[
C_{ijk}=\mathbf S_i\cdot(\mathbf S_j\times\mathbf S_k)
=\frac{i}{2}\sum_{(a,b,c)=(i,j,k),(j,k,i),(k,i,j)}
(S_a^+S_b^- - S_a^-S_b^+)S_c^z.
\]

`scalar_chirality_mpo(sites,(i,j,k);angles=(0,0,0))` はこの演算子を
U(1)保存の複素MPOで返す。巡回置換は同じ符号、向きの反転は逆符号になる。
任意の局所角 `alpha` を渡すと `U_alpha*C*U_alpha†` を返し、
`S+_a S-_b` の係数に `exp(i*(alpha_a−alpha_b))` が掛かる。

`triangle_chiralities(psi,lattice,theta;gauge)` は三角形の順に正規化期待値を返す。
入力MPSは変更しない。既存のDMRG結果には次のように適用できる。

```julia
triangles = oriented_triangles(lattice)
values = triangle_chiralities(result.psi, lattice, result.theta; gauge=:seam)
# values[a] corresponds to triangles[a].sites and triangles[a].image_y.
```

非零fluxでは次の局所移送角で定義した **dressed chirality** を測る。
周期画像を `m_a`、既存の `gauge_angles` を `chi_a` とすると、

\[
\alpha_a^{\rm seam}=-m_a\theta,\qquad
\alpha_a^{\rm uniform}=-m_a\theta+\chi_a.
\]

seamでは `alpha_a−alpha_b=(m_b−m_a)theta` なので、辺の交換位相と一致する。
uniformでは状態と演算子の双方を同じ `U_chi` で変換する。
この定義は `theta=0` で通常のscalar chiralityに一致するが、非零fluxで
継ぎ目を跨ぐ三角形では、位相を付けない裸の三スピン積とは異なる。
`theta` はそのまま使い、有限fluxで異なる定義の観測値を混ぜない。

## 独立校正の基準

三スピンの密行列をPauli行列のLevi-Civita和から作り、productionの
ladder展開・MPO構築を参照せずに照合する。固有値は
`−sqrt(3)/4` と `+sqrt(3)/4` が各2個、0が4個。
`omega=exp(2pi*i/3)` として

\[
|+\rangle=\frac{1}{\sqrt3}
(|\downarrow\uparrow\uparrow\rangle+
\omega|\uparrow\downarrow\uparrow\rangle+
\omega^2|\uparrow\uparrow\downarrow\rangle)
\]

は `Q=1` で期待値 `+sqrt(3)/4`、その複素共役は逆符号になる。
これは測定符号の解析的対照であり、kagome模型のCSL基底状態ではない。

## 実行結果

`test/test_chirality.jl` の独立参照で上記の三spin行列・符号・固有値、
U(1)保存、Hermiticity、向きの反転、任意の局所回転を確認した。
周期画像を持つ平面距離graphから三角形集合を独立に列挙し、
`(Lx,Ly)=(1,3),(2,4)` の個数・CCW方向・windingを照合した。
`(2,3)` のNNと異なるJ1/J2/J3模型でも三角形列は一致した。

有限fluxの対照は `N=18,Q=2,theta=0.61` の**合成状態**で、
継ぎ目を跨ぐ `(10,2,18)` に上記の正chirality状態を埋め込んだ。
3個の積状態の和なので保持次元は最大3、N18の密Hilbertテンソルは作らない。
入力normを2に変えた場合も、正規化値と全入力テンソルの保存を確認した。

| 校正量 | 結果 |
| --- | ---: |
| 独立Pauli行列とMPOの最大要素差 | 約 `8.33e-17` |
| 正の三spin対照の期待値 | `+sqrt(3)/4 ≈ 0.4330127`、誤差 `<2e-14` |
| seam/uniformの全三角形期待値の最大差 | 約 `2.22e-16` |
| 非正規化合成状態の期待値と解析値との差 | `0.0` |
| 測定前後の入力テンソルの最大差 | `0.0` |

誤差表示はJulia loggerの丸め精度。設定・検証値・180個のsource/manifest hash・
13個のtest/helper hashは[機械可読記録](data/p3_chirality_validation.toml)に保存した。
理論・実装・独立参照を別担当で照合した後、次の重点検証を1回実行した。

```sh
julia --project=. --startup-file=no --threads=1 -e 'using Pkg; Pkg.test(; test_args=["lattice", "observables", "checkpoint"])'
```

Julia 1.13.0、Julia/BLAS各1 threadで **736 / 736件が成功**。
Test計測94.0秒、外側Julia起動を除く `Pkg.test` 呼出103.058秒、
package precompile約3秒（呼出時間の内数）。初回MPO/MPS・Schmidt等のJITを含む
3群の実行で60秒の重点目安を超えた。chirality単独の時間は未分離であり、
大系の測定性能へ外挿しない。追加の全体suiteは実行せず、Hamiltonian・
DMRG・切断・continuationのアルゴリズムを変更しない今回の影響範囲を確認した。

## 保存状態と残る範囲

新ソースをcheckpointのallowlistに登録する。既存の厳密なsource/runtime一致を
維持するため、以前のN27 trialを新ソースの結果として読替え・再保存しない。
変更前に、N27記録にある179個のproduction/manifestファイルと6個のdriverを
`outputs/source-snapshots/n27-refine-before-chirality-20260919/` にコピーし、
全SHA-256の一致を確認した。既存のtrial payloadは変更していない。
コミット時には旧実装を `48d15e2b566b4d4b50a6c5295b84348cc5f30f8d` に分離し、
そのGit treeの179個のproduction/manifestと6個のdriverも、元の数値記録のSHAと全一致した。
このコミットは後から実装を保存したものであり、計算時に記録したGit revisionは書き換えない。
これは旧ソースの保全であり、異なる版のcheckpoint移行を実装したものではない。

初版は三角形ごとに小さいMPOを構築・収縮する。研究サイズでの費用は未計測。
非零chiralityだけでCSLとはしない。次のP3a作業は拡張模型
`J1=1,J2=J3=0.5,N=18,Q=0` の残差・射影・移送の照合であり、
その後に複数cutを持つ状態準備と非零ポンプへ進む。
