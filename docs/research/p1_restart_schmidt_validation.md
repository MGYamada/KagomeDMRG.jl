# P1 の保存・再開と Schmidt 電荷診断の検証

実施日: 2026-09-19。区分: **実装と有限系の数値検証**。
今回の goal は、同一環境での checkpoint 再開と、flux 追跡に必要な
絶対 Schmidt 電荷の読み出しを整備することである。
1/9 plateau、量子化ポンプ、熱力学的な相についての新しい判定は行わない。

## 実装と保存契約

`save_checkpoint`、`load_checkpoint`、`resume_dmrg` を追加した。
保存対象は完了した DMRG 点の MPS と、実測零 flux profile を与える基準 MPS。
格子の全 site・bond・winding・couplings、縦磁場、gauge、`Q=2Sz_total`、
site identity、未折返し theta、呼出側が指定した theta 列、solver 設定、
実測 energy・variance・sweep ごとの切断誤差を allowlist で記録する。
variance 未測定とゼロは区別する。接続設定や環境変数の一括保存は行わない。

試行と受理済み snapshot は別ディレクトリに置く。各 snapshot は新規作成し、
metadata、payload、checksum を閉じた後、同一の親への rename で公開する。
既存の受理済み状態を上書きする経路はない。
これは読者から見た原子的公開であり、電源断に対する耐久性の保証ではない。

読込は schema、サイズ、SHA-256、Julia・package・architecture、ソースと
dependency manifest の一致を確認してから Julia payload を deserialize する。
MPS の norm、各 tensor の charge conservation、全 charge、site indices、
局所 Sz、基準 MPS と基準密度を検査し、保存した模型から Hamiltonian を
再構築して energy を照合する。元の Sz profile が一致していても、
局所位相変更で energy が変わった状態は拒否する。

形式は信頼済みローカル計算の同一環境 snapshot に限定する。
Julia の型・package 実装に依存するため、長期交換形式として扱わない。
[Serialization の公式契約](https://docs.julialang.org/en/v1/stdlib/Serialization/)を参照。
HDF5 のような可搬形式への移行と、基準 MPS の重複を省く I/O 最適化は未実装である。

再開は保存した solver 設定と正確な site indices から、新しい DMRG batch を開始する。
Hamiltonian と environments は新しい theta で作り直す。
途中 sweep の再開、実行途中の Krylov 状態の保存、自動の枝判定は行わない。
`status=:accepted` と theta 列は呼出側の指定であり、断熱連続性を意味しない。

## Schmidt 電荷の独立校正

`schmidt_diagnostics(psi,b)` は MPS prefix `1:b` の全 Schmidt 確率を求める。
入力のコピーを正準化し、cutoff・maxdim を指定せず block-sparse SVD を行う。
各確率に対応する整数電荷は、左 tensor 群の flux の和から
境界 link の向き付き QN を引いて求める。物理 Sz はこの電荷の半分である。
密度を使って電荷の offset を調整せず、両者を独立に照合する。

次の校正が通過した。

- 左磁化の異なる積状態、全 Up・全 Down、非零の全 Q。
- 複素二成分状態の既知確率 `0.8/0.2`、縮退 `0.5/0.5`、積状態極限。
- 自然対数 entropy、物理 Sz の平均・分散、規格化前 norm²。
- 9 サイト、Q=1 の独立スピン基底からの charge ごとの SVD との一致。
  theta=0 と 0.37 の ED 固有ベクトルの複素重ね合わせを用い、全 8 cut を確認した。
- canonical center の移動、link QN の `+7` shift と arrow 反転に対する不変性。
- 入力 MPS の tensor と正準中心情報が変更されないこと。

追加の幾何学的校正は `Lx=3,Ly=3,N=27,Q=3` の低 bond dimension 状態で行った。
端の Up/Dn を交換する 2 つの積状態を、移動側の重み `p=0,0.2,1` で重ね合わせた。
2 つの切断 `c=1,2`、すなわち MPS bond `9,18` の両方で、
`Delta Sz_right = p = -Delta Schmidt Sz_left` を `1e-13` の許容誤差で確認した。
全磁化保存、左右相殺、切断間一致、逆向き比較、零変化対照も通過した。
これは既知の合成状態による読み出し校正であり、Hamiltonian の flux 応答でも
分数量子化の証拠でもない。dense 27-site 化や DMRG は使っていない。

## 統合テストと環境

Julia 1.13.0、ITensors 0.9.31、ITensorMPS 0.4.1、NDTensors 0.4.31、
KrylovKit 0.10.4、Julia/BLAS 各 1 thread、block-sparse threading 無効。
`Manifest-v1.13.toml` を使用した。

```sh
julia --project=. --startup-file=no --threads=1 -e 'using Pkg; Pkg.test()'
```

全 **3,782 assertions が成功**した。内訳は既存検証 2,205、Schmidt 検証 1,477、
checkpoint 検証 100。geometry・磁場・gauge・site ID・solver・flux の不一致、
破損ファイル、schema/runtime/source の不一致、norm・charge・profile・energy と
基準状態の不整合、試行 snapshot の通常再開拒否を含む。
破損 metadata の checksum を意図的に更新したケースでも、意味的な検査で拒否した。

独立レビューでは、最初に energy と MPS の関連、および MPS のない baseline を
受理する穴を発見した。Hamiltonian 期待値検査と基準 MPS 保存で修正し、
DMRG を使わない別の 12 assertions でも拒否動作と正常な再保存を確認した。
これらの独立 probe は上記の package test 件数には含めない。

SHA と Serialization の標準ライブラリ依存を追加したため、Julia 1.12.7 用の
`Manifest.toml` もその runtime の Pkg で解決し直した。外部 package の版は変更していない。
拡張後の suite を Julia 1.12.7 では実行していない。

## 別 process からの再開

再現コマンド:

```sh
julia +1.13 --project=. --startup-file=no --threads=1 examples/validate_restart.jl
```

9 サイトの最近接等方 Heisenberg 模型、Q=1、seam gauge、seed=11 を使う。
手動 fixture の flux 列は `0 → 0.37 → 2pi+0.71`。
中間点を保存し、同じ Julia 実行ファイルを使う別 process で最後の点に進め、
保存せず続けた結果と照合する。この大きな step は未折返し theta の保存と
Hamiltonian 更新を試すためのもので、断熱的な走行ではない。

零 flux の ED との比較は、任意の単一固有ベクトルでなく二重縮退した基底空間への
射影を用いる。非零 flux の点も独立 ED の energy・residual・密度と比較する。
局所 solver の energy と最終 MPS の energy、設定 cutoff と実測切断誤差は区別する。

実行は成功し、[機械可読な検証記録](data/p1_restart_validation.toml)を保存した。
コードと manifest の SHA-256 は計算前後で一致し、記録の `code.status` は `verified`。
作業ツリーの変更は `worktree_dirty=true` と記録されており、git revision だけで
今回のソースを表しているとはしない。

| 比較量 | 実測値 | 受入上限 |
| --- | ---: | ---: |
| 再開と直接継続の energy 差 | 0 | 1e-10 |
| 再開と直接継続の最大局所 Sz 差 | 0 | 1e-7 |
| `abs(1-abs(overlap))` | 4.44e-16 | 1e-10 |
| 再開と直接継続の Schmidt 左 Sz 差 | 0 | 1e-7 |
| ED energy 誤差/site（4 記録点の最大） | 1.98e-16 未満 | 1e-8 |
| 独立 ED residual norm（最大） | 1.06e-13 未満 | 2e-6 |
| ED 基底空間射影と局所 Sz の差（最大） | 5.28e-14 未満 | 1e-6 |
| Schmidt 左 Sz と直接和の差（最大） | 3.06e-16 未満 | 2e-12 |

ここで 4 記録点は baseline、中間保存点、直接継続した終点、別 process で再開した
同じ終点であり、4 つの独立な flux 点ではない。site indices、全 charge、基準 MPS・
profile、未折返し theta 列の保存と、元の受理 snapshot の不変も確認した。
物理点の受入誤差は期待する Hall 応答から決めていない。

## 残る研究課題

この段階で整備したのは保存・再開と診断基盤である。一般的な縮退 QN 切断境界の
修正・誤差校正は未達であり、既存 observer による空状態の検出は kernel 修正ではない。
次はこの制限を扱い、P2 の複数診断による受理・棄却、刻みの適応、rollback を
独立 ED と零応答対照で検証する。

既知 CSL の pump、長さ・幅・bond dimension・step 依存の収束、1/9 plateau の
ポンプと相の判定は未実施。9 サイト cylinder は内部の軸方向幾何 cut を持たず、
今回の高精度 ED 一致から軸方向 pump を推論できない。
