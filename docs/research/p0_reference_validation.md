# P0 と小系 ITensor 基準実装の検証

実施日: 2026-09-19。区分: **実装と有限系の数値検証**。
対象は最近接スピン 1/2 kagome Heisenberg 模型、`Jxy=Jz=1`、固定 `M/Msat=1/9`。
この記録は plateau、量子化ポンプ、熱力学的な相を主張しない。

## 達成した goal

格子・磁化・flux の規約をテスト可能な Julia パッケージにし、
独立な ED で複素 U(1) DMRG の最小基準実装を照合した。
P0 は完了、P1 は小系照合まで完了。checkpoint と再開は未実装である。

実行環境は Julia 1.12.7、ITensors 0.9.31、ITensorMPS 0.4.1、KrylovKit 0.10.4。
依存関係はリポジトリの `Manifest.toml` に固定した。数値計算は BLAS 1 thread、
block-sparse threading 無効。GPU・MPI・大規模 DMRG scan は使用していない。

## 独立性と検証範囲

本体の格子は設計した 6 種類の結合テンプレートから生成する。
ED 側は Cartesian 座標と周期像における距離 1/2 から最近接結合を探す。
ED の Hamiltonian は Up を bit 1 とするスピン基底で構築し、本体の MPO、
bond generator、ゲージ関数を呼ばない。

- `Ly=3,4,6` を含む格子で site ordering、終端の配位数、結合数 `(6Lx−2)Ly`、二重計数を検証。
- seam を横切る三角形と六角形で位相和ゼロ、円周を巡る経路で位相和 `±θ` を検証。
- 9 サイト全 512 状態の MPO 行列を ED と照合。
  θ=`0,0.37,−0.71` の両ゲージで Hermiticity と `[H,Q]=0` を検証。
- seam gauge の 2π 周期、uniform gauge の `Huniform=U Hseam U†`、
  異方的な `Jxy=0.7,Jz=1.3` の規約を検証。
- N=18、Q=2、Nup=10 の疎行列を 43,758 次元で構築し Hermiticity を確認。
  18 サイトの Heisenberg 基底状態の対角化・DMRG 照合は未実施。
- 18 サイトの `Jxy=Jz=0` と符号を選んだ局所縦磁場で一意な積状態を固定し、
  非零 flux の DMRG と零 flux のエネルギー・密度・移送が一致することを確認。
  元の模型とは異なる読み出し用対照である。
- 3 列の合成密度で複数 cut の左右相殺と累積移送 `±2` の保存を確認。
  modulo 1 に丸めず、実測の基準 profile を引く。

## 9 サイト DMRG と ED

格子は `Lx=1,Ly=3`、wrap vector `3a2`、完全な単位胞の端を持つ。
`N=9,Q=1,M=1/2,Nup=5,Ndown=4`、sector 次元は 126、結合は 12 本。
MPS は複素数、初期状態は固定電荷で seed を指定する。
8 sweeps、maxdim schedule `[16,32,64]`、cutoff `1e-12`、noise `0`、
Krylov tolerance `1e-12`、Krylov dimension `20`、maxiter `10` を使用した。

| θ | seed | ED E₀ | 基底空間の次元 | 基底空間より上の有限系 gap |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 11 | −3.561967560726586 | 2 | 0.416357921998 |
| 0.37 | 11, 29 | −3.569278674796576 | 1 | 0.0161371964532 |

零 flux は二重縮退している。ED が返す最初の固有ベクトルとの密度比較は適切でない。
DMRG ベクトルを ED の基底空間に射影し、その漏れと residual を検証した上で、
規格化した射影状態と密度・`SzSz`・`S+S−` 相関を比較した。
非零 flux の一意な基底状態でも同じ手順を使った。
これらの gap はこの有限系の値であり、bulk neutral gap の推定ではない。

package test と保存用実行例の照合で、エネルギー差/site は `7e-16` 以下、局所 Sz 差は `6e-14` 以下、
相関差は `6e-14` 以下、独立 ED による residual norm は `3e-13` 以下だった。
受入条件はそれぞれ `1e-8`、`1e-6`、`1e-6`、`2e-6` であり、
期待する Hall 応答値とは無関係に設定した。

各 sweep の実測切断誤差は今回の小系ではゼロだった。9 サイトを表現できる十分な bond dimension
を使った結果であり、研究サイズで切断誤差が消えるという意味ではない。
variance の大きさは `3e-14` 以下で、差し引きの丸めにより負の微小値もそのまま保持する。
局所 Krylov の convergence flag は upstream が公開せず、確認済みとはしない。

## 切断誤差の校正で見つかった制限

非縮退の二サイト対照 `H=S₁·S₂−(3/8)Sz₁+(3/8)Sz₂` は、厳密 energy `−7/8`、
Schmidt 確率 `(4/5,1/5)` を持つ。`maxdim=1,cutoff=0,noise=0` の切断誤差 `1/5` と、
切断後の積状態の energy `−5/8` を検証する。設定 cutoff と測定誤差、
局所 eigensolver の返す energy と最終 MPS の energy を区別する校正である。

同じ設定で等方的二サイト singlet の等しい Schmidt 確率 `(1/2,1/2)` を切ると、
検証した upstream QN backend は両方のセクターを落とし、norm ゼロの状態を作った。
それでも切断誤差は `1/2` と報告された。`maxdim=2` では正常に保持される。
[独立診断の記録](data/qn_truncation_boundary.toml)を保存した。
observer は各更新後の正準中心テンソルの norm を検査し、ゼロまたは非有限なら即座に失敗とする。
異常な状態を次の計算や成功した結果に渡さない。

これは upstream の切断カーネルそのものの修正ではない。
一部の状態が残る縮退境界でも報告切断誤差が正確とは限らず、その一般的な校正は未達である。
研究サイズへ進む際は、十分な maxdim、独立 residual/variance、bond dimension 依存を調べ、
縮退境界の修正・検証を P1 の残課題として扱う。今回の 9 サイト基準値は切断を要しない幅で検証した。

## 再現と保存

```sh
julia --project=. --startup-file=no -e 'using Pkg; Pkg.instantiate()'
julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
julia --project=. --startup-file=no examples/validate_small_system.jl
```

標準 package test で **2,205 assertions が成功**した。
テスト用の全 Hilbert 空間への変換は小系だけに制限し、実際の DMRG solver は MPS/MPO を用いる。
最小実行例は θ 点を独立に最適化し、`outputs/p0-reference-*/validation.toml` に保存する。
連続的な flux trajectory ではない。出力先の引数を指定する場合、既存ファイルへの上書きを拒否する。
失敗した判定・例外も成功した点とともに記録し、判定失敗時は異常終了する。

保存項目は模型・格子・電荷・ゲージ・seed・solver 設定・数値診断・依存バージョン・
git revision・dirty flag・使用したソースと Manifest の SHA-256 に限定する。
ソースは計算前後に照合し、変更や取得失敗があれば数値が一致しても記録全体を失敗とする。
環境変数や接続情報は保存しない。MPS checkpoint の保存・再開機能ではない。

初回の保存処理はコマンド構築の不具合で数値照合後に停止した。
[保存失敗の部分記録](data/p0_initial_recording_failure.toml)を残し、処理を修正した。
[ガード追加前の保存成功記録](data/p0_reference_validation.toml)と
[最終コードでの記録](data/p0_reference_validation_guarded.toml)を区別して保持する。

## 未達事項と次の研究単位

1. 縮退した QN 切断境界の修正・校正。原子的な checkpoint 保存、site indices・格子・Q・設定を照合した再開。
2. Schmidt charge の積状態校正と、実空間移送との照合。
3. 受理・棄却・rollback を伴う逐次 flux continuation と、枝の診断。
4. 原論文の結合図を照合した既知 CSL 正対照、step/bond dimension/長さ依存。
5. その後に 1/9 plateau の有限サイズ測定と競合状態の比較。

N=9 の系には内部軸方向 cut がないため、今回の精密な ED 一致から量子化ポンプは結論できない。
SU(3)₁、Hall-active/inactive D(Z₃)、VBC、磁気秩序、gapless 状態の区別は未検証。

公開 API、docstrings、保存キーと [導入文書](../getting_started.md)は英語に揃えた。
研究記録と設計は当面日本語とし、英語リリース前に規約・検証・制限を含めて翻訳する。
