# dmrg_mps  
一次元量子系の基底エネルギーを密度行列くりこみ群 (DMRG) で計算するプログラムです。
データ構造として行列積状態 (MPS)とMatrix-product operator (MPO)を使っています。MATLABコードです。  
ハミルトニアンをmatrix-product operator (MPO)で記述することでエネルギーを高速に計算します。  

例として臨界点直上の横磁場イジングの基底エネルギーを求めるスクリプトを書きました。   他の系に使う時はtransverseIsingMPO.mを参考にしてハミルトニアンをMPO表示してください。次にminimizeE.mの23行目を書き換えて下さい。

## ドキュメント
手法・記法は[Schollwoeck(2011)](http://arxiv.org/abs/1008.3477)に従った。詳しくは[論文](http://arxiv.org/abs/1008.3477)を参照せよ。

トップレベルの呼び出しは  
`E_GS = minimizeE(J, h, N, D, isGS);`  
`isGS=1`なら基底状態を求める。`isGS=0`なら第一励起状態を求める（まだ実装していない。）

## 関数
[`minimizeE(J, h, N, D, isGS)`](https://github.com/YoshiHotta/DMRG/wiki/DMRG_Hotta#function-e_gs--minimizeej-h-n-d-isgs)  
[`contractTensors(X, rankX, indXContr, Y, rankY, indYContr)`](https://github.com/YoshiHotta/DMRG/wiki/DMRG_Hotta#function-z--contracttensorsx-rankx-indxcontr-y-ranky-indycontr)  
[`rightNormalize( mps )`](https://github.com/YoshiHotta/DMRG/wiki/DMRG_Hotta#bs--rightnormalize-mps-
)  
[`transverseIsingMPO(N, J, h)`](https://github.com/YoshiHotta/DMRG/wiki/DMRG_Hotta#hamiltonianmpo--transverseisingmpon-j-h)  
[`exactTransverseIsing(J, h, N)`](https://github.com/YoshiHotta/DMRG/wiki/DMRG_Hotta#function-e_gs--exacttransverseisingj-h-n)

## スクリプト
onCriticalPoint.m  
energyGap.m  

## テスト
test_contractTensors.m  
test_rightNormalize.m    


***

##### `function E_GS = minimizeE(J, h, N, D, isGS)` 
横磁場イジングモデルの基底エネルギーを求める関数。  
引数  
       `J, h` : モデルのパラメータ  
       `N` : サイト数  
       `D` : MPSのボンド次元  
返り値   
       `E_GS` : `isGS=1`なら基底エネルギー。`isGS=0`なら第一励起準位。  
**[バグ](https://github.com/YoshiHotta/DMRG/issues/3)**
**[未解決課題](https://github.com/YoshiHotta/DMRG/issues/4)**

MPS/MPOを使って実装している。 
本関数内で`transverseIsingMPO()`を呼び出してMPO形式のHamiltonianを取得している。  
![MPO](https://dl.dropboxusercontent.com/u/16614273/GitHubFigures/def_MPO.jpg)  
MPOで書けるHamiltonianであれば本関数を若干修正することで基底エネルギーを計算できる。MPOで書けないHamiltonianは大幅な修正をしなければ基底エネルギーを計算できない。
  
初期状態はランダムに選んでいる。

本関数内で多次元配列（テンソル）の縮約を何度も行っている。テンソルの添え字の定義は以下の通りである：  
![MPS_Index](https://dl.dropboxusercontent.com/u/16614273/GitHubFigures/defOfIndex.jpg)  
![LR_Index](https://dl.dropboxusercontent.com/u/16614273/GitHubFigures/defOfRIndex.jpg)  
図中で添え字に付いている丸nは、多次元配列のn次元目がそのボンドに対応することを表す。
例えば`M{l}(i,j,k)`は①が`i`に、②が`j`に、③が`k`に対応する。
***

##### `function Z = contractTensors(X, rankX, indXContr, Y, rankY, indYContr)`  
テンソルXの添え字indXContrとテンソルYの添え字indYContrを縮約する。  
Parameters  
`X` = A tensor  
    `indXContr`     = The index of X to contract  
    `rankX`         = Xのrank  
 `Y` = A tensor  
    `indYContr`     = The index of Y to contract  
    `rankY`         = Yのrank  
 Returns  
    `Z`        = The contracted tensor of X and Y   

Zの添え字の順番はXの縮約していない添え字が先に来て、Yの縮約していない添え字が後に来る。  
例１)  
```matlab    
X(i_1, ..., i_p, j_1, ..., j_q)  
Y(k_1, ..., k_r, j_1, ..., j_q)  
Z = contractTensors(X, p+q, [j_1, ..., j_q], Y, r+q, [j_1, ..., j_q]) 
-> sum_j X Y = Z(i_1, ..., i_p, k_1, ..., k_r)
```   

例２）
```
X(i,p,k)  
Y(r,p,s)
Z = contractTensors(X, 3, 2, Y, 3, 2)
-> Z(i,k,r,s)
```

例３）
```
X(i,p,q,j)
Y(k,l,m,p,q)
Z = contractTensors(X, 4, [2 3], Y, 5, [4 5])
-> Z(i,j,k,l,m)
```

`test_testTensors`で動作をテストした。

＊注意
MATLABは`size(T)=[n m 1]`のテンソルを勝手に`n*m`の行列に変えてしまうので、この関数は正常に動かないことがある。  
MATLABの多次元配列はリトルエンディアンである。すなわち`D*D*D`次元の多次元配列`A`は`A(i,j,k)=A(i + D*j + D^2*k)`の関係を満たす。

***

##### `function Bs = rightNormalize( mps )`  
任意のMPSをnormalizeされたright-canonical formに変換する。
入力 : `mps = {M_1, M_2, ..., M_N}`
出力 : `Bs  = {B_1, B_2, ..., B_N}`
ここでNは粒子数であり、M_iはrank-3のテンソルである。
B_iはSchollwoeckのB-matrixである。

テストは`test_rightNormalized`で行った。

***

##### `function HamiltonianMPO = transverseIsingMPO(N, J, h)`  
横磁場イジングのハミルトニアンをMatrix product operator形式で返す
出力：`HamiltonianMPO` : `cell(1,N)`で各セルは4次元配列、下図のWである。
![transverseIsing](https://dl.dropboxusercontent.com/u/16614273/GitHubFigures/TransverseIsing.png)
***

##### `function E_GS = exactTransverseIsing(J, h, N)`
横磁場イジングモデルの基底エネルギーの厳密解。
