# demo from the package
install.packages('fastICA')
# Algorithm
# First, the data are centered by subtracting the mean of each column of the data matrix X.
# The data matrix is then ‘whitened’ by projecting the data onto its principal component directions
# i.e. X -> XK where K is a pre-whitening matrix. The number of components can be specified by
# the user.
# The ICA algorithm then estimates a matrix W s.t XKW = S . W is chosen to maximize the negentropy approximation under the constraints that W is an orthonormal matrix. This constraint ensures that the estimated components are uncorrelated. The algorithm is based on a fixed-point
# iteration scheme for maximizing the neg-entropy.
# Projection Pursuit
# In the absence of a generative model for the data the algorithm can be used to find the projection
# pursuit directions. Projection pursuit is a technique for finding ‘interesting’ directions in multidimensional datasets. These projections and are useful for visualizing the dataset and in density
# estimation and regression. Interesting directions are those which show the least Gaussian distribution, which is what the FastICA algorithm does.
# dev.off()
#---------------------------------------------------
#Example 1: un-mixing two mixed independent uniforms
#---------------------------------------------------
S <- matrix(runif(10000), 5000, 2) # sourse
A <- matrix(c(1, 1, -1, 3), 2, 2, byrow = TRUE) # mixing matrix
X <- S %*% A
a <- fastICA::fastICA(X, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1,
             method = "C", row.norm = FALSE, maxit = 200,
             tol = 0.0001, verbose = TRUE)
par(mfrow = c(1, 3))
plot(a$X, main = "Pre-processed data")
plot(a$X %*% a$K, main = "PCA components")
plot(a$S, main = "ICA components")
#--------------------------------------------
#Example 2: un-mixing two independent signals
#--------------------------------------------
S <- cbind(sin((1:1000)/20), rep((((1:200)-100)/100), 5)) 
dim(S) # 1000 times 2
A <- matrix(c(0.291, 0.6557, -0.5439, 0.5572), 2, 2)
X <- S %*% A
a <- fastICA::fastICA(X, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1,
             method = "R", row.norm = FALSE, maxit = 200,
             tol = 0.0001, verbose = TRUE)
par(mfcol = c(2, 3))
plot(1:1000, S[,1 ], type = "l", main = "Original Signals",
     xlab = "", ylab = "")
plot(1:1000, S[,2 ], type = "l", xlab = "", ylab = "")
plot(1:1000, X[,1 ], type = "l", main = "Mixed Signals",
     xlab = "", ylab = "")
plot(1:1000, X[,2 ], type = "l", xlab = "", ylab = "")
plot(1:1000, a$S[,1 ], type = "l", main = "ICA source estimates",
     xlab = "", ylab = "")
plot(1:1000, a$S[, 2], type = "l", xlab = "", ylab = "")
#-----------------------------------------------------------
#Example 3: using FastICA to perform projection pursuit on a
# mixture of bivariate normal distributions
#-----------------------------------------------------------
if(require(MASS)){
  x <- mvrnorm(n = 1000, mu = c(0, 0), Sigma = matrix(c(10, 3, 3, 1), 2, 2))
  x1 <- mvrnorm(n = 1000, mu = c(-1, 2), Sigma = matrix(c(10, 3, 3, 1), 2, 2))
  X <- rbind(x, x1)
  a <- fastICA::fastICA(X, 2, alg.typ = "deflation", fun = "logcosh", alpha = 1,method = "R", row.norm = FALSE, maxit = 200,
               tol = 0.0001, verbose = TRUE)
  par(mfrow = c(1, 3))
  plot(a$X, main = "Pre-processed data")
  plot(a$X %*% a$K, main = "PCA components")
  plot(a$S, main = "ICA components")
}

a$X
a$K
a$W^2
a$A
a$S


##### EEG Dataset
# パッケージの読み込み
install.packages("PEIP")
library(eegkit) # EEGデータの処理に便利なパッケージ
library(fastICA) # ICAの実行に必要なパッケージ

##########   EXAMPLE   ##########
# https://rdrr.io/cran/eegkit/man/eegica.html
# get "c" subjects of "eegdata" data
data(eegdata)
eegdata[,1]
# > eegdata[1,]
# subject group condition trial channel time voltage
# 1 co2a0000364     a        S1     0     FP1    0  -8.921
# unique(eegdata$group)
# [1] a c
# Levels: a c
eegdata[2,]

# URLのデモ
idx <- which(eegdata$group=="c")
eegdata <- eegdata[idx,]
# get average data (across subjects)
eegmean <- tapply(eegdata$voltage,list(eegdata$channel,eegdata$time),mean)
dim(eegmean)
# > dim(eegmean)
# [1]  64 256 feature sample size

# remove ears and nose
acnames <- rownames(eegmean)
# rownames(eegmean)
# [1] "AF1" "AF2" "AF7" "AF8" "AFZ" "C1"  "C2"  "C3"  "C4"  "C5"  "C6"  "CP1" "CP2" "CP3" "CP4" "CP5" "CP6" "CPZ" "CZ"  "F1"  "F2" 
# [22] "F3"  "F4"  "F5"  "F6"  "F7"  "F8"  "FC1" "FC2" "FC3" "FC4" "FC5" "FC6" "FCZ" "FP1" "FP2" "FPZ" "FT7" "FT8" "FZ"  "nd"  "O1" 
# [43] "O2"  "OZ"  "P1"  "P2"  "P3"  "P4"  "P5"  "P6"  "P7"  "P8"  "PO1" "PO2" "PO7" "PO8" "POZ" "PZ"  "T7"  "T8"  "TP7" "TP8" "X"  
# [64] "Y" 
idx <- c(which(acnames=="X"),which(acnames=="Y"),which(acnames=="nd"))
eegmean <- eegmean[-idx,]

# get spatial coordinates (for plotting)
data(eegcoord)
cidx <- match(rownames(eegmean),rownames(eegcoord))

# temporal ICA with 4 components
eegica
# function (X, nc, center = TRUE, maxit = 100, tol = 1e-06, Rmat = diag(nc), 
#           type = c("time", "space"), method = c("imax", "fast", "jade"),
icatime <- eegica(eegmean,4)
icatime$vafs
# icatime$vafs Variance-accounted-for by each component.
# [1] 0.58986476 0.19985796 0.11850832 0.04495305
# EEGデータにおけるVariance-accounted-for by each componentとは、
# EEGデータの分散を各成分（独立成分分析（ICA）で抽出された信号源）がどれだけ説明しているかを示す指標。
# EEGデータから各成分を引いたときの分散と、元のEEGデータの分散との比率で計算する。
# 例えば、ある成分がEEGデータの全てを説明している場合、EEGデータからその成分を引くと分散は0になり、
# Variance-accounted-for by each componentは100％になる。逆に、ある成分がEEGデータに全く影響を与えていない場合、
# EEGデータからその成分を引いても分散は変わらず、Variance-accounted-for by each componentは0%になる。
# この指標は、各成分がEEGデータにどれだけ寄与しているかを評価するために用いられます。
# quartz()
par(mfrow=c(4,2))
tseq <- (0:255)*1000/255
for(j in 1:4){
  par(mar=c(5.1,4.6,4.1,2.1))
  # グラフの周りの余白を指定するパラメータです。`mar` は、下、左、上、右の順にグラフと図枠の間の空白を表す数値ベクトルです¹。
  sptitle <- bquote("VAF:  "*.(round(icatime$vafs[j],4))) 
  # 変数 `sptitle` に、独立成分分析（ICA）で得られた成分の分散説明率（VAF）を表す文字列を代入するコードです。`bquote` は、式中に変数の値を埋め込む関数です²。`round(icatime$vafs[j],4)` は、`icatime$vafs[j]` の値を小数点以下4桁に丸める関数です³。
  eegtime(tseq,icatime$S[,j],main=bquote("Component  "*.(j)),cex.main=1.5)
  # EEGの時間経過をプロットする関数です⁴。引数として、時間軸のベクトル `tseq` 、ICAで得られた成分の電圧ベクトル `icatime$S[,j]` 、プロットのタイトルとして成分番号を表す文字列 `bquote("Component  "*.(j))` 、タイトルの文字サイズ `cex.main=1.5` を与えています。
  eegspace(eegcoord[cidx,4:5],icatime$M[,j],main=sptitle)
  # EEGの空間分布をプロットする関数です⁵。引数として、電極座標の行列 `eegcoord[cidx,4:5]` 、ICAで得られた成分の空間マップベクトル `icatime$M[,j]` 、プロットのタイトルとしてVAFを表す文字列 `sptitle` を与えています。
}

# spatial ICA with 4 components
icaspace <- eegica(eegmean,4,type="space")
icaspace$vafs
# quartz()
par(mfrow=c(4,2))
tseq <- (0:255)*1000/255
for(j in 1:4){
  par(mar=c(5.1,4.6,4.1,2.1))
  sptitle <- bquote("VAF:  "*.(round(icaspace$vafs[j],4)))
  eegtime(tseq,icaspace$M[,j],main=bquote("Component  "*.(j)),cex.main=1.5)
  eegspace(eegcoord[cidx,4:5],icaspace$S[,j],main=sptitle)
}

