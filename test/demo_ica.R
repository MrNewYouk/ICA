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



# パッケージの読み込み
install.packages("PEIP")
library(eegkit) # EEGデータの処理に便利なパッケージ
library(fastICA) # ICAの実行に必要なパッケージ
library(PEIP)

# データの読み込み
# ここでは、OpenBCI Communityで紹介されているPublicly Available EEG Datasetsの中から、
# Left/Right Hand MI: Includes 52 subjects (38 validated subjects with discriminative features), results of physiological and psychological questionnares, EMG Datasets, location of 3D EEG electrodes, and EEGs for non-task related states
# というデータセット¹を使用します。
# データセットには、S01〜S52までの被験者ごとに、EEGデータ（A01T.mat〜A03T.mat）、EMGデータ（B01T.mat〜B03T.mat）、アンケート結果（questionnaire.mat）、電極位置（electrodes.mat）が含まれています。
# ここでは、被験者S01のEEGデータ（A01T.mat）を読み込んでみます。
data <- loadmat("A01T.mat") # データの読み込み
eeg <- data$X # EEGデータの抽出
eeg <- eeg[,,1] # 最初の試行のみ使用
eeg <- t(eeg) # 行列の転置（fastICAの入力として適切な形にするため）

# ICAの実行
ica <- fastICA(eeg, n.comp = 22) # 22個の独立成分を求める
S <- ica$S # 独立成分信号
A <- ica$A # 混合行列

# 結果の可視化
par(mfrow = c(2,1)) # 2行1列にプロットを並べる
plot(eeg[1,], type = "l", main = "Original EEG signal (channel 1)") # 元のEEG信号（チャンネル1）をプロット
plot(S[,1], type = "l", main = "Independent component signal (component 1)") # 独立成分信号（成分1）をプロット









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
idx <- which(eegdata$group=="c")
eegdata <- eegdata[idx,]
# get average data (across subjects)
eegmean <- tapply(eegdata$voltage,list(eegdata$channel,eegdata$time),mean)

# remove ears and nose
acnames <- rownames(eegmean)
idx <- c(which(acnames=="X"),which(acnames=="Y"),which(acnames=="nd"))
eegmean <- eegmean[-idx,]

# get spatial coordinates (for plotting)
data(eegcoord)
cidx <- match(rownames(eegmean),rownames(eegcoord))

# temporal ICA with 4 components
icatime <- eegica(eegmean,4)
icatime$vafs
# quartz()
par(mfrow=c(4,2))
tseq <- (0:255)*1000/255
for(j in 1:4){
  par(mar=c(5.1,4.6,4.1,2.1))
  sptitle <- bquote("VAF:  "*.(round(icatime$vafs[j],4)))
  eegtime(tseq,icatime$S[,j],main=bquote("Component  "*.(j)),cex.main=1.5)
  eegspace(eegcoord[cidx,4:5],icatime$M[,j],main=sptitle)
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

