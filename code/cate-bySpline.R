library(tidyverse)
library(ggplot2)
library(MASS)
library(glmnet)
library(mgcv)

library(gt)

## シュミレーション設定----
# 初期化
rm(list = ls(all.names = TRUE))　

# set.seed(12345)
# サンプルサイズ
n = 1000

# シミュレーションの繰り返し数
kk_T <- 50
# kk_T <- 500

pb <- txtProgressBar(min = 1, max = kk_T, style = 3)

# 共変量行列の次元
# 100以上でglmはエラーになる
nval <- 2

# results.beta_hat <- data.frame(matrix(vector(), kk_T, nval))
results.corr <- data.frame(matrix(NaN,nrow = kk_T, ncol=2))
colnames(results.corr) <- c("lm", "gam")

results.x0 <- data.frame(matrix(NaN, nrow = kk_T, ncol = 2))
colnames(results.x0) <- c("lm", "gam")

results.x05 <- data.frame(matrix(NaN, nrow = kk_T, ncol = 2))
colnames(results.x05) <- c("lm", "gam")

results.x1 <- data.frame(matrix(NaN, nrow = kk_T, ncol = 2))
colnames(results.x1) <- c("lm", "gam")


for (i in 1:kk_T) {
  setTxtProgressBar(pb,i)
  
  ## 共変量の分布の設定、乱数の発生----
  mu<- rep(0,nval)
  
  # 相関パラメータ
  rho <- 0
  
  Sigma <- diag(nval)
  
  XX <- mvrnorm(n, mu, Sigma)
  
  ## 処置の分布の設定,生成----
  
  # T_lm <- 0.2+XX[,1]+XX[,2]+XX[,3]+XX[,4]+XX[,5]+XX[,6]+XX[,7]+XX[,8]+XX[,9]+XX[,10]
  # p_T <- exp(T_lm)/(1+exp(T_lm))
  
  TT <- rep(0,n)
  for (j in 1:n) {
    # TT[j] <- rbinom(1,size=1,prob=p_T[j])
    TT[j] <- rbinom(1,size=1,prob=0.5)
    if(TT[j]==0){
      TT[j]=-1
    }
  }
  
  
  ## アウトカムの回帰パラメータの設定----
  
  # 交互作用項
  # beta_0 <- replace(rep(0,nval + 1), c(1,2,3,4,5), c(0.4, 0.8, -0.8, 0.8, -0.8))
  beta_0 <- c(0.4, 0.8, 0.8)
  # beta_0 <- c(0.4, 0.8)
  
  # 主効果項
  # alpha_0 <- replace(rep(0,nval + 1), 
  #                    c(1,4,5,6,7,8,9,10,11),
  #                    c(sqrt(6)^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1)
  # )
  
  alpha_0 <- replace(rep(0,nval + 1),
                     c(1,2,3),
                     c(sqrt(6)^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1)
  )
  
  # alpha_0 <- c(sqrt(6)^-1, (2*sqrt(6))^-1)
  
  # 連続の場合
  # Y_lm <- (rep(alpha_0[1], n) + (XX %*% alpha_0[-1]))^2 + (rep(beta_0[1], n) + (XX %*% beta_0[-1]) + (XX[,1]*XX[,2] *0.8))*TT + sqrt(2)*rnorm(n)
  # p_Y <- exp(Y_lm)/(1+exp(Y_lm))
  # 
  # score_true <- 1.6 * (0.5 + XX[,1] - XX[,2] + XX[,3] - XX[,4] + XX[,1]*XX[,2])
  
  # ２値の場合
  Y_lm <- (rep(alpha_0[1], n) + (XX %*% alpha_0[-1]))^2 + (rep(beta_0[1], n) + (XX %*% beta_0[-1]) + (XX[,1]*XX[,2] *0.8))*TT
  # Y_lm <- (rep(alpha_0[1], n) + (XX %*% alpha_0[-1]))^2 + (rep(beta_0[1], n) + (XX %*% beta_0[-1]))*TT
  # Y_lm <- (rep(beta_0[1],n) + beta_0[2]*XX[,1])*TT
  
  # Y <- ifelse(Y_lm >=0, 1,0)
  
  p_Y <- exp(Y_lm)/(1+exp(Y_lm))
  
  Y <- rep(0, n)
  for(j in 1:n){
    Y[j] <-rbinom(1, size=1, prob=p_Y[j])
  }
  
  
  # lm <- 1.6 * (0.5 + XX[,1] - XX[,2] + XX[,3] - XX[,4] + XX[,1]*XX[,2])
  # bx <- 2*(beta_0[1] + beta_0[2]*XX[,1])
  
  true_score_fn <- function(XX){
    bx <- 2*(beta_0[1] + beta_0[2]*XX[,1] + beta_0[3]*XX[,2] + 0.8*XX[,1]*XX[,2])
    # bx <- 2*(beta_0[1] + beta_0[2]*XX[,1] )
    score_true <- (exp(bx/2)-1) / (exp(bx/2)+1)
    return(score_true)
  }
  
  score_true <- true_score_fn(XX)
  
  
  
  # score_true <- (1-exp(bx)) / (1+exp(bx))
  
  
  
  ##　誤判別を含むアウトカムデータの生成----
  # Misclassified outcomesの割合の設定
  # p_Y <- exp(Y_lm)/(1+exp(Y_lm))
  # pp10 <- 0; pp01 <- 0
  # 
  # Y <- rep(0,n)
  # for(j in 1:n){
  #   p_Yc <- pp10*(1-p_Y[j])+(1-pp01)*p_Y[j]
  #   Y[j] <-rbinom(1, size=1, prob=p_Yc)
  # }

  
  # 感度・特異度に関する設定
  # thet0 <- -2.5; thet1 <- 2; thet2 <- 1; thet3 <- -1; thet4 <- -0.2;
  # thet0 <- -2.5; thet1 <- 2; thet2 <- 1;
  
  # 修正共変量
  
  
  W <- cbind(rep(1, nrow(XX)), XX)
  # W_s <- naturalSpline(W, intercept = FALSE)
  W_star <- W* TT/2
  
  # Y_starの生成
  # Y_star_lm <- thet0 + thet1 * Y + thet2 * W_star[,2] + thet3 * W_star[,3] + thet4 * W_star[,4]
  # Y_star_lm <- thet0 + thet1 * Y + thet2 * W_star[,2]
  # p_Y_star <- exp(Y_star_lm)/(1+exp(Y_star_lm))
  # 
  # Y_star <- rep(0,n)
  # for(j in 1:n){
  #   Y_star[j] <-rbinom(1, size=1, prob=p_Y_star[j])
  # }
  # Y_star <- ifelse(Y_star_lm >=0, 1,0)
  
  
  # 簡単なケース
  # SEx <- rep(0.85, n); SPx <- rep(0.8, n);
  # # 
  # 
  # p_Y_star <- SEx*p_Y + (1-SPx)*(1-p_Y)
  # Y_star <- rep(0,n)
  # for(j in 1:n){
  #   Y_star[j] <-rbinom(1, size=1, prob=p_Y_star[j])
  # }
  
  ## バリデーションデータの分割----
  
  # Y_star_v <- Y_star[1:20]; Y_star_m <- Y_star[21:100];
  # Y_v <- Y[1:20]; Y_m <- Y[21:100]; # Y_mは観測できない
  # XX_v <- XX[1:20, ]; XX_m <- XX[21:100, ];
  # TT_v <- TT[1:20]; TT_m <- TT[21:100];
  # W_star_v <- W_star[1:20,]; W_star_m <- W_star[21:100,];
  # W_v <- W[1:20,]; W_m <- W[21:100,];

  ## 感度・特異度の推定----
  
  # YW_v <- cbind(Y_v, W_star_v)
  # 
  # # glm(Y_star ~ Y + W_star[,2:4], family=binomial)
  # model <- glm(Y_star_v ~ YW_v, family=binomial) #推定がおかしい、完全分離している
  # 
  # YW_v1 <- YW_v[YW_v[,1]==1,] %>% data.frame()
  # YW_v0 <- YW_v[YW_v[,1]==0,] %>% data.frame()
  # 
  # SEx = predict(model, newdata=YW_v1, type="response")
  # SPx = 1 - predict(model, newdata=YW_v0, type="response")
  
  
  ## 尤度関数の定義----
  
  # like_fun <- function(par_list){
  #   # W_star_lm <- 0
  #   # index <- 1
  #   # for (par in par_list) {
  #   #   W_star_lm <- W_star_lm + W_star[,index]%*%par
  #   #   index <- index +1
  #   # }
  #   
  #   W_star_lm <- W_star %*% par_list
  # 
  #   py1 <- exp(W_star_lm)/(1+exp(W_star_lm))
  #   
  #   # 尤度
  #   term1mn <- ((1-SPx)*(1-py1) + SEx*py1)^Y_star
  #   term2mn <- (SPx*(1-py1) + (1-SEx)*py1)^(1-Y_star)
  #   
  #   term1v <- (SEx*py1)^(Y[1:20]*Y_v)
  #   term2v <- ((1 - SPx)*(1-py1))^(Y_star_v*(1 - Y_v))
  #   term3v <- ((1 - SEx)*py1)^((1 - Y_star_v)*Y_v)
  #   term4v <- (SPx*(1-py1))^((1 - Y_star_v)*(1 - Y_v))
  #   
  #   intval <- c(rep(1,20),rep(0,80))
  #   
  #   like = ((term1mn*term2mn)^(1 - intval))*((term1v*term2v*term3v*term4v)^intval)
  #   return(-sum(like))
  #   
  # }
  
  
  ## パラメータ推定----
  
  # dfの準備
  par_list <- c()
  for ( j in 1:ncol(W_star)) {
    par_list <- c(par_list, paste("W_star",j,sep = ""))
  }

  colnames(W)<-par_list
  colnames(W_star)<-par_list
  
  W_df <- data.frame(W)

  YW_star_df <- data.frame(Y,W_star)
  YW_df <- data.frame(Y,W)
  
  # モデルの作成
  # ml <- glm(Y~W_star-1, family=binomial)
  lm.model <- gam(Y ~ W_star1 + W_star2 + W_star3 + W_star2*W_star3, data=YW_star_df, family=binomial(link="logit"), drop.intercept=TRUE)        #　線形回帰
  # gam.model <- gam(Y ~ W_star1 + s(W_star2) +s(W_star3) + s(W_star2,W_star3), data = YW_star_df, family=binomial(link="logit"), drop.intercept = TRUE)    #　平滑化スプライン
  gam.model <- gam(Y ~ W_star1 + poly(W_star2, degree = 10, raw = TRUE)+ poly(W_star3, degree = 10, raw = TRUE)+poly(W_star2*W_star3,degree = 10,row=TRUE), data = YW_star_df, family=binomial(link="logit"), drop.intercept = TRUE)    #　平滑化スプライン
  
  # スコア計算
  lm.pred <- predict(lm.model, newdata=W_df, type="link")
  gam.pred <- predict.gam(gam.model, newdata = W_df, type="link")
  
  
  score_fn <- function(pred){
    score <- (exp(pred/2)-1) / (exp(pred/2)+1)
    return(score)
  }
  
  lm.score <- score_fn(lm.pred)
  gam.score <- score_fn(gam.pred)
  
  lm.corr <- cor(score_true, lm.score)
  gam.corr <- cor(score_true, gam.score)
  
  
  # x=(0,0)
  w0 <- data.frame(c(1,0,0)) %>% t %>% data.frame()
  colnames(w0)<-par_list
  lm.pred_x0 <- predict(lm.model, newdata=w0, type="link")
  gam.pred_x0 <- predict.gam(gam.model, newdata = w0, type="link")
  lm.score_x0 <- score_fn(lm.pred_x0)
  gam.score_x0 <- score_fn(gam.pred_x0)
  
  # x=(0.5,0.5)
  w05 <- data.frame(c(1,0.5,0.5)) %>% t %>% data.frame()
  colnames(w05)<-par_list
  lm.pred_x05 <- predict(lm.model, newdata=w05, type="link")
  gam.pred_x05 <- predict.gam(gam.model, newdata = w05, type="link")
  lm.score_x05 <- score_fn(lm.pred_x05)
  gam.score_x05 <- score_fn(gam.pred_x05)
  
  # x=(1,1)
  w1 <- data.frame(c(1,1,1)) %>% t %>% data.frame()
  colnames(w1)<-par_list
  lm.pred_x1 <- predict(lm.model, newdata=w1, type="link")
  gam.pred_x1 <- predict.gam(gam.model, newdata = w1, type="link")
  lm.score_x1 <- score_fn(lm.pred_x1)
  gam.score_x1 <- score_fn(gam.pred_x1)
  
  
  
  # 結果保存
  results.corr[i,1] <- lm.corr
  results.corr[i,2] <- gam.corr
  
  results.x0[i,1] <- lm.score_x0
  results.x0[i,2] <- gam.score_x0
  
  results.x05[i,1] <- lm.score_x05
  results.x05[i,2] <- gam.score_x05
  
  results.x1[i,1] <- lm.score_x1
  results.x1[i,2] <- gam.score_x1
  
  
  
  
  
}



# true score
beta_0 <- c(0.4, 0.8, 0.8)
true_score_fn <- function(XX){
  bx <- 2*(beta_0[1] + beta_0[2]*XX[,1] + beta_0[3]*XX[,2] + 0.8*XX[,1]*XX[,2])
  # bx <- 2*(beta_0[1] + beta_0[2]*XX[,1])
  score_true <- (exp(bx/2)-1) / (exp(bx/2)+1)
  return(score_true)
}

#x=0,0
x0 <- data.frame(c(0,0)) %>% t 
true_x0 <- true_score_fn(x0)
#x=1,1
x05 <- data.frame(c(0.5,0.5)) %>% t
true_x05 <- true_score_fn(x05)
#x=2,2
x1 <- data.frame(c(1,1)) %>% t
true_x1 <- true_score_fn(x1)




# ロジスティック回帰の差と同じ結果に

# alpha_0 <- replace(rep(0,nval + 1),
#                    c(1,2,3),
#                    c(sqrt(6)^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1)
# )

# score_fn <- function(XX){
#   n <- dim(XX)[1]
#   Yp1_lm <- (rep(alpha_0[1], n) + (XX %*% alpha_0[-1]))^2 + (rep(beta_0[1], n) + (XX %*% beta_0[-1]) + (XX[,1]*XX[,2] *0.8)) 
#   Ym1_lm <- (rep(alpha_0[1], n) + (XX %*% alpha_0[-1]))^2 - (rep(beta_0[1], n) + (XX %*% beta_0[-1]) + (XX[,1]*XX[,2] *0.8)) 
#   p_Yp1 <- exp(Yp1_lm)/(1+exp(Yp1_lm))
#   p_Ym1 <- exp(Ym1_lm)/(1+exp(Ym1_lm))
#   score <- p_Yp1 - p_Ym1
#   return(score)
# }
# 
# score_fn(x05)







theme_gtsummary_mean_sd()
result_df <- results.x0 %>% tbl_summary(digits = everything()~3) %>% 
  modify_header(label="手法") %>%
  modify_footnote(label ='真値:(0.197)') %>%
  modify_caption("x=0の時")
  
result_df %>% as_gt() %>% gtsave("./results/table/0813/res_2d_0.png")



# gtsave(export_tb, "./results/fig/res_table_ex2.png")


## アウトプット
# export_data <- results.corr

# write_csv2(export_data, file = "~/Projects/gamma_DR_ver0.1/results/ite_sesp/0709.csv")

# 
# df <- read_csv2("~/Projects/gamma_DR_ver0.1/results/0423/g0_m01.csv")
# df <- read_csv2("~/gamma_DR_ver0.1/results/tian/misbin_nval50.csv")
# df %>% summary()
# df %>% boxplot()
