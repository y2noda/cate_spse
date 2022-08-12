library(tidyverse)
library(ggplot2)
library(MASS)
library(glmnet)
library(mgcv)
library(gee)
library(gt)
library(gtsummary)
library(geepack)
library(nleqslv)
library(webshot)
  

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

miss_type_list <- c('min','mid','max')

# resultの箱
# results.corrected <- data.frame(matrix(NaN,nrow = kk_T, ncol=4))
# colnames(results.corrected) <- c("beta1", "beta2", "beta3", "miss_type")
# 
# results.naive <- data.frame(matrix(NaN,nrow = kk_T, ncol=4))
# colnames(results.naive) <- c("beta1", "beta2", "beta3","miss_type")

results.corrected <- data.frame(matrix(NaN,nrow = kk_T, ncol=3))
colnames(results.corrected) <- c("beta1", "beta2", "beta3")

results.naive <- data.frame(matrix(NaN,nrow = kk_T, ncol=3))
colnames(results.naive) <- c("beta1", "beta2", "beta3")


sigmoid <- function(x){
  return( 1 / (1 + exp(-x)) )
}

# for (miss_type in miss_type_list) {
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
  
  # 主効果項
  # alpha_0 <- replace(rep(0,nval + 1), 
  #                    c(1,4,5,6,7,8,9,10,11),
  #                    c(sqrt(6)^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1)
  # )
  
  # alpha_0 <- replace(rep(0,nval + 1),
  #                    c(1,2,3),
  #                    c(sqrt(6)^-1, (2*sqrt(6))^-1, (2*sqrt(6))^-1)
  # )
  

  alpha_0 <- c(0.2, 0.7, 0.7)
  
  # 連続の場合
  # Y_lm <- (rep(alpha_0[1], n) + (XX %*% alpha_0[-1]))^2 + (rep(beta_0[1], n) + (XX %*% beta_0[-1]) + (XX[,1]*XX[,2] *0.8))*TT + sqrt(2)*rnorm(n)
  # p_Y <- exp(Y_lm)/(1+exp(Y_lm))
  # 
  # score_true <- 1.6 * (0.5 + XX[,1] - XX[,2] + XX[,3] - XX[,4] + XX[,1]*XX[,2])
  
  # ２値の場合
  # Y_lm <- (rep(alpha_0[1], n) + (XX %*% alpha_0[-1]))^2 + (rep(beta_0[1], n) + (XX %*% beta_0[-1]) + (XX[,1]*XX[,2] *0.8))*TT + sqrt(2)*rnorm(n)
  # Y_lm <- (rep(beta_0[1],n) + beta_0[2]*XX[,1])*TT
  
  # 主効果のみ
  Y_lm <- rep(alpha_0[1], n) + (XX %*% alpha_0[-1])
  
  
  # Y <- ifelse(Y_lm >=0, 1,0)
  
  p_Y <- exp(Y_lm)/(1+exp(Y_lm))
  
  Y <- rep(0, n)
  for(j in 1:n){
    Y[j] <-rbinom(1, size=1, prob=p_Y[j])
  }
  
  
  ##　誤判別を含むアウトカムデータの生成----
  
  # 感度・特異度に関する設定
  # thet0 <- -2.5; thet1 <- 2; thet2 <- 1; thet3 <- -1; thet4 <- -0.2;
  # thet0 <- -2.5; thet1 <- 2; thet2 <- 1;
  
  # 修正共変量
  # W <- cbind(rep(1, nrow(XX)), XX)
  # W_star <- W * TT/2
  
  # Y_starの生成
  
  # 簡単なケース
  # SEx <- rep(1.0, n); SPx <- rep(1.0, n);
  SEx <- rep(0.5, n); SPx <- rep(0.55, n);
  
  # SEx <- switch(miss_type,
  #              "min" = rep(0.7, n),
  #              "mid" = rep(0.8, n),
  #              "max" = rep(0.9, n),
  #              stop("error")
  # )
  # 
  # SPx <- switch(miss_type,
  #               "min" = rep(0.75, n),
  #               "mid" = rep(0.85, n),
  #               "max" = rep(0.95, n),
  #               stop("error")
  # )
  
  # p_Y_star <- SEx*Y + (1-SPx)*(1-Y)
  p_Y_star <- SEx*p_Y + (1-SPx)*(1-p_Y)
  Y_star <- rep(0,n)
  for(j in 1:n){
    Y_star[j] <-rbinom(1, size=1, prob=p_Y_star[j])
  }
  
  
  ## パラメータ推定----
  
  # Naive estimate
  ml <- glm(Y_star~XX, family=binomial)
  
  # Corrected Outcome
  CO <- (Y_star - (1-SPx)) / (SEx + SPx -1)
  
  # 推定方程式
  fn <- function(X, co, beta){
    # t(X) %*% (co - sigmoid(X %*% beta))
    t(X) %*% (sigmoid(X %*% beta)*(1-sigmoid(X %*% beta))*(co - sigmoid(X %*% beta)))
  }
  
  XX_1 <- cbind(rep(1,n), XX)
  
  fn1 <- function(beta) fn(XX_1, CO, beta)
  ans <- nleqslv(c(0.1, 0.1, 0.1), fn1)
  
  ### 結果の保存
  # j <- switch(miss_type,
  #               "min" = 0,
  #               "mid" = kk_T,
  #               "max" = 2*kk_T,
  #               stop("error")
  # )
  
  # results.corrected[i+j,] <- c(ans$x, miss_type)
  # results.naive[i+j,] <- c(ml$coefficients, miss_type)
  
  results.corrected[i,] <- c(ans$x)
  results.naive[i,] <- c(ml$coefficients)

}
# }


##　ボックスプロット----
# results.corr %>% boxplot()
# results.bias1 %>% boxplot()
# results.bias10 %>% boxplot()
# results.bias100 %>% boxplot()
# results.bias1000 %>% boxplot()

## result_dfの作成----
# result_df <- data.frame(
cor_df <- data.frame(results.corrected, method='corrected')
naive_df <- data.frame(results.naive, method='naive')
temp_df<- rbind(cor_df,naive_df) 

# min_df <- temp_df %>% filter(miss_type=='min') %>% select(-miss_type)


theme_gtsummary_mean_sd()
result_df <- temp_df %>% tbl_summary(by=method, digits = everything()~2) %>% 
  modify_header(label="Odds ratio") %>% 
  modify_footnote(label ='真値:(0.2, 0.7, 0.7)') %>% 
  # modify_caption("SE, SP = ( 0.8, 0.85 )の時")
  modify_caption("SE, SP = ( 0.5, 0.55 )の時")

## テーブル出力----
result_df %>% as_gt() %>% gtsave("./results/table/0808/res_table05055.png")


## アウトプット
# export_data <- results.corr

# write_csv2(export_data, file = "~/Projects/gamma_DR_ver0.1/results/ite_sesp/0709.csv")

# 
# df <- read_csv2("~/Projects/gamma_DR_ver0.1/results/0423/g0_m01.csv")
# df <- read_csv2("~/gamma_DR_ver0.1/results/tian/misbin_nval50.csv")
# df %>% summary()
# df %>% boxplot()

