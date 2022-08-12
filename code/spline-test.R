set.seed(200)
n <- 500
x <- seq(-2*pi, 2*pi, length.out=n)
y <- sin(x) + rnorm(n, 0, 0.5)

# fit_gam_cr <- gam(y~s(x, bs="cr",k=5), knots = list(knots))
fit_gam_cr <- gam(y~s(x, bs="cr"), knots = list(knots))
a <- predict(fit_gam_cr)


knots <- c(0.3, 0.5, 0.6)
nsMat <- naturalSpline(x, knots = knots, intercept = FALSE)
x_s <-predict(nsMat,x)
gam_model <- lm(y~nsMat)

b <- gam_model %>% predict(x_s)

par(mfrow = c(1, 1))
plot(x,y)
par(new=T)

plot(x,a)
par(new=T)
plot(x,b)
            
            
par(mfrow=c(3,3))
for(i in 1:9){
  if(i==10){
    plot(fit_gam_cc, xlab="", ylab="", main="Estimated and True function")
    lines(x, sin(x), col="red")
    legend("bottomleft", legend=c("Estimated", "True"), col=c("black", "red"), lty=1)
  }else{
    plot(model.matrix(fit_gam_cc) %*% diag(fit_gam_cc$coefficients)[,i], type="l",
         xlab="", ylab="", main="basis function")
  }
}

library(splines2)

knots <- c(0.3, 0.5, 0.6)
x <- seq(0, 1, 0.01)
bsMat <- bSpline(x, knots = knots, degree = 1, intercept = TRUE)
matplot(x, bsMat, type = "l", ylab = "y")
abline(v = knots, lty = 2, col = "gray")

nsMat <- naturalSpline(x, knots = knots, intercept = TRUE)
insMat <- naturalSpline(x, knots = knots, intercept = TRUE, integral = TRUE)
par(mfrow = c(1, 2))
matplot(x, nsMat, type = "l", ylab = "Basis")
matplot(x, insMat, type = "l", ylab = "Integrals")
