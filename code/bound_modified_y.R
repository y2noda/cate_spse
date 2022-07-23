library(plotly)
# volcano is a numeric matrix that ships with R


modified_y <- function(se,sp){
  res <- rep(0.5,length(se))
  for (i in 1:length(se)) {
    res[i] <- (0.5-(1-sp[i]))/(se[i]+sp[i]-1)
  }
  return(res)
}


se <- seq(0.501,0.999, by=0.001)
sp <- seq(0.999,0.501, by=-0.001)

z <- modified_y(se=se, sp=sp)

plot_ly(x=se, y=sp, z=modified_y(se,sp))
