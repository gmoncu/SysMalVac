get.ellipse <- function( fit, n=360 ) 
{
  # Calculate points on an ellipse described by 
  # the fit argument as returned by fit.ellipse 
  # 
  # n is the number of points to render 
  
  tt <- seq(0, 2*pi, length=n) 
  sa <- sin(fit$angle) 
  ca <- cos(fit$angle) 
  ct <- cos(tt) 
  st <- sin(tt) 
  
  x <- fit$center[1] + fit$maj * ct * ca - fit$min * st * sa 
  y <- fit$center[2] + fit$maj * ct * sa + fit$min * st * ca 
  
  cbind(x=x, y=y) 
}