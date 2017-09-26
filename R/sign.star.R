sign.star <- function(y)
{
  if(y <= 0.001)
    y <- noquote(paste0(y, " ***"))
  if(y > 0.001 && y <= 0.01)
    y <- noquote(paste0(y, " **"))
  if(y > 0.01 && y <= 0.05)
    y <- noquote(paste0(y, " *"))
  if(y > 0.05 && y <= 0.1)
    y <- noquote(paste0(y, " ."))
  return(y)
}
