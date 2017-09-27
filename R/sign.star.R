sign.star <- function(y)
{
  if(y <= 0.001)
    y <- "***"
  if(y > 0.001 && y <= 0.01)
    y <- "**"
  if(y > 0.01 && y <= 0.05)
    y <- "*"
  if(y > 0.05 && y <= 0.1)
    y <- "."
  if(y > 0.1)
    y <- " "
  return(y)
}
