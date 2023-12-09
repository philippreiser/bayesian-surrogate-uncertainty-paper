monotonic_1 <- function(x1, x2, x3, a=0, b=0) {
  x1+x2**3+x3**5
}

radial <- function(x1, x2) {
  (x1)**2 + (x2)**2
}

ishigami <- function(x1, x2, x3, a, b) {
  a * sin(x2)^2 + (1+b*(x3^4)) * sin(x1)
}

ishigami_2 <- function(x1, x2, a=7, b=0.1) {
  a * sin(pi*x2)^2 + sin(pi*x1)
}

sin_fct <- function(x1) {
  sin(0.5*pi*x1)
}

logistic_fct <- function(x1, a=10, b=0) {
  2/(1+exp(-a*(x1-b)))-1
}
