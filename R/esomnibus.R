#' Epps-Singleton two-sample test for equality of distribution
#'
#' @aliases esomnibus.default esomnibus.formula
#' @param x A vector or a formula
#' @param ...
#'
#' @return An S3 object of type
#' @export
#'
#' @examples
#' x = runif(50)
#' y = rnorm(50)
#' esomnibus(x,y)
esomnibus <- function(x, ...) {
  UseMethod("esomnibus")
}

#' Calculate Epps-Singleton two-sample test
#'
#' @param y1 Vector of data from first sample
#' @param y2 Vector of data from second sample
#' @param t Vector of points to calculate empirical characteristic function at. Default c(0.4,0.8)
#' @param iqr_type How to calculate semi-inter-quartile-range. Reference the type argument in the IQR function. Defaults to 7. Special value: 99, as referenced in the appendix to the original publication.
#' @param small_correction Should the small sample correction be applied? (TRUE or FALSE)
#'
#' @return A list (which is also of class 'esomnibustest') with the following elements:
#'
#' @export
#'

esomnibus.default <- function(y1 = numeric(), y2, t = c(0.4, 0.8), iqr_type = 7, small_correction = FALSE) {
  if(!missing(t) && (length(t) < 2 || is.na(t)))
    stop("'t' must be a vector of length >1")
  DNAME <- paste(deparse1(substitute(y1)),"vs",
                 deparse1(substitute(y2)))
  y1ok <- !is.na(y1)
  y2ok <- !is.na(y2)
  y1 <- y1[y1ok]
  y2 <- y2[y2ok]
  y <- c(y1, y2)
  n1 <- length(y1)
  n2 <- length(y2)
  n <- n1+n2

  if (n1<5 || n2 <5) stop("Should be at least 5 observations in each group")
  if (iqr_type == 99) {
    # This is the case of the invalid IQR measure taken from the Appendix
    # of the original publication.
    ys <- sort(y)
    l1 <- floor(n/4)
    l2 <- ceiling(n/4)
    u1 <- floor(3*n/4)
    u2 <- ceiling(3*n/4)
    sigma <- 0.25 * (ys[l1] + ys[l2] + ys[u1] + ys[u2])
  } else {
    sigma <- stats::IQR(y, type = iqr_type)/2
  }
  ts =t/sigma
  my1 <- as.matrix(y1)
  my2 <- as.matrix(y2)
  mts <- t(as.matrix(ts))
  g1 <- cbind( cos(my1 %*% mts), sin(my1 %*% mts))
  g2 <- cbind( cos(my2 %*% mts), sin(my2 %*% mts))
  cov1 <- stats::cov(g1)
  cov2 <- stats::cov(g2)
  est_cov <- (n/n1)*cov1 + (n/n2)*cov2
  stopifnot("Covariance matrix has not got full rank!" = qr(est_cov)$rank == nrow(est_cov))
  est_cov_inv <- solve(est_cov)
  r <- nrow(est_cov_inv)
  g_diff <- colMeans(g1) - colMeans(g2)
  w <- n * t(g_diff) %*% est_cov_inv %*% g_diff
  w <- as.numeric(w)

  correcting <- FALSE
  if (small_correction == TRUE) {
    if (max(n1,n2)<25) {
      correcting <- TRUE
      correction_factor <- 1 / (1 + n^(-0.45) + 10.1*(n1^(-1.7) + n2^(-1.7)))
      w <- w * correction_factor
    }
  }

  p <- 1 - stats::pchisq(w,r)
  output <- list("method"="Epps-Singleton two-sample test for equality of distribution",
                      "p.value"=p,
                      "parameter"=w,
                      "t"=t,
                      "n1"=n1,
                      "n2"=n2,
                      "outcome"=DNAME,
                      "Small-sample correction"=correcting,
                      "iqr_type"=iqr_type,
                      "df"=r)
  class(output) <- "esomnibustest"
  output
}

esomnibus.formula <- function(formula, data, subset, na.action, ...) {
  if(missing(formula)
     || (length(formula) != 3L)
     || (length(attr(stats::terms(formula[-2L]), "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if(is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if(nlevels(g) != 2L)
    stop("grouping factor must have exactly 2 levels")
  DATA <- split(mf[[response]], g)
  y <- esomnibus(y1 = DATA[[1L]], y2 = DATA[[2L]], ...)

  y$outcome = DNAME
  y
}

print.esomnibustest <- function(x, digits = getOption("digits"), prefix = "\t", ...) {

  cat("\n")
  cat(strwrap(x$method, prefix = prefix), sep = "\n")
  cat("\n")
  cat("H0: Distributions are equal\n")
  cat("data: ", x$outcome, "\n", sep="")
  cat(paste0("n1: ", x$n1, ", n2: ", x$n2,
             ifelse(x$`Small-sample correction`, ", small-sample correction applied\n","\n")))
  cat(paste0("W(",x$df,"): "), format(x$parameter),"\n")
  cat("p-value: ", format(x$p.value, digits=digits),"\n")
  if (x$iqr_type==99) {
    cat("WARNING: Original, but invalid IQR-method has been applied!")
  }
}
