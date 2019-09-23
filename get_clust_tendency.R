#' get_clust_tendency from factoextra modified
#'
#' Made possible to include NAs, use naomit=F
#'
#' @param data
#' @param n
#' @param gradient
#' @param seed
#'
#' @return
#' @export
#'
#'
get_clust_tendency<-function (data, n,  gradient = list(low = "red",
                                                 mid = "white", high = "blue"), seed = 123, naomit=T)
{
  set.seed(seed)
  if (is.data.frame(data))
    data <- as.matrix(data)
  if (!(is.matrix(data)))
    stop("data must be data.frame or matrix")
  if (n >= nrow(data))
    stop("n must be no larger than num of samples")
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("reshape2 package needed for this function to work. Please install it.")
  }
  if(naomit==T){
  data <- na.omit(data)
  }
  rownames(data) <- paste0("r", 1:nrow(data))

  p <- apply(data, 2, function(x, n) {
    suppressWarnings(runif(n, min(x, na.rm=T), max(x, na.rm=T)))
  }, n)

k <- round(suppressWarnings(runif(n, 1, nrow(data))))

  q <- as.matrix(data[k, ])
  distp = rep(0, nrow(data))
  distq = 0
  minp = rep(0, n)
  minq = rep(0, n)
  for (i in 1:n) {
    distp[1] <- stats::dist(rbind(p[i, ], data[1, ]))
    minqi <- stats::dist(rbind(q[i, ], data[1, ]))
    for (j in 2:nrow(data)) {
      distp[j] <- stats::dist(rbind(p[i, ], data[j, ]))
      error <- q[i, ] - data[j, ]
      if (sum(abs(error), na.rm = T) != 0) {
        distq <- stats::dist(rbind(q[i, ], data[j, ]))
        if (distq < minqi)
          minqi <- distq
      }
    }
    minp[i] <- min(distp)
    minq[i] <- minqi
  }
  list(hopkins_stat = sum(minq)/(sum(minp) + sum(minq)))
}