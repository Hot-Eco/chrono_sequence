rmsep <- function(y.obs, y.pred) {
  stopifnot(length(y.obs) == length(y.pred))
  sqrt( sum((y.obs - y.pred)^2) / length(y.obs) )
}

dat.wild <- dat %>%
  filter(fireType == "wildfire") %>%
  select(percentC, tsf.years)

dat.pb <- dat %>%
  filter(fireType == "prescribed") %>%
  select(percentC, tsf.years)

partition <- function(dat, p.train = 0.7) {
  n.test <- round((1 - p.train) * nrow(dat))
  i.test <- 1:nrow(dat)
  
  while (length(i.test) > n.test) {
    i <- sample(length(i.test), 1)
    i.test <- i.test[-i]
  }
  
  list(train = dat[-i.test, ], test = dat[i.test, ])
}

do.cv.span <- function(dat, p.train, spans) {
  d <- partition(dat, p.train)
  perf <- lapply(spans, 
                 function(s) {
                   fit <- loess(percentC ~ tsf.years, 
                                data=d$train, family="symmetric", span=s)
                   p <- predict(fit, newdata = d$test)
                   c(span = s, error = rmsep(d$test$percentC, p))
                 })
  
  do.call(rbind, perf)
}

cv.span <- function(dat, p.train, spans = seq(0.5, 0.9, 0.1), niter = 1000) {
  res <- lapply(1:niter, function(i) do.cv.span(dat, p.train, spans))
  as.data.frame( do.call(rbind, res) )
}

