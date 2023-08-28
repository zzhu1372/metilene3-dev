## check for input characters
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3){
  stop("3 input parameters required: beta.shape1 beta.shape2 chromatin_annotation_chr10.txt")
}

## global options
coverage.mean <- 30
coverage.sd <- 5
coverage.min <- 15
beta.shape1 <- as.numeric(args[1])
beta.shape2 <- as.numeric(args[2])
smoothing <- TRUE
smoothing.span <- 0.2 ## default: 0.75 but is too smooth for my taste
smoothing.range.min <- 25
smoothing.range.max <- 50
smoothing.weight <- 0.5

## read data
data <- read.table(as.character(args[3]), stringsAsFactors=FALSE, col.names=c("chr", "pos", "strand", "dinucleotide", "chromatin.state"))

## init columns
data$meth <- 0
data$cov <- 0
data$meth.rates <- -1
data$smooth <- -1

## simulate coverages
data$cov <- round(rnorm(nrow(data), coverage.mean, coverage.sd))
data$cov[data$cov < coverage.min] <- coverage.min

## categorize into unmethylated and methylated positions
## depending on the chromatin state (promotor or not)
data$promoter <- as.integer(grepl("Promoter", data$chromatin.state))

## simulate rates at methylated positions (non-promoter)
idx <- which(data$promoter == 0)
data$meth.rates[idx] <- rbeta(length(idx), beta.shape1, beta.shape2)
rm(idx)

## simulate rates at unmethylated positions (promoter)
idx <- which(data$promoter == 1)
data$meth.rates[idx] <- rbeta(length(idx), beta.shape2, beta.shape1)
rm(idx)

## optional: smoothing
if (smoothing){
  ## calculate loess smoothing
  data$smooth <- loess(data$meth.rates ~ data$pos, span=smoothing.span)$fitted

  ## get boundaries between promoter/non-promoter
  ## and identify CpGs in range of a sampled
  ## smoothing range between min and max for
  ## smoothing by the weighted mean between
  ## raw and fitted methylation rate
  ##
  ## identify CpG position before change
  breaks.idx <- which(diff(data$promoter) != 0)
  #breaks.margin <- 50
  breaks.margin <- runif(length(breaks.idx), smoothing.range.min, smoothing.range.max)
  
  range <- breaks.idx
  ## search for CpG positions after change within margin
  i <- 1  
  breaks.tmp <- abs(diff(data$pos, lag=i))[breaks.idx]
  while(any(breaks.tmp <= breaks.margin)){    
    range <- c(range, breaks.idx[breaks.tmp <= breaks.margin] + i)
    i <- i + 1
    breaks.tmp <- abs(diff(data$pos, lag=i))[breaks.idx]
    stopifnot(i <= 1000)
  }
  ## search for CpG positions before change within margin
  i <- 1
  breaks.tmp <- abs(diff(rev(data$pos), lag=i))[nrow(data)-breaks.idx+1]
  while(any(breaks.tmp <= breaks.margin)){    
    range <- c(range, breaks.idx[breaks.tmp <= breaks.margin] - i)
    i <- i + 1
    breaks.tmp <- abs(diff(rev(data$pos), lag=i))[nrow(data)-breaks.idx+1]
    stopifnot(i <= 1000)
  }
  ## calculate weighted mean of raw and smoothing value
  data$meth.rates[range] <- data$meth.rates[range] * 0.5 + data$smooth[range] * 0.5
}

## simulate number of methylated reads
data$meth <- rbinom(nrow(data), data$cov, data$meth.rates)

## remove not required columns
data$meth.rates <- NULL
data$chromatin.state <- NULL
data$smooth <- NULL

## pretty formatting
data$pos <- format(data$pos, scientific=F, trim=T)
data$meth <- format(data$meth, scientific=F, trim=T)
data$cov <- format(data$cov, scientific=F, trim=T)

## output file
write.table(file=stdout(), data, quote=F, sep="\t", row.names=F, col.names=F)

## cleanup
rm(data)
