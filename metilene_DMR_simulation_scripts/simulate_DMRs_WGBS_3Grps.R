args <- commandArgs(trailingOnly = TRUE)

#### GENERAL PARAMETERS #####
require(plyr)
set.seed(1337)
min.boundary <- 5
max.boundary <- 10
min.length <- 15
max.length <- 50
region.min.length <- 70

shape.1 <- as.numeric(args[1])
shape.2 <- as.numeric(args[2])
shape.ratio <- as.numeric(args[3])

d.bg <- 500
d.fg <- 500

#### simulate chr10 ####
BL01 <- read.table(paste0(args[4], "sample_1.txt"), col.names=c('chr','pos','strand','context','BL01_m','BL01_c','promoter'))
BL02 <- read.table(paste0(args[4], "sample_2.txt"), col.names=c('chr','pos','strand','context','BL02_m','BL02_c','promoter'))
BL03 <- read.table(paste0(args[4], "sample_3.txt"), col.names=c('chr','pos','strand','context','BL03_m','BL03_c','promoter'))
BL04 <- read.table(paste0(args[4], "sample_4.txt"), col.names=c('chr','pos','strand','context','BL04_m','BL04_c','promoter'))
BL05 <- read.table(paste0(args[4], "sample_5.txt"), col.names=c('chr','pos','strand','context','BL05_m','BL05_c','promoter'))
BL06 <- read.table(paste0(args[4], "sample_6.txt"), col.names=c('chr','pos','strand','context','BL06_m','BL06_c','promoter'))
BL07 <- read.table(paste0(args[4], "sample_7.txt"), col.names=c('chr','pos','strand','context','BL07_m','BL07_c','promoter'))
BL08 <- read.table(paste0(args[4], "sample_8.txt"), col.names=c('chr','pos','strand','context','BL08_m','BL08_c','promoter'))
BL09 <- read.table(paste0(args[4], "sample_9.txt"), col.names=c('chr','pos','strand','context','BL09_m','BL09_c','promoter'))
BL10 <- read.table(paste0(args[4], "sample_10.txt"), col.names=c('chr','pos','strand','context','BL10_m','BL10_c','promoter'))

FL01 <- read.table(paste0(args[4], "sample_11.txt"), col.names=c('chr','pos','strand','context','FL01_m','FL01_c','promoter'))
FL02 <- read.table(paste0(args[4], "sample_12.txt"), col.names=c('chr','pos','strand','context','FL02_m','FL02_c','promoter'))
FL03 <- read.table(paste0(args[4], "sample_13.txt"), col.names=c('chr','pos','strand','context','FL03_m','FL03_c','promoter'))
FL04 <- read.table(paste0(args[4], "sample_14.txt"), col.names=c('chr','pos','strand','context','FL04_m','FL04_c','promoter'))
FL05 <- read.table(paste0(args[4], "sample_15.txt"), col.names=c('chr','pos','strand','context','FL05_m','FL05_c','promoter'))
FL06 <- read.table(paste0(args[4], "sample_16.txt"), col.names=c('chr','pos','strand','context','FL06_m','FL06_c','promoter'))
FL07 <- read.table(paste0(args[4], "sample_17.txt"), col.names=c('chr','pos','strand','context','FL07_m','FL07_c','promoter'))
FL08 <- read.table(paste0(args[4], "sample_18.txt"), col.names=c('chr','pos','strand','context','FL08_m','FL08_c','promoter'))
FL09 <- read.table(paste0(args[4], "sample_19.txt"), col.names=c('chr','pos','strand','context','FL09_m','FL09_c','promoter'))
FL10 <- read.table(paste0(args[4], "sample_20.txt"), col.names=c('chr','pos','strand','context','FL10_m','FL10_c','promoter'))

TL01 <- read.table(paste0(args[4], "sample_21.txt"), col.names=c('chr','pos','strand','context','TL01_m','TL01_c','promoter'))
TL02 <- read.table(paste0(args[4], "sample_22.txt"), col.names=c('chr','pos','strand','context','TL02_m','TL02_c','promoter'))
TL03 <- read.table(paste0(args[4], "sample_23.txt"), col.names=c('chr','pos','strand','context','TL03_m','TL03_c','promoter'))
TL04 <- read.table(paste0(args[4], "sample_24.txt"), col.names=c('chr','pos','strand','context','TL04_m','TL04_c','promoter'))
TL05 <- read.table(paste0(args[4], "sample_25.txt"), col.names=c('chr','pos','strand','context','TL05_m','TL05_c','promoter'))
TL06 <- read.table(paste0(args[4], "sample_26.txt"), col.names=c('chr','pos','strand','context','TL06_m','TL06_c','promoter'))
TL07 <- read.table(paste0(args[4], "sample_27.txt"), col.names=c('chr','pos','strand','context','TL07_m','TL07_c','promoter'))
TL08 <- read.table(paste0(args[4], "sample_28.txt"), col.names=c('chr','pos','strand','context','TL08_m','TL08_c','promoter'))
TL09 <- read.table(paste0(args[4], "sample_29.txt"), col.names=c('chr','pos','strand','context','TL09_m','TL09_c','promoter'))
TL10 <- read.table(paste0(args[4], "sample_30.txt"), col.names=c('chr','pos','strand','context','TL10_m','TL10_c','promoter'))


data <- cbind(BL01[,-c(5,6)],BL01[,c(5,6)],BL02[,c(5,6)],BL03[,c(5,6)],BL04[,c(5,6)],BL05[,c(5,6)],
                             BL06[,c(5,6)],BL07[,c(5,6)],BL08[,c(5,6)],BL09[,c(5,6)],BL10[,c(5,6)],
                             FL01[,c(5,6)],FL02[,c(5,6)],FL03[,c(5,6)],FL04[,c(5,6)],FL05[,c(5,6)],
                             FL06[,c(5,6)],FL07[,c(5,6)],FL08[,c(5,6)],FL09[,c(5,6)],FL10[,c(5,6)],
                             TL01[,c(5,6)],TL02[,c(5,6)],TL03[,c(5,6)],TL04[,c(5,6)],TL05[,c(5,6)],
                             TL06[,c(5,6)],TL07[,c(5,6)],TL08[,c(5,6)],TL09[,c(5,6)],TL10[,c(5,6)])

## additional data
data$idx <- 1:nrow(data)
data$dist <- c(0,diff(data$pos))
data$breaks <- cumsum(c(0,abs(diff(data$promoter))))
data$group <- cumsum(as.integer(data$dist>=300))+1
data$region <- data$breaks + data$group
region.bg <- ddply(data, .(region), summarise, st=min(idx), end=max(idx), len=length(pos), promoter=0)
region.bg <- subset(region.bg, promoter==0)
region.fg <- ddply(data, .(group), function(i){
  breaks <- which(diff(i$promoter) != 0)
  if (length(breaks) == 0){
    NULL
  }
  else {
    pos <- unique(sort(unlist(lapply(breaks, function(j){max(1,j-max.length):min(nrow(i),j+max.length)}))))
    pos.breaks <- which(diff(c(-1,pos)) > 1)
    ival.st <- pos[pos.breaks]
    ival.end <- c(pos[pos.breaks[-1]-1],pos[length(pos)])
    data.frame(st=i$idx[ival.st], end=i$idx[ival.end], len=ival.end-ival.st, promoter=1)
  }
})

data.methyl.BL <- c(6,8,10,12,14,16,18,20,22,24)
data.methyl.FL <- c(26,28,30,32,34,36,38,40,42,44)
data.methyl.TL <- c(46,48,50,52,54,56,58,60,62,64)

## DMR background data simulation
## (i.e., in non-promotor regions)
sample.bg.len <- sample(min.length:max.length, d.bg, replace=TRUE)
	
j <- 1
DMR.bg <- do.call(rbind.fill, lapply(1:d.bg, function(i){
  ## init data line
  entry <- data.frame(chr=0, st=0, end=0, diff=0, st.idx=0, end.idx=0, comparison=0)  
  
  ## sample start and calculate end
  len <- sample.bg.len[i]
  idx <- sample(1:nrow(region.bg), size=1, prob=pmax.int(0, region.bg$len-len))
  st <- sample(region.bg$st[idx]:(region.bg$end[idx]-len+1), size=1)
  end <- st+len-1
	
  ## store DMR boundaries
  entry$chr <- data$chr[st]
  entry$st <- data$pos[st]-1
  entry$end <- data$pos[end]
  entry$st.idx <- st
  entry$end.idx <- end

  ## alternating variable for
  ## direction of change
  j <<- j+1
  
  ## use simulated coverages
  cov.BL <- as.vector(apply(data[st:end,c(data.methyl.BL)+1], 1, unlist))
  cov.FL <- as.vector(apply(data[st:end,c(data.methyl.FL)+1], 1, unlist))
  cov.TL <- as.vector(apply(data[st:end,c(data.methyl.TL)+1], 1, unlist))

  ## down-regulation
  if (j%%2 == 1){
    if (j%%5 == 0) {
      mean.BL <- c(rbeta(length(data.methyl.BL), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl.BL), shape.2, shape.1)*(1-shape.ratio))
      mean.FL <- c(rbeta(length(data.methyl.FL), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl.FL), shape.1, shape.2)*(1-shape.ratio))
      mean.TL <- c(rbeta(length(data.methyl.TL), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl.TL), shape.2, shape.1)*(1-shape.ratio))
      ## store true methylation difference
      entry$diff <- mean(mean.FL)-((mean(mean.BL)+mean(mean.TL))/2)
      entry$comparison <- "FL-(BL,TL)"
    }
    else {
      mean.BL <- c(rbeta(length(data.methyl.BL), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl.BL), shape.2, shape.1)*(1-shape.ratio))
      mean.FL <- c(rbeta(length(data.methyl.FL), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl.FL), shape.1, shape.2)*(1-shape.ratio))
      mean.TL <- c(rbeta(length(data.methyl.TL), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl.TL), shape.1, shape.2)*(1-shape.ratio))
      ## store true methylation difference
      entry$diff <- mean(mean.BL)-((mean(mean.FL)+mean(mean.TL))/2)
      entry$comparison <- "BL-(FL,TL)"
    }

  }
  ## up-regulation
  else {
    if (j%%10 == 0) {
      mean.BL <- c(rbeta(length(data.methyl.BL), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl.BL), shape.1, shape.2)*(1-shape.ratio))
      mean.FL <- c(rbeta(length(data.methyl.FL), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl.FL), shape.1, shape.2)*(1-shape.ratio))
      mean.TL <- c(rbeta(length(data.methyl.TL), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl.TL), shape.2, shape.1)*(1-shape.ratio)) 
      ## store true methylation difference
      entry$diff <- mean(mean.TL)-((mean(mean.BL)+mean(mean.FL))/2)
      entry$comparison <- "TL-(BL,FL)"
    } 
    else {
      mean.BL <- c(rbeta(length(data.methyl.BL), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl.BL), shape.1, shape.2)*(1-shape.ratio))
      mean.FL <- c(rbeta(length(data.methyl.FL), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl.FL), shape.2, shape.1)*(1-shape.ratio))
      mean.TL <- c(rbeta(length(data.methyl.TL), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl.TL), shape.2, shape.1)*(1-shape.ratio)) 
      ## store true methylation difference
      entry$diff <- mean(mean.BL)-((mean(mean.FL)+mean(mean.TL))/2)
      entry$comparison <- "BL-(FL,TL)"
    }
       
  }

  ## draw methylated read counts
  data[st:end, data.methyl.BL] <<- matrix(rbinom(length(cov.BL), cov.BL, mean.BL), nrow=len, byrow=T)
  data[st:end, data.methyl.FL] <<- matrix(rbinom(length(cov.FL), cov.FL, mean.FL), nrow=len, byrow=T)
  data[st:end, data.methyl.TL] <<- matrix(rbinom(length(cov.TL), cov.TL, mean.TL), nrow=len, byrow=T)

  # ## store true methylation difference
  # entry$diff <- mean(mean.BL)-((mean(mean.FL)+mean(mean.TL))/2)
  # entry$comparison <- "BL-(FL,TL)"

  ## remove current region
  ## and add remainder to
  ## the end of region.bg
  region.bg$len[idx] <<- 0
  region.bg <<- rbind(region.bg, data.frame(region=region.bg$region[idx], st=region.bg$st[idx], end=st-1, len=st-region.bg$st[idx], promoter=region.bg$promoter[idx]))
  region.bg <<- rbind(region.bg, data.frame(region=region.bg$region[idx], st=end+1, end=region.bg$end[idx], len=region.bg$end[idx]-end, promoter=region.bg$promoter[idx]))

  ## return
  entry
}))

## DMR foreground data simulation
## (i.e., in promotor regions)
sample.fg.len <- sample(min.length:max.length, d.fg, replace=TRUE)
	
j <- 1
DMR.fg <- do.call(rbind.fill, lapply(1:d.fg, function(i){
  ## init data line
  entry <- data.frame(chr=0, st=0, end=0, diff=0, st.idx=0, end.idx=0)  
  
  ## sample start and calculate end
  len <- sample.fg.len[i]
  idx <- sample(1:nrow(region.fg), size=1, prob=pmax.int(0, region.fg$len-len))
  st <- sample(region.fg$st[idx]:(region.fg$end[idx]-len+1), size=1)
  end <- st+len-1
	
  ## store DMR boundaries
  entry$chr <- data$chr[st]
  entry$st <- data$pos[st]-1
  entry$end <- data$pos[end]
  entry$st.idx <- st
  entry$end.idx <- end

  ## alternating variable for
  ## direction of change
  j <<- j+1
  
  ## use simulated coverages
  cov.BL <- as.vector(apply(data[st:end,c(data.methyl.BL)+1], 1, unlist))
  cov.FL <- as.vector(apply(data[st:end,c(data.methyl.FL)+1], 1, unlist))
  cov.TL <- as.vector(apply(data[st:end,c(data.methyl.TL)+1], 1, unlist))

  ## down-regulation
  if (j%%2 == 1){
    if (j%%5 == 0) {
      mean.BL <- c(rbeta(length(data.methyl.BL), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl.BL), shape.2, shape.1)*(1-shape.ratio))
      mean.FL <- c(rbeta(length(data.methyl.FL), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl.FL), shape.1, shape.2)*(1-shape.ratio))
      mean.TL <- c(rbeta(length(data.methyl.TL), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl.TL), shape.2, shape.1)*(1-shape.ratio))
      ## store true methylation difference
      entry$diff <- mean(mean.FL)-((mean(mean.BL)+mean(mean.TL))/2)
      entry$comparison <- "FL-(BL,TL)"
    }
    else {
      mean.BL <- c(rbeta(length(data.methyl.BL), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl.BL), shape.2, shape.1)*(1-shape.ratio))
      mean.FL <- c(rbeta(length(data.methyl.FL), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl.FL), shape.1, shape.2)*(1-shape.ratio))
      mean.TL <- c(rbeta(length(data.methyl.TL), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl.TL), shape.1, shape.2)*(1-shape.ratio))    
      ## store true methylation difference
      entry$diff <- mean(mean.BL)-((mean(mean.FL)+mean(mean.TL))/2)
      entry$comparison <- "BL-(FL,TL)"  
    }

  }
  ## up-regulation
  else {
    if (j%%10 == 0) {
      mean.BL <- c(rbeta(length(data.methyl.BL), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl.BL), shape.1, shape.2)*(1-shape.ratio))
      mean.FL <- c(rbeta(length(data.methyl.FL), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl.FL), shape.1, shape.2)*(1-shape.ratio)) 
      mean.TL <- c(rbeta(length(data.methyl.TL), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl.TL), shape.2, shape.1)*(1-shape.ratio))
      ## store true methylation difference
      entry$diff <- mean(mean.TL)-((mean(mean.BL)+mean(mean.FL))/2)
      entry$comparison <- "TL-(BL,FL)"
    }
    else {
      mean.BL <- c(rbeta(length(data.methyl.BL), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl.BL), shape.1, shape.2)*(1-shape.ratio))
      mean.FL <- c(rbeta(length(data.methyl.FL), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl.FL), shape.2, shape.1)*(1-shape.ratio)) 
      mean.TL <- c(rbeta(length(data.methyl.TL), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl.TL), shape.2, shape.1)*(1-shape.ratio))  
      ## store true methylation difference
      entry$diff <- mean(mean.BL)-((mean(mean.FL)+mean(mean.TL))/2)
      entry$comparison <- "BL-(FL,TL)"      
    }
    
  }

  ## draw methylated read counts
  data[st:end, data.methyl.BL] <<- matrix(rbinom(length(cov.BL), cov.BL, mean.BL), nrow=len, byrow=T)
  data[st:end, data.methyl.FL] <<- matrix(rbinom(length(cov.FL), cov.FL, mean.FL), nrow=len, byrow=T)
  data[st:end, data.methyl.TL] <<- matrix(rbinom(length(cov.TL), cov.TL, mean.TL), nrow=len, byrow=T)

  # ## store true methylation difference
  # entry$diff <- mean(mean.BL)-mean(mean.FL)

  ## remove current region
  ## and add remainder to
  ## the end of region.fg
  region.fg$len[idx] <<- 0
  region.fg <<- rbind(region.fg, data.frame(group=region.fg$group[idx], st=region.fg$st[idx], end=st-1, len=st-region.fg$st[idx], promoter=region.fg$promoter[idx]))
  region.fg <<- rbind(region.fg, data.frame(group=region.fg$group[idx], st=end+1, end=region.fg$end[idx], len=region.fg$end[idx]-end, promoter=region.fg$promoter[idx]))
  
  ## return
  entry
}))


names <- c("0_BL01","0_BL02","0_BL03","0_BL04","0_BL05","0_BL06","0_BL07","0_BL08","0_BL09","0_BL10",
           "1_FL01","1_FL02","1_FL03","1_FL04","1_FL05","1_FL06","1_FL07","1_FL08","1_FL09","1_FL10",
           "2_TL01","2_TL02","2_TL03","2_TL04","2_TL05","2_TL06","2_TL07","2_TL08","2_TL09","2_TL10")

# BSmooth
j <- 1
for (i in c(data.methyl.BL, data.methyl.FL, data.methyl.TL)){
	sample <- data[,c(1:4,i,(i+1))]
	sample <- na.omit(sample)
	write.table(file=paste0(as.character(args[5]), names[j], "_beta.", shape.1, ".", shape.2, "_ratio.", shape.ratio, ".bsmooth"), format(sample, scientific=F, trim=T), col.names=F, row.names=F, sep="\t", quote=F)
	j <- j+1
}

# metilene
j <- (ncol(data) + 1)
k <- (ncol(data) + 1)
for (i in c(data.methyl.BL, data.methyl.FL, data.methyl.TL)){
	data[,j] <- data[,i]/data[,(i+1)]
	j <- j+1
}
data <- data[,c(1,2,k:(j-1))]
colnames(data)[3:ncol(data)] <- names
data$chr <- paste("chr", data$chr, sep="")
write.table(file=paste0(as.character(args[5]), "beta.", shape.1, ".", shape.2, "_ratio.", shape.ratio, ".metilene"), format(data, scientific=F, trim=T), col.names=T, row.names=F, sep="\t", quote=F)

# annotation DMRs
DMR.bg$promoter <- 0
DMR.fg$promoter <- 1
DMR.bg$chr <- paste("chr", DMR.bg$chr, sep="")
DMR.fg$chr <- paste("chr", DMR.fg$chr, sep="")
write.table(file=paste0(as.character(args[5]), "DMR_annotation_beta.", shape.1, ".", shape.2, "_ratio.", shape.ratio, ".bed"), format(DMR.bg, scientific=F, trim=T), col.names=F, row.names=F, sep="\t", quote=F)
write.table(file=paste0(as.character(args[5]), "DMR_annotation_beta.", shape.1, ".", shape.2, "_ratio.", shape.ratio, ".bed"), format(DMR.fg, scientific=F, trim=T), col.names=F, row.names=F, sep="\t", quote=F, append=T)







