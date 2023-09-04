#args <- commandArgs(trailingOnly = TRUE)

#### GENERAL PARAMETERS #####
require(plyr)
require(hash)

set.seed(1337)
min.boundary <- 5
max.boundary <- 10
min.length <- 15
max.length <- 50
region.min.length <- 70

shape.1 <- 40 #as.numeric(args[1])
shape.2 <- 3 #as.numeric(args[2])
shape.3 <- 1.2 #as.numeric(args[1])
shape.4 <- 0.8 #as.numeric(args[2])
shape.ratio <- 1 #as.numeric(args[3])

d.bg <- 1000
d.fg <- 1000

groupNames = c('0','1','2','3','4')

groupNames.L1A = c('0','1','2')
groupNames.L1B = c('3','4')

groupSizes = hash()
groupIDs = hash()
groupSizes[['0']] = 10
groupSizes[['1']] = 10
groupSizes[['2']] = 10
groupSizes[['3']] = 10
groupSizes[['4']] = 10
for (gn in c(1:length(groupNames))) {
  groupIDs[[groupNames[gn]]] = c()
}

#### simulate chr10 ####
setwd("/Users/zzhu/Documents/research/metilene/metilene2_vscode/metilene_DMR_simulation_scripts/bg_N50/")
files = list.files(pattern = "^sample")
files
i = 1
for (gn in groupNames) {
  for (j in c(1:groupSizes[[gn]])) {
    groupIDs[[gn]][j] = files[i]
    if (i==1) {
      tmp <- read.table(files[i], 
                        col.names=c('chr','pos','strand','context',paste0(paste(gn,files[i],sep='_'),'_m'),paste0(paste(gn,files[i],sep='_'),'_c'),'promoter'))
      data <- cbind(tmp[,-c(5,6)],tmp[,c(5,6)])
      print(paste0(paste(gn,files[i],sep='_'),'_m'))
      rm(tmp)
    } else {
      tmp <- read.table(files[i], 
                        col.names=c('chr','pos','strand','context',paste0(paste(gn,files[i],sep='_'),'_m'),paste0(paste(gn,files[i],sep='_'),'_c'),'promoter'))
      data <- cbind(data,tmp[,c(5,6)])
      rm(tmp)
      print(paste0(paste(gn,files[i],sep='_'),'_m'))
    }
    i = i+1
  }
}

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

data.methyl = hash()
i=1
for (gn in groupNames) {
  data.methyl[[gn]] = c()
  for (j in c(1:groupSizes[[gn]])) {
    data.methyl[[gn]][j] = 4+i*2 # idx of methyl read counts in the data
    i = i+1
  }
}

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
  
  ## use simulated coverages
  cov.groups = hash()
  for (gn in groupNames) {
    cov.groups[[gn]] <- as.vector(apply(data[st:end,c(data.methyl[[gn]])+1], 1, unlist))
  }
  
  mean.groups = hash()
  if (j <= 600){ # major differences
    ## down-regulation
    if (j%%2 == 1){
      for (gn in groupNames.L1A) {
        mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.2, shape.1)*(1-shape.ratio))
      }
      for (gn in groupNames.L1B) {
        mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.1, shape.2)*(1-shape.ratio))
      }
    }
    ## up-regulation
    else {
      for (gn in groupNames.L1A) {
        mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.1, shape.2)*(1-shape.ratio))
      }
      for (gn in groupNames.L1B) {
        mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.2, shape.1)*(1-shape.ratio))
      }    
    }
    ## store true methylation difference
    mean.L1A = mean(mean.groups[[groupNames[1]]])*0
    gn.L1A = ''
    for (gn in groupNames.L1A) {
      mean.L1A = mean.L1A+mean(mean.groups[[gn]])
      gn.L1A = paste(gn.L1A,gn,sep=',')
    }
    mean.L1B = mean(mean.groups[[groupNames[1]]])*0
    gn.L1B = ''
    for (gn in groupNames.L1B) {
      mean.L1B = mean.L1B+mean(mean.groups[[gn]])
      gn.L1B = paste(gn.L1B,gn,sep=',')
    }
    entry$diff <- mean.L1A/length(groupNames.L1A)-mean(mean.L1B)/length(groupNames.L1B)
    entry$comparison <- paste(substr(gn.L1A,2,999),substr(gn.L1B,2,999),sep='|')
  }
  
  if ((j > 600) && (j<=700)){ # one vs others
    gn_cmp = groupNames[((j-1)/2)%%length(groupNames)+1]
    ## down-regulation
    if (j%%2 == 1){
      for (gn in groupNames) {
        mean.groups[[gn_cmp]] <- c(rbeta(length(data.methyl[[gn_cmp]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn_cmp]]), shape.2, shape.1)*(1-shape.ratio))
        if (gn!=gn_cmp){
          mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.1, shape.2)*(1-shape.ratio))
        }
      }
    }
    ## up-regulation
    else {
      for (gn in groupNames) {
        mean.groups[[gn_cmp]] <- c(rbeta(length(data.methyl[[gn_cmp]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn_cmp]]), shape.1, shape.2)*(1-shape.ratio))
        if (gn!=gn_cmp){
          mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.2, shape.1)*(1-shape.ratio))
        }
      }   
    }
    ## store true methylation difference
    mean.L2 = mean(mean.groups[[groupNames[1]]])*0
    for (gn in groupNames) {
      if (gn!=gn_cmp){
        mean.L2 = mean.L2+mean(mean.groups[[gn]])
      }
    }
    entry$diff <- mean(mean.groups[[gn_cmp]])-mean.L2/(length(groupNames)-1)
    entry$comparison <- paste0(gn_cmp,sep='|others')
  }
  
  if ((j > 700) && (j<=850)){ # one vs two
    gn_cmp1 = '0'
    gn_cmp21 = '1'
    gn_cmp22 = '2'
    ## down-regulation
    if (j%%2 == 1){
      for (gn in groupNames) {
        mean.groups[[gn_cmp1]] <- c(rbeta(length(data.methyl[[gn_cmp1]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn_cmp1]]), shape.2, shape.1)*(1-shape.ratio))
        mean.groups[[gn_cmp21]] <- c(rbeta(length(data.methyl[[gn_cmp21]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn_cmp21]]), shape.1, shape.2)*(1-shape.ratio))
        mean.groups[[gn_cmp22]] <- c(rbeta(length(data.methyl[[gn_cmp22]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn_cmp22]]), shape.1, shape.2)*(1-shape.ratio))
        if ((gn!=gn_cmp1)&&(gn!=gn_cmp21)&&(gn!=gn_cmp22)){
          if ((j%%4)/2 >= 1){
            mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.4, shape.3)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.3, shape.4)*(1-shape.ratio))
          } else {
            mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.3, shape.4)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.4, shape.3)*(1-shape.ratio))
          }
        }
      }
    }
    ## up-regulation
    else {
      for (gn in groupNames) {
        mean.groups[[gn_cmp1]] <- c(rbeta(length(data.methyl[[gn_cmp1]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn_cmp1]]), shape.1, shape.2)*(1-shape.ratio))
        mean.groups[[gn_cmp21]] <- c(rbeta(length(data.methyl[[gn_cmp21]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn_cmp21]]), shape.2, shape.1)*(1-shape.ratio))
        mean.groups[[gn_cmp22]] <- c(rbeta(length(data.methyl[[gn_cmp22]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn_cmp22]]), shape.2, shape.1)*(1-shape.ratio))
        if ((gn!=gn_cmp1)&&(gn!=gn_cmp21)&&(gn!=gn_cmp22)){
          if ((j%%4)/2 >= 1){
            mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.4, shape.3)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.3, shape.4)*(1-shape.ratio))
          } else {
            mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.3, shape.4)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.4, shape.3)*(1-shape.ratio))
          }
        }
      }
    }
    ## store true methylation difference
    entry$diff <- mean(mean.groups[[gn_cmp1]])-(mean(mean.groups[[gn_cmp21]])+mean(mean.groups[[gn_cmp22]]))/2
    entry$comparison <- '0|1,2'
  }
  
  if ((j > 850) && (j<=1000)){ # one vs one
    gn_cmp1 = '3'
    gn_cmp2 = '4'
    ## down-regulation
    if (j%%2 == 1){
      for (gn in groupNames) {
        mean.groups[[gn_cmp1]] <- c(rbeta(length(data.methyl[[gn_cmp1]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn_cmp1]]), shape.2, shape.1)*(1-shape.ratio))
        mean.groups[[gn_cmp2]] <- c(rbeta(length(data.methyl[[gn_cmp2]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn_cmp2]]), shape.1, shape.2)*(1-shape.ratio))
        if ((gn!=gn_cmp1)&&(gn!=gn_cmp2)){
          if ((j%%4)/2 >= 1){
            mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.4, shape.3)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.3, shape.4)*(1-shape.ratio))
          } else {
            mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.3, shape.4)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.4, shape.3)*(1-shape.ratio))
          }
        }
      }
    }
    ## up-regulation
    else {
      for (gn in groupNames) {
        mean.groups[[gn_cmp1]] <- c(rbeta(length(data.methyl[[gn_cmp1]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn_cmp1]]), shape.1, shape.2)*(1-shape.ratio))
        mean.groups[[gn_cmp2]] <- c(rbeta(length(data.methyl[[gn_cmp2]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn_cmp2]]), shape.2, shape.1)*(1-shape.ratio))
        if ((gn!=gn_cmp1)&&(gn!=gn_cmp2)){
          if ((j%%4)/2 >= 1){
            mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.4, shape.3)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.3, shape.4)*(1-shape.ratio))
          } else {
            mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.3, shape.4)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.4, shape.3)*(1-shape.ratio))
          }
        }
      }
    }
    ## store true methylation difference
    entry$diff <- mean(mean.groups[[gn_cmp1]])-mean(mean.groups[[gn_cmp2]])
    entry$comparison <- '3|4'
  }
  
  
  ## draw methylated read counts
  for (gn in groupNames){
    data[st:end, data.methyl[[gn]]] <<- matrix(rbinom(length(cov.groups[[gn]]), cov.groups[[gn]], mean.groups[[gn]]), nrow=len, byrow=T)
  }
  print(j)
  ## alternating variable for
  ## direction of change
  j <<- j+1
  
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
  
  ## use simulated coverages
  cov.groups = hash()
  for (gn in groupNames) {
    cov.groups[[gn]] <- as.vector(apply(data[st:end,c(data.methyl[[gn]])+1], 1, unlist))
  }
  
  mean.groups = hash()
  if (j <= 600){ # major differences
    ## down-regulation
    if (j%%2 == 1){
      for (gn in groupNames.L1A) {
        mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.2, shape.1)*(1-shape.ratio))
      }
      for (gn in groupNames.L1B) {
        mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.1, shape.2)*(1-shape.ratio))
      }
    }
    ## up-regulation
    else {
      for (gn in groupNames.L1A) {
        mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.1, shape.2)*(1-shape.ratio))
      }
      for (gn in groupNames.L1B) {
        mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.2, shape.1)*(1-shape.ratio))
      }    
    }
    ## store true methylation difference
    mean.L1A = mean(mean.groups[[groupNames[1]]])*0
    gn.L1A = ''
    for (gn in groupNames.L1A) {
      mean.L1A = mean.L1A+mean(mean.groups[[gn]])
      gn.L1A = paste(gn.L1A,gn,sep=',')
    }
    mean.L1B = mean(mean.groups[[groupNames[1]]])*0
    gn.L1B = ''
    for (gn in groupNames.L1B) {
      mean.L1B = mean.L1B+mean(mean.groups[[gn]])
      gn.L1B = paste(gn.L1B,gn,sep=',')
    }
    entry$diff <- mean.L1A/length(groupNames.L1A)-mean(mean.L1B)/length(groupNames.L1B)
    entry$comparison <- paste(substr(gn.L1A,2,999),substr(gn.L1B,2,999),sep='|')
  }
  
  if ((j > 600) && (j<=700)){ # one vs others
    gn_cmp = groupNames[((j-1)/2)%%length(groupNames)+1]
    ## down-regulation
    if (j%%2 == 1){
      for (gn in groupNames) {
        mean.groups[[gn_cmp]] <- c(rbeta(length(data.methyl[[gn_cmp]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn_cmp]]), shape.2, shape.1)*(1-shape.ratio))
        if (gn!=gn_cmp){
          mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.1, shape.2)*(1-shape.ratio))
        }
      }
    }
    ## up-regulation
    else {
      for (gn in groupNames) {
        mean.groups[[gn_cmp]] <- c(rbeta(length(data.methyl[[gn_cmp]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn_cmp]]), shape.1, shape.2)*(1-shape.ratio))
        if (gn!=gn_cmp){
          mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.2, shape.1)*(1-shape.ratio))
        }
      }   
    }
    ## store true methylation difference
    mean.L2 = mean(mean.groups[[groupNames[1]]])*0
    for (gn in groupNames) {
      if (gn!=gn_cmp){
        mean.L2 = mean.L2+mean(mean.groups[[gn]])
      }
    }
    entry$diff <- mean(mean.groups[[gn_cmp]])-mean.L2/(length(groupNames)-1)
    entry$comparison <- paste0(gn_cmp,sep='|others')
  }
  
  if ((j > 700) && (j<=850)){ # one vs two
    gn_cmp1 = '0'
    gn_cmp21 = '1'
    gn_cmp22 = '2'
    ## down-regulation
    if (j%%2 == 1){
      for (gn in groupNames) {
        mean.groups[[gn_cmp1]] <- c(rbeta(length(data.methyl[[gn_cmp1]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn_cmp1]]), shape.2, shape.1)*(1-shape.ratio))
        mean.groups[[gn_cmp21]] <- c(rbeta(length(data.methyl[[gn_cmp21]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn_cmp21]]), shape.1, shape.2)*(1-shape.ratio))
        mean.groups[[gn_cmp22]] <- c(rbeta(length(data.methyl[[gn_cmp22]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn_cmp22]]), shape.1, shape.2)*(1-shape.ratio))
        if ((gn!=gn_cmp1)&&(gn!=gn_cmp21)&&(gn!=gn_cmp22)){
          if ((j%%4)/2 >= 1){
            mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.4, shape.3)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.3, shape.4)*(1-shape.ratio))
          } else {
            mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.3, shape.4)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.4, shape.3)*(1-shape.ratio))
          }
        }
      }
    }
    ## up-regulation
    else {
      for (gn in groupNames) {
        mean.groups[[gn_cmp1]] <- c(rbeta(length(data.methyl[[gn_cmp1]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn_cmp1]]), shape.1, shape.2)*(1-shape.ratio))
        mean.groups[[gn_cmp21]] <- c(rbeta(length(data.methyl[[gn_cmp21]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn_cmp21]]), shape.2, shape.1)*(1-shape.ratio))
        mean.groups[[gn_cmp22]] <- c(rbeta(length(data.methyl[[gn_cmp22]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn_cmp22]]), shape.2, shape.1)*(1-shape.ratio))
        if ((gn!=gn_cmp1)&&(gn!=gn_cmp21)&&(gn!=gn_cmp22)){
          if ((j%%4)/2 >= 1){
            mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.4, shape.3)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.3, shape.4)*(1-shape.ratio))
          } else {
            mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.3, shape.4)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.4, shape.3)*(1-shape.ratio))
          }
        }
      }
    }
    ## store true methylation difference
    entry$diff <- mean(mean.groups[[gn_cmp1]])-(mean(mean.groups[[gn_cmp21]])+mean(mean.groups[[gn_cmp22]]))/2
    entry$comparison <- '0|1,2'
  }
  
  if ((j > 850) && (j<=1000)){ # one vs one
    gn_cmp1 = '3'
    gn_cmp2 = '4'
    ## down-regulation
    if (j%%2 == 1){
      for (gn in groupNames) {
        mean.groups[[gn_cmp1]] <- c(rbeta(length(data.methyl[[gn_cmp1]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn_cmp1]]), shape.2, shape.1)*(1-shape.ratio))
        mean.groups[[gn_cmp2]] <- c(rbeta(length(data.methyl[[gn_cmp2]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn_cmp2]]), shape.1, shape.2)*(1-shape.ratio))
        if ((gn!=gn_cmp1)&&(gn!=gn_cmp2)){
          if ((j%%4)/2 >= 1){
            mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.4, shape.3)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.3, shape.4)*(1-shape.ratio))
          } else {
            mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.3, shape.4)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.4, shape.3)*(1-shape.ratio))
          }
        }
      }
    }
    ## up-regulation
    else {
      for (gn in groupNames) {
        mean.groups[[gn_cmp1]] <- c(rbeta(length(data.methyl[[gn_cmp1]]), shape.2, shape.1)*shape.ratio+rbeta(length(data.methyl[[gn_cmp1]]), shape.1, shape.2)*(1-shape.ratio))
        mean.groups[[gn_cmp2]] <- c(rbeta(length(data.methyl[[gn_cmp2]]), shape.1, shape.2)*shape.ratio+rbeta(length(data.methyl[[gn_cmp2]]), shape.2, shape.1)*(1-shape.ratio))
        if ((gn!=gn_cmp1)&&(gn!=gn_cmp2)){
          if ((j%%4)/2 >= 1){
            mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.4, shape.3)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.3, shape.4)*(1-shape.ratio))
          } else {
            mean.groups[[gn]] <- c(rbeta(length(data.methyl[[gn]]), shape.3, shape.4)*shape.ratio+rbeta(length(data.methyl[[gn]]), shape.4, shape.3)*(1-shape.ratio))
          }
        }
      }
    }
    ## store true methylation difference
    entry$diff <- mean(mean.groups[[gn_cmp1]])-mean(mean.groups[[gn_cmp2]])
    entry$comparison <- '3|4'
  }
  
  
  ## draw methylated read counts
  for (gn in groupNames){
    data[st:end, data.methyl[[gn]]] <<- matrix(rbinom(length(cov.groups[[gn]]), cov.groups[[gn]], mean.groups[[gn]]), nrow=len, byrow=T)
  }
  print(j)
  ## alternating variable for
  ## direction of change
  j <<- j+1
  
  ## remove current region
  ## and add remainder to
  ## the end of region.fg
  region.fg$len[idx] <<- 0
  region.fg <<- rbind(region.fg, data.frame(group=region.fg$group[idx], st=region.fg$st[idx], end=st-1, len=st-region.fg$st[idx], promoter=region.fg$promoter[idx]))
  region.fg <<- rbind(region.fg, data.frame(group=region.fg$group[idx], st=end+1, end=region.fg$end[idx], len=region.fg$end[idx]-end, promoter=region.fg$promoter[idx]))
  
  ## return
  entry
}))


names <- c()
for (gn in groupNames) {
  for (j in c(1:groupSizes[[gn]])) {
    names <- c(names,paste(gn,groupIDs[[gn]][j],sep='_'))
  }
}
names
    
# # BSmooth
# j <- 1
# for (gn in groupNames) {
#   for (i in data.methyl[[gn]]){
#     sample <- data[,c(1:4,i,(i+1))]
#     sample <- na.omit(sample)
#     write.table(file=paste0(names[j], "_beta.", shape.1, ".", shape.2, "_ratio.", shape.ratio, ".bsmooth"), format(sample, scientific=F, trim=T), col.names=F, row.names=F, sep="\t", quote=F)
#     j <- j+1
#   }
# }

# metilene
j <- (ncol(data) + 1)
k <- (ncol(data) + 1)

for (gn in groupNames) {
  for (i in data.methyl[[gn]]){
    data[,j] <- data[,i]/data[,(i+1)]
    j <- j+1
  }
}
data <- data[,c(1,2,k:(j-1))]
colnames(data)[3:ncol(data)] <- names
data$chr <- paste("chr", data$chr, sep="")
write.table(file=paste0("beta.", shape.1, ".", shape.2, "_ratio.", shape.ratio, ".metilene"), format(data, scientific=F, trim=T), col.names=T, row.names=F, sep="\t", quote=F)

# annotation DMRs
DMR.bg$promoter <- 0
DMR.fg$promoter <- 1
DMR.bg$chr <- paste("chr", DMR.bg$chr, sep="")
DMR.fg$chr <- paste("chr", DMR.fg$chr, sep="")
write.table(file=paste0("DMR_annotation_beta.", shape.1, ".", shape.2, "_ratio.", shape.ratio, ".bed"), format(DMR.bg, scientific=F, trim=T), col.names=F, row.names=F, sep="\t", quote=F)
write.table(file=paste0("DMR_annotation_beta.", shape.1, ".", shape.2, "_ratio.", shape.ratio, ".bed"), format(DMR.fg, scientific=F, trim=T), col.names=F, row.names=F, sep="\t", quote=F, append=T)




# Apply beta function         
y_beta <- rbeta(1000, shape1 = shape.1, shape2 = shape.2) 
y_beta_r <- rbeta(1000, shape1 = shape.2, shape2 = shape.1) 
# Plot beta values
p1 <- hist(y_beta)                     # centered at 4
p2 <- hist(y_beta_r)                     # centered at 6
plot( p2, col=rgb(0,0,1,1/4), xlim=c(0,1), ylim=c(0,250))  # first histogram
plot( p1, col=rgb(1,0,0,1/4), xlim=c(0,1), ylim=c(0,250), add=T)  # second
# Apply beta function         

# Apply beta function         
y_beta <- rbeta(1000, shape1 = shape.3, shape2 = shape.4) 
y_beta_r <- rbeta(1000, shape1 = shape.4, shape2 = shape.3) 
# Plot beta values
p1 <- hist(y_beta)                     # centered at 4
p2 <- hist(y_beta_r)                     # centered at 6
plot( p2, col=rgb(0,0,1,1/4), xlim=c(0,1), ylim=c(0,250))  # first histogram
plot( p1, col=rgb(1,0,0,1/4), xlim=c(0,1), ylim=c(0,250), add=T)  # second
# Apply beta function                                        
