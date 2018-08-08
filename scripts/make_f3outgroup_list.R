args <- commandArgs(T)

target <- args[1]
known <- args[2]
fname <- args[3]
a <- read.table("data/europe/europe.ind", as.is=T)
n <- unique(a$V3)
n <- n[!n %in%c(target, known)]
write.table(cbind(known, t(combn(n, 1)), target),  fname,
            row.names=F, col.names=F, quote=F)

