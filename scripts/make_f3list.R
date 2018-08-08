args <- commandArgs(T)

target <- args[1]
a <- read.table("data/world/world.ind", as.is=T)
n <- unique(a$V3)
n <- n[!n %in%target]
write.table(cbind(t(combn(n, 2)), target), args[2], 
            row.names=F, col.names=F, quote=F)

