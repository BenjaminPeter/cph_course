library(tidyverse)
library(data.table)

#estimate het
rm_missing <- function(data){
    d2  <- data %>% group_by(SNP) %>% 
    summarize(no_missing=all(NCHROBS>0)) %>%
    filter(no_missing) %>% 
    select(-no_missing) %>%
    left_join(data)

}

fst_pair <- function(data, pop1, pop2){
    data <- rm_missing(data)
    pop <- c(pop1, pop2)
     n <- data %>% select(SNP, NCHROBS, CLST) %>% filter(CLST %in% pop) %>% spread(key=CLST, value=NCHROBS)  
     mac <- data %>% select(SNP, MAC, CLST) %>% filter(CLST %in% pop) %>% spread(key=CLST, value=MAC)        
    if(length(unique(pop)) == 2){
     data <- left_join(mac, n, by="SNP")
    } else{
        data <- cbind(mac, mac[,2], n[,2], n[,2])
    }
     names(data)[-1] <- c("a1", "a2", "n1", "n2")
     data <- data %>% 
         filter(n1>0, n2>0) %>%
         mutate(h1=a1 * (n1-a1) / n1 / (n1-1)) %>%
         mutate(h2=a2 * (n2-a2) / n2 / (n2-1)) %>%
         mutate(f2= (a1/n1 - a2/n2)^2 - h1 / n1 - h2 / n2) %>%
         mutate(fst = f2 / (f2 + h1 + h2)) %>%
         mutate(fst = ifelse(is.nan(fst), 0, fst))
     return(c(mean(data$f2, na.rm=T), mean(data$fst, na.rm=T), 
              cov(data$a1/data$n1, data$a2/data$n2)))
}


fst_mat <- function(data){
    #d <- hets(data)
    pops <- unique(data$CLST)
    n <- length(pops)
    fst <- matrix(NA, nrow=n, ncol=n)
    rownames(fst) <- pops
    colnames(fst) <- pops
    f2 <- fst
    cov_mat <- fst
    diag(fst) <- 0
    diag(f2) <- 0

    
    for(i in 1:(n-1)){
        for(j in (i+1):n){
            x <- fst_pair(data, pops[i], pops[j])
            print(c(pops[i], pops[j], x))
            f2[i, j] <- x[1]
            f2[j, i] <- x[1]
            fst[i, j] <- x[2]
            fst[j, i] <- x[2]
            cov_mat[i,j] <- x[3]
            cov_mat[j,i] <- x[3]
        }
        x <- fst_pair(data, pops[i], pops[i])
        cov_mat[i,i] <- x[3]
    }
    cov_mat[n, n] <- fst_pair(data, pops[n], pops[n])[3]

    return(list(fst=fst, f2=f2, cov=cov_mat, pops=pops))


}


read_data <- function(fname, ...){
    as_tibble(fread(fname, ...))
}
