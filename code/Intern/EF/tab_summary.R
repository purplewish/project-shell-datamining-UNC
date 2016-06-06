variables1 <- c("S2","Tmax","Romeasured")
tab_kernel <- read.csv("docs/tabs/univariate/tab_kernel.csv")
tab_uni <- read.csv("docs/tabs/univariate/tab_uni.csv")
tab_uni_all <- cbind(tab_uni,tab_kernel[,2])
colnames(tab_uni_all) <- c("variable","old","new","kernel")
write.csv(tab_uni_all,"tabs/univariate/tab_uni_all.csv")

name_comb <- unlist(lapply(1:(length(variables1)-1),function(x){paste(variables1[x],variables1[(x+1):length(variables1)],sep="_")}))

tab_all <- data.frame(method = c(rep(c("LMC_vg","Matern","LMC_mle"),length(name_comb)),"univariate"),S2 = 0, Tmax = 0, Romeasured = 0)
tab_lmc_mle <- read.csv("docs/tabs/bivariate/tab_lmc_mle_S2_Tmax_Ro.csv")

for(j in 1:3)
{
  path <- paste("tabs/bivariate/","tab_bivariate_",name_comb[j],".csv",sep="")
  tab1 <- read.csv(path,stringsAsFactors = FALSE)
  temp <- tab_lmc_mle[rowSums(tab_lmc_mle[,tab1[,1]] ==0)==0,tab1[,1]]
  tab_all[(1+(j-1)*3):(3*j),tab1[,1]] <- rbind(as.matrix(t(tab1[,-1])),as.numeric(temp))
}

value_uni <- apply(tab_uni_all[tab_uni_all$variable %in% variables1,-1],1,min)

tab_all[10,-1] <- value_uni

write.csv(tab_all,"docs/tabs/tab_all_S2_Tmax_Ro.csv",row.names = FALSE)


