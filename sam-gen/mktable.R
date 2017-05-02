# Make the table that is used for EUT in Chapter 3

# Save the final Dataset
samdat <- "../data/sam-dat/"
filename <- paste0(samdat,"2.5Mil-EUT.Rda")
load(file=filename)

PP <- D$Pattern
PP <- strsplit(PP, split = ",")
PP <- lapply(PP, as.numeric)
PP <- do.call(rbind, PP)
PP <- PP[1:10,]

PP

o <- 1:10

TT <- cbind(o, PP, D$PC[o], D$EE[o], D$WP[o], D$WC[o], D$E.0[o], D$E.1[o])
TT <- data.frame(TT)
colnames(TT) <- c("RANK","1","2","3","4","5","6","7","8","9","10","PC","EE","WP","WC","E.0","E.1")
TT

write.csv(TT, file=paste0(samdat,"TopTenEUT.csv"), row.names=F, quote=F)
write.csv(TT, file=paste0("~/Dropbox/Thesis/Chapter 3/tables/","TopTenEUT.csv"), row.names=F, quote=F)
