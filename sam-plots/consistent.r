# Clear All
rm(list=ls())

library(plyr)
library(dplyr)

# Grab our data - A dataframe called "D"
load("../data/sam-dat/2.5Mil.Rda")

# Transform it into a table dataframe for use with dplyr
D <- tbl_df(D)

# Need to make at least 4 groups
#	- Consistent
#	- FOSD Only
#	- Light MSB
#	- Light MSB + FOSD

D <- arrange(D,Pattern)

# Recreate the pattern variable in D, this actually is backwards from the way it is in D,
# but it is still the most efficient way that I'm aware of.
patmat <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))
pats <- data.frame(t(patmat))
pat2 <- patmat

# Just flip the thing around
for(i in 10:1){
	pat2[,i] <- patmat[,11-i]
}
patmat <- pat2
rm(pat2)


getCategory <- function(pat){

	cat <- c()

	# FOSD is easy, 0 in row 10 is FOSD
	FOSD <- ifelse(pat[10]==0,T,F)

	# Consistent means 1 switch point with 0 in row 1 or 0 switch points
	switches <- 0

	for(i in 2:9){
		if(pat[i] != pat[i-1]) switches <- switches + 1
	}

	consistent <- (pat[1] == 0 & switches ==1 ) | (pat[1] ==0 & switches ==0 ) | (pat[1] ==1 & switches ==0 )

	LightMSB <- (pat[1]==0 & switches == 3 ) | (pat[1]==0 & switches == 2 )

	if( consistent & !FOSD ){
		cat <- c(0)
	}
	else if( consistent & FOSD ){
		cat <- c(1)
	}
	else if( LightMSB & !FOSD ){
		cat <- c(2)
	}
	else if( LightMSB & FOSD ){
		cat <- c(3)
	}
	else{
		cat <- c(4)
	}

	cat

}

D$Type <- apply(patmat,1,getCategory)

E <- filter(D,Type < 4)

E$Type <- factor(E$Type, labels=c("Consistent","FOSD Only","Light MSB","Light MSB + FOSD"))

cons <- E %>%
	filter(Type=="Consistent") %>%
	select(Pattern,PC,EE,WP,WC,E.0,E.1) %>%
	arrange(desc(PC))

dir <- "../data/sam-plots/"
file <- "Consistent.csv"
fname <- paste(dir,file,sep="")

write.csv(cons, file = fname)
