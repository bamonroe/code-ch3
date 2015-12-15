# Clear All
rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)

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

E$ll <- log(E$PC)

# There are three points to make, A, B and C. A is the most likely consistent point
# B is the most likely FOSD only point, C is the consistent point that is closest in
# liklihood to B that has greater likelihood

cons <- E %>%
	select(Type,WP,ll) %>%
	filter(Type=="Consistent") %>%
	arrange(desc(ll))

FOSD <- E %>%
	select(Type,WP,ll) %>%
	filter(Type=="FOSD Only") %>%
	arrange(desc(ll))

# A and B are straight forward
A <- c( cons$ll[1] , cons$WP[1] )
B <- c( FOSD$ll[1] , FOSD$WP[1] )

# Need to find the closes point

diff <- cons %>%
	mutate(diff =(ll - B[1])) %>%
	filter(diff>0) %>%
	arrange(diff)

C <- c( diff$ll[1] , diff$WP[1] )

# Use arrange so that Consistent is plotted last
E <- arrange(E,desc(Type),ll)


# Now to configure the graph aesthetically

# Configure the points
## Shape of the points
shape <- 20
## Size of the points
point.size <- 5
## alpha
alph <- .8

# Configure x-axis title
## title size
x.title.size <- 14
## title vertical justification
x.title.vjust <- 1

# Configure x-axis text
## text size
x.text.size <- 11
## angle of the text
x.text.angle <- 45

# Configure y-axis title
## title size
y.title.size <- 14
## title horizontal justification
y.title.hjust <- .5
## title vertical justification
y.title.vjust <- 1

# Configure y-axis text
## text size
y.text.size <- 11
## angle of the text
y.text.angle <- 45

# Configure the annotate for A, B, and C stuff
## A needs to be pushed down and to the right
A.label <- c(A[1]+.1,A[2]-.007)
B.label <- c(B[1]+.1,B[2]-.007)
C.label <- c(C[1]+.1,C[2]-.007)

p <- ggplot(data=E, aes_string(x="ll",y="WP"))
p <- p + annotate("rect", xmin=B[1],xmax=Inf,ymin=-Inf,ymax=B[2],alpha=.2)
p <- p + geom_point(shape=shape,size=point.size, alpha=alph, aes_string(color="Type") )
p <- p + labs(x="Log of Simulated Likelihood",y="Expected Ratio of Obtained to Optimal Welfare")
p <- p + theme( axis.title.x=element_text(size=x.title.size),
				axis.text.x=element_text(size=x.text.size,angle=x.text.angle),
				axis.title.y=element_text(size=y.title.size,hjust=y.title.hjust,vjust=y.title.vjust),
				axis.text.y=element_text(size=y.text.size,angle=y.text.angle)
				)
p <- p + annotate("text", x=A.label[1],y=A.label[2], label="A")
p <- p + annotate("text", x=B.label[1],y=B.label[2], label="B")
p <- p + annotate("text", x=C.label[1],y=C.label[2], label="C")

p

#stop("here")

# Configurations for saving the file
## Width of the plot
w <- 6.9	# Width of paragraph area
## Height of the plot
h <- 5		# About half a page
## Scale of the plot
save.scale <- 1.5

dir <- "../data/sam-plots/"
file <- "Figure1.jpg"
fname <- paste(dir,file,sep="")

ggsave(filename=fname, plot=p, width=w, height=h, units="in",scale=save.scale )

