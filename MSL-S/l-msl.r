
source("/home/woodape/git/thesis/sim/sim.r", echo = TRUE, print.eval = TRUE, chdir = TRUE)

E.0 <- data.frame(p0=D$pA0 , p1=D$pA1 , pz0=D$A0 , pz1=D$A1 , min=D$min , max=D$max,  r=D$r , mu=D$mu , U=D$UA , p=D$pA )
E.1 <- data.frame(p0=D$pB0 , p1=D$pB1 , pz0=D$B0 , pz1=D$B1 , min=D$min , max=D$max,  r=D$r , mu=D$mu , U=D$UB , p=D$pB )

E.0$c <- ifelse(D$c == 0 , 1 , 0)
E.1$c <- ifelse(D$c == 1 , 1 , 0)

E.0$id=D$ID
E.1$id=D$ID

E.0$gid <- 1:nrow(E.0)
E.1$gid <- 1:nrow(E.1)

E.0$alt <- 0
E.1$alt <- 1

E <- rbind(E.0,E.1)

E <- E[order(E$id,E$gid,E$alt),]

write.dta(E,"choice.dta")
