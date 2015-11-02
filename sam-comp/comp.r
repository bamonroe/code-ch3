library(np)

#load("Sam1500000-sim1000--0.21.3.Rda")
load("../data/aggdat/Test.Rda")
# Multiple Linear Regression  on the row 6-7 in depth


plot(MM$rm,MM$M.WP)
plot(MM$rm,MM$V.WP)

plot(MM$rm,MM$M.EE)
plot(MM$rm,MM$V.EE)

plot(MM$rm,MM$M.WE0)
plot(MM$rm,MM$V.WE0)


stop("here")

M1 <- MM[which(MM$rm > .14 & MM$rm < .45),]
#M1 <- MM[which(MM$rm > .3 & MM$rm < .55),]
M1 <- MM


rm2 <- M1$rm^2
rs2 <- M1$rs^2
um2 <- M1$um^2
us2 <- M1$us^2

fit <- lm(M.WP ~ rm + rs + rm2 + rs2 + um + us +um2 +us2, data=M1)
sf <- summary(fit) # show results


rm2 <- MM$rm^2
rs2 <- MM$rs^2
um2 <- MM$um^2
us2 <- MM$us^2


fit1 <- lm(M.WP ~ rm + rs + rm2 + rs2 + um + us +um2 +us2, data=MM)
sf1 <- summary(fit1) # show results

# Some interesting plots ?
ll <- list(
	rm =seq(from=(min(M1$rm)-.2),to=(max(M1$rm)+.2),by=.001),
	rs =seq(from=0,to=3,by=.01),
	um =seq(from=0,to=3,by=.01),
	us =seq(from=0,to=3,by=.01)
)

fc <- fit$coefficients

line <- list(
	rm=fc["rm"]*ll$rm + fc["rm2"]*ll$rm^2 + fc[1],
	rs=fc["rs"]*ll$rs + fc["rs2"]*ll$rs^2 + fc[1],

	um=fc["um"]*ll$um + fc["um2"]*ll$um^2 + fc[1],
	us=fc["us"]*ll$us + fc["us2"]*ll$us^2 + fc[1]

)

plot(ll$rm,line$rm)

print(sf)
print(sf1)

M2 <- MM[which(MM$rm<.41 & MM$rm >.14),]
M2 <- MM
U <- runif(nrow(M2))
#M2 <- M2[which(U<.001),]

#M2$M.WP <- M2$M.WP *2


model.par <- lm(V.EE ~ rm + I(rm^2), data = M2)
summary(model.par)

model.np <- npreg(V.EE ~ rm,
                 regtype = "ll",
                 bwmethod = "cv.aic",
                 gradients = TRUE,
                 data = M2)

summary(model.np)


plot(M2$rm, M2$V.EE, xlab = "r mean", ylab = "Weighted Var WP", cex=.1)
lines(M2$rm, fitted(model.np), lty = 1, col = "blue")
lines(M2$rm, fitted(model.par), lty = 2, col = " red")

stop("here")

plot(model.np, plot.errors.method = "asymptotic")


plot(model.np, gradients = TRUE)
lines(M2$rm, coef(model.par)[2]+2*M2$rm*coef(model.par)[3],
 	lty = 2,
 	col = "red")

plot(model.np, gradients = TRUE, plot.errors.method = "asymptotic")
