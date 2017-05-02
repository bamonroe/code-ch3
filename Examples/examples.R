# Functions
crra <- function(x){
	x^(1-r) / (1-r)
}

Pr <- function(x,y){
	exp(x) / ( exp(x) + exp(y) )
}

CE <- function(u){
	((1 - r) * u) ^ (1 / (1-r))
}

# Set up the HL-MPL
A0 <- 1.6
A1 <- 2
B0 <- 0.1
B1 <- 3.85

pA0 <- seq(from = 0.9, to = 0, by = -0.1)
pB0 <- seq(from = 0.9, to = 0, by = -0.1)
pA1 <- seq(from = 0.1, to = 1, by = 0.1)
pB1 <- seq(from = 0.1, to = 1, by = 0.1)

# The example choice pattern is 5 A's followed by 5 B's
choice <- c(0,0,0,0,0,1,1,1,1,1)

# Our example subject's preferences
r <- 0.65
mu <- 0.35

# The utilities of the outcomes
ua0 <- crra(A0)
ua1 <- crra(A1)
ub0 <- crra(B0)
ub1 <- crra(B1)

# The utilities of the options
UA <- ua0*pA0 + ua1*pA1
UB <- ub0*pB0 + ub1*pB1

# The certainty equivalents of the options
CEA <- CE(UA)
CEB <- CE(UB)
CE_chosen <- ifelse(choice == 1, CEB, CEA)
CE_unchosen <- ifelse(choice == 0, CEB, CEA)
CE_diff <- CE_chosen - CE_unchosen
CEM <- ifelse(CEB > CEA, CEB, CEA)
CEMin <- ifelse(CEB > CEA, CEA, CEB)

# Adjust the utilities of the options for Contextual utility
UA[1:9] <- UA[1:9] / (ub1 - ub0) / mu
UA[10]  <- UA[10] / (ub1 - ua1) / mu

UB[1:9] <- UB[1:9] / (ub1 - ub0) / mu
UB[10]  <- UB[10] / (ub1 - ua1) / mu

# The choice probabilities of the options
pA <- Pr(UA,UB)
pB <- Pr(UB,UA)

c_prob <- ifelse(choice == 1, pB, pA)
# The likelihood of the choice pattern 
L <- prod(c_prob)
# Save the results to a csv
write.csv(rbind(c_prob, L), file = "Example_Prob.csv", row.names=FALSE)

# Now the CE stuff
out <- cbind(1:10, CEA, CEB, CE_chosen, CE_unchosen, CE_diff, CEM)
out <- data.frame(out)
colnames(out) <- c("Task","CEA","CEB","CE_chosen","CE_unchosen","CE_diff","CEM")
out

# Save the results to a csv
write.csv(out, file = "Example_CE.csv", row.names=FALSE)


csums <- colSums(out)
csums
csums[4] / csums[7]
# Just wondering here if this would be a useful metric...
sum(CE_diff) / sum(CEM - CEMin)


