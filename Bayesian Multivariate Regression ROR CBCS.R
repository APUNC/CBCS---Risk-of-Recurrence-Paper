
#Black (follow same pattern for White women and for subtype centroid correlations as outcome)
grex = as.data.frame(fread('analyzeOnThis.csv'))
grex = subset(grex,Race=='Black')
exp = as.data.frame(fread('CBCS_exp_2189samples_166genes.txt'))
grex = subset(grex,BCAC_ID %in% colnames(exp))
exp = exp[,c('Gene',grex$BCAC_ID)]


grex_genes <- c("MMP1")
grex.mat = as.matrix(as.data.frame(grex)[,grex_genes])
rownames(grex.mat) = grex$BCAC_ID
#colnames(grex.mat)= "MMP1"
# Not needed, but double check - colnames(grex.mat) = colnames(grex)[10:15]
for (i in 1:ncol(grex.mat)){
  grex.mat[,i] = (grex.mat[,i] - mean(grex.mat[,i]))/sd(grex.mat[,i])
}

covs.mat <- grex[,c("BCAC_ID","agesel","er","stage","phaseCode")]

covs.mat$stage2 <- ifelse(covs.mat$stage==2,1,0)
covs.mat$stage3 <- ifelse(covs.mat$stage==3,1,0)

covs.mat$phaseCode_new[covs.mat$phaseCode==3] <- 1
covs.mat$phaseCode_new[covs.mat$phaseCode==12] <- 0

covs.mat$er_new[covs.mat$er==2] <- 1
covs.mat$er_new[covs.mat$er==1] <- 0

covs.mat$phaseCode_new <- as.numeric(covs.mat$phaseCode_new)
covs.mat$er_new <- as.numeric(covs.mat$er_new)

table(covs.mat$er_new)
table(covs.mat$er)

agesel_v <- covs.mat$agesel
er_new_v <- covs.mat$er_new
phaseCode_new_v <- covs.mat$phaseCode_new
stage2_v <- covs.mat$stage2
stage3_v <- covs.mat$stage3

covs.mat2 <- cbind(agesel_v, er_new_v, phaseCode_new_v,stage2_v,stage3_v)

outcome_genes <- c("NAT1","MAPT","KIF2C","GRB7","KNTC2","CCNB1","UBE2C","EXO1","SLC39A6",
                   "PTTG1","RRM2","BIRC5","ACTR3B","TMEM45B","FOXC1","MELK","MYBL2",
                   "ERBB2","UBE2T","CDC6","TYMS","CDC20","ORC6L","FOXA1","BAG1",
                   "CEP55","PGR","ESR1","MMP11","MKI67","CXXC5",
                   "ANLN","BLVRA","GPR160","CENPF","CCNE1","CDCA1")
outcome <- subset(exp, Gene %in% outcome_genes)
outcome.mat = t(as.matrix(outcome[,-1]))
colnames(outcome.mat) <- outcome_genes
require(missForest)
covs.mat2 = missForest(covs.mat2)$ximp

# Check the X= code line in relation to above edits
X = cbind(grex.mat,covs.mat2)
Y = outcome.mat

XY = as.data.frame(cbind(X,Y))

runConjugate <- function(Y,X,reps=1e5,p.cred=.975){
  require(MCMCpack)
  require(matrixNormal)
  
  n = nrow(Y)
  k = ncol(X)
  m = ncol(Y)
  
  V.0 = diag(1,m)
  nu.0 = m
  B.0 = matrix(0,nrow=k,
               ncol=m)
  Lambda.0 = diag(1,k)
  nu.n = nu.0 + n
  B.n = ginv(t(X) %*% X + Lambda.0) %*% 
    (t(X) %*% Y + Lambda.0 %*% B.0)
  Lambda.n = t(X) %*% X + Lambda.0
  V.n = V.0 + 
    t(Y-X%*%B.n)%*%(Y-X%*%B.n) + 
    t(B.n - B.0)%*%Lambda.0%*%(B.n-B.0)
  
  getPosteriors <- function(nu.n,V.n,B.n,Lambda.n){
    Sigma.post = riwish(nu.n,V.n)
    B.post = rmatnorm(s=1,B.n,ginv(Lambda.n),Sigma.post)
    return(B.post)
  }
  
  total = pbapply::pbreplicate(reps,
                               getPosteriors(nu.n = nu.n,V.n = V.n,B.n = B.n,
                                             Lambda.n = Lambda.n))
  
  posteriorMeans = apply(total,1:2,mean)
  posteriorSD = apply(total,1:2,sd)
  posteriorLower = apply(total,1:2,quantile,probs=0.025)
  posteriorUpper = apply(total,1:2,quantile,probs=0.975)
  
  return(list(Distributions = total,
              Means = posteriorMeans,
              SD = posteriorSD,
              Lower = posteriorLower,
              Upper = posteriorUpper))
}

TestAA <- runConjugate(Y,X,1e5,0.975)

TestAA[["Means"]]
TestAA[["SD"]]
TestAA[["Lower"]]
TestAA[["Upper"]]