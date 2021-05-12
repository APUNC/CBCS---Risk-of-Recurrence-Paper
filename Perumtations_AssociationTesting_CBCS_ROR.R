
analyzeOnThis = fread('analyzeOnThis_herit.csv')
analyzeOnThis_White = subset(analyzeOnThis,Race=="White")
analyzeOnThis_Black = subset(analyzeOnThis,Race=="Black")

aaCV = fread('African_cvH2.csv')
wwCV = fread('Euro_CVH2.csv')

eGenes = eGenes[eGenes %in% wwCV$Gene[round(wwCV$CV,3) >= 0.01]]
aGenes = aGenes[aGenes %in% aaCV$Gene[round(aaCV$CV,3) >= 0.01]]


### RORS R2 - WW

gg2 = c("MCM10","FAM64A","CCNB2","MMP1","VAV3","PCSK6")
gg2

thisGene = as.data.frame(analyzeOnThis[,c ('Race',gg2,'RORS','Proliferation','RORP','agesel','er','stage','phaseCode'),with=F])
thisGenesubset = thisGene[which(thisGene$Race=='White'),]
regRORSWW = lm(RORS ~ MCM10 + FAM64A + CCNB2 + MMP1 + VAV3 + PCSK6 + agesel + er + as.factor(stage) + phaseCode ,data = thisGenesubset)
summary(regRORSWW)

#Permute
thisGenesubset$ResidRORS_WW = resid(lm(RORS ~ agesel+as.factor(er)+as.factor(stage)+as.factor(phaseCode),data = thisGenesubset))
regResidRORS_WW = lm(ResidRORS_WW ~ MCM10 + FAM64A + CCNB2 + MMP1 + VAV3 + PCSK6 ,data=thisGenesubset)
summary(regResidRORS_WW)


forPermute = thisGenesubset[,c('ResidRORS_WW',gg2)]

permuteR2Outcome <- function(df,outcome){
  
  outcome = sample(outcome)
  new.df = as.data.frame(cbind(outcome,df))
  colnames(new.df)[1] = 'Outcome'
  return(summary(lm(Outcome ~ .-outcome,data=new.df))$adj.r.squared)
  
}

null.rorsww.permute = replicate(100000,permuteR2Outcome(df=forPermute,outcome=forPermute$ResidRORS_WW))
hist(null.rorsww.permute, main="ROR-S actual R2 = 0.092",
     xlab="Permuted R2",xlim=c(0,0.1))
null.rorsww.permute_df = as.data.frame(null.rorsww.permute)


### RORP R2 - WW

gg2 = c("MCM10","FAM64A","CCNB2","MMP1","VAV3","PCSK6","GNG11")
gg2
thisGene = as.data.frame(analyzeOnThis[,c ('Race',gg2,'RORS','Proliferation','RORP','agesel','er','stage','phaseCode'),with=F])
thisGenesubset = thisGene[which(thisGene$Race=='White'),]
regRORPWW = lm(RORP ~ MCM10 + FAM64A + CCNB2 + MMP1 + VAV3 + PCSK6 + GNG11 + agesel + er + as.factor(stage) + phaseCode ,data = thisGenesubset)
summary(regRORPWW)

#Permute
thisGenesubset$ResidRORP_WW = resid(lm(RORP ~ agesel+er+as.factor(stage)+phaseCode,data = thisGenesubset))
regResidRORP_WW = lm(ResidRORP_WW ~ MCM10 + FAM64A + CCNB2 + MMP1 + VAV3 + GNG11, data = thisGenesubset)
summary(regResidRORP_WW)


forPermute = thisGenesubset[,c('ResidRORP_WW',gg2)]

permuteR2Outcome <- function(df,outcome){
  
  outcome = sample(outcome)
  new.df = as.data.frame(cbind(outcome,df))
  colnames(new.df)[1] = 'Outcome'
  return(summary(lm(Outcome ~ .-outcome,data=new.df))$adj.r.squared)
  
}

null.rorpww.permute = replicate(100000,permuteR2Outcome(df=forPermute,outcome=forPermute$ResidRORP_WW))
hist(null.rorpww.permute, main="ROR-P actual R2 = 0.079",
     xlab="Permuted R2",xlim=c(0,0.1))
warnings()



### Proliferation R2 - WW

gg2 = c("MCM10","FAM64A","CCNB2","MMP1","VAV3")
gg2
thisGene = as.data.frame(analyzeOnThis[,c ('Race',gg2,'RORS','Proliferation','RORP','agesel','er','stage','phaseCode'),with=F])
thisGenesubset = thisGene[which(thisGene$Race=='White'),]
regProlifWW = lm(Proliferation ~ MCM10 + FAM64A + CCNB2 + MMP1 + VAV3 + agesel + er + as.factor(stage) + phaseCode ,data = thisGenesubset)
summary(regProlifWW)

#Permute
thisGenesubset$ResidProlif_WW = resid(lm(Proliferation ~ agesel+er+as.factor(stage)+phaseCode,data = thisGenesubset))
regResidProlif_WW = lm(ResidProlif_WW ~ MCM10 + FAM64A + CCNB2 + MMP1 + VAV3 ,data=thisGenesubset)
summary(regResidProlif_WW)


forPermute = thisGenesubset[,c('ResidProlif_WW',gg2)]

permuteR2Outcome <- function(df,outcome){
  
  outcome = sample(outcome)
  new.df = as.data.frame(cbind(outcome,df))
  colnames(new.df)[1] = 'Outcome'
  return(summary(lm(Outcome ~ .-outcome,data=new.df))$adj.r.squared)
  
}

null.prolifww.permute = replicate(100000,permuteR2Outcome(df=forPermute,outcome=forPermute$ResidProlif_WW))
hist(null.prolifww.permute, main="Proliferation actual R2 = 0.071",
     xlab="Permuted R2",xlim=c(0,0.1))




### RORP R2 - AA

gg2 = 'MMP1'

thisGene = as.data.frame(analyzeOnThis[,c ('Race',gg2,'RORS','Proliferation','RORP','agesel','er','stage','phaseCode'),with=F])
thisGenesubset = thisGene[which(thisGene$Race=='Black'),]
regRORPAA = lm(RORP ~ MMP1 + agesel + er + as.factor(stage) + phaseCode ,data = thisGenesubset)
summary(regRORPAA)

#Permute
thisGenesubset$ResidRORP_AA = resid(lm(RORP ~ agesel+er+as.factor(stage)+phaseCode,data = thisGenesubset))
regResidRORP_AA = lm(ResidRORP_AA ~ MMP1 ,data=thisGenesubset)
summary(regResidRORP_AA)


forPermute = thisGenesubset[,c('ResidRORP_AA',gg2)]

permuteR2Outcome <- function(df,outcome){
  
  outcome = sample(outcome)
  new.df = as.data.frame(cbind(outcome,df))
  colnames(new.df)[1] = 'Outcome'
  return(summary(lm(Outcome ~ .-outcome,data=new.df))$adj.r.squared)
  
}

null.rorpaa.permute = replicate(100000,permuteR2Outcome(df=forPermute,outcome=forPermute$ResidRORP_AA))
hist(null.rorpaa.permute, main="ROR-P actual R2 = 0.009",
     xlab="Permuted R2",xlim=c(0,0.05))

null.rorpaa.permute_df = as.data.frame(null.rorpaa.permute)
fwrite(null.rorpaa.permute_df,"null.rorp.aa.p.csv")


### Prolif R2 - AA

gg2 = 'MMP1'

thisGene = as.data.frame(analyzeOnThis[,c ('Race',gg2,'RORS','Proliferation','RORP','agesel','er','stage','phaseCode'),with=F])
thisGenesubset = thisGene[which(thisGene$Race=='Black'),]
regProlifAA = lm(Proliferation ~ MMP1 + agesel + er + as.factor(stage) + phaseCode ,data = thisGenesubset)
summary(regProlifAA)

#Permute
thisGenesubset$ResidProlif_AA = resid(lm(Proliferation ~ agesel+er+as.factor(stage)+phaseCode,data = thisGenesubset))
regResidProlif_AA = lm(ResidProlif_AA ~ MMP1 ,data=thisGenesubset)
summary(regResidProlif_AA)


forPermute = thisGenesubset[,c('ResidProlif_AA',gg2)]

permuteR2Outcome <- function(df,outcome){
  
  outcome = sample(outcome)
  new.df = as.data.frame(cbind(outcome,df))
  colnames(new.df)[1] = 'Outcome'
  return(summary(lm(Outcome ~ .-outcome,data=new.df))$adj.r.squared)
  
}

null.prolifaa.permute = replicate(100000,permuteR2Outcome(df=forPermute,outcome=forPermute$ResidProlif_AA))
hist(null.prolifaa.permute, main="Proliferation actual R2 = 0.009",
     xlab="Permuted R2",xlim=c(0,0.05))


analyzeOnThis_White = as.data.frame(analyzeOnThis_White)
analyzeOnThis_White_random7 = analyzeOnThis_White[,colnames(analyzeOnThis_White) %in% eGenes]
analyzeOnThis_White_random7 = cbind(analyzeOnThis_White_random7,analyzeOnThis_White[,153:155])


## Random set of 7 GReX of genes (WW)

nullR2 = function(df,genes,outcome,num=7){
  
  ggset = sample(genes,num)
  test = df[,c(outcome,ggset)]
  colnames(test)[1] = 'Outcome'
  reg = lm(Outcome~.,data=test)
  sum = summary(reg)
  return(sum$adj.r.squared)
  
}


ww_nullrandom7.rors.r2 = pbapply::pbreplicate(n = 1e5, expr = nullR2(df = as.data.frame(analyzeOnThis_White_random7),
                                                                  genes = colnames(analyzeOnThis_White_random7)[1:59], outcome = 'RORS'))
ww_nullrandom7.rorp.r2 = pbapply::pbreplicate(n = 1e5, expr = nullR2(df = as.data.frame(analyzeOnThis_White_random7),
                                                                  genes = colnames(analyzeOnThis_White_random7)[1:59], outcome = 'RORP'))
ww_nullrandom7.prolif.r2 = pbapply::pbreplicate(n = 1e5, expr = nullR2(df = as.data.frame(analyzeOnThis_random7),
                                                                    genes = colnames(analyzeOnThis_White_random7)[1:59], outcome = 'Proliferation'))


analyzeOnThis_Black = as.data.frame(analyzeOnThis_Black)
analyzeOnThis_Black_random1 = analyzeOnThis_Black[,colnames(analyzeOnThis_Black) %in% aGenes]
analyzeOnThis_Black_random1 = cbind(analyzeOnThis_Black_random1,analyzeOnThis_Black[,153:155])

## Random set of 1 GReX of genes (BW)

nullR2 = function(df,genes,outcome,num=1){
  
  ggset = sample(genes,num)
  test = df[,c(outcome,ggset)]
  colnames(test)[1] = 'Outcome'
  reg = lm(Outcome~.,data=test)
  sum = summary(reg)
  return(sum$adj.r.squared)
  
}


aa_nullrandom1.rors.r2 = pbapply::pbreplicate(n = 1e5, expr = nullR2(df = as.data.frame(analyzeOnThis_Black_random1),
                                                                  genes = colnames(analyzeOnThis_Black_random1)[1:45], outcome = 'RORS'))
aa_nullrandom1.rorp.r2 = pbapply::pbreplicate(n = 1e5, expr = nullR2(df = as.data.frame(analyzeOnThis_Black_random1),
                                                                  genes = colnames(analyzeOnThis_Black_random1)[1:45], outcome = 'RORP'))
aa_nullrandom1.prolif.r2 = pbapply::pbreplicate(n = 1e5, expr = nullR2(df = as.data.frame(analyzeOnThis_Black_random1),
                                                                    genes = colnames(analyzeOnThis_Black_random1)[1:45], outcome = 'Proliferation'))


k2euroexp = fread("k2euroexp.txt")
k2africanexp = fread("k2africanexp.txt")

k2euroexp_sub = subset(k2euroexp,Gene %in% colnames(analyzeOnThis_White_random7))
k2africanexp_sub = subset(k2africanexp,Gene %in% colnames(analyzeOnThis_Black_random1))

# Transposing

aaa = as.data.frame(k2euroexp_sub)[-1]
aaa = as.matrix(aaa)
rownames(aaa) = k2euroexp_sub$Gene
colnames(aaa) = colnames(k2euroexp_sub)[-1]
k2euroexp_sub = as.data.frame(t(aaa))

aaa = as.data.frame(k2africanexp_sub)[-1]
aaa = as.matrix(aaa)
rownames(aaa) = k2africanexp_sub$Gene
colnames(aaa) = colnames(k2africanexp_sub)[-1]
k2africanexp_sub = as.data.frame(t(aaa))


# Adding covariates

k2euroexp_sub$BCAC_ID = row.names(k2euroexp_sub)
k2africanexp_sub$BCAC_ID = row.names(k2africanexp_sub)


# analyzeOnThis with BCAC_IDs for outcome merging
test = as.data.frame(fread("analyzeOnThis.csv"))
test_euro= subset(test,BCAC_ID %in% k2euroexp_sub$BCAC_ID)
test_euro = test_euro[,c(1,433:435)]
test_african= subset(test,BCAC_ID %in% k2africanexp_sub$BCAC_ID)
test_african = test_african[,c(1,433:435)]

k2euroexp_final = merge(k2euroexp_sub,test_euro,by="BCAC_ID",all=F)
k2africanexp_final = merge(k2africanexp_sub,test_african,by="BCAC_ID",all=F)

k2euroexp_final = as.data.frame(k2euroexp_final)[-1]
k2africanexp_final = as.data.frame(k2africanexp_final)[-1]



## Random set of tumor expression of TWAS test set genes (WW) 

nullR2 = function(df,genes,outcome,num=7){
  
  ggset = sample(genes,num)
  test = df[,c(outcome,ggset)]
  colnames(test)[1] = 'Outcome'
  reg = lm(Outcome~.,data=test)
  sum = summary(reg)
  return(sum$adj.r.squared)
  
}


tumor_ww_nullrandom7.rors.r2 = pbapply::pbreplicate(n = 1e5, expr = nullR2(df = as.data.frame(k2expeuro_final),
                                                                     genes = colnames(k2euroexp_final)[1:59], outcome = 'ROR-S (Subtype Only'))
tumor_ww_nullrandom7.rorp.r2 = pbapply::pbreplicate(n = 1e5, expr = nullR2(df = as.data.frame(k2euroexp_final),
                                                                     genes = colnames(k2euroexp_final)[1:59], outcome = 'ROR-P (Subtype + Proliferation'))
tumor_ww_nullrandom7.prolif.r2 = pbapply::pbreplicate(n = 1e5, expr = nullR2(df = as.data.frame(k2euroexp_final),
                                                                       genes = colnames(k2euroexp_final)[1:59], outcome = 'Proliferation Score'))



## Random set of tumor expression of TWAS test set genes (BW) 

nullR2 = function(df,genes,outcome,num=1){
  
  ggset = sample(genes,num)
  test = df[,c(outcome,ggset)]
  colnames(test)[1] = 'Outcome'
  reg = lm(Outcome~.,data=test)
  sum = summary(reg)
  return(sum$adj.r.squared)
  
}


tumor_aa_nullrandom1.rors.r2 = pbapply::pbreplicate(n = 1e5, expr = nullR2(df = as.data.frame(k2africanexp_final),
                                                                     genes = colnames(k2africanexp_final)[1:45], outcome = 'ROR-S (Subtype Only)'))
tumor_aa_nullrandom1.rorp.r2 = pbapply::pbreplicate(n = 1e5, expr = nullR2(df = as.data.frame(k2africanexp_final),
                                                                     genes = colnames(k2africanexp_final)[1:45], outcome = 'ROR-P (Subtype + Proliferation)'))
tumor_aa_nullrandom1.prolif.r2 = pbapply::pbreplicate(n = 1e5, expr = nullR2(df = as.data.frame(k2africanexp_final),
                                                                       genes = colnames(k2africanexp_final)[1:45], outcome = 'Proliferation Score'))






















