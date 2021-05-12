
### Identifying significant GReX against ROR measures by race;
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("bacon")

install.packages("data.table")
install.packages("xlsx")
install.packages("curl")
require(bacon)
require(data.table)
require(curl)

analyzeOnThis = fread('analyzeOnThis_herit.csv')
ll = readRDS('prioritizedGenes.RDS')
eGenes = ll$eGenes
aGenes = ll$aGenes

aaCV = fread('African_cvH2.csv')
wwCV = fread('Euro_CVH2.csv')

eGenes = eGenes[eGenes %in% wwCV$Gene[round(wwCV$CV,3) >= 0.01]]
aGenes = aGenes[aGenes %in% aaCV$Gene[round(aaCV$CV,3) >= 0.01]]

eGenes

intGenes = unique(c(aGenes,eGenes))
resRORS_AA = resProlif_AA = resRORP_AA = as.data.frame(matrix(ncol = 5,nrow = length(aGenes)))
resRORS_WW = resProlif_WW = resRORP_WW = as.data.frame(matrix(ncol = 5,nrow = length(eGenes)))
colnames(resRORS_AA) = colnames(resProlif_AA) = colnames(resRORP_AA) = c('Gene','Beta','Z','SE','P')
colnames(resRORS_WW) = colnames(resProlif_WW) = colnames(resRORP_WW) = c('Gene','Beta','Z','SE','P')
resRORS_AA$Gene = resProlif_AA$Gene = resRORP_AA$Gene = aGenes
resRORS_WW$Gene = resProlif_WW$Gene = resRORP_WW$Gene = eGenes
analyzeOnThis = analyzeOnThis[,c('Race',intGenes,'RORS','Proliferation','RORP','agesel','er','stage','phaseCode'),with=F]

for (i in 1:nrow(resRORS_AA)){
  
  Gene = aGenes[i]
  print(Gene)
  thisGene = analyzeOnThis[,c('Race',Gene,'RORS','Proliferation','RORP','agesel','er','stage','phaseCode'),with=F]
  colnames(thisGene)[2] = 'Gene'
  thisGene$Gene = (thisGene$Gene - mean(thisGene$Gene))/sd(thisGene$Gene)
  
  reg = lm(RORS ~ Gene + agesel + as.factor(er) + as.factor(stage) + as.factor(phaseCode),data = subset(thisGene,Race == 'Black'))
  sss = coef(summary(reg))
  resRORS_AA$Beta[i] = sss[2,1]
  resRORS_AA$Z[i] = sss[2,1]/sss[2,2]
  resRORS_AA$SE[i] = sss[2,2]
  resRORS_AA$P[i] = sss[2,4]
  
  reg = lm(RORP ~ Gene + agesel + as.factor(er) + as.factor(stage) + as.factor(phaseCode),data = subset(thisGene,Race == 'Black'))
  sss = coef(summary(reg))
  resRORP_AA$Beta[i] = sss[2,1]
  resRORP_AA$Z[i] = sss[2,1]/sss[2,2]
  resRORP_AA$SE[i] = sss[2,2]
  resRORP_AA$P[i] = sss[2,4]
  
  reg = lm(Proliferation ~ Gene + agesel + as.factor(er) + as.factor(stage) + as.factor(phaseCode),data = subset(thisGene,Race == 'Black'))
  sss = coef(summary(reg))
  resProlif_AA$Beta[i] = sss[2,1]
  resProlif_AA$Z[i] = sss[2,1]/sss[2,2]
  resProlif_AA$SE[i] = sss[2,2]
  resProlif_AA$P[i] = sss[2,4]
  
  
}

for (i in 1:nrow(resRORS_WW)){
  
  Gene = eGenes[i]
  print(Gene)
  thisGene = analyzeOnThis[,c('Race',Gene,'RORS','Proliferation','RORP','agesel','er','stage','phaseCode'),with=F]
  colnames(thisGene)[2] = 'Gene'
  thisGene$Gene = (thisGene$Gene - mean(thisGene$Gene))/sd(thisGene$Gene)
  
  reg = lm(RORS ~ Gene + agesel + as.factor(er) + as.factor(stage) + as.factor(phaseCode),data = subset(thisGene,Race == 'White'))
  sss = coef(summary(reg))
  resRORS_WW$Beta[i] = sss[2,1]
  resRORS_WW$Z[i] = sss[2,1]/sss[2,2]
  resRORS_WW$SE[i] = sss[2,2]
  resRORS_WW$P[i] = sss[2,4]
  
  reg = lm(RORP ~ Gene + agesel + as.factor(er) + as.factor(stage) + as.factor(phaseCode),data = subset(thisGene,Race == 'White'))
  sss = coef(summary(reg))
  resRORP_WW$Beta[i] = sss[2,1]
  resRORP_WW$Z[i] = sss[2,1]/sss[2,2]
  resRORP_WW$SE[i] = sss[2,2]
  resRORP_WW$P[i] = sss[2,4]
  
  reg = lm(Proliferation ~ Gene + agesel + as.factor(er) + as.factor(stage) + as.factor(phaseCode),data = subset(thisGene,Race == 'White'))
  sss = coef(summary(reg))
  resProlif_WW$Beta[i] = sss[2,1]
  resProlif_WW$Z[i] = sss[2,1]/sss[2,2]
  resProlif_WW$SE[i] = sss[2,2]
  resProlif_WW$P[i] = sss[2,4]
  
  
}

#### FOR AA
require(bacon)
bc <- bacon(teststatistics = resRORS_AA$Z,verbose=T)
resRORS_AA$Bacon <- pval(bc)
resRORS_AA$FDR <- p.adjust(resRORS_AA$Bacon,method='BH')
resRORS_AA = resRORS_AA[order(resRORS_AA$FDR),]

bc <- bacon(teststatistics = resRORP_AA$Z,verbose=T)
resRORP_AA$Bacon <- pval(bc)
resRORP_AA$FDR <- p.adjust(resRORP_AA$Bacon,method='BH')
resRORP_AA = resRORP_AA[order(resRORP_AA$FDR),]

bc <- bacon(teststatistics = resProlif_AA$Z,verbose=T)
resProlif_AA$Bacon <- pval(bc)
resProlif_AA$FDR <- p.adjust(resProlif_AA$Bacon,method='BH')
resProlif_AA = resProlif_AA[order(resProlif_AA$FDR),]

require(xlsx)
write.xlsx(resRORS_AA,'AA_TWAS_multigenescores.xlsx',sheetName = 'RORS',row.names=F,col.names=T)
write.xlsx(resRORP_AA,'AA_TWAS_multigenescores.xlsx',sheetName = 'RORP',append=T,row.names=F,col.names=T)
write.xlsx(resProlif_AA,'AA_TWAS_multigenescores.xlsx',sheetName = 'Proliferation Score',append=T,row.names=F,col.names=T)

#### FOR WW
require(bacon)
bc <- bacon(teststatistics = resRORS_WW$Z,verbose=T)
resRORS_WW$Bacon <- pval(bc)
resRORS_WW$FDR <- p.adjust(resRORS_WW$Bacon,method='BH')
resRORS_WW = resRORS_WW[order(resRORS_WW$FDR),]

bc <- bacon(teststatistics = resRORP_WW$Z,verbose=T)
resRORP_WW$Bacon <- pval(bc)
resRORP_WW$FDR <- p.adjust(resRORP_WW$Bacon,method='BH')
resRORP_WW = resRORP_WW[order(resRORP_WW$FDR),]

bc <- bacon(teststatistics = resProlif_WW$Z,verbose=T)
resProlif_WW$Bacon <- pval(bc)
resProlif_WW$FDR <- p.adjust(resProlif_WW$Bacon,method='BH')
resProlif_WW = resProlif_WW[order(resProlif_WW$FDR),]

require(xlsx)
write.xlsx(resRORS_WW,'WW_TWAS_multigenescores.xlsx',sheetName = 'RORS',row.names=F,col.names=T)
write.xlsx(resRORP_WW,'WW_TWAS_multigenescores.xlsx',sheetName = 'RORP',append=T,row.names=F,col.names=T)
write.xlsx(resProlif_WW,'WW_TWAS_multigenescores.xlsx',sheetName = 'Proliferation Score',append=T,row.names=F,col.names=T)

