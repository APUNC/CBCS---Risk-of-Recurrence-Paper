require(data.table)
require(boot)
setwd('/proj/milovelab/patel')
geneList = c(rep('MCM10',3),
          rep('FAM64A',3),
          rep('CCNB2',3),
          rep('MMP1',6),
          rep('VAV3',3),
          rep('PCSK6',3),
          rep('GNG11',3))
popList = c(rep('WW',3),
            rep('WW',3),
            rep('WW',3),
            rep('WW',3),
            rep('AA',3),
            rep('WW',3),
            rep('WW',3),
            rep('WW',3))
outcomeList = rep(c('ROR-S (Subtype Only)',
                    "Proliferation Score",
                    "ROR-P (Subtype + Proliferation)"),
                  8)

calculateTWAS <- function(effects,
                          snps,
                          analyze,
                          outcome,
                          indices){
    effects = effects[indices]
    grex = as.numeric(effects %*% as.matrix(snps[,-1]))
    colnames(analyze)[which(colnames(analyze) == outcome)] = 'Outcome'
    reg = lm(analyze$Outcome ~ 
                 scale(grex) + 
                 analyze$agesel + 
                 as.factor(analyze$er) + 
                 as.factor(analyze$stage) + 
                 as.factor(analyze$phaseCode))
    colnames(analyze)[which(colnames(analyze) == 'Outcome')] = outcome
    return(coef(reg)[2])
}

for (i in 1:length(geneList)){
    modelFolder = '/proj/milovelab/bhattacharya/CBCS_TWAS_Paper-master'
    gene = geneList[i]
    pop = popList[i]
    population = ifelse(pop == 'AA','african','euro')
    race = ifelse(pop == 'AA','Black','White')
    outcome = outcomeList[i]
    outFile = 'TWAS_results_ROR_permutation.tsv'
    
    analyze = fread('analyzeOnThis.csv')
    analyze = subset(analyze,Race == race)
    geneModel = fread(file.path(modelFolder,
                                paste0(gene,'_betaMatrix_',pop,'.csv')))
    geneModel = geneModel[-1,]
    chr = unique(geneModel$Chromosome)[1]
    
    
    snps = fread(paste0('ncbcs_',
                           population,
                           '_oncoarray_imputed_dosages_',
                           chr,
                           '.txt.gz_dosages.txt'))
    snpList = intersect(geneModel$Feature,snps$SNP)
    
    snps = subset(snps,SNP %in% snpList)
    geneModel = subset(geneModel, Feature %in% snpList)
    snps = snps[match(snps$SNP,geneModel$Feature),]
    analyze = analyze[match(analyze$BCAC_ID,
                            colnames(snps)[-1]),]
    
    colnames(analyze)[which(colnames(analyze) == outcome)] = 'Outcome'
    colnames(analyze)[which(colnames(analyze) == gene)] = 'Gene_Int'
    original_reg = lm(Outcome ~ 
                          scale(Gene_Int) +
                          agesel +
                          as.factor(er) + 
                          as.factor(stage) +
                          as.factor(phaseCode),
                      data = analyze)
    colnames(analyze)[which(colnames(analyze) == 'Outcome')] = outcome
    colnames(analyze)[which(colnames(analyze) == 'Gene_Int')] = gene
    effect_size = as.numeric(round(coef(original_reg)[2],2))
    conf_int = paste0('(',paste(round(confint(original_reg)[2,],2),collapse = ', '),
                      ')')
    snps = as.data.frame(snps)
    
   
    
    permutation = boot(data = geneModel$Beta,
                       statistic = calculateTWAS,
                       R = 5e3,
                       sim = 'permutation',
                       snps = snps,
                       analyze = analyze,
                       outcome = outcome)
    perm_pval = (5e3 * mean(abs(permutation$t) >
                                   abs(coef(original_reg)[2])) + 1)/(5e3+1)
    df = data.frame(Population = pop,
                    Outcome = outcome,
                    Gene = gene,
                    Effect = effect_size,
                    CI = conf_int,
                    Permutation_P = perm_pval)
    fwrite(df,
           outFile,
           append=T,
           quote=F,
           row.names=F)
}

a = fread(outFile)
a$Permutation_FDR = p.adjust(a$Permutation_P,method='BH')
a$Permutation_Bonf = p.adjust(a$Permutation_P,method='bonf')
fwrite(a,outFile,col.names=T,row.names=F,quote=F,sep='\t')