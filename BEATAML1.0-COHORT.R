library(SummarizedExperiment)
library(TCGAbiolinks)
require(EDASeq)
require(dplyr)
require(NOISeq)
library(DESeq2)

qry.rna <- GDCquery(project = "BEATAML1.0-COHORT",
                    data.category= "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts",)

GDCdownload(qry.rna)

d <- qry.rna[[1]][[1]]
table(as.factor(d$sample_type))

rnas <- GDCprepare(qry.rna, summarizedExperiment = TRUE)

data <- assay(rnas)
rownames(data) <- rowData(rnas)$external_gene_name
rownames(rnas) <- rowData(rnas)$external_gene_name
head(rownames(data))

dim(data)
dataFilt <- TCGAanalyze_Filtering(tabDF = data,
                                  method = "quantile",
                                  qnt.cut = 0.25)
threshold <- round(dim(data)[2]/2)
ridx <- rowSums(dataFilt == 0) <= threshold
dataFilt <- dataFilt[ridx, ]
print(dim(dataFilt))    
ridx <- rowMeans(dataFilt) >= 10
dataFilt <- dataFilt[ridx, ]
print(dim(dataFilt))
rnas <- rnas[rownames(rnas) %in% rownames(dataFilt), ]

dataNorm <- TCGAanalyze_Normalization(tabDF = assay(rnas), geneInfo = geneInfo)
dim(dataNorm)
rnas <- rnas[rownames(rnas) %in% rownames(dataNorm), ]
rnas <- rnas[!duplicated(rownames(rnas)),]

#ANOTACION CON BIOMART
data.biomart <- assay(rnas)
rownames(data.biomart) <- rowData(rnas)$ensembl_gene_id # pongo el EnsemblID en lugar del external name
rownames(data.biomart)
data.bio <- cbind(data.biomart, c(rownames(data.biomart))) # agregué una columna para que merge dejara de darme error
View(data.bio)
annot<-read.delim(file="mart_export.txt", sep="\t")
names(annot)<-c("Gene.name", "Chr", "Start", "End", "GC", "Type", "ensembl_gene_id")
annot$Length<-abs(annot$End - annot$Start)
rownames(annot) <- annot$ensembl_gene_id
dim(annot)
dim(data.bio)
annot <- annot[rownames(annot) %in% rownames(data.bio), ] # mantengo solo con la información de Biomart de los genes que ya estan filtrados y normalizados

uniq.annot <- length(unique(annot$ensembl_gene_id)) == nrow(annot)
if(uniq.annot) {
  cat('Unique EnsemblIDs in annotation file\n')
} else {
  cat('Repeated EnsemblIDs in annotation file. \n')
  stop()
}

rnas.Biomart<-merge(x=annot, y=data.bio, by.x="ensembl_gene_id", by.y="V511",
                       all=FALSE, all.x=TRUE, all.y=FALSE, sort=FALSE)
extra.rows <- nrow(rnas.Biomart)-nrow(data.bio) 
cat('There are ', extra.rows, ' extra rows in the counts matrix.\n') # perdí 28 genes

rnas <- rnas[rownames(rnas) %in% rnas.Biomart$ensembl_gene_id, ] #aseguro que no hayan duplicados de EnsemblID
dim(rnas)

ln.data <- withinLaneNormalization(assay(rnas), 
                                   rnas.Biomart$Length, which = "full")
gcn.data <- withinLaneNormalization(ln.data , rnas.Biomart$GC,
                                    which = "full")
norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
noiseqData <- NOISeq::readData( norm.counts , factors = as.data.frame(rnas$primary_diagnosis), as.data.frame(rnas$sample_type)) # ¿puedo añadir aquí el diagnosis y sample type a la vez?
mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE)
assay(rnas) <- exprs(mydata2corr1)

# Filtro Primary Diagnosis

AML.NPM1 <- rnas[ ,rnas$primary_diagnosis == "Acute myeloid leukemia with mutated NPM1"]
AML.mrc <- rnas[ ,rnas$primary_diagnosis == "Acute myeloid leukemia with myelodysplasia-related changes"]

# Filtro Sample Type 

PB.normal <- rnas[ ,rnas$sample_type == "Blood Derived Normal"]

AML.NPM1.PB <- AML.NPM1[ ,AML.NPM1$sample_type == "Primary Blood Derived Cancer - Peripheral Blood"]
AML.NPM1.rPB <- AML.NPM1[ ,AML.NPM1$sample_type == "Recurrent Blood Derived Cancer - Peripheral Blood"]
AML.NPM1.BM <- AML.NPM1[ ,AML.NPM1$sample_type == "Primary Blood Derived Cancer - Bone Marrow"]
AML.NPM1.rBM <- AML.NPM1[ ,AML.NPM1$sample_type == "Recurrent Blood Derived Cancer - Bone Marrow"]

AML.mrc.PB <- AML.mrc[ ,AML.mrc$sample_type == "Primary Blood Derived Cancer - Peripheral Blood"] 
AML.mrc.rPB <- AML.mrc[ ,AML.mrc$sample_type == "Recurrent Blood Derived Cancer - Peripheral Blood"]
AML.mrc.BM <- AML.mrc[ ,AML.mrc$sample_type == "Primary Blood Derived Cancer - Bone Marrow"]
AML.mrc.rBM <- AML.mrc[ ,AML.mrc$sample_type == "Recurrent Blood Derived Cancer - Bone Marrow"]
