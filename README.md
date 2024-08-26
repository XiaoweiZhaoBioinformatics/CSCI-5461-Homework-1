# CSCI 5461 Homework 1
 Identification of differentially expressed genes in patients who are diagnosed with Glioblastoma Multiforme (GBM)
1. Acquiring the gene expression data: Go to the GSE6944 repository page and read the
descriptions of the dataset
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944). We have already
downloaded the count and clinical data in the file getRNASeq_data.R and included it for
you in HW1-GSE62944-count.csv and HW1-GSE62944-clinical.csv
respectively. For the count data, the rows are the genes and the columns are the patient
identifiers. For the clinical data, there are two columns, the “sampleName” and “Group” that
identifies patients as “short” or “long” survivors. The final piece of data is the analysis from
the DESeq2 method. We have completed that analysis and saved it to the file HW1-
DESEQ2.csv. If you’re interested in the details of the pipeline that was used to process the
raw RNAseq data, you can read about those details in this paper
(https://www.ncbi.nlm.nih.gov/pubmed/26209429).
