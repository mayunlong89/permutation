#Rscript
#Author: Yunlong Ma
#100,000 times in silico permutation analysis for genes identified from Sherlock and MAGMA analysis


#Set the work directory
setwd("C:\\Users\\Administrator\\Desktop\\06-Simulation_analysis")
set.seed(12345)

#Part I Read data on significant genes and background genes

#Read significant genes of Geneset #1
Sig_1 <- read.table("Zeller_sig_CAD.txt", header=T)
Sig_Zeller <- Sig_1$Gene_name

#Read significant genes from Geneset #2
Sig_2 <- read.table("Dixon_sig_CAD.txt", header=T)
Sig_Dixon <- Sig_2$Gene_name

#Read background genes of Dixon eQTL data
Backgroud_2<- read.table("Dixon_all_CAD.txt", header=T)
Backgroud_Dixon <- Backgroud_2$Gene_name

#Read significant genes from Geneset #3
Sig_3 <- read.table("MAGMA_sig_CAD.txt", header=T)
Sig_MAGMA <- Sig_3$Gene_name

#Read background genes of MAGMA analysis on CAD GWAS summary dataset
Background_3<- read.table("MAGMA_all_CAD.txt", header=T)
Background_MAGMA <- Background_3$Gene_name

#Calculate the numebr of genes in each gene set
len_Sig_Zeller <- length(Sig_Zeller)
len_Sig_Dixon <- length(Sig_Dixon)
len_Backgroud_Dixon <- length(Backgroud_Dixon)
len_Sig_MAGMA <- length(Sig_MAGMA)
len_Background_MAGMA <- length(Background_MAGMA)


#Part II establish a function for permutation analysis

#Permutation Function
Permut_analysis <- function(x,y,z){
       
   random_selected_genes <- sample(x,y)
   
   temp <- match(random_selected_genes,z)
   
   random_overlaped_genes <- na.omit(temp)
    
    num<-length(random_overlaped_genes)
   
  return(num)
  
}


#100000 times permutation analysis for Sherlock analysis of Zeller VS. Dixon
results_Dixon <- replicate(100000,Permut_analysis(Backgroud_Dixon,len_Sig_Dixon,Sig_Zeller))

#100000 times permutation analysis for Sherlock analysis of Zeller VS. MAGMA on CAD GWAS
results_MAGMA <- replicate(100000,Permut_analysis(Background_MAGMA,len_Sig_MAGMA,Sig_Zeller))


#Ploting function
Fig_random <- function(x,y,z){
  
  hist(x, col="red",xlab="Counts of overlapped genes",main=NULL)
  temp1 <- match(y,z)
  Observed_overlaped_genes <- na.omit(temp1)
  Observed_gene_num <- length(Observed_overlaped_genes)
  abline(v=Observed_gene_num,col="darkblue",lty="longdash")
  P_value=length(x[x>Observed_gene_num])/length(x)
  x1= Observed_gene_num
  freq <- table(x)
  y1 = max(freq)
  text(x1,y1,P_value)
  
}


#Visulization for Sherlock analysis of Zeller VS. Dixon
Fig_random(results_Dixon,Sig_Dixon,Sig_Zeller)

#Visulization for Sherlock analysis of Zeller VS. MAGMA on CAD GWAS
Fig_random(results_MAGMA,Sig_MAGMA,Sig_Zeller)


#End

