
#Installing the packages
setRepositories() #no. with the repositories will show up 
1 2 3 4 5 6 7 8 #This will enforce all the options given in the package -> install.

#This is how we load the packages
install.packages(
  "datapasta", 
  repos = c(mm = "https://milesmcbain.r-universe.dev", getOption("repos")))

install.packages('devtools')
install.packages("gplots")
devtools::install_github("aljrico/gameofthrones")
devtools::install_github("ropensci/plotly") # you will probably benefit from the latest version of plotly
devtools::install_github('talgalili/heatmaply')
devtools::install_github("rstudio/d3heatmap")
#GSEABase
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSEABase")
#GSVA
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSVA")

install.packages("gprofiler2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

install.packages("msigdbr")


library(rhdf5)
library(tximport)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(biomaRt)
library(beepr)
library(datapasta)
library(tidyverse)
library(edgeR)
library(matrixStats)
library(cowplot)
library(limma)
library(RColorBrewer)
library(gplots)
library(gameofthrones)
library(heatmaply)
library(d3heatmap)
library(DT)
library(GSEABase)
library(Biobase)
library(GSVA)
library(gprofiler2)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(dbplyr)
library(dplyr)
#More info in Lecture 4 part 1 
#studydesign should contain as much as info as possible - like experimental variable/notes - so its a text file which can be uploaded to R. 

#  targets <- read_tsv("studydesign.txt") ?read_tsv ??read_tsv
#Target has column of sample, accession no. and group which consists of hralthy or unhealthy individuals
targets_Sample_Labels <- read_csv("sample_labels_PCA.csv")
targets_geneID <- read_csv("sample_labels_PCA_geneID.csv")
targets <- read_tsv("sample_labels_R.txt")
getwd()
setwd("C:/Users/asim/Documents/ReviewPaper_antioxidants/Obesity_RNA_Seq_Data")
targets
targets_geneID
targets_Sample_Labels

path <- file.path(targets$Sample_labels, "")

#path <- file.path(targets$sample, "abundance.tsv")
#Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
#Tx <- as_tibble(Tx)
#Tx <- dplyr::rename(Tx, target_id = tx_id)
#Tx <- dplyr::select(Tx, "target_id", "gene_name")
#Txi_gene <- tximport(path, 
#                    type = "kallisto", 
#                     tx2gene = Tx,
#                     txOut = FALSE, 
#                     countsFromAbundance = "lengthScaleTPM", 
#                     ignoreTxVersion = TRUE
#                     )

myTPM <- read.csv("FPKM_B6_LVR_R.csv", row.names = "geneID")#row.names defines which particular column okwe need to put as teh label in R. 
my_count <- read.csv("readcount_B6_LVR_R.csv", row.names = "geneID")
colSums(myTPM)
colSums(my_count)
matrix_myTPM <- as.matrix(myTPM) #this converts the data frame or anything into a matrix form


myTPM.stats <- transform(matrix_myTPM, 
                         SD=rowSds(matrix_myTPM), 
                         AVG=rowMeans(matrix_myTPM), 
                         MED=rowMedians(matrix_myTPM))

class(matrix_myTPM) #this is to check the status of the table - is it a data frame or is it a matrix. 
str(myTPM) #to check what is the string for the elemenets in the table.

head(myTPM.stats) #to check the top 6 rows of any data matrix

#my_count <- read.delim(file.choose(), row.names=1)
myDGEList <- DGEList(my_count)#DGEList is from EdgeR #what we get is a data matrix also we have 2 parts - one count and the other is the sample
save(myDGEList, file = "myDGEList") #cehck the file in the directory 
load(file = "myDGEList")#R object is taken from the harddrive into the R


# data_stat <- data.matrix(FPKM_B6_LVR_R, rownames.force = NA)#Since it was not accepting in the fromat that was been uploaded directly from the Excel, we need to convert this into a matrix. This code helps to convert it into the matrix.
# my_count <- data.matrix(readcount_B6_LVR_R, rownames.force = NA) #When we convert the data into a matrix, it needs to be in the numeric order, it doesnt accept the character vector.

#Highest std will show the highest expression of the gene.

ggplot(myDEG.stats) +
  aes(x = SD, y = MED) + 
  geom_hex(shape=25, size=0.2) + 
  geom_smooth(method=lm) +
  geom_hex(show.legend = FALSE) + 
  labs(y="Median", x = "Standard deviation", 
       title= "Log10(FPKM+1)", 
       subtitle= "Unfiltered data - Non DEG", 
       caption= "RNA Seq B6 LVR") + 
  theme_classic() + 
  theme_dark() + 
  theme_bw()
?file.choose()



cpm <- cpm(myDGEList) #Counts per Million or Reads per Kilobase per Million
colSums(cpm) #Now they all tally up to a million, cpm is from EdgeR
log2.cpm <- cpm(myDGEList, log = TRUE) #Log2 (is default so when we say True it goes for Log2) is important foe dealing with data of heteroscedasticity
??cpm

#Coerce is useful becuase the data we have currently is not in the form of dataframe. Its just a numeric data where the genes are rownames and smaples are column names. 
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID") #Dataframe does not understand rownames, it believes that it should be just a part of the data. so, we will drop the row names.   
log2.cpm.df #df stands for data frame

colnames(log2.cpm.df) <- c("geneID",sampleLabels)
log2.cpm.df

#Start correcting from here#
#Now we will add our sample names to the dataframe, 
#colnames(log2.cpm.df) <- as_tibble(log2.cpm, rownames = "geneID")
#sampleLabels <- targets$Sample_labels #Sample labels is a column coming from the file sample_labels_R
#sampleLabels_gene_ID <- targets$geneID
#colnames(log2.cpm.df) <- c("geneID", sampleLabels) #concatanate is done by c

#To convert the data in a tidy form - having all the columns designated by their category
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, 
                                 cols = B1_C:Sac5A_C4, 
                                 names_to = "samples", 
                                 values_to = "expression")

Sac1A_C1vsB1_C_DEG_all$X=NULL #This is to remove the unnecessary columns
Sac1A_C1vsB1_C_DEG_all$X.1=NULL
Sac1A_C1vsB1_C_DEG_all$X.2=NULL
Sac1A_C1vsB1_C_DEG_all
?pivot_longer #to know the applicaiton of the fucntion 



p1 <- ggplot(log2.cpm.df.pivot) + 
  aes(x=samples, y=expression, fill=samples ) + 
  geom_violin(trim = FALSE, show.legend = FALSE) + 
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) + 
  labs(y="Log2 expression", x = "Sample", 
       title = "Log2 count per million (CPM)",
       subtitle = "Unfiltered, non-normalized", 
       caption = paste0("Produced on ", Sys.time())) + 
  theme_bw() + 
  coord_flip()

table(rowSums(myDGEList$counts==0)==17) #This shows the no. of genes or transcripts which have no reads at all in 17 samples. 
#22879 (false) were not 0 count but 12396 of the samples were 0 count in 17 samples. 

keepers <- rowSums(cpm>1)>=8 #We are filtering now. We can filter using this function in any case. But we are doing for cpm. 
#How many cpm are more than 1 count in 8 samples (because we have 8 samples for both HFD and STD group - we can take 17 samples as well but we want to avoid the genes with the low count to be omitted)


myDGEList.filtered <- myDGEList[keepers,]#Subsetting (shown by square brackets) - whcih rows adn columns we need to remove or keep
dim(myDGEList.filtered) #To check the dimension of the table, so now we have 12413 rows in 17 samples

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels) 
class(log2.cpm.filtered.df)

log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, 
                                           cols = B1_C:Sac5A_C4,
                                             names_to = "samples", 
                                           values_to = "expression")

p2 <- ggplot(log2.cpm.filtered.df.pivot) + 
  aes(x=samples, y=expression, fill=samples ) + 
  geom_violin(trim = FALSE, show.legend = FALSE) + 
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) + 
  labs(y="Log2 expression", x = "Sample", 
       title = "Log2 count per million (CPM)",
       subtitle = "Filtered, non-normalized", 
       caption = paste0("Produced on ", Sys.time())) + 
  theme_bw() + 
  coord_flip()

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM") #Normalization 
#In result Normalization factor is not 1.

log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)

log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, 
                                           cols = B1_C:Sac5A_C4,
                                           names_to = "samples", 
                                           values_to = "expression")

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) + 
  aes(x=samples, y=expression, fill=samples ) + 
  geom_violin(trim = FALSE, show.legend = FALSE) + 
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) + 
  labs(y="Log2 expression", x = "Sample", 
       title = "Log2 count per million (CPM)",
       subtitle = "Filtered, Normalized", 
       caption = paste0("Produced on ", Sys.time())) + 
  theme_bw() + 
  coord_flip()

#There will be a subtle change in the normalized vs non-normalized graph. 

plot_grid(p1, p2, p3, labels = c('A','B','C'), label_size = 12)
print("Step 2 complete")


#PCA plot 
install.packages("DT")
install.packages("plotly")
install.packages("gt")

library(tidyverse)
library(DT)
library(plotly)
library(gt)

targets1 <- read_csv("sample_labels_PCA.csv")
getwd()
setwd("C:/Users/asim/Documents/ReviewPaper_antioxidants/Obesity_RNA_Seq_Data")
targets

sampleLabels <- targets$Sample_labels #Sample labels is a column coming from the file sample_labels_R
sampleLabels

group <- targets$Treatment 
group <- factor(group) #This will provide the level that means lets say if there are 2 groups - STD and HFD, levels will distinguish STD and HFD

log2.cpm.filtered.norm.df

distance <- dist(t(log2.cpm.filtered.norm), method = "maximum")
clusters <- hclust(distance, method = "complete")
plot(clusters, labels=sampleLabels) #something wrong with labels, since its giving an error. Tried running without labels it shows the desired output though. However, if i just go with plot(clusters) its giving me the desired output.
?plot()
clusters <- hclust(distance, method = "average") #The outliers will prominently show up away from the other population

pca.res <- prcomp(t(log2.cpm.filtered.norm), scale=F, retx=T)
ls(pca.res)
summary(pca.res) #there are p17 because there are 17 samples
pca.res$rotation
pca.res$x

screeplot(pca.res)
pc.var <- pca.res$sdev^2
pc.per <- round(pc.var/sum(pc.var)*100,1)
pc.per

pca.res.df <- as_tibble(pca.res$x)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2,label=sampleLabels,color=group) +
  geom_point(size=4) + 
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + #we can twitch the PCs by playing around different combinations like PC2 vs PC3 to check the best variation observed 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())
       ) +
  coord_fixed () + #regardless of the windows size the figure would be of the same size - imp to export the image.
  theme_bw()

pca.res.df <- pca.res$x[,1:4] %>%
  as_tibble() %>%
  add_column(sample = sampleLabels, 
             group = group)

pca.pivot <- pivot_longer(pca.res.df, 
                          cols = PC1:PC4, 
                          names_to = "PC", 
                          values_to = "loadings")

#instead of going through different combination of PCs variation/combination we can go fro plotting the small multiples plot - which has been used since 19th century.
ggplot(pca.pivot) + 
  aes(x=sample, y=loadings, fill=group) + 
  geom_bar(stat="identity") + 
  facet_wrap(~PC) + 
  labs(title = "PCA 'small multiples' plot", 
       caption = paste0("produced on", Sys.time())) + 
  theme_bw() + 
  coord_flip()

mydata.df <- log2.cpm.filtered.norm.df %>%
  mutate(STD.AVG = (Sac2A_C1 + Sac2A_C2 + Sac2A_C3 + Sac2A_C4)/4,
         HFD.AVG = (Sac4A_C1 + Sac4A_C2 + Sac4A_C3 + Sac4A_C4)/4,
         LogFC = (STD.AVG - HFD.AVG)) %>%
mutate_if(is.numeric, round, 2)
mydata.df
#colnames(mydata.df) 

#not useful for what we will do next, but just for the practice
mydata.sort <- mydata.df %>%
  dplyr::arrange(desc(LogFC)) %>%
  dplyr::select(geneID, LogFC)

#This way we can select individual gene and check their expressions.
mydata.filter <- mydata.df %>%
  dplyr::filter(geneID=="ENSMUSG00000054417"| geneID=="ENSMUSG00000038917") %>%#we can give any genes we want, | is a boolean operator meaning or. %>% is a pipe fucntion. 
  #ALso inorder to check the exact text matches, we can seearch, follow the codes in the lecture.
  dplyr::select(geneID, STD.AVG, HFD.AVG, LogFC) %>%
  dplyr::arrange(desc(LogFC)) 

#For publication quality tables - we can use the following gt code
gt(mydata.filter)

#With few more options to it - there is some error, but check lecture whenever there is time - also interactive tables can be made like drop down list etc. 
gt()%>%
  fmt_number(columns=2:4, decimals = 1) %>%
  tab_header(title = md("**blah blah**"), 
             subtitle = md("*something something*"))%>%
  tab_footnote(
    footnote = "this is a foonote", 
    locations = cells_body(
      columns = vars(geneID),
      rows = c(6,7))) %>%
  tab_footnote(
    footnote = "another footnote here", 
    location = cells_body(
      columns = vats(geneID), 
      rows = c(2:10))) %>%
  tab_source_note(
    source_note = md("Reference: Dave *et al*., (2021)")
  )
  
#My Fav - making the interactive plots

myplot <- ggplot(mydata.df) + 
  aes(x=STD.AVG, y=HFD.AVG) + 
  geom_point(shape=16, size=1) + 
  ggtitle("STD vs HFD") + 
  theme_bw()
  
#now we will convert this to make the interactive plot, but we cannot see which gene it is.
ggplotly(myplot) 
  
#Now we will be able to see the gene symbol.. We can also add the p value and FDRs - coming up in the next module
myplot <- ggplot(mydata.df) + 
  aes(x=STD.AVG, y= HFD.AVG,
      text = paste("geneID", geneID)) + 
  geom_point(shape=16, size=1) + 
  ggtitle("STD vs HFD") + 
  theme_bw()

ggplotly(myplot) 

#how to access publicaly accessible RNA Seq data - check out the module fro more infromation - also that way we can compare the gene that we have against the data on the portal.

#Kallisto bootstraps - is a way of doing the quantification - uncertainity of the count and we can analyze DGE. Its not in the module but its a great way - still its developeing. 

# Making the DGE list

library(limma)

design <- model.matrix(~0 + group) #Gives the smaples and group designation with 0 and 1 - this way we have a linear model.
colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE) #voom - variance stabalize data and model the mean varaince realtionship. 
#Here lowly expressed genes have higher varriance and we see a skewed line. 
#This changed our DGE list -  we have a new elemebnt - weights that helps apply adn improve our ability to choose DEG list. 
fit <- lmFit(v.DEGList.filtered.norm, design) #The are always built around count and not the fpkm etc. so we show th ecounts to them
#gives us all the stats

#From here only HFD and STD is considered
contrast.matrix <- makeContrasts(Disease = HFD - STD, levels=design) #this is to make a comparison, the pairwise comparison can be done here. Put column headers as input. 
#HFD-STD will give the genes whcih are highly upregulated in the HFD state (with a positive sign), if it is STD - HFD, it will give genes which are upregulated in STD. 

fits <- contrasts.fit(fit, contrast.matrix)
#again we get the stats

ebFit <- eBayes(fits)
  
#Now we are done with creating the DEG list!!!!
#Lets try looking at it..

#Look at the top DEG by pairwise comparison
myTopHits <- topTable(ebFit, adjust = "BH", coef=1, number=40000, sort.by = "logFC") #Coef is contrast disease (designated as 1, no. 2 would be the next condition and so on).
#Number is how many genes should it return
#We have a p value but we are doing multiple test - which can give us type 1 error. 5% of the time there will be a false positive. 
#Therefore, we are doing multiple test correction method which is Benjamini-Hochberg method represented by BH...
#So with BH method we get FDR or adjusted P vlaue. 

myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

gt(myTopHits.df)
  

#Now we will plot the volcano plot

vplot <- ggplot(myTopHits.df) + 
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("GeneID:", geneID)) + 
  geom_point(size = 2) + 
  geom_hline(yintercept = -log10(0.1), linetype="longdash", colour="grey", linewidth=1) + 
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", linewidth=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A",linewidth=1) +
  #annotate("rect", xmin = 1, xmax = 6, ymin = -log10(0.1), ymax = 3, alpha = .2, fill = "#BE684D") + 
  #annotate("rect", xmin = -1, xmax = -6, ymin = -log10(0.1), ymax = 3, alpha = .2, fill = "#2C467A") + 
  labs(title = "Volcano plot", 
       subtitle = "Cutaneous leishmaniasis", 
       caption = paste0("produced on ", Sys.time())) + 
    theme_bw()

#genes on right is upregulated in disease and left is less regulated in disease. 

ggplotly(vplot)

#Now we will pull out the genes with a cut offs
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=0.5)
#It will go row by row and if those transcripts meet  that cut off of HFD - STD
# 0 means it did not, but if it is -1 or 1, that means it passed (sign like negative shows the direction, like downregulation)
#Now we can make a venn diagram - 1 and 0 in 2 conditions that means not intersecting the venn diagram.  

head(results)
summary(results) #down means downregulated 
vennDiagram(results, include="both") 
#venndiagram can be misleading if there are too many circles - since there can be a condition where one made through -1 or +1 but it did not in another condition. But still it will show up in the venndiagram, people think its a DEG but its not. Since it passes through the cut off in only one condition but not the other. 

#Retrieve expression data for DEG
head(v.DEGList.filtered.norm$E) #$E is where the expression data resides
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,] #We want the rows in which the result ([,1 is showing the column]) is not equal to 0
#This code can be written as:
#diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0 | results[,2] !=0,]
#Meaning of the code: [,1] is column 1, [,2] is column 2. !=0 is not equal to 0, | is or
#This is a condition where either column 1 is not equal to 0 or column 2 is not equal to 0. 

head(diffGenes)
dim(diffGenes)

diffGenes.df <- as_tibble(diffGenes, rownames="geneID")

#now we will create an interactive table

datatable(diffGenes.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: DEGs in HFD vs STD', 
          options = list(keys=TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns = c(2:18), digits = 2)
#Now we get an interactive table where we can search

#I have skipped DTU analysis

#We will now chcek the same functional genes - Lec 10, part 1 - we will do clustering now. 

myheatcolors1 <- bluered(75) #from gplots package - going to ramp 2 colors together - 75 is the gradiation
display.brewer.all()
display.brewer.all(colorblindFriendly = TRUE)

myheatcolors2 <- colorRampPalette(colors =c("blue", "white", "red"))(100)

myheatcolors3 <- brewer.pal(name="RdBu", n=-11) #use sip ot search online with the code

myheatcolors3 <- c("#fed976", "#268f9c")

got_palette <- got(75, option = "Arya") #GOT pallete 

#Now we will make the heat map - Module 10 - part 2 


results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=0.5)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
dim(diffGenes)

clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete")

clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")

module.assign <- cutree(clustRows, k=2) #we will get 2 clusters - HFD and STD, if there is more than 2 clusters then just try changing overcut or under cut what is more appropriate..
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9)
module.color <- module.color[as.vector(module.assign)]

heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color, 
          col=myheatcolors1, scale='row', labRow=NA, 
          density.info="none", trace="none", 
          cexRow=1, cexCol=1, margins=c(8,20))

#HTML widget version of the heatmap

d3heatmap(diffGenes, 
          colors = myheatcolors2, 
          Rowv=as.dendrogram(clustRows), 
          row_side_colors = module.color, 
          scale='row')

colnames(diffGenes) <- targets$Treatment

diffGenes.AVG <- avearrays(diffGenes)
heatmap.2(diffGenes.AVG, 
          Rowv=as.dendrogram(clustRows), 
          #Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color, 
          col=myheatcolors1, scale='row', labRow=NA, 
          density.info="none", trace="none", 
          cexRow=1, cexCol=1, margins=c(8,20))

#now we are going to slice out 
names(module.color) <- names(module.assign)
?as_tibble
module.assign.df <- as_tibble(as.list(module.assign))
module.assign.pivot <- pivot_longer(module.assign.df, 
                                    cols = 1:17, 
                                    names_to = "geneID", 
                                    values_to = "module")


module.assign.pivot <- module.assign.pivot %>%
  mutate(moduleColor = case_when(
    module == 1 ~"#FF0099", 
    module == 2 ~"#FF9900"
  ))

ggplot(module.assign.pivot) + 
  aes(module) + 
  geom_bar(aes(fill=moduleColor)) + 
  theme_bw()

modulePick <- 2 # we would like to pick the 2 

myModule <- diffGenes[names(module.assign[module.assign %in% modulePick]),]
hrsub <- hclust(as.dist(1-cor(t(myModule), method="pearson")), method="complete")

#somehow this coed is not working as expected or shown in the module
heatmap.2(myModule, 
          Rowv=as.dendrogram(hrsub), 
          Colv=NA,
          labRow=NA, 
          col=myheatcolors2, scale="row",
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20))

#there is more to heat map but I did not go over it all, go through more in the module..

#now we will cluster the DEG - Module 10 - part 3 - Bi-CoPaM approach - Its in immunity paper by Andrea. 
#its not describe much but woul db einteresting to know how its done - find more


#Functional enrichment analysis - Module 11 - part 1 - GSEA 
#the enrichment score in GSEA can be explained by running sum staistics - eg. black jack - low cards +1 and high cards -1.. Mid cards like 7, 8, 9 - nothing.
#There are 2 GSEA - 1. Self-contained - Are any o fthe genes in my set differnetially expressed? Does not consider genes outside your set. 
#and 2. Competitive (most popular) - Are genes in my set more commonly differnetially exporessed than genes outside the set? Does consider genes outside your set..

#Module 11 - part 2 - GSEA 
#We will do the analysis using Gprofiler - the same that Dr. Park mentioend
myTopHits <- topTable(ebFit, adjust="BH", coef=1, number=50, sort.by="logFC") #Still working with HFD vs STD
#gost fucntion - helps to run from gprofiler - R is easier, faster and more customizable - but we can also do using the webpage. 
gost.res <- gost(rownames(myTopHits), organism = "mmusculus", correction_method = "fdr")
myModule
mygostplot <- gostplot(gost.res, interactive = F, capped = F) #we can put capped as F as well. THat way nothing will be capped and the dots will spread according to the axis .

publish_gostplot(
  mygostplot, 
  highlight_term = c("GO:0033993"), 
  filename = NULL,
  width = NA,
  height = NA)

gost.res <- gost(rownames(myModule), organism = "mmusculus", correction_method = "fdr")
mygostplot <- gostplot(gost.res, interactive = F, capped = F) #we can put capped as F as well. THat way nothing will be capped and the dots will spread according to the axis .
publish_gostplot(
  mygostplot, 
  highlight_term = c("GO:0033993"), 
  filename = NULL,
  width = NA,
  height = NA)









