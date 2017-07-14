Key Learning Outcomes
---------------------

After completing this module the trainee should be able to:

-   Find gene ontology enrichment in a list of differentially expressed
    genes using R-based packages.

-   Running GO enrichment analysis using the web tool DAVID and the web
    tool Gorilla

-   To run webtools such as REVIGO and STRING

Resources You’ll be Using
-------------------------

### Tools Used

Goana from Limma:  
https://bioconductor.org/packages/release/bioc/html/limma.html

DAVID:  
http://david.abcc.ncifcrf.gov

GOrilla:  
http://cbl-gorilla.cs.technion.ac.il

REVIGO:  
http://revigo.irb.hr

STRING:  
http://string-db.org

## Author Information

*Primary Author(s):*
    Susan M Corley s.corley@unsw.edu.au
  
*Contributor(s):*
    Sonika Tyagi sonika.tyagi@agrf.org.au   
    Nandan Deshpande n.deshpane@unsw.edu.au

Introduction
------------

The goal of this hands-on session is to allow you to develop some
familiarity with commonly used, freely available R based packages and
web tools which can be used to gain biological insight from a
differential expression experiment. We will use the differentially
expressed genes (DEGs) identified in the last session. First, we will
look at whether these genes are enriched for gene ontology terms which
gives us some insight as to whether the DEGs are involved in particular
functions. Then we will use a tool that constructs an interaction
network from these genes. This will allow us to identify clusters of
DEGs that are known to interact.

In using any database tools it is always advisable to check whether they
are regularly updated. We suggest that you experiment with more than one
tool.

Gene ontology analysis with GOana
---------------------------------

First we will go back to the R environment and use the function GOana
associated with the limma package. To use this function we need to have
our DEGs annotated with the entrez gene identifier for each gene.We did
this early on in our data processing. We use the fit object (fit_v) generated
using the voom function in limma for this analysis. To obtain more
information regarding the goana function type ?goana within your R
session.

Open the Terminal and go to the `rnaseq/edgeR` working directory:

    cd /home/trainee/rnaseq/edgeR

!!! failure "STOP"
    All commands entered into the terminal for this tutorial should be from
    within the **`/home/trainee/rnaseq/edgeR`** directory.

R (press enter)

Check that the directory you are in contains the above-mentioned fit_v file by
typing:

    ls()
    library(limma)
    library(RColorBrewer)
    library(gplots)
    library(org.Hs.eg.db)

We will use the goana function to obtain the gene ontology terms associated with the DEGs.

    DE_GOana<-goana(fit_v, coef=2, geneid=fit_v$genes$Entrez, FDR=0.05, species = "Hs", trend=F, plot=F )

Now we will look at the most significant biological process (BP)ontology
terms

    DE_GOana_top_BP<- topGO(DE_GOana, ontology=c("BP"), number=150L, truncate.term=50)
    head(DE_GOana_top_BP, 20)      
    DE_GOana_top_BP_down<- topGO(DE_GOana, ontology=c("BP"), sort = "down", number=150L, truncate.term=50)
    head(DE_GOana_top_BP_down, 10)
    DE_GOana_top_BP_up<- topGO(DE_GOana, ontology=c("BP"), sort = "up", number=150L, truncate.term=50)
    head(DE_GOana_top_BP_up, 10)

Rather than looking at biological process (BP) let’s now look at
molecular function (MF) terms

      
    DE_GOana_top_MF_down<- topGO(DE_GOana, ontology=c("MF"), sort = "down", number=150L, truncate.term=50)
    head(DE_GOana_top_MF_down, 10)
    DE_GOana_top_MF_up<- topGO(DE_GOana, ontology=c("MF"), sort = "up", number=150L, truncate.term=50)
    head(DE_GOana_top_MF_up, 10)

### Questions and Answers

Based on the above section:

!!! note "Question"
    What is the general theme emerging when we look at biological process in the down-regulated genes and the up-regulated genes?


!!! success ""
    ??? "**Answer**"
        We see terms involving the cell cycle and cell division are enriched in the down-regulated genes and terms involving ER stress and protein folding are enriched in the up-regulated genes. 

!!! note "Question"
    Looking at the biological process (BP) term “cell division: GO:0051301”, how many down-regulated and up-regulated DEGs are annotated with this term, and what statistical significance is associated with this enrichment?


!!! success ""
    ??? "**Answer**"
        Look for the row for GO:0051301 in the top 20 BP ontology terms. 


There are a number of tools and packages available with the
R-bioconductor repositories that you can use with your R code to run
ontologies and pathway analysis. 

Gene ontology analysis with DAVID
---------------------------------

Click on your Firefox web browser. Go to the DAVID website:
<http://david.abcc.ncifcrf.gov>. Go to your edgeR folder and open the
file `voom_res_sig_lfc.txt` using LibreOffice Calc. For the separator
options in LibreOffice Calc choose *Separated by tab*. Once you have
opened this file copy the *Ensembl Gene Ids* (Column A). This list can
then be pasted into DAVID.  

Screenshots of the DAVID website and the steps to move through the
website are provided in the presentation prepared for this session. Use
that material to work through this exercise.  

We will use the Functional Annotation Clustering tool in DAVID. First we
will uncheck all the defaults and look only at the GO terms involving
biological process. After unchecking all the defaults, expand Gene
Ontology and select GOTERM\_BP\_5. Then select the button Functional
Annotation Clustering.  

This will bring up a screen where GO terms are clustered. Statistical
testing is performed to assess whether the GO terms are more enriched in
the list of DEGs than would be expected by chance. You will see a column
of P\_Value and also adjusted P values. Have a look at the brief
description of the statistical test used in DAVID
(<https://david.ncifcrf.gov/helps/functional_annotation.html>).

### Questions and Answers

Based on the above section:

!!! note "Question"
    What functional themes emerge in Cluster 1 and Cluster 2?


!!! success ""
    ??? "**Answer**"
        Cluster 1 is enriched for cell singling and Cluster 2 is enriched for cell cycle. Note that answers may vary based on the genes entered and the options selected. 

!!! note "Question"
    How many of the DEGs are annotated with the term “intracellular signal transduction” and name five genes?


!!! success ""
    ??? "**Answer**"
        Looking at the first row of the DAVID results we see that 181 DEGs are annotated with “Intracellular signal transduction”. The first 5 genes listed are DHCR24, HMGCR, ADAP12, ARL5A, ARFGEF3.

!!! note "Question"
    Does this seem to be sensible in an experiment that looks at the response of cancer cells to a stimulant?


### Gene ontology analysis with GOrilla

Click on your Firefox web browser. Go to the GOrilla website:
<http://cbl-gorilla.cs.technion.ac.il>  

For this tool we will use a background list of genes. Open the file
`voom_res.txt` using LibreOffice Calc. Copy the *Ensembl Gene Ids*
(Column A) and paste this into GOrilla as the background set.


As in the previous exercise we will use the DEGs found by voom
(`voom_res_sig_lfc.txt`) as the Target set.  

We will firstly look for enriched GO process terms. Screenshots of the
website showing the steps you need to follow are in the presentation for
this session.  

GOrilla will display the GO term hierarchy. This shows you which terms
are parent and child terms and how the terms are related.Under this you
will find a table of the most significantly enriched GO terms. Have a
look at the DEGs associated with the most enriched clusters.

### Questions and Answers

Based on the above section:

!!! note "Question"
    Find the GO term *regulation of cell proliferation* how far can you trace this back to the parent terms. 

!!! success ""
    ??? "**Answer**"
        regulation of cell proliferation – regulation of cellular process – regulation of biological process – biological regulation – biological process

!!! note "Question"
    What are the direct child terms of *regulation of cell proliferation*?. 

!!! success ""
    ??? "**Answer**"
        regulation of stem cell proliferation, regulation of sooth muscle cell proliferation, positive regulation of cell proliferation, regulation of mesenchymal cell proliferation

!!! note "Question"
    What is the enrichment score for *regulation of cell proliferation* ? How is this calculated (Hint: scroll down the page for the heading Enrichment). 

!!! success ""
    ??? "**Answer**"
        1.58 ((104/870)/(919/12137))


### REVIGO to reduce redundancy and visualise

We can use the results generated by the GOrilla web tool as input to
REVIGO which will summarise the GO data and allow us to visualize the
simplified data. Click on the link Visualize output in REViGO. Follow the
screen shots in the presentation.Go to the treemap view.

!!! note "Question"
    What are the main functional categories emerging in this analysis? 


STRING
------

Using STRING to look at networks that may be formed by the DEGsClick on
your Firefox web browser. Go to the STRING website:
http://string-db.org. For this exercise we will only use the top 500
DEGs. Go to the file (`voom_res_sig_lfc.txt`) copy only the top 500.
These will be pasted as input to the STRING website. Follow the screen
shots in the presentation.  

You will see a large interaction network being built from the 500
DEGs.We will refine this by clicking on the Data settings tab and
selecting *high confidence* (see the screen shot in the presentation).
Look at the gene clusters that are generated.

Find the gene CDK1. Look at the cluster generated around this gene. What
other DEGs are interaction partners of CDK1.  

Click on CDK1 and find what functions it is involved in. Click on the
interaction partners of CDK1 and find their functions. 

!!! note "Question"
    Does this help in further explaining some of the gene ontology results? 


Explore other clusters that are formed in this analysis.
