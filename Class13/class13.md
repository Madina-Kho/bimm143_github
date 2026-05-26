# Lab 13: RNA Seq
Madina Khorami (A18555185)

- [Background](#background)
- [Data Import](#data-import)
- [DESeq Analysis](#deseq-analysis)
- [Volcano Plot](#volcano-plot)
- [Save out results to date](#save-out-results-to-date)
- [Adding Annotaion Data](#adding-annotaion-data)
- [Pathway Analysis](#pathway-analysis)
- [Save our annotated results](#save-our-annotated-results)

## Background

Today we are going to do an RNA-sew analysis of a data set on the
RNA-seq experiment where airway smooth muscle cells were treated with
dexamethasone, a synthetic glucocorticoid steroid with anti-inflammatory
effects.

## Data Import

Let’s read the \`count data and metadata about this eperiment setup
fromt he supplied csv file.

``` r
counts <- read.csv("airway_scaledcounts.csv", row.names=1)
metadata <-  read.csv("airway_metadata.csv")
```

``` r
head(counts)
```

                    SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
    ENSG00000000003        723        486        904        445       1170
    ENSG00000000005          0          0          0          0          0
    ENSG00000000419        467        523        616        371        582
    ENSG00000000457        347        258        364        237        318
    ENSG00000000460         96         81         73         66        118
    ENSG00000000938          0          0          1          0          2
                    SRR1039517 SRR1039520 SRR1039521
    ENSG00000000003       1097        806        604
    ENSG00000000005          0          0          0
    ENSG00000000419        781        417        509
    ENSG00000000457        447        330        324
    ENSG00000000460         94        102         74
    ENSG00000000938          0          0          0

``` r
ncol(counts)
```

    [1] 8

``` r
colnames(counts) == metadata$id
```

    [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

and the metadata that tells us what is actually in the column of our
`count` object.

``` r
metadata
```

              id     dex celltype     geo_id
    1 SRR1039508 control   N61311 GSM1275862
    2 SRR1039509 treated   N61311 GSM1275863
    3 SRR1039512 control  N052611 GSM1275866
    4 SRR1039513 treated  N052611 GSM1275867
    5 SRR1039516 control  N080611 GSM1275870
    6 SRR1039517 treated  N080611 GSM1275871
    7 SRR1039520 control  N061011 GSM1275874
    8 SRR1039521 treated  N061011 GSM1275875

> Q1. How many genes are in this dataset?

There are 38694 genes int this dataset.

> Q2. How many control cell lines do we have?

``` r
sum(metadata$dex == "control")
```

    [1] 4

OR

``` r
table(metadata$dex)
```


    control treated 
          4       4 

> Q3. How would you make the above code in either approach more robust?
> Is there a function that could help here?

- Find the control columns in our `counts` object.
- Extracts just the “control” column values for all genes
- Calculate the average value per gene in these “control” columns.

``` r
## Find the contol columns
control.inds <- metadata$dex == "control"
## Extract the control columns
control.counts <- counts[ , control.inds]
## take the avergae of control column over all genes
control.mean <- rowMeans(control.counts)
```

``` r
head(control.mean)
```

    ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
             900.75            0.00          520.50          339.75           97.25 
    ENSG00000000938 
               0.75 

> Q4. Follow the same procedure for the treated samples (i.e. calculate
> the mean per gene across drug treated samples and assign to a labeled
> vector called treated.mean)

- Find the treated columns in our `counts` object.
- Extracts just the “treated” column values for all genes
- Calculate the average value per gene in these “treated” columns.

``` r
treated.inds <- metadata$dex =="treated"

treated.counts <- counts[ , treated.inds]

treated.mean <- rowMeans(treated.counts)
```

``` r
head(treated.mean)
```

    ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
             658.00            0.00          546.00          316.50           78.75 
    ENSG00000000938 
               0.00 

First we combine both mean for that dataset, to store them together as a
new object called `meancounts`.

``` r
meancounts <- data.frame(control.mean, treated.mean)
head(meancounts)
```

                    control.mean treated.mean
    ENSG00000000003       900.75       658.00
    ENSG00000000005         0.00         0.00
    ENSG00000000419       520.50       546.00
    ENSG00000000457       339.75       316.50
    ENSG00000000460        97.25        78.75
    ENSG00000000938         0.75         0.00

> Q5 (a). Create a scatter plot showing the mean of the treated samples
> against the mean of the control samples. Your plot should look
> something like the following.

Now we make the plot of the `control.mean` vs `treated.mean`

``` r
plot(meancounts[,1],meancounts[,2], xlab="Control", ylab="Treated")
```

![](class13_files/figure-commonmark/unnamed-chunk-13-1.png)

> Q5 (b).You could also use the ggplot2 package to make this figure
> producing the plot below. What geom\_?() function would you use for
> this plot?

``` r
library(ggplot2)

ggplot(meancounts, aes( x= control.mean, y=treated.mean)) + 
  geom_point(alpha = 0.4) + 
  labs(x = "Control", y = "Treated")
```

![](class13_files/figure-commonmark/unnamed-chunk-14-1.png)

> Q6. Try plotting both axes on a log scale. What is the argument to
> plot() that allows you to do this?

Our counts data is highly skewed where a huge amount of data falls in
one side and very few on the other side, so we see a apttern like this
plot it SCREEMS log transfrom this.

``` r
plot(meancounts, log="xy")
```

    Warning in xy.coords(x, y, xlabel, ylabel, log): 15032 x values <= 0 omitted
    from logarithmic plot

    Warning in xy.coords(x, y, xlabel, ylabel, log): 15281 y values <= 0 omitted
    from logarithmic plot

![](class13_files/figure-commonmark/unnamed-chunk-15-1.png)

We most often use log2 transform for this kind of data cause it has an
easier interpretation in bioinformatics area.

``` r
## Treated / Control
log2(20/20)
```

    [1] 0

``` r
log2(10/20)
```

    [1] -1

``` r
log2(40/20)
```

    [1] 1

``` r
log2(80/20)
```

    [1] 2

We call this fraction the “log2 fold change” as it tells us how much
more or less gene expression we have in units of doubling, etc.

Let’s calculate the log2 fold change for our `treated.mean` vs
`control.mean` counts adn call this `log2fc`

``` r
meancounts$log2fc <- log2(meancounts$treated.mean / meancounts$control.mean)
head(meancounts)
```

                    control.mean treated.mean      log2fc
    ENSG00000000003       900.75       658.00 -0.45303916
    ENSG00000000005         0.00         0.00         NaN
    ENSG00000000419       520.50       546.00  0.06900279
    ENSG00000000457       339.75       316.50 -0.10226805
    ENSG00000000460        97.25        78.75 -0.30441833
    ENSG00000000938         0.75         0.00        -Inf

A common “rule of thumb” threshold for calling a gene “up regulated” or
“down regulated” is a log2 fold change value of +2 or -2 (or greater).

``` r
zero.vals <- which(meancounts[,1:2]==0, arr.ind=TRUE)

to.rm <- unique(zero.vals[,1])
mycounts <- meancounts[-to.rm,]
head(mycounts)
```

                    control.mean treated.mean      log2fc
    ENSG00000000003       900.75       658.00 -0.45303916
    ENSG00000000419       520.50       546.00  0.06900279
    ENSG00000000457       339.75       316.50 -0.10226805
    ENSG00000000460        97.25        78.75 -0.30441833
    ENSG00000000971      5219.00      6687.50  0.35769358
    ENSG00000001036      2327.00      1785.75 -0.38194109

> Q7. What is the purpose of the arr.ind argument in the which()
> function call above? Why would we then take the first column of the
> output and need to call the unique() function?

The arr.ind=TRUE argument makes which() return the row and column
positions of the zero values as a matrix, instead of just returning a
single index number.

We then take the first column because it contains the row numbers where
zeros occur. Since a row could contain more than one zero, the same row
number might appear multiple times, so unique() is used to keep only one
copy of each row number before removing those rows from the dataset.

> Q8. Using the up.ind vector above can you determine how many up
> regulated genes we have at the greater than 2 fc level?

``` r
up.ind <- mycounts$log2fc > 2
sum(up.ind)
```

    [1] 250

There are 250 up regulated gene that are greater than 2.

> Q9. Using the down.ind vector above can you determine how many down
> regulated genes we have at the greater than 2 fc level?

``` r
down.ind <- mycounts$log2fc < (-2)
sum(down.ind)
```

    [1] 367

There are 367 down regulated genes

> Q10. Do you trust these results? Why or why not?

Our analysis so far has only focused on fold change values. However, a
gene can show a large increase or decrease in expression without the
result actually being statistically significant based on p-values.
Because of this, the current results could be misleading.

## DESeq Analysis

Let’s do thi analyze properly and not forget about the signifcant of the
difference.

For this we will use the **DESeq2** package.

``` r
library(DESeq2)
```

To run a DESeq analysis we need at least two inputs

- `countData` i.e. are gene counts across different experiments.
- `colData` i.e. our metadata about those count columns.

``` r
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~dex)
```

    converting counts to integer mode

``` r
dds
```

    class: DESeqDataSet 
    dim: 38694 8 
    metadata(1): version
    assays(1): counts
    rownames(38694): ENSG00000000003 ENSG00000000005 ... ENSG00000283120
      ENSG00000283123
    rowData names(0):
    colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
    colData names(4): id dex celltype geo_id

Now we can run the DESeq analysis pipeline using the `dds` object that
has all the inputs we need.

``` r
dds <- DESeq(dds)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

``` r
res <- results(dds)
head(res)
```

    log2 fold change (MLE): dex treated vs control 
    Wald test p-value: dex treated vs control 
    DataFrame with 6 rows and 6 columns
                      baseMean log2FoldChange     lfcSE      stat    pvalue
                     <numeric>      <numeric> <numeric> <numeric> <numeric>
    ENSG00000000003 747.194195      -0.350703  0.168242 -2.084514 0.0371134
    ENSG00000000005   0.000000             NA        NA        NA        NA
    ENSG00000000419 520.134160       0.206107  0.101042  2.039828 0.0413675
    ENSG00000000457 322.664844       0.024527  0.145134  0.168996 0.8658000
    ENSG00000000460  87.682625      -0.147143  0.256995 -0.572550 0.5669497
    ENSG00000000938   0.319167      -1.732289  3.493601 -0.495846 0.6200029
                         padj
                    <numeric>
    ENSG00000000003  0.163017
    ENSG00000000005        NA
    ENSG00000000419  0.175937
    ENSG00000000457  0.961682
    ENSG00000000460  0.815805
    ENSG00000000938        NA

## Volcano Plot

This is ubiquitous and common visualizaiton for this type of data that
puts the log2 fold change and the adjusted p-value together in one plot
that people can get insight for what is going on in the whole dataset
results.

``` r
library(ggplot2)
```

``` r
ggplot(res) +
  aes(log2FoldChange, padj) + 
  geom_point(alpha = 0.4)
```

    Warning: Removed 23549 rows containing missing values or values outside the scale range
    (`geom_point()`).

![](class13_files/figure-commonmark/unnamed-chunk-25-1.png)

That plot is not very useful because we don’t need graph with very high
p-values, we want the very low values below our alpha threshold
(e.g. 0.01).

Let’s log the y axis so we can see thses genes/points more clearily.

``` r
ggplot(res) +
  aes(log2FoldChange, log(padj)) + 
  geom_point(alpha = 0.4)
```

    Warning: Removed 23549 rows containing missing values or values outside the scale range
    (`geom_point()`).

![](class13_files/figure-commonmark/unnamed-chunk-26-1.png)

We need to flip the y-axis so out Volcano is not upside down.

``` r
ggplot(res) +
  aes(log2FoldChange, -log(padj)) + 
  geom_point()
```

    Warning: Removed 23549 rows containing missing values or values outside the scale range
    (`geom_point()`).

![](class13_files/figure-commonmark/unnamed-chunk-27-1.png)

Now to make the two side devided we can draw lines.

``` r
# Setup custom colors
mycols <- rep("gray", nrow(res))
mycols[abs(res$log2FoldChange) > 2] <- "red"

inds <- (res$padj < 0.01) & (abs(res$log2FoldChange) > 2)
mycols[inds] <- "blue"

# Volcano plot with ggplot
ggplot(res) +
  aes(x = log2FoldChange, y = -log(padj)) +
  geom_point(color = mycols, alpha= 0.4) +
  geom_vline(xintercept = c(-2, 2),
             color = "gray", linetype = "dashed") +
  geom_hline(yintercept = -log(0.1),
             color = "gray", linetype = "dashed") +
  xlab("Log2(FoldChange)") +
  ylab("-Log(P-value)") +
  theme_classic()
```

    Warning: Removed 23549 rows containing missing values or values outside the scale range
    (`geom_point()`).

![](class13_files/figure-commonmark/unnamed-chunk-28-1.png)

> Add annotaion to this volcano plot including the log2 Fold Change
> threshold of +2 and -2 and the p-value threshold of 0.05. Also color
> up just the genes that meet both these thresholds. these are the ones
> we will focus on next day!

## Save out results to date

``` r
write.csv(res, file = "myresults.csv")
```

## Adding Annotaion Data

We need to “map” or translate our ENSEMBLE gene identifiers in our
results object the date to the identifiers used in the different
databases we want to learn from.

For this we will use a couple of BioCoductor packages, with
`BiocManager::install("AnnotationDbi")` and
`BiocManager::install("org.Hs.eg.db")`

``` r
library(AnnotationDbi)
library(org.Hs.eg.db)
```

We can see the columns in `org.Hs.eg.db` that list the different
databases

``` r
columns(org.Hs.eg.db)
```

     [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
     [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
    [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"         
    [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
    [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
    [26] "UNIPROT"     

> Q11. Run the mapIds() function two more times to add the Entrez ID and
> UniProt accession and GENENAME as new columns called
> res$entrez, res$uniprot and res\$genename.

We can use the `mapIDs()` function to map between these different
database identifiers formats.

``` r
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res), 
                     keytype="ENSEMBL",       
                     column="SYMBOL")
```

    'select()' returned 1:many mapping between keys and columns

> Q. Can you map to “GENENAME” and add asa new col to our `res` object?

``` r
res$genename <- mapIds(org.Hs.eg.db,
                     keys=row.names(res), 
                     keytype="ENSEMBL",       
                     column="GENENAME")
```

    'select()' returned 1:many mapping between keys and columns

> Q Add “ENTERZID” as `res$entrez`

``` r
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res), 
                     keytype="ENSEMBL",       
                     column="ENTREZID")
```

    'select()' returned 1:many mapping between keys and columns

``` r
res$uniprot <- mapIds(org.Hs.eg.db,
                     keys=row.names(res), 
                     keytype="ENSEMBL",       
                     column="UNIPROT")
```

    'select()' returned 1:many mapping between keys and columns

``` r
head(res)
```

    log2 fold change (MLE): dex treated vs control 
    Wald test p-value: dex treated vs control 
    DataFrame with 6 rows and 10 columns
                      baseMean log2FoldChange     lfcSE      stat    pvalue
                     <numeric>      <numeric> <numeric> <numeric> <numeric>
    ENSG00000000003 747.194195      -0.350703  0.168242 -2.084514 0.0371134
    ENSG00000000005   0.000000             NA        NA        NA        NA
    ENSG00000000419 520.134160       0.206107  0.101042  2.039828 0.0413675
    ENSG00000000457 322.664844       0.024527  0.145134  0.168996 0.8658000
    ENSG00000000460  87.682625      -0.147143  0.256995 -0.572550 0.5669497
    ENSG00000000938   0.319167      -1.732289  3.493601 -0.495846 0.6200029
                         padj      symbol               genename      entrez
                    <numeric> <character>            <character> <character>
    ENSG00000000003  0.163017      TSPAN6          tetraspanin 6        7105
    ENSG00000000005        NA        TNMD            tenomodulin       64102
    ENSG00000000419  0.175937        DPM1 dolichyl-phosphate m..        8813
    ENSG00000000457  0.961682       SCYL3 SCY1 like pseudokina..       57147
    ENSG00000000460  0.815805       FIRRM FIGNL1 interacting r..       55732
    ENSG00000000938        NA         FGR FGR proto-oncogene, ..        2268
                        uniprot
                    <character>
    ENSG00000000003  A0A087WYV6
    ENSG00000000005      Q9H2S6
    ENSG00000000419      H0Y368
    ENSG00000000457      X6RHX1
    ENSG00000000460      A6NFP1
    ENSG00000000938      B7Z6W7

## Pathway Analysis

Now we have our annotated results with. their log2 fold-change and
p-values we can figure out which biologoical pathways and process these
genes are involved with.

We will use the **gage** and **pathview** packages for this step and we
can install them with,
`BiocManager::install( c("pathview", "gage", "gageData") )`.

``` r
library(gage)
library(gageData)
library(pathview)
```

Let’s have a peek at gageData,

``` r
data(kegg.sets.hs)

# Examine the first 2 pathways in this kegg set for humans
head(kegg.sets.hs, 2)
```

    $`hsa00232 Caffeine metabolism`
    [1] "10"   "1544" "1548" "1549" "1553" "7498" "9"   

    $`hsa00983 Drug metabolism - other enzymes`
     [1] "10"     "1066"   "10720"  "10941"  "151531" "1548"   "1549"   "1551"  
     [9] "1553"   "1576"   "1577"   "1806"   "1807"   "1890"   "221223" "2990"  
    [17] "3251"   "3614"   "3615"   "3704"   "51733"  "54490"  "54575"  "54576" 
    [25] "54577"  "54578"  "54579"  "54600"  "54657"  "54658"  "54659"  "54963" 
    [33] "574537" "64816"  "7083"   "7084"   "7172"   "7363"   "7364"   "7365"  
    [41] "7366"   "7367"   "7371"   "7372"   "7378"   "7498"   "79799"  "83549" 
    [49] "8824"   "8833"   "9"      "978"   

We need a named vector of importance (e.g. fold-change values) that has
gene ids as names. These name need to be in the correct format (using
the correct database format for the IDs)

``` r
x <- c(10, 9, 7)
names(x) <- c("Alice", "Maddy", "Barry")
x
```

    Alice Maddy Barry 
       10     9     7 

``` r
names(x)
```

    [1] "Alice" "Maddy" "Barry"

Here we will make a wee input vector called `foldchanges` that has
“entrez” ids as names

``` r
foldchanges <- res$log2FoldChange
names(foldchanges) <- res$entrez
```

Now we can run `gage()` to do our pathway analysis.

``` r
# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs)
```

``` r
attributes(keggres)
```

    $names
    [1] "greater" "less"    "stats"  

The top 3 overlapping pathways from KEGG

``` r
head(keggres$less, 3)
```

                                          p.geomean stat.mean        p.val
    hsa05332 Graft-versus-host disease 0.0004250607 -3.473335 0.0004250607
    hsa04940 Type I diabetes mellitus  0.0017820379 -3.002350 0.0017820379
    hsa05310 Asthma                    0.0020046180 -3.009045 0.0020046180
                                            q.val set.size         exp1
    hsa05332 Graft-versus-host disease 0.09053792       40 0.0004250607
    hsa04940 Type I diabetes mellitus  0.14232788       42 0.0017820379
    hsa05310 Asthma                    0.14232788       29 0.0020046180

Now we can use the **pathview** package with the found KEGG pathway IDs
(e.g. “hsa05310” for the Asthma pathway) to make a pathway figure
showing our Differential Expressed Genes (DEGs)

``` r
pathview(gene.data = foldchanges, pathway.id="hsa05310")
```

    Info: Downloading xml files for hsa05310, 1/1 pathways..

    Info: Downloading png files for hsa05310, 1/1 pathways..

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/madinakho/Desktop/UCSD/Classes/Spring 26/BIMM 143/Class R/Lab 15/bimm143_github/Class13

    Info: Writing image file hsa05310.pathview.png

![](hsa05310.pathway.png) \> Q12. Can you do the same procedure as above
to plot the pathview figures for the top 2 down-regulated pathways?

## Save our annotated results

``` r
write.csv(res, file = "myresults_annotated.csv")
```
