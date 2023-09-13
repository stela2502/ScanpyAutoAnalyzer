```R
library(SingleR)
library(scuttle)
library(zellkonverter)

system('ls /opt/celldex/dataset/', intern=T)
```

```R
ifile = "Raw.h5ad"
load("/opt/celldex/dataset/ImmGenData.RData")
ls()
```

```R
ownDataC = as_cellexalvrR( ifile, specie='mouse')
ownDataC
```

```R
ownData = 
SingleCellExperiment (
    list(counts=ownDataC@data ),
    colData=DataFrame( ownDataC@usedObj$samples),
    rowData = DataFrame( ownDataC@meta.gene)
)       
ownData
```

```R
ownData <- logNormCounts(ownData, )
pred.ownData <- SingleR(test=ownData, ref=hpca.se, labels=hpca.se$label.main, de.method="wilcox")
table(pred.ownData$labels)
```

```R
write.table( pred.ownData, file="SingleR_predicted_Raw.csv", sep="\t", quote=F )
```