# Processing Root/Shoot Pi treatment libraries
Analysis of Eutrema root and shoot RNA-seq libraries under Pi deficiency treatments, these libraries are from Yukon cabinet grown plants of single-seed descent from a Yukon field plant. Fertilization occurs weekly for 4 weeks with the respective Pi treatments, then the roots and shoots were harvested at week 4 for RNA-seq analysis. 

## Raw read location
```cd PiLibs```

## Trimmed reads
```cd trimOut```

## STAR outputs
``` cd STAR```

## QoRTs read counts
```cd QoRTs```

## Scipts
```cd scripts```

## Workflow
* Trim and pair reads together (`trim.sh`) as well as mapping them with STAR (`mapping.sh`), the read counting with QoRTs (`QoRTs.sh`)

* Meta analysis done on the properties of reads, i.e. % trimmed, % uniquely mapped (`meta.sh`, `metaAnalysis.R`)

* Retain reasonably uniquely mapped reads, i.e. mapQ ≥ 10

* Differential gene expression analysis and comparison with the Pi/S libraries (`DEGs.R`)
