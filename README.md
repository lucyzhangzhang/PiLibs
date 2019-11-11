# Processing Root/Shoot Pi treatment libraries

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
Trim and pair reads together (`trim.sh`) as well as mapping them with STAR (`STAR.sh`), the read counting with QoRTs (`QoRTs.sh`)

Meta analysis done on the properties of reads, i.e. % trimmed, % uniquely mapped (`meta.sh`, `metaAnalysis.R`)

Differential gene expression analysis and comparison with the Pi/S libraries (`DEGs.R`)
