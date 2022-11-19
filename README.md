# DGE-with-microarray-data

The aim of this study is to find the role of BRD4 isoforms in both breast and ovarian cancer. For that we performed differential gene expression analysis with microarray data downloaded from GEO database (https://www.ncbi.nlm.nih.gov/geo/). 

Accession numbers for breast cancer data - GSE29431, GSE42568, and GSE61304 

| Accession       | # tumor tissues | # non-tumor tissues |
| --------------- | --------------- | ------------------- |
| GSE29431        |       54        |         54          |
| GSE42568        |      104        |         17          |
| GSE61304        |       58        |          4          |

Accession numbers for ovarian cancer data - GSE14407, GSE36668, and GSE38666

| Accession       | # tumor tissues | # non-tumor tissues |
| --------------- | --------------- | ------------------- |
| GSE14407        |       12        |         12          |
| GSE36668        |        4        |          4          |
| GSE38666        |       18        |         12          |


All datasets are in CEL file format created by Affymetrix DNA microarray image analysis software. They contain information about the intensity values taken from probes on an Affymetrix GeneChip. 

.CEL files were analyzed by the Affy package in Bioconductor using the Robust Multiarray Average (RMA) algorithm. 

Used the Combat algorithm from the Bioconductor’s Surrogate Variable Analysis (SVA) package for batch normalization and performed class-specific quantile normalization.

Limma package used for DEG identification.

Converted the probe sets to relevant genes were using the Annotation package and hgu133plus2.db (Affymetrix HG-U133 Plus 2 Array) annotation data package in Bioconductor.

Note: the CEL files of the relevant cancer GEO datasets should be separated into tumor and non-tumor groups inorder to excute the codes.

Order to run:

<ol>
  <li>read_CELFiles.r</li>
  <li>batchEffect_quantileNorm.r</li>
  <li>pca.r (to find the patterns of the datapoints - optional)</li>
  <li>limma_DGE_analysis.r</li>
  <li>volcano_plot.r (to visualize limma output)</li>
  <li>probe_annotation.r</li>
</ol>

