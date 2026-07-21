# Real-data inputs

Raw third-party datasets are not redistributed in this repository. The formal
real-data scripts expect the following local layout under the directory supplied
with `--data-root`:

```text
data-root/
  pbmcs/
    filtered_feature_bc_matrix/
      matrix.mtx
      features.tsv
      barcodes.tsv
    clusters.csv
  pancreas_GSM2230759_human3_umifm_counts.csv
  bhattacherjee/
    GSE124952_expression_matrix.csv.gz
    GSE124952_meta_data.csv.gz
```

## Provenance summary

| Dataset | Official source | Cells | Reference labels used |
|---|---|---:|---|
| PBMCs | [10x `pbmc_1k_protein_v3`](https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_web_summary.html), Cell Ranger 3.0.0 | 713 | `analysis/clustering/kmeans_6_clusters/clusters.csv`; automated six-cluster assignments |
| Pancreas | [GSM2230759](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2230759) from GSE84133 | 3,605 | `assigned_cluster` (14 classes) |
| Bhattacherjee | [GSE124952](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124952) processed expression matrix and metadata | 35,360 | `CellType` (8 classes) |

Official download files:

- PBMCs: [filtered feature-barcode matrix](https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.tar.gz) and [Cell Ranger analysis archive](https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_analysis.tar.gz).
- Pancreas: [GSM2230759 processed UMI table](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2230759&file=GSM2230759_human3_umifm_counts.csv.gz&format=file).
- Bhattacherjee: [expression matrix](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE124952&file=GSE124952_expression_matrix.csv.gz&format=file) and [metadata](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE124952&file=GSE124952_meta_data.csv.gz&format=file).

## Verified input checksums

All local manuscript inputs were verified byte-for-byte against the official
downloads on 12 July 2026. The Pancreas input is the decompressed GEO file.

| Dataset | Local input | MD5 |
|---|---|---|
| PBMCs | `filtered_feature_bc_matrix/matrix.mtx` | `472d822337ac3a11a6a2f11200e7c3da` |
| PBMCs | `filtered_feature_bc_matrix/features.tsv` | `1e0949c8bc793752a476b074bab793f0` |
| PBMCs | `filtered_feature_bc_matrix/barcodes.tsv` | `745f599c9dce159437e7022969bf9505` |
| PBMCs | `clusters.csv` | `63190cfe9972f90413b833186e4ca68e` |
| Pancreas | `pancreas_GSM2230759_human3_umifm_counts.csv` | `fee29f30b6a6ab121f027e8436667eaa` |
| Bhattacherjee | `GSE124952_expression_matrix.csv.gz` | `81de6c10e73e41e06ec91128bc2c16f7` |
| Bhattacherjee | `GSE124952_meta_data.csv.gz` | `cde3dc48b47926ba753f0204c973a290` |

The Bhattacherjee analysis uses all 35,360 matched cells in the official GEO
processed files, not the 24,822-cell derivative reported in the SCENA
benchmark. Raw third-party data are not redistributed in this repository.

The same checksum records are provided in `DATA_MD5SUMS.txt`. Run the check from
the repository root after preparing the files. The Pancreas checksum applies to
the decompressed `pancreas_GSM2230759_human3_umifm_counts.csv` file.
