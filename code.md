### 数据准备

#### 创建环境

使用的工具：
- [parallel-fastq-dump](https://github.com/rvalieris/parallel-fastq-dump)：`conda create -n fastq-download parallel-fastq-dump 'sra-tools>=3.0.0' -y`
- [STAR](https://github.com/alexdobin/STAR)：`conda create -n STAR bioconda::star -y`，我使用的版本为2.7.11b
- [stellarscope](https://github.com/nixonlab/stellarscope)：`conda create -n stellarscope stellarscope samtools -y`，我使用的版本为1.5
- [scTE](https://github.com/JiekaiLab/scTE)：`conda create -n scTE python scte samtools -y`，我使用的版本为1.0.0
- [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks)：`conda create -n pygenometracks pygenometracks -y`，我使用的版本为3.9

R环境：本次分析使用的包除了Seurat以外，未见明显的版本差异（因版本冲突而导致的报错）
- 服务器（构建Seurat对象/数据预处理）：4.2.0
  - Seurat：5.3.0
  - glmGamPoi：1.8.0
  - MAST：1.22.0
  - harmony：1.2.4
  - SingleCellExperiment：1.18.0
  - scDblFinder：1.10.0
  - tidyverse：2.0.0
  - readxl：1.4.2
  - ggpubr：0.4.0
  - ggrepel：0.9.1
  - ggVennDiagram：1.5.7
  - ggbreak：0.1.7
  - ggplotify：0.1.3
  - patchwork：1.3.0
  - clustree：0.5.1
  - pheatmap：1.0.13
  - ComplexHeatmap：2.14.0
  - circlize：0.4.17
  - Matrix：1.6-5
  - data.table：1.14.10
  - WGCNA：1.74
  - hdWGCNA：0.4.11
- 本地（具体细胞类型内分析）：4.4.1
  - Seurat：5.1.0
  - glmGamPoi：1.16.0
  - MAST：1.30.0
  - SingleCellExperiment：1.26.0
  - tidyverse：2.0.0
  - tidytext：0.4.3
  - readxl：1.4.3
  - ggpubr：0.6.0
  - ggrepel：0.9.5
  - patchwork：1.2.0
  - cowplot：1.1.3
  - magick：2.9.1
  - pdftools：3.7.0
  - grid：4.4.1
  - clustree：0.5.1
  - pheatmap：1.0.12
  - enrichplot：1.24.4
  - Matrix：1.7-0
  - data.table：1.14.8
  - clusterProfiler：4.12.6
  - org.Hs.eg.db：3.19.1
  - WGCNA：1.72-5
  - hdWGCNA：0.4.11

#### 数据下载

以下几个文件如果下载后是gz的，都用`gunzip`命令解压
- 人类基因组[GRCh38.p14.genome.fa.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.p14.genome.fa.gz)
- 人类基因组注释[gencode.v49.annotation.gtf.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz)
- stellarscope使用的hERV专用注释[transcripts.gtf](https://media.githubusercontent.com/media/mlbendall/telescope_annotation_db/refs/heads/master/builds/HERV_rmsk.hg38.v2/transcripts.gtf?download=true)（来自stellarscope作者）
- scTE使用的RepeatMasker注释[rmsk.txt.gz](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz)
- [10X Genomics的whitelist文件](https://github.com/noamteyssier/10x_whitelist_mirror)：包含多个版本的whitelist，我这里用的两个数据用的都是[3M-february-2018.txt](https://github.com/noamteyssier/10x_whitelist_mirror/blob/main/3M-february-2018.txt.gz)
- 其它所需的metadata文件（包括两个数据集的SRR+GSM编号、样本信息）均在我的GitHub仓库[hERV_AD_snRNA-seq_analyze](https://github.com/lwstkhyl/hERV_AD_snRNA-seq_analyze)中获取
  - `SRR_Acc_List.txt`：从ncbi上下载的该数据集包含的所有SRR号
  - `srr2gsm.csv`：SRR与GSM编号的对应关系
  - `family2group.xlsx`：根据stellarscope作者提供的hERV的subfamily-group映射[families.tsv](https://github.com/mlbendall/telescope_annotation_db/blob/master/builds/HERV_rmsk.hg38.v2/families.tsv)整理得到
  - `LTR2HERV.tsv`：根据families.tsv整理得到，看stellarscope作者把哪些LTR归入到了哪个HERV中
  - `LTR2HERV_overlap.tsv`：LTR2HERV.tsv中被归入多个HERV的LTR
  - `LTR2Class.tsv`：根据rmsk.txt，把repName映射到repFamily（ERVL/ERVK/ERV1）
  - 其它文件都是样本信息：从SRA Run Selector上下载的metadata和论文中附带的表格进行精简而来

**构建参考基因组**：

```sh
module load miniconda3/base
# STAR
conda activate STAR
mkdir -p /public/home/wangtianhao/Desktop/STAR_ref/
cd /public/home/wangtianhao/Desktop/STAR_ref/
STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir hg38 \
  --genomeFastaFiles GRCh38.p14.genome.fa \
  --sjdbGTFfile gencode.v49.annotation.gtf
# scTE
conda activate scTE
mkdir -p /public/home/wangtianhao/Desktop/scTE_ref
cd /public/home/wangtianhao/Desktop/scTE_ref
scTE_build -g hg38  # 需要联网
```

**fastq原始reads下载**：不需要解压
- [GSE157827](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157827)：可以从上面我的GitHub仓库中直接下载`SRR_Acc_List.txt`，也可以从[SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA662923)上下载，下同
  - 共92个，直接全选，下载Metadata和Accession List
- [GSE174367](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367)：[SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA729525)
  - Assay Type选`rna-seq`
  - isolation选`unbiased total nuclei isolation`
  - 共152个，全选，在Selected栏里面下载Metadata和Accession List

```sh
module load miniconda3/base
conda activate fastq-download
cd /public/home/wangtianhao/Desktop/GSE157827/fastq  # 用于存放原始fastq文件的文件夹，把SRR_Acc_List.txt放到这里面
mkdir fastq_temp
for SRR in $(cat SRR_Acc_List.txt); do
    parallel-fastq-dump --sra-id ${SRR} --tmpdir ./fastq_temp --threads 8 --split-files --gzip
done
# cd /public/home/wangtianhao/Desktop/GSE174367/fastq
# mkdir fastq_temp
# for SRR in $(cat SRR_Acc_List.txt); do
#     parallel-fastq-dump --sra-id ${SRR} --tmpdir ./fastq_temp --threads 8 --split-files --gzip
# done
```
- `parallel-fastq-dump`这个下载工具不是很稳定，有时会发生下载失败的情况（可能发生“即使下载失败，fastq.gz文件也存在，只是大小不对”的情况），务必检查log输出，重新下载中间报错的SRR
- 如果频繁报错，可尝试更换较稳定的SRR下载工具
- 如果想放在后台执行（使用`nohup`命令），就把上面的代码放到写到一个`.sh`文件中，然后`nohup sh download.sh &`，不要在`parallel-fastq-dump`前加`nohup`，这样会同时开很多个`parallel-fastq-dump`容易炸服务器
- GSE157827数据集下载后约1.1T，GSE174367约733G

最后工作目录下应有两个文件夹`GSE157827`和`GSE174367`，每个文件夹内都有一个fastq文件夹保存`SRRxxx.fastq.gz`，为了方便运行后续代码，我把`srr2gsm.csv`也放到了fastq文件夹内

#### 比对和计数

##### STARsolo

- GSE157827：`_2`是cDNA，`_1`是barcode+UMI（16位bc+10位UMI）
- GSE174367：`_3`是cDNA，`_2`是barcode+UMI（16位bc+12位UMI）
  - 上面两个样本在运行STARsolo时，唯一差别是`--soloUMIlen`参数一个设置为10，另一个为12
- GEO数据库里面说用的是10X Genomics v3，所以whitelist用`3M-february-2018.txt`（其实用v2的那个差别不是很大）
- 每个样本（GSM）大约需要1-2h，推荐16线程+256G内存

```sh
workDir=/public/home/GENE_proc/wth/GSE157827/  # 工作目录
# workDir=/public/home/GENE_proc/wth/GSE174367/
genomeDir=/public/home/wangtianhao/Desktop/STAR_ref/hg38/  # STAR索引
fqDir=${workDir}/fastq  # fastq.gz
mapFile=${fqDir}/srr2gsm.csv  # SRR-GSM映射
whitelist=/public/home/wangtianhao/Desktop/STAR_ref/whitelist/3M-february-2018.txt  # whitelist文件
res_barcodes=barcodes  # 计数结果命名
res_features=features
res_counts=counts
module load miniconda3/base
conda activate STAR
cd ${workDir}
dos2unix ${mapFile}
GSM_LIST=$(cut -d',' -f2 "${mapFile}" | sort -u)
for GSM in ${GSM_LIST}; do
    GSM=${GSM//$'\r'/}
    echo "=============================="
    echo "[GSM] ${GSM}"
    SRR_LIST=$(awk -F',' -v g="${GSM}" '$2==g{print $1}' "${mapFile}")
    nSRR=$(echo "${SRR_LIST}" | wc -l | awk '{print $1}')
    echo "  SRR count: ${nSRR}"
    echo "${SRR_LIST}" | sed 's/^/  - /'
    cdna_files=$(echo "${SRR_LIST}" | sed "s|^|${fqDir}/|; s|$|_2.fastq.gz|" | paste -sd, -)
    bcumi_files=$(echo "${SRR_LIST}" | sed "s|^|${fqDir}/|; s|$|_1.fastq.gz|" | paste -sd, -)
    for f in $(echo "${cdna_files}" | tr ',' ' ') $(echo "${bcumi_files}" | tr ',' ' '); do
        if [[ ! -f "${f}" ]]; then
        echo "[ERROR] Missing fastq: ${f}" >&2
        exit 1
        fi
    done
    mkdir -p star
    STAR \
        --runMode alignReads \
        --runThreadN 16 \
        --genomeDir ${genomeDir} \
        --readFilesIn ${cdna_files} ${bcumi_files} \
        --readFilesCommand zcat \
        --outFileNamePrefix star/${GSM} \
        --soloType CB_UMI_Simple \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 10 \
        # --soloUMIlen 12 \
        --soloBarcodeReadLength 0 \
        --soloCBwhitelist ${whitelist} \
        --soloFeatures Gene GeneFull \
        --clipAdapterType CellRanger4 \
        --soloCellFilter EmptyDrops_CR \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloUMIdedup 1MM_CR \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
        --outSAMunmapped Within \
        --limitOutSJcollapsed 5000000 \
        --outFilterScoreMin 30 \
        --outFilterMultimapNmax 500 \
        --outFilterMultimapScoreRange 5
    mkdir -p mtx/${GSM}/gene
    cp star/${GSM}Solo.out/GeneFull/filtered/barcodes.tsv mtx/${GSM}/gene/${res_barcodes}.tsv
    cp star/${GSM}Solo.out/GeneFull/filtered/features.tsv mtx/${GSM}/gene/${res_features}.tsv
    cp star/${GSM}Solo.out/GeneFull/filtered/matrix.mtx mtx/${GSM}/gene/${res_counts}.mtx
done
```

##### stellarscope

- 使用`best_average`模式：如果有多个最佳比对就取平均值，默认`best_exclude`是直接丢弃，因为想让hERV计数多一些，所以改成用平均值
- 每个样本需要4-13h不等（依据bam文件大小不同，序列越多需要时间越长），推荐16线程+128G内存
- 这两组数据的代码相同

```sh
workDir=/public/home/GENE_proc/wth/GSE157827/
# workDir=/public/home/GENE_proc/wth/GSE174367/
fqDir=${workDir}/fastq
mapFile=${fqDir}/srr2gsm.csv
hERV_gtf=/public/home/wangtianhao/Desktop/STAR_ref/transcripts.gtf  # hERV的gtf文件
res_barcodes=barcodes
res_features=features
res_counts=counts
module load miniconda3/base
conda activate stellarscope
cd ${workDir}
GSM_LIST=$(cut -d',' -f2 "${mapFile}" | sort -u)
for GSM in ${GSM_LIST}; do
    GSM=${GSM//$'\r'/}
    echo "=============================="
    echo "[GSM] ${GSM}"
    mkdir -p stellarscope/${GSM}
    mkdir -p tmp/${GSM}
    samtools view -@4 -u -F 4 \
        -D CB:<(tail -n+1 star/${GSM}Solo.out/GeneFull/filtered/barcodes.tsv) \
        star/${GSM}Aligned.sortedByCoord.out.bam \
        | samtools sort -@16 -n -t CB -T ./tmp/${GSM} \
        > stellarscope/${GSM}/Aligned.sortedByCB.bam
    stellarscope assign \
        --outdir stellarscope/${GSM} \
        --nproc 16 \
        --stranded_mode F \
        --whitelist star/${GSM}Solo.out/GeneFull/filtered/barcodes.tsv \
        --pooling_mode individual \
        --reassign_mode best_average \
        --max_iter 500 \
        --updated_sam \
        stellarscope/${GSM}/Aligned.sortedByCB.bam \
        ${hERV_gtf}
    mkdir -p mtx/${GSM}/hERV
    cp stellarscope/${GSM}/stellarscope-barcodes.tsv mtx/${GSM}/hERV/${res_barcodes}.tsv
    cp stellarscope/${GSM}/stellarscope-features.tsv mtx/${GSM}/hERV/${res_features}.tsv
    cp stellarscope/${GSM}/stellarscope-TE_counts.mtx mtx/${GSM}/hERV/${res_counts}.mtx
done
```

##### scTE

```sh
workDir=/public/home/GENE_proc/wth/GSE157827/
# workDir=/public/home/GENE_proc/wth/GSE174367/
fqDir=${workDir}/fastq
mapFile=${fqDir}/srr2gsm.csv
scte_idx=/public/home/wangtianhao/Desktop/scTE_ref/hg38.exclusive.idx
module load miniconda3/base
conda activate scTE
GSM_LIST=$(cut -d',' -f2 "${mapFile}" | sort -u)
for GSM in ${GSM_LIST}; do
    GSM=${GSM//$'\r'/}
    echo "=============================="
    echo "[GSM] ${GSM}"
    mkdir -p scTE/${GSM}/
    samtools view -@4 -h -F 4 \
        -D CB:star/${GSM}Solo.out/GeneFull/filtered/barcodes.tsv \
        star/${GSM}Aligned.sortedByCoord.out.bam \
        | samtools view -@4 -b -o scTE/${GSM}/${GSM}.filtered_by_STARsolo.bam
    scTE \
        -i scTE/${GSM}/${GSM}.filtered_by_STARsolo.bam \
        -o ${GSM} \
        -x ${scte_idx} \
        -g hg38 \
        -CB True \
        -UMI True \
        --expect-cells $(wc -l < star/${GSM}Solo.out/GeneFull/filtered/barcodes.tsv) \
        --min_genes 0 \
        --hdf5 False \
        -p 8
    mkdir -p scTE_res
    mv ${GSM}.csv scTE_res
done
```

---

注：作者的snRNA-seq数据中虽然有GSM5292839这个样本，但在论文的补充材料中，却说这个样本并未用于snRNA分析，后来计数时发现这个样本的平均UMI和基因数确实偏低，所以这里我也直接删掉了这个样本的计数结果

```sh
rm -rf /public/home/GENE_proc/wth/GSE174367/mtx/GSM5292839/
rm -rf /public/home/GENE_proc/wth/GSE174367/scTE_res/GSM5292839.csv
```

### 构建Seurat对象和数据预处理

```r
# Seurat
library(Seurat)
library(glmGamPoi)
library(MAST)
# harmony
library(harmony)
# scDblFinder
library(SingleCellExperiment)
library(scDblFinder)
# tidyverse
library(plyr)
library(readxl)
library(tidyverse)
library(tidytext)
# 画图相关
library(patchwork)
library(cowplot)
library(grid)
library(scales)
library(clustree)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(enrichplot)
library(ggpubr)
library(ggrepel)
library(ggVennDiagram)
library(ggbreak)
library(ggplotify)
library(magick)
library(pdftools)
# 矩阵/df处理
library(Matrix)
library(data.table)
# GO
library(clusterProfiler)
library(org.Hs.eg.db)
# WGCNA
library(WGCNA)
library(hdWGCNA)
enableWGCNAThreads(nThreads = 8)  # hdWGCNA开启8线程，否则ConstructNetwork步骤可能需要1h以上
# 并行运算
library(future)
plan("multiprocess", workers = 4)  # 启用并行处理（FindIntegrationAnchors步骤可能需要数小时），但需考虑服务器兼容性（例如内存不足可能导致Rsession卡死）
options(future.globals.maxSize = 128*1024^3)  # 如果出现有关future的内存报错就执行这句（设置最大内存为128G）
# 结果存放目录
res_dir <- "/public/home/GENE_proc/wth/res/"
set.seed(520)
```

#### 构建Seurat对象

```r
# 读取计数矩阵
normalize_barcode <- function(x) {  # 如果已经有`-数字`后缀（比如`-1`）就保留，否则补`-1`
  ifelse(grepl("-\\d+$", x), x, paste0(x, "-1"))
}
read_one_gsm <- function(gsm, gsm_info, data_root) {
  gene_dir <- file.path(data_root, gsm, "gene")
  herv_dir <- file.path(data_root, gsm, "hERV")
  # 读取矩阵
  gene_counts <- ReadMtx(
    mtx = file.path(gene_dir, "counts.mtx"),
    features = file.path(gene_dir, "features.tsv"),
    cells = file.path(gene_dir, "barcodes.tsv"),
    unique.features = TRUE
  )
  herv_counts <- ReadMtx(
    mtx = file.path(herv_dir, "counts.mtx"),
    features = file.path(herv_dir, "features.tsv"),
    cells = file.path(herv_dir, "barcodes.tsv"),
    feature.column = 1,
    unique.features = TRUE
  )
  # 取该GSM的样本信息
  info <- gsm_info %>% filter(GSM == gsm)
  indv <- info$orig.ident
  # 使用原论文作者的细胞名格式
  gene_bc <- normalize_barcode(colnames(gene_counts))
  herv_bc <- normalize_barcode(colnames(herv_counts))
  colnames(gene_counts) <- paste(gsm, indv, gene_bc, sep = "_")
  colnames(herv_counts) <- paste(gsm, indv, herv_bc, sep = "_")
  # 对齐两套assay的细胞集合
  common_cells <- intersect(colnames(gene_counts), colnames(herv_counts))
  gene_counts <- gene_counts[, common_cells, drop = FALSE]
  herv_counts <- herv_counts[, common_cells, drop = FALSE]
  # 构建Seurat对象
  seu <- CreateSeuratObject(
    counts = gene_counts,
    assay = "RNA",
    project = "seu"
  )
  seu[["HERV"]] <- CreateAssayObject(counts = herv_counts)
  # 添加metadata
  seu$GSM <- gsm
  seu$orig.ident <- indv
  seu$group <- info$group
  seu$tissue <- info$tissue
  extra_cols <- setdiff(colnames(info), c("GSM", "group", "orig.ident", "tissue"))
  for (cc in extra_cols) {
    seu[[cc]] <- info[[cc]][1]
  }
  return(seu)
}
# 读取metadata-GSE157827
meta_map_GSE157827 <- read_xlsx("/public/home/GENE_proc/wth/metadata/GSE157827/metadata.xlsx")
sample_info_GSE157827 <- read_xlsx("/public/home/GENE_proc/wth/metadata/GSE157827/sample_info.xlsx")
names(meta_map_GSE157827) <- make.names(names(meta_map_GSE157827))
names(sample_info_GSE157827) <- make.names(names(sample_info_GSE157827))
gsm_meta_GSE157827 <- meta_map_GSE157827 %>%  # GSM级别信息
  transmute(
    GSM = as.character(sample_id),
    group = as.character(diagnosis),
    orig.ident = as.character(individuals),
    tissue = as.character(tissue)
  ) %>%
  distinct(GSM, .keep_all = TRUE)
donor_info_GSE157827 <- sample_info_GSE157827 %>%  # 样本信息
  mutate(orig.ident = as.character(sample_id)) %>%
  select(-sample_id)
gsm_info_GSE157827 <- gsm_meta_GSE157827 %>%  # 合并
  left_join(donor_info_GSE157827, by = "orig.ident") %>% 
  mutate(dataset="GSE157827")
# 读取metadata-GSE174367
meta_map_GSE174367 <- read_xlsx("/public/home/GENE_proc/wth/metadata/GSE174367/metadata.xlsx")
sample_info_GSE174367 <- read_xlsx("/public/home/GENE_proc/wth/metadata/GSE174367/sample_info.xlsx")
names(meta_map_GSE174367) <- make.names(names(meta_map_GSE174367))
names(sample_info_GSE174367) <- make.names(names(sample_info_GSE174367))
gsm_meta_GSE174367 <- meta_map_GSE174367 %>%
  transmute(
    GSM = as.character(sample_id),
    group = as.character(diagnosis),
    orig.ident = as.character(individuals),
    tissue = as.character(tissue)
  ) %>%
  distinct(GSM, .keep_all = TRUE)
donor_info_GSE174367 <- sample_info_GSE174367 %>%
  mutate(orig.ident = as.character(sample_id)) %>%
  select(-sample_id)
gsm_info_GSE174367 <- gsm_meta_GSE174367 %>%
  left_join(donor_info_GSE174367, by = "orig.ident") %>% 
  mutate(dataset="GSE174367")
# 合并metadata
gsm_info <- rbind.fill(gsm_info_GSE157827, gsm_info_GSE174367) %>% 
  select(-diagnosis)
gsm_info$orig.ident <- gsub("Sample-", "", gsm_info$orig.ident)
gsm_info$sample_uid <- paste(gsm_info$dataset, gsm_info$GSM, sep = "_")
gsm_info$donor_uid  <- paste(gsm_info$dataset, gsm_info$orig.ident, sep = "_")
write.csv(gsm_info, "/public/home/GENE_proc/wth/my_data/gsm_info.csv", row.names=F)
gsm_info[] <- lapply(gsm_info, function(x) {  # 将字符串NA转为真实NA
  x[as.character(x) %in% c("NA", "")] <- NA
  x
})mapFile
# 检查有无NA
any(is.na(gsm_info$orig.ident))  # FALSE
nrow(gsm_info)  # 40
# 构建Seurat对象
seu_list <- list()
data_root <- "/public/home/GENE_proc/wth/GSE157827/mtx/"
gsm_ids <- list.dirs(data_root, full.names = FALSE, recursive = FALSE)
for (gsm in gsm_ids) {
  tmp <- read_one_gsm(gsm, gsm_info = gsm_info, data_root = data_root)
  if (!is.null(tmp)) {
    seu_list[[gsm]] <- tmp
  }
}
data_root <- "/public/home/GENE_proc/wth/GSE174367/mtx/"
gsm_ids <- list.dirs(data_root, full.names = FALSE, recursive = FALSE)
for (gsm in gsm_ids) {
  tmp <- read_one_gsm(gsm, gsm_info = gsm_info, data_root = data_root)
  if (!is.null(tmp)) {
    seu_list[[gsm]] <- tmp
  }
}
# 合并每个样本
seu <- Reduce(function(x, y) merge(x, y), seu_list)
rm(seu_list)
seu <- JoinLayers(seu, assay = "RNA")
Layers(seu, assay = "RNA")  # 检查一下layers："counts"
# 线粒体/核糖体比例
seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent_mito")
seu <- PercentageFeatureSet(seu, pattern = "^RP[SL]", col.name = "percent_ribo")
seu$HERV_fraction <- seu$nCount_HERV / (seu$nCount_HERV + seu$nCount_RNA)
# 查看结果情况
seu
colnames(seu@meta.data)
table(seu$dataset)
table(seu$group)
table(seu$orig.ident)
seu$group <- factor(seu$group, levels = c("NC", "AD"))
# 保存数据
saveRDS(
  seu, 
  file = "/public/home/GENE_proc/wth/my_data/all.rds", 
  compress = "xz"
)
```

#### 质控

**过滤前画图**：

```r
seu <- readRDS("/public/home/GENE_proc/wth/my_data/all.rds")
p1 <- VlnPlot(seu, features = "nFeature_RNA",
              pt.size = 0, ncol = 1, group.by = "dataset") + NoLegend() + theme(axis.title.x = element_blank())
p1_1 <- VlnPlot(seu, features = "nCount_RNA",
              pt.size = 0, ncol = 1, group.by = "dataset") + NoLegend() + theme(axis.title.x = element_blank())
p2 <- VlnPlot(seu, features = "nCount_RNA",
              pt.size = 0, ncol = 1, group.by = "orig.ident") + NoLegend() + theme(axis.title.x = element_blank())
p3 <- VlnPlot(seu, features = "nFeature_HERV",
              pt.size = 0, ncol = 1, group.by = "dataset") + NoLegend() + theme(axis.title.x = element_blank())
p3_1 <- VlnPlot(seu, features = "nCount_HERV",
              pt.size = 0, ncol = 1, group.by = "dataset") + NoLegend() + theme(axis.title.x = element_blank())
p5 <- VlnPlot(seu, features = "percent_mito",
              pt.size = 0, ncol = 1, group.by = "dataset") + NoLegend() + theme(axis.title.x = element_blank())
p5_1 <- VlnPlot(seu, features = "percent_ribo",
              pt.size = 0, ncol = 1, group.by = "dataset") + NoLegend() + theme(axis.title.x = element_blank())
p6 <- VlnPlot(seu, features = "percent_mito",
              pt.size = 0, ncol = 1, group.by = "orig.ident") + NoLegend() + theme(axis.title.x = element_blank())
p7 <- VlnPlot(seu, features = "HERV_fraction",
              pt.size = 0, ncol = 1, group.by = "dataset") + NoLegend() + theme(axis.title.x = element_blank())
p8 <- VlnPlot(seu, features = "HERV_fraction",
              pt.size = 0, ncol = 1, group.by = "orig.ident") + NoLegend() + theme(axis.title.x = element_blank())
pdf(file.path(res_dir, "before_filter.pdf"), width = 24, height = 20)
p <- (p1 | p1_1 | p5 | p5_1 | p3 | p3_1 | p7) / (p2 | p6 | p8)
print(p)
dev.off()
fivenum(filter(seu@meta.data, dataset=="GSE174367")$nFeature_RNA)
fivenum(filter(seu@meta.data, dataset=="GSE174367")$nCount_RNA)
fivenum(filter(seu@meta.data, dataset=="GSE174367")$percent_mito)
fivenum(filter(seu@meta.data, dataset=="GSE174367")$percent_ribo)
fivenum(filter(seu@meta.data, dataset=="GSE157827")$nFeature_RNA)
fivenum(filter(seu@meta.data, dataset=="GSE157827")$nCount_RNA)
fivenum(filter(seu@meta.data, dataset=="GSE157827")$percent_mito)
fivenum(filter(seu@meta.data, dataset=="GSE157827")$percent_ribo)
```

**设置阈值**：

```r
# 细胞阈值
cutoff_tbl <- tibble::tribble(
  ~dataset,     ~min_features, ~max_features, ~min_counts, ~max_counts, ~max_mito, ~max_ribo,
  "GSE157827",  200,           8000,          500,        40000,       14,        3,
  "GSE174367",  300,          12000,         1000,        50000,       10,        2
)
keep <- rep(FALSE, ncol(seu))
for (i in seq_len(nrow(cutoff_tbl))) {
  ds <- cutoff_tbl$dataset[i]
  idx <- which(seu$dataset == ds)
  if (length(idx) == 0) next
  x <- seu@meta.data[idx, , drop = FALSE]
  keep[idx] <-
    x$nFeature_RNA > cutoff_tbl$min_features[i] &
    x$nFeature_RNA < cutoff_tbl$max_features[i] &
    x$nCount_RNA   > cutoff_tbl$min_counts[i] &
    x$nCount_RNA   < cutoff_tbl$max_counts[i] &
    x$percent_mito < cutoff_tbl$max_mito[i] &
    x$percent_ribo < cutoff_tbl$max_ribo[i]
}
scRNA_cells <- subset(seu, cells = colnames(seu)[keep])
# 基因阈值
counts_rna <- GetAssayData(scRNA_cells, assay = "RNA", slot = "counts")
feature_rowsum <- Matrix::rowSums(counts_rna > 0)
retained_f_low <- feature_rowsum > ncol(scRNA_cells) * 0.0005
genes_keep <- names(feature_rowsum)[retained_f_low]
scRNA_cells[["RNA"]] <- subset(scRNA_cells[["RNA"]], features = genes_keep)
# hERV位点：去掉--no-feature
counts_herv <- GetAssayData(scRNA_cells, assay = "HERV", slot = "counts")
herv_keep <- rownames(counts_herv)
herv_keep <- herv_keep[! herv_keep %in% c("--no-feature", "__no_feature")]
scRNA_cells[["HERV"]] <- subset(scRNA_cells[["HERV"]], features = herv_keep)
# 查看过滤后结果
scRNA_cells
table(scRNA_cells$dataset)
table(scRNA_cells$group)
table(scRNA_cells$orig.ident)
# 更新nCount/nFeature
rna_counts <- GetAssayData(scRNA_cells, assay = "RNA", slot = "counts")
scRNA_cells$nCount_RNA <- Matrix::colSums(rna_counts)
scRNA_cells$nFeature_RNA <- Matrix::colSums(rna_counts > 0)
herv_counts <- GetAssayData(scRNA_cells, assay = "HERV", slot = "counts")
scRNA_cells$nCount_HERV <- Matrix::colSums(herv_counts)
scRNA_cells$nFeature_HERV <- Matrix::colSums(herv_counts > 0)
scRNA_cells$HERV_fraction <- scRNA_cells$nCount_HERV / (scRNA_cells$nCount_HERV + scRNA_cells$nCount_RNA)
```

看看哪些指标过滤了细胞+基因阈值图

```r
datasets <- unique(seu$dataset)
for (ds in datasets) {
  cut1 <- cutoff_tbl %>% filter(dataset == ds)
  if (nrow(cut1) != 1) next
  sub <- subset(seu, subset = dataset == ds)
  retained_c_feature_low  <- sub$nFeature_RNA >= cut1$min_features
  retained_c_feature_high <- sub$nFeature_RNA <= cut1$max_features
  retained_c_count_low    <- sub$nCount_RNA   >= cut1$min_counts
  retained_c_count_high   <- sub$nCount_RNA   <= cut1$max_counts
  retained_c_mito         <- sub$percent_mito <= cut1$max_mito
  retained_c_ribo         <- sub$percent_ribo <= cut1$max_ribo
  keep_cells <- retained_c_feature_low &
                retained_c_feature_high &
                retained_c_count_low &
                retained_c_count_high &
                retained_c_mito &
                retained_c_ribo
  # Venn图
  venn_eles <- setNames(
    list(
      Cells(sub)[!retained_c_feature_low],
      Cells(sub)[!retained_c_feature_high],
      Cells(sub)[!retained_c_count_low],
      Cells(sub)[!retained_c_count_high],
      Cells(sub)[!retained_c_mito],
      Cells(sub)[!retained_c_ribo]
    ),
    c(
      paste0("Filt_low_feature\n(", cut1$min_features, ")"),
      paste0("Filt_high_feature\n(", cut1$max_features, ")"),
      paste0("Filt_low_count\n(", cut1$min_counts, ")"),
      paste0("Filt_high_count\n(", cut1$max_counts, ")"),
      paste0("Filt_high_mito\n(", cut1$max_mito, "%)"),
      paste0("Filt_high_ribo\n(", cut1$max_ribo, "%)")
    )
  )
  pdf(file.path(res_dir, paste0(ds, "_cell_filter_venn.pdf")), width = 14, height = 10)
  p1 <- ggVennDiagram(venn_eles, set_size = 4) +
    scale_fill_gradient(low = "#ffffb2", high = "#f03b20") +
    ggtitle(paste0(ds, " cell filtering overview"))
  print(p1)
  dev.off()
  #gene filter rankplot
  sub_keep <- subset(sub, cells = Cells(sub)[keep_cells])
  counts <- GetAssayData(sub_keep, assay = "RNA", slot = "counts")
  feature_rowsum <- Matrix::rowSums(counts > 0)
  rankplot <- data.frame(
    feature_count = sort(feature_rowsum),
    gene = names(sort(feature_rowsum)),
    Rank = seq_along(feature_rowsum)
  )
  gene_cutoff <- ncol(sub_keep) * 0.0005
  pdf(file.path(res_dir, paste0(ds, "_gene_filter_rankplot.pdf")), width = 12, height = 8)
  p2 <- ggplot(rankplot, aes(x = Rank, y = feature_count)) +
    geom_point(size = 0.4) +
    scale_y_break(c(10000, 100000)) +
    geom_hline(yintercept = gene_cutoff, color = "red") +
    annotate(
      "text",
      x = max(rankplot$Rank) * 0.15,
      y = gene_cutoff * 1.2,
      label = "Feature cutoff : ncol(filtered_cells) * 0.05%",
      size = 5,
      hjust = 0
    ) +
    ggtitle(paste0(ds, " gene filter rankplot")) +
    theme_classic()
  print(p2)
  dev.off()
}
```

**去双细胞**：

```r
sample_col <- "sample_uid"
dbl_class <- rep(NA_character_, ncol(scRNA_cells))
dbl_score <- rep(NA_real_, ncol(scRNA_cells))
sample_ids <- unique(scRNA_cells[[sample_col]][, 1])
for (sid in sample_ids) {
  message("Running scDblFinder on ", sid)
  cells <- colnames(scRNA_cells)[scRNA_cells[[sample_col]][, 1] == sid]
  sub <- subset(scRNA_cells, cells = cells)
  sce <- SingleCellExperiment(
    assays = list(counts = GetAssayData(sub, assay = "RNA", layer = "counts"))
  )
  sce <- scDblFinder(sce)
  dbl_class[match(cells, colnames(scRNA_cells))] <- as.character(colData(sce)$scDblFinder.class)
  dbl_score[match(cells, colnames(scRNA_cells))] <- as.numeric(colData(sce)$scDblFinder.score)
}
scRNA_cells$doublet_class <- dbl_class
scRNA_cells$doublet_score <- dbl_score
scDblFinder_res <- as.data.frame.array(table(scRNA_cells$GSM, scRNA_cells$doublet_class)) %>%
  mutate(doublet_freq=doublet/(doublet+singlet))
write.csv(scDblFinder_res, file = file.path(res_dir, "scDblFinder_res.csv"))
# 过滤
seu_filt <- subset(scRNA_cells, subset = doublet_class == "singlet")
# 查看过滤后结果
seu_filt
table(seu_filt$dataset)
table(seu_filt$group)
table(seu_filt$orig.ident)
```

**过滤后画图**：

```r
p1 <- VlnPlot(seu_filt, features = "nFeature_RNA",
              pt.size = 0, ncol = 1, group.by = "dataset") + NoLegend() + theme(axis.title.x = element_blank())
p1_1 <- VlnPlot(seu_filt, features = "nCount_RNA",
              pt.size = 0, ncol = 1, group.by = "dataset") + NoLegend() + theme(axis.title.x = element_blank())
p2 <- VlnPlot(seu_filt, features = "nCount_RNA",
              pt.size = 0, ncol = 1, group.by = "orig.ident") + NoLegend() + theme(axis.title.x = element_blank())
p3 <- VlnPlot(seu_filt, features = "nFeature_HERV",
              pt.size = 0, ncol = 1, group.by = "dataset") + NoLegend() + theme(axis.title.x = element_blank())
p3_1 <- VlnPlot(seu_filt, features = "nCount_HERV",
              pt.size = 0, ncol = 1, group.by = "dataset") + NoLegend() + theme(axis.title.x = element_blank())
p5 <- VlnPlot(seu_filt, features = "percent_mito",
              pt.size = 0, ncol = 1, group.by = "dataset") + NoLegend() + theme(axis.title.x = element_blank())
p5_1 <- VlnPlot(seu_filt, features = "percent_ribo",
              pt.size = 0, ncol = 1, group.by = "dataset") + NoLegend() + theme(axis.title.x = element_blank())
p6 <- VlnPlot(seu_filt, features = "percent_mito",
              pt.size = 0, ncol = 1, group.by = "orig.ident") + NoLegend() + theme(axis.title.x = element_blank())
p7 <- VlnPlot(seu_filt, features = "HERV_fraction",
              pt.size = 0, ncol = 1, group.by = "dataset") + NoLegend() + theme(axis.title.x = element_blank())
p8 <- VlnPlot(seu_filt, features = "HERV_fraction",
              pt.size = 0, ncol = 1, group.by = "orig.ident") + NoLegend() + theme(axis.title.x = element_blank())
pdf(file.path(res_dir, "after_filter.pdf"), width = 24, height = 20)
p <- (p1 | p1_1 | p5 | p5_1 | p3 | p3_1 | p7) / (p2 | p6 | p8)
print(p)
dev.off()
# 保存数据
saveRDS(
    seu_filt, 
    file = "/public/home/GENE_proc/wth/my_data/all_filt.rds", 
    compress = "xz"
)
```

#### 数据预处理+去批次效应

```r
seu <- readRDS("/public/home/GENE_proc/wth/my_data/all_filt.rds")
DefaultAssay(seu) <- "RNA"
```

**integration流程**（和下面的harmony流程二选一）：

```r
seu <- NormalizeData(seu)  # 计算RNA assay的data
obj.list <- SplitObject(seu, split.by = "orig.ident")  # 按样本拆分
obj.list <- lapply(obj.list, function(x) {  # 每个样本分别做SCTransform
  SCTransform(
    x,
    assay = "RNA",
    new.assay.name = "SCT",
    vars.to.regress = "percent_mito",
    return.only.var.genes = TRUE,
    conserve.memory = TRUE,
    method = "glmGamPoi"
  )
})
features <- SelectIntegrationFeatures(  # 选integration features
  object.list = obj.list,
  nfeatures = 2000
)
obj.list <- PrepSCTIntegration(
  object.list = obj.list,
  anchor.features = features,
)
obj.list <- lapply(obj.list, function(x) {
  RunPCA(
    x,
    assay = "SCT",
    features = features,
    npcs = 30,
  )
})
anchors <- FindIntegrationAnchors(  # 找anchors，为了节省时间和内存，使用rpca
  object.list = obj.list,
  normalization.method = "SCT",
  anchor.features = features,
  dims = 1:20,
  reduction = "rpca",  # 使用原方法需要大约5h，rpca仅需30min
)
seu <- IntegrateData(  # 整合，可能需要128G+内存，运行时间需要30min+
  anchorset = anchors,
  normalization.method = "SCT",
  dims = 1:20,
)
DefaultAssay(seu) <- "integrated"  # 用integrated assay做降维和聚类
seu <- RunPCA(seu, npcs = 50)
pdf(file.path(res_dir, "ElbowPlot_all.pdf"))
p <- ElbowPlot(seu, ndims = 50)
print(p)
dev.off()
pc.use <- 1:30  # 选主成分数量
seu <- RunUMAP(seu, dims = pc.use)
p_dim_group <- DimPlot(seu, group.by = "group")
p_dim_dataset <- DimPlot(seu, group.by = "dataset", raster = TRUE)
# 看一下效果
DimPlot(seu, group.by = "orig.ident", raster = TRUE)
p_dim_group
p_dim_dataset
# 再画clustree
seu <- FindNeighbors(seu, dims = pc.use)
seu <- FindClusters(seu, resolution = c(0.01, 0.05, 0.1, 0.2, 0.5))
p_cutree <- clustree(seu@meta.data, prefix = "RNA_snn_res.")
p_cutree
```

**harmony流程**：

```r
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = 3000)
seu <- ScaleData(  # 大约30min-1h
  seu,
  vars.to.regress = c("nCount_RNA", "percent_mito"),
  features = VariableFeatures(seu),
)
seu <- RunPCA(
  seu, 
  npcs = 50, 
  features = VariableFeatures(seu),
)
pdf(file.path(res_dir, "ElbowPlot_all.pdf"))
p <- ElbowPlot(seu, ndims = 50)
print(p)
dev.off()
pc.use <- 1:30  # 选主成分数量
seu <- RunHarmony(
  object = seu,
  group.by.vars = "sample_uid",
  reduction.use = "pca",
  dims.use = pc.use,
  assay.use = "RNA",
  verbose = TRUE
)
seu <- RunUMAP(seu, reduction = "harmony", dims = pc.use)
p_dim_group <- DimPlot(seu, group.by = "group")
p_dim_dataset <- DimPlot(seu, group.by = "dataset", raster = TRUE)
# 看一下效果
DimPlot(seu, group.by = "orig.ident", raster = TRUE)
p_dim_group
p_dim_dataset
# 再画clustree
seu <- FindNeighbors(seu, reduction = "harmony", dims = pc.use)
seu <- FindClusters(seu, resolution = c(0.01, 0.05, 0.1, 0.2, 0.5))
p_cutree <- clustree(seu@meta.data, prefix = "RNA_snn_res.")
p_cutree
```

**选分辨率**：

```r
Idents(seu) <- seu$RNA_snn_res.0.1
p_dim_cluster <- DimPlot(seu, label = TRUE)
pdf(file.path(res_dir, "data_preprocess.pdf"), width = 16, height = 16)
p <- (p_dim_group | p_dim_dataset) / (p_cutree | p_dim_cluster)
print(p)
dev.off()
```

#### 细胞类型注释

以下是接续harmony流程的，如果用SCT流程，除了设置`DefaultAssay(seu) <- "SCT"`，还可能需要先运行`PrepSCTFindMarkers()`

```r
# marker基因
gene_list <- list(
  Astro = c("AQP4", "GFAP", "SLC1A2", "ADGRV1", "GPC5", "RYR3", "CHI3L1"),
  Endo  = c("CLDN5", "FLT1", "ABCB1", "EBF1", "VWF", "PECAM1"),
  Excit = c("CAMK2A", "CBLN2", "LDB2", "SNAP25", "SYT1", "SLC17A7", "SATB2"),
  Inhib = c("GAD1", "GAD2", "LHFPL3", "PCDH15", "SV2C"),
  Mic   = c("C3", "LRMDA", "DOCK8", "P2RY12", "CSF1R", "CX3CR1", "CD74", "SPP1", "CD163"),
  Oligo = c("MBP", "PLP1", "ST18", "MOBP", "MOG", "OPALIN", "CNP"),
  OPC   = c("PDGFRA", "CSPG4", "VCAN")
)
# 在各个聚类中的表达量
pdf(file.path(res_dir, "cluster_marker_gene.pdf"), width = 24, height = 12)
p <- DotPlot(seu, features = gene_list) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust = 0.8))
print(p)
dev.off()
# 计算每个cluster对每个细胞类型的marker-set的平均表达
score_by_cluster <- function(seu, gene_list, group_col) {
  out <- lapply(names(gene_list), function(ct) {
    feats <- intersect(gene_list[[ct]], rownames(seu))
    if (length(feats) == 0) return(NULL)
    avg <- AverageExpression(
      seu,
      assays = "RNA",
      features = feats,
      group.by = group_col,
      slot = "data",
      verbose = FALSE
    )$RNA  # features x clusters
    data.frame(
      cluster = colnames(avg),
      celltype = ct,
      score = colMeans(avg, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }) |> bind_rows()
  tab <- out |>
    pivot_wider(names_from = celltype, values_from = score) |>
    arrange(as.numeric(as.character(cluster)))
  mat <- as.matrix(tab[, intersect(names(gene_list), colnames(tab)), drop = FALSE])
  tab$pred <- colnames(mat)[max.col(mat, ties.method = "first")]
  tab
}
score_tab <- score_by_cluster(seu, gene_list, "RNA_snn_res.0.1")
View(score_tab)
```

**根据表达量图指定细胞类型**：

```r
# 指定细胞类型
cl_Astro <- c("2")
cl_Endo <- c("12")
cl_Mic <- c("6")
cl_Oligo <- c("0", "16")
cl_OPC <- c("3")
cl_Inhib <- c("4", "5", "10", "13") 
# 其它细胞归为Excit
cl <- as.character(seu$RNA_snn_res.0.1)
seu$celltype <- ifelse(cl %in% cl_Astro, "Astro",
                  ifelse(cl %in% cl_Endo, "Endo",
                    ifelse(cl %in% cl_Mic, "Mic",
                      ifelse(cl %in% cl_Oligo, "Oligo",
                        ifelse(cl %in% cl_OPC, "OPC",
                          ifelse(cl %in% cl_Inhib, "Inhib", "Excit"))))))
seu$celltype <- factor(seu$celltype, levels = c("Astro","Endo","Excit","Inhib","Mic","Oligo", "OPC"))
Idents(seu) <- seu$celltype
table(seu$celltype, seu$group)
table(seu$celltype)
# 保存数据
saveRDS(
    seu, 
    file = "/public/home/GENE_proc/wth/my_data/all_celltype.rds",
    compress = "xz"
)
```

**画图展示**（模仿原论文）：

```r
DefaultAssay(seu) <- "RNA"
plot_tag_fz <- 22
ct_cols <- c(
  Astro = "#1b9e77",
  Endo  = "#d95f02",
  Excit = "#7570b3",
  Inhib = "#e7298a",
  Mic   = "#66a61e",
  Oligo = "#e6ab02",
  OPC   = "#a6761d"
)
dim_plot_raw <- DimPlot(  # 细胞类型在umap上的分布
  seu, 
  reduction = "umap", 
  cols = ct_cols, 
  label = FALSE
)
dim_plot <- dim_plot_raw +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.text = element_text(size = 16),
    plot.tag = element_text(size = plot_tag_fz, face = "bold"),
    plot.tag.position = c(-0.03, 0.995),
  )
ct_stat <- as.data.frame(table(seu$celltype))  # 统计各细胞类型的细胞数和比例
colnames(ct_stat) <- c("celltype", "total")
ct_stat$celltype <- factor(ct_stat$celltype, levels = levels(seu@active.ident))
ct_stat$percentage <- 100 * ct_stat$total / sum(ct_stat$total)
celltype_bar <- ggplot(ct_stat, aes(x = celltype, y = percentage, fill = celltype)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = ct_cols) +
  theme_classic() +
  ylab("Proportion of total nuclei (%)") +
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 18, face = "bold"), 
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    legend.position = "none",
    plot.tag = element_text(size = plot_tag_fz, face = "bold"),
  )
genes_paper <- c(  # 各细胞类型marker热图
  "AQP4", "GFAP", "SLC1A2", "GPC5", "ADGRV1",
  "CLDN5", "FLT1", "EBF1",
  "SLC17A7", "CAMK2A", "SATB2", "SNAP25", "CBLN2",
  "GAD1", "GAD2", "LHFPL3", "PCDH15",
  "P2RY12", "CSF1R", "CX3CR1", "CD74", "C3",
  "MBP", "PLP1", "MOBP", "OPALIN", "MOG",
  "PDGFRA", "CSPG4", "VCAN"
)
genes_use <- intersect(genes_paper, rownames(seu))
# 如果做出来的热图不好看可以试试↓，注意先保存seu
seu <- ScaleData(seu, features = genes_use)
heatmap_raw <- DoHeatmap(
  seu,
  features = genes_use,
  disp.min = -1.5,
  disp.max = 1.5,
  angle = 0,
  group.colors = unname(ct_cols),
  label = FALSE
)
heatmap <- heatmap_raw + 
  scale_fill_gradientn(
    colors = c("#f100ef", "black", "#f0e442"),
    limits = c(-1.5, 1.5)
  ) +
  theme(
    axis.text.y = element_text(size = 14, face = "italic"),
    legend.title = element_text(size = 16.5),
    legend.text = element_text(size = 16),
    plot.tag = element_text(size = plot_tag_fz, face = "bold"),
  )
right_plot <- celltype_bar / heatmap +
  plot_layout(heights = c(1, 1.5))
# 拼图
final_plot <- dim_plot + right_plot +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(
    theme = theme(
      plot.tag.position = c(0, 1),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(20, 10, 20, 50) 
    )
  )
ggsave(
  file.path(res_dir, "celltype.pdf"),
  final_plot, width = 24, height = 16
)
```

#### 计算差异基因

```r
scRNA <- readRDS("/public/home/GENE_proc/wth/my_data/all_celltype.rds")
DefaultAssay(scRNA) <- "RNA"
scRNA$group <- factor(scRNA$group, levels = c("NC","AD"))
scRNA$celltype <- factor(
  scRNA$celltype,
  levels = c("Astro","Endo","Excit","Inhib","Mic","Oligo","OPC")
)
Idents(scRNA) <- scRNA$celltype
balance_cells <- function(obj, ct){
  sub <- subset(obj, idents = ct)
  ad <- WhichCells(sub, expression = group == "AD")
  nc <- WhichCells(sub, expression = group == "NC")
  if(length(ad) > 2000){
    n <- min(length(ad), length(nc))
    sub[, c(sample(ad, n), sample(nc, n))]
  } else {
    sub
  }
}
cts <- levels(scRNA$celltype)
obj_bal <- merge(x = balance_cells(scRNA, cts[1]),
                 y = lapply(cts[-1], \(ct) balance_cells(scRNA, ct)))
obj_bal$celltype <- factor(
  obj_bal$celltype,
  levels = c("Astro","Endo","Excit","Inhib","Mic","Oligo","OPC")
)
obj_bal <- JoinLayers(obj_bal, assay = "RNA")
obj_bal$compare <- paste(obj_bal$group, obj_bal$celltype, sep = "_")
scRNA$compare <- paste(scRNA$group, scRNA$celltype, sep = "_")
ct <- levels(obj_bal$celltype)
all_markers <- lapply(ct, function(x){
  message("Running: ", x)
  FindMarkers(
    obj_bal,
    group.by = "compare",
    ident.1 = paste0("AD_", x),
    ident.2 = paste0("NC_", x),
    test.use = "MAST",
    # latent.vars = c("nCount_RNA", "percent_mito"),
    logfc.threshold = 0.25,
    min.pct = 0.1
  )
})
names(all_markers) <- ct
lapply(all_markers, nrow)
all_markers_sig <- lapply(all_markers, function(df){
  df2 <- subset(
    df, 
    p_val_adj < 0.05
  )
  return(df2)
})
View(all_markers_sig[["Astro"]])
save(all_markers_sig, file = "/public/home/GENE_proc/wth/my_data/celltype_markers_sig.RData")
```

#### hERV按subfamily/class聚合+标准化

**构建位点-家族映射**：

```r
get_attr <- function(x, key) {
  pat <- paste0(key, ' "([^"]+)"')
  m <- regexec(pat, x, perl = TRUE)
  rr <- regmatches(x, m)
  out <- rep(NA_character_, length(x))
  ok <- lengths(rr) >= 2
  out[ok] <- vapply(rr[ok], `[`, character(1), 2)
  out
}
gtf_path <- "/public/home/wangtianhao/Desktop/STAR_ref/transcripts.gtf"
lines <- readLines(gtf_path, warn = FALSE)
lines <- lines[!grepl("^\\s*#", lines)]
lines <- lines[nzchar(lines)]
gtf <- data.table::fread(
  text = lines,
  sep = "\t",
  header = FALSE,
  quote = "",
  fill = TRUE,
  data.table = TRUE
)
gtf <- gtf[, 1:9]
setnames(gtf, c("seqname","source","feature","start","end","score","strand","frame","attribute"))
gtf[, gene_id    := get_attr(attribute, "gene_id")]
gtf[, locus      := get_attr(attribute, "locus")]
gtf[, intModel   := get_attr(attribute, "intModel")]
gtf[, repName    := get_attr(attribute, "repName")]
gtf[, repFamily  := get_attr(attribute, "repFamily")]
gtf[, geneRegion := get_attr(attribute, "geneRegion")]
gtf[, category := get_attr(attribute, "category")]
# gene行：intModel（subfamily）
gene_dt <- gtf[feature == "gene", .(locus_id = gene_id, intModel, category)][!is.na(locus_id)]
# exon行：repFamily（class）
exon_dt <- gtf[feature == "exon", .(locus_id = gene_id, repFamily)][!is.na(locus_id)]
# 每个位点的class：取出现最多的repFamily
repfam_dt <- exon_dt[!is.na(repFamily),
                     .(repFamily = {
                       x <- repFamily
                       ux <- unique(x); ux[which.max(tabulate(match(x, ux)))]
                     }),
                     by = locus_id]
map <- merge(gene_dt, repfam_dt, by = "locus_id", all = TRUE)
map[, subfamily := intModel]
map[, class := repFamily]
map$subfamily_raw <- map$subfamily
map$subfamily <- gsub("-int", "", map$subfamily)  # 去掉结尾的-int标记
herv_map <- unique(map[, .(locus_id, subfamily_raw, subfamily, class, category)])
group2family <- readxl::read_xlsx("/public/home/GENE_proc/wth/metadata/hERV/family2group.xlsx")
herv_group_map <- merge(herv_map, group2family, by = "subfamily_raw") %>% 
  select(locus_id, category, subfamily_raw, subfamily, subfamily_alias, group, class, include)
write.csv(herv_group_map, file = "/public/home/GENE_proc/wth/metadata/hERV/hERV_locus2family.csv", row.names = F)
```

**聚合并标准化**：

```r
seu <- readRDS("/public/home/GENE_proc/wth/my_data/all_celltype.rds")
herv_map <- read.csv("/public/home/GENE_proc/wth/metadata/hERV/hERV_locus2family.csv")
herv_map$locus_id <- gsub("_", "-", herv_map$locus_id)
herv_counts <- GetAssayData(seu, assay = "HERV", slot = "counts")
aggregate_herv_counts <- function(herv_counts, map, group_by = c("subfamily", "subfamily_raw", "class", "group")) {
  group_by <- match.arg(group_by)
  keep_loci <- intersect(rownames(herv_counts), map$locus_id)
  map2 <- map[match(keep_loci, map$locus_id), ]
  map2 <- map2[!is.na(map2[[group_by]]) & nzchar(map2[[group_by]]), ]
  mat <- herv_counts[map2$locus_id, , drop = FALSE]
  grp <- map2[[group_by]]
  f <- factor(grp, levels = unique(grp))
  G <- sparseMatrix(
    i = seq_along(f),
    j = as.integer(f),
    x = 1,
    dims = c(length(f), nlevels(f)),
    dimnames = list(rownames(mat), levels(f))
  )
  out <- Matrix::t(G) %*% mat
  out
}
herv_subfamily_counts <- aggregate_herv_counts(herv_counts, herv_map, group_by = "subfamily_raw")
herv_class_counts <- aggregate_herv_counts(herv_counts, herv_map, group_by = "class")
herv_group_counts <- aggregate_herv_counts(herv_counts, herv_map, group_by = "group")
# 将新计数矩阵的细胞与原来的对齐
align_to_cells <- function(mat, cells){
  mat <- as(mat, "dgCMatrix")
  extra_cells <- setdiff(colnames(mat), cells)
  if (length(extra_cells) > 0) {  # 如果herv计数结果中有多余细胞：删掉
    mat <- mat[, setdiff(colnames(mat), extra_cells), drop = FALSE]
  }
  miss_cells <- setdiff(cells, colnames(mat))
  if (length(miss_cells) > 0) {  # 如果herv计数结果中缺少某些细胞：补0
    zero <- Matrix(
      0,
      nrow = nrow(mat),
      ncol = length(miss_cells),
      sparse = TRUE,
      dimnames = list(rownames(mat), miss_cells)
    )
    mat <- cbind(mat, zero)
  }
  mat[, cells, drop = FALSE]
}
cells <- colnames(seu)
herv_counts <- align_to_cells(herv_counts, cells)
herv_subfamily_counts <- align_to_cells(herv_subfamily_counts, cells)
herv_class_counts <- align_to_cells(herv_class_counts, cells)
herv_group_counts <- align_to_cells(herv_group_counts, cells)
# 加到seu上
replace_assay <- function(seu, assay_name, counts){
  counts <- as(counts, "dgCMatrix")
  seu[[assay_name]] <- CreateAssayObject(counts = counts)
  seu
}
seu <- replace_assay(seu, "HERV_family", herv_subfamily_counts)
seu <- replace_assay(seu, "HERV_class", herv_class_counts)
seu <- replace_assay(seu, "HERV_group", herv_group_counts)
# 标准化
normalize_herv_by_rna_depth <- function(seu, assay, scale.factor = 10000) {
  rna_counts <- GetAssayData(seu, assay = "RNA",  slot = "counts")
  herv_counts <- GetAssayData(seu, assay = assay, slot = "counts")
  lib_rna <- Matrix::colSums(rna_counts)
  herv_norm <- t(t(herv_counts) / lib_rna) * scale.factor
  herv_norm@x <- log1p(herv_norm@x)
  seu <- SetAssayData(seu, assay = assay, slot = "data", new.data = herv_norm)
  return(seu)
}
seu <- normalize_herv_by_rna_depth(seu, "HERV")
seu <- normalize_herv_by_rna_depth(seu, "HERV_family")
seu <- normalize_herv_by_rna_depth(seu, "HERV_class")
seu <- normalize_herv_by_rna_depth(seu, "HERV_group")
# 保存数据
saveRDS(
    seu, 
    file = "/public/home/GENE_proc/wth/my_data/all_final.rds",
    compress = "xz"
)
```

**简单验证一下**：

```r
herv_subfamily_counts[, 1]
sum(herv_subfamily_counts[, 1])  # 都是30
herv_class_counts[, 1]
sum(herv_class_counts[, 1])
herv_group_counts[, 1]
sum(herv_group_counts[, 1])
herv_counts[, 1]
sum(herv_counts[, 1])
herv_subfamily_counts[, 2]
sum(herv_subfamily_counts[, 2])  # 都是14
herv_class_counts[, 2]
sum(herv_class_counts[, 2])
herv_group_counts[, 2]
sum(herv_group_counts[, 2])
herv_counts[, 2]
sum(herv_counts[, 2])
herv_data <- GetAssayData(seu, assay = "HERV", slot = "data")
herv_data[, 1]
HERV_family_data <- GetAssayData(seu, assay = "HERV_family", slot = "data")
HERV_family_data[, 1]
HERV_class_data <- GetAssayData(seu, assay = "HERV_class", slot = "data")
HERV_class_data[, 1]
HERV_group_data <- GetAssayData(seu, assay = "HERV_group", slot = "data")
HERV_group_data[, 1]
```

### scTE补充分析

使用之前读取stellarscope计数矩阵时的函数`normalize_barcode`/`align_to_cells`/`replace_assay`/`normalize_herv_by_rna_depth`

#### 读取scTE计数矩阵

```r
gsm_info <- read.csv("/public/home/GENE_proc/wth/my_data/gsm_info.csv")
roots <- c(
  "/public/home/GENE_proc/wth/GSE157827/scTE_res",
  "/public/home/GENE_proc/wth/GSE174367/scTE_res"
)
ltr2herv <- fread("/public/home/GENE_proc/wth/metadata/hERV/LTR2HERV.tsv") %>%
  distinct(repName, subfamily)
ltr2class <- fread("/public/home/GENE_proc/wth/metadata/hERV/LTR2Class.tsv") %>%
  distinct(repName, repClass)
rmsk_repname <- sort(unique(c(ltr2herv$repName, ltr2class$repName)))
csv_files <- c(
  list.files(roots[1], pattern = "^GSM.*\\.csv$", full.names = TRUE),
  list.files(roots[2], pattern = "^GSM.*\\.csv$", full.names = TRUE)
)
csv_files <- csv_files[order(basename(csv_files))]
all_features <- unique(unlist(lapply(csv_files, function(f) {
  names(fread(f, nrows = 0, check.names = FALSE))[-1]
})))
# HERV_subfamily_scTE：只保留以HERV/ERV开头，且属于rmsk repName的feature
feat_subfamily <- intersect(
  grep("^(HERV|ERV)", rmsk_repname, value = TRUE),
  all_features
)
feat_subfamily <- sort(feat_subfamily)
# HERV_subfamily_LTR_scTE：按LTR2HERV.tsv聚合
ltr2herv <- ltr2herv %>% filter(repName %in% all_features)
feat_subfamily_ltr <- sort(unique(ltr2herv$repName))
subfamily_levels <- sort(unique(ltr2herv$subfamily))
# HERV_class_scTE：按LTR2Class.tsv聚合
ltr2class <- ltr2class %>% filter(repName %in% all_features)
feat_class <- sort(unique(ltr2class$repName))
class_levels <- sort(unique(ltr2class$repClass))
# 开始读取
feat_read_all <- unique(c(feat_subfamily, feat_subfamily_ltr, feat_class))
mat_subfamily_list <- list()
mat_subfamily_ltr_list <- list()
mat_class_list <- list()
for (f in csv_files) {
  gsm <- sub("\\.csv$", "", basename(f))
  indv <- gsm_info %>% filter(GSM == gsm) %>% pull(orig.ident) %>% unique()
  hdr <- names(fread(f, nrows = 0, check.names = FALSE))
  feat_read <- intersect(feat_read_all, hdr[-1])
  dt <- fread(
    f,
    select = c(hdr[1], feat_read),
    check.names = FALSE
  )
  bc <- normalize_barcode(dt[[1]])
  cells_this <- paste(gsm, indv, bc, sep = "_")
  mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- cells_this
  colnames(mat) <- feat_read
  # HERV_subfamily_scTE
  m1 <- t(mat[, intersect(feat_subfamily, feat_read), drop = FALSE])
  miss1 <- setdiff(feat_subfamily, rownames(m1))
  if (length(miss1) > 0) {
    zero1 <- Matrix(
      0,
      nrow = length(miss1),
      ncol = ncol(m1),
      sparse = TRUE,
      dimnames = list(miss1, colnames(m1))
    )
    m1 <- rbind(m1, zero1)
  }
  m1 <- Matrix(m1[feat_subfamily, , drop = FALSE], sparse = TRUE)
  # HERV_subfamily_LTR_scTE
  m2 <- t(mat[, intersect(feat_subfamily_ltr, feat_read), drop = FALSE])
  grp2 <- ltr2herv$subfamily[match(rownames(m2), ltr2herv$repName)]
  m2 <- rowsum(m2, group = grp2)
  miss2 <- setdiff(subfamily_levels, rownames(m2))
  if (length(miss2) > 0) {
    zero2 <- Matrix(
      0,
      nrow = length(miss2),
      ncol = ncol(m2),
      sparse = TRUE,
      dimnames = list(miss2, colnames(m2))
    )
    m2 <- rbind(m2, zero2)
  }
  m2 <- Matrix(m2[subfamily_levels, , drop = FALSE], sparse = TRUE)
  # HERV_class_scTE
  m3 <- t(mat[, intersect(feat_class, feat_read), drop = FALSE])
  grp3 <- ltr2class$repClass[match(rownames(m3), ltr2class$repName)]
  m3 <- rowsum(m3, group = grp3)
  miss3 <- setdiff(class_levels, rownames(m3))
  if (length(miss3) > 0) {
    zero3 <- Matrix(
      0,
      nrow = length(miss3),
      ncol = ncol(m3),
      sparse = TRUE,
      dimnames = list(miss3, colnames(m3))
    )
    m3 <- rbind(m3, zero3)
  }
  m3 <- Matrix(m3[class_levels, , drop = FALSE], sparse = TRUE)
  mat_subfamily_list[[gsm]] <- m1
  mat_subfamily_ltr_list[[gsm]] <- m2
  mat_class_list[[gsm]] <- m3
}
herv_subfamily_scTE_counts <- do.call(cbind, mat_subfamily_list)
herv_subfamily_ltr_scTE_counts <- do.call(cbind, mat_subfamily_ltr_list)
herv_class_scTE_counts <- do.call(cbind, mat_class_list)
meta_scTE <- data.frame(
  GSM        = sub("_.*$", "", colnames(herv_subfamily_scTE_counts)),
  orig.ident = sub("^[^_]+_([^_]+)_.*$", "\\1", colnames(herv_subfamily_scTE_counts)),
  row.names  = colnames(herv_subfamily_scTE_counts),
  stringsAsFactors = FALSE
)
scTE_seu <- CreateSeuratObject(
  counts = herv_subfamily_scTE_counts,
  assay = "HERV_subfamily_scTE",
  meta.data = meta_scTE
)
scTE_seu <- replace_assay(scTE_seu, "HERV_subfamily_LTR_scTE", herv_subfamily_ltr_scTE_counts)
scTE_seu <- replace_assay(scTE_seu, "HERV_class_scTE", herv_class_scTE_counts)
save(herv_subfamily_scTE_counts, herv_subfamily_ltr_scTE_counts, herv_class_scTE_counts, file = "/public/home/GENE_proc/wth/my_data/scTE_counts.RData")
saveRDS(
  scTE_seu, 
  file = "/public/home/GENE_proc/wth/my_data/scTE.rds", 
)
```

#### 加到已有seu对象上

```r
add_scTE <- function(seu_name){
  my_print(seu_name, ":")
  my_print("- reading rds")
  seu <- readRDS(paste0("/public/home/GENE_proc/wth/my_data/", seu_name, ".rds"))
  ct <- tolower(identity_seu(seu))
  cells <- colnames(seu)
  my_print("- mapping")
  herv_subfamily_scTE_counts <- align_to_cells(herv_subfamily_scTE_counts, cells)
  herv_subfamily_ltr_scTE_counts <- align_to_cells(herv_subfamily_ltr_scTE_counts, cells)
  herv_class_scTE_counts <- align_to_cells(herv_class_scTE_counts, cells)
  seu <- replace_assay(seu, "HERV_subfamily_scTE", herv_subfamily_scTE_counts)
  seu <- replace_assay(seu, "HERV_subfamily_LTR_scTE", herv_subfamily_ltr_scTE_counts)
  seu <- replace_assay(seu, "HERV_class_scTE", herv_class_scTE_counts)
  seu <- normalize_herv_by_rna_depth(seu, "HERV_subfamily_scTE")
  seu <- normalize_herv_by_rna_depth(seu, "HERV_subfamily_LTR_scTE")
  seu <- normalize_herv_by_rna_depth(seu, "HERV_class_scTE")
  my_print("- saving")
  assign(ct, seu)
  save(list = ct, file = paste0("/public/home/GENE_proc/wth/my_data/", seu_name, ".RData"))
}
lapply(c("astro", "oligo", "excit", "inhib", "OPC", "endo", "mic", "all_final"), add_scTE)
```

### hERV重调情况

#### hERV总量/class重调

**获得总hERV计数**：

```r
# hERV计数矩阵
herv_counts <- GetAssayData(seu, assay = "HERV_class", slot = "counts")
herv_total_counts <- Matrix::colSums(herv_counts)
# RNA库深度
rna_counts <- GetAssayData(seu, assay = "RNA", slot = "counts")
rna_lib <- Matrix::colSums(rna_counts)
# 标准化
sf <- 1e4
herv_total_data <- log1p((herv_total_counts / rna_lib) * sf)
# 计算总hERV
new_counts_row <- Matrix::Matrix(herv_total_counts, nrow = 1, sparse = TRUE)
rownames(new_counts_row) <- "HERV_total"
colnames(new_counts_row) <- colnames(seu)
new_data_row <- Matrix::Matrix(herv_total_data, nrow = 1, sparse = TRUE)
rownames(new_data_row) <- "HERV_total"
colnames(new_data_row) <- colnames(seu)
# 添加到HERV_class的assay中
mat_counts <- GetAssayData(seu, assay = "HERV_class", slot = "counts")
mat_counts2 <- rbind(mat_counts, new_counts_row)
mat_data <- GetAssayData(seu, assay = "HERV_class", slot = "data")
mat_data2 <- rbind(mat_data, new_data_row)
new_assay <- CreateAssayObject(counts = mat_counts2)
new_assay <- SetAssayData(new_assay, slot="data", new.data = mat_data2)
seu[["HERV_class"]] <- new_assay
# 准备差异分析
Idents(seu) <- factor(seu$group, levels = c("NC","AD"))
DefaultAssay(seu) <- "HERV_class"
```

**全体细胞**：

```r
res_all <- FindMarkers(
  object = seu,
  ident.1 = "AD", ident.2 = "NC",
  test.use = "MAST",
  # latent.vars = c("nCount_RNA", "percent_mito"),
  logfc.threshold = 0, min.pct = 0, return.thresh = 1
)
write.csv(res_all, file.path(res_dir, "hERV_class_all.csv"))
View(res_all)
```

**各细胞类型**：

```r
celltypes <- sort(unique(seu$celltype))
run_one_ct <- function(ct) {
  obj <- subset(seu, subset = celltype == ct)
  n_tab <- table(obj$group)
  print(paste("running:", ct))
  fm <- FindMarkers(
    object = obj,
    ident.1 = "AD", ident.2 = "NC",
    test.use = "MAST",
    # latent.vars = c("nCount_RNA", "percent_mito"),
    logfc.threshold = 0, min.pct = 0, return.thresh = 1
  )
  out <- fm %>%
    tibble::rownames_to_column("feature") %>%
    mutate(celltype = ct,
           n_NC = unname(n_tab["NC"]),
           n_AD = unname(n_tab["AD"]))
  out
}
res_by_ct <- do.call(rbind, lapply(celltypes, run_one_ct))
write.csv(res_by_ct, file.path(res_dir, "hERV_class_celltype.csv"))
# 只看p值<0.05的
View(res_by_ct %>%
  dplyr::filter(p_val_adj<0.05) %>% 
  arrange(p_val_adj) %>%
  select(celltype, n_NC, n_AD, avg_log2FC, p_val, p_val_adj, everything()))
# 全部数据的logFC
View(
  res_by_ct %>%
    select(celltype, avg_log2FC, feature) %>% 
    tidyr::pivot_wider(id_cols = celltype, names_from = feature, values_from = avg_log2FC)
)
# 过滤p值后的logFC
View(
  res_by_ct %>%
    dplyr::filter(p_val_adj<0.05) %>% 
    select(celltype, avg_log2FC, feature) %>% 
    tidyr::pivot_wider(id_cols = celltype, names_from = feature, values_from = avg_log2FC,)
)
# 画图展示
plot_class_violin <- function(seu, class_ct, assay = "HERV_class", celltype_order = c("Astro","Endo","Excit","Inhib","Mic","Oligo","OPC"), features = c("ERV1","ERVK","ERVL","HERV-total"), class_order = c("ERVL","ERVK","ERV1","HERV-total"), is_plot = T){
  p_to_star <- function(p){
    dplyr::case_when(
      is.na(p) ~ "",
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    )
  }
  herv_counts <- GetAssayData(seu, assay = assay,  layer = "counts")
  herv_total_counts <- Matrix::colSums(herv_counts)
  rna_counts <- GetAssayData(seu, assay = "RNA",  layer = "counts")
  rna_lib <- Matrix::colSums(rna_counts)
  herv_total_data <- log1p((herv_total_counts / rna_lib) * 1e4)
  new_counts_row <- Matrix::Matrix(herv_total_counts, nrow = 1, sparse = TRUE)
  rownames(new_counts_row) <- "HERV_total"
  colnames(new_counts_row) <- colnames(seu)
  new_data_row <- Matrix::Matrix(herv_total_data, nrow = 1, sparse = TRUE)
  rownames(new_data_row) <- "HERV_total"
  colnames(new_data_row) <- colnames(seu)
  mat_counts <- GetAssayData(seu, assay = assay,  layer = "counts")
  mat_counts2 <- rbind(mat_counts, new_counts_row)
  mat_data <- GetAssayData(seu, assay = assay,  layer = "data")
  mat_data2 <- rbind(mat_data, new_data_row)
  new_assay <- CreateAssayObject(counts = mat_counts2)
  new_assay <- SetAssayData(new_assay,  layer="data", new.data = mat_data2)
  seu[[assay]] <- new_assay
  Idents(seu) <- factor(seu$group, levels = c("NC","AD"))
  DefaultAssay(seu) <- assay
  mat <- GetAssayData(seu, assay = assay, layer = "data")[features, , drop = FALSE]
  meta <- seu@meta.data %>%
    transmute(
      cell = rownames(seu@meta.data),
      group = factor(group, levels = c("NC", "AD")),
      celltype = factor(celltype, levels = celltype_order)
    )
  df_expr <- as.data.frame(as.matrix(t(mat))) %>%
    rownames_to_column("cell") %>%
    left_join(meta, by = "cell") %>%
    pivot_longer(cols = all_of(features), names_to = "feature", values_to = "expr") %>%
    mutate(feature = factor(feature, levels = class_order))
  res_stars <- class_ct %>%
    transmute(
      celltype,
      feature,
      star = p_to_star(p_val_adj)
    ) %>%
    filter(star != "")
  ymax_df <- df_expr %>%
    group_by(feature, celltype) %>%
    summarise(ymax = max(expr, na.rm = TRUE), .groups = "drop")
  anno_df <- res_stars %>%
    inner_join(ymax_df, by = c("feature", "celltype")) %>%
    mutate(y = ymax + 0.15)
  cols <- c("NC" = "#F8766D", "AD" = "#00BFC4")
  if(!is_plot) return(list(df_expr, anno_df))
  ggplot(df_expr, aes(x = celltype, y = expr, fill = group)) +
    geom_violin(
      position = position_dodge(width = 0.85),
      trim = TRUE,
      scale = "width",
      linewidth = 0.25
    ) +
    facet_wrap(~ feature, ncol = 2, scales = "free_y") +
    geom_text(
      data = anno_df,
      aes(x = celltype, y = y, label = star),
      inherit.aes = FALSE,
      size = 5,
      vjust = 0
    ) +
    scale_fill_manual(values = cols, name = "Group") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(x = "Cell type", y = "Expression level") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
    )
}
pdf(file.path(res_dir, "hERV_class_celltype.pdf"), width = 16, height = 10)
print(plot_class_violin(seu, res_by_ct))
dev.off()
```

#### 细胞大类中的hERV家族/位点重调

```r
run_de_one_celltype <- function(seu, ct, assay, filter = FALSE, keep_top_percent = 0.1, expr_pct_threshold = 0.05) {
  if(ct!="All"){
    obj <- subset(seu, subset = celltype == ct)
  }else{
    obj <- seu
  }
  DefaultAssay(obj) <- assay
  if(filter){
    mat <- GetAssayData(seu, assay = assay, slot = "data")
    avg_expr <- Matrix::rowMeans(mat)
    pct_expr <- Matrix::rowMeans(mat > 0)
    cutoff_avg <- quantile(avg_expr, 1 - keep_top_percent, na.rm = TRUE)  # 保留表达量前10%的位点
    keep_loci <- names(avg_expr)[avg_expr >= cutoff_avg & pct_expr >= expr_pct_threshold]  # 至少在5%的细胞表达
    obj[[assay]] <- subset(obj[[assay]], features = keep_loci)
  }
  n_tab <- table(obj$group)
  n_NC <- unname(n_tab["NC"])
  n_AD <- unname(n_tab["AD"])
  Idents(obj) <- factor(obj$group, levels = c("NC","AD"))
  res <- FindMarkers(
    object = obj,
    ident.1 = "AD",
    ident.2 = "NC",
    test.use = "MAST",
    # latent.vars = c("nCount_RNA", "percent_mito"),
    min.pct = 0.05,
    logfc.threshold = 0
  )
  res <- res %>%
    tibble::rownames_to_column("feature") %>%
    mutate(
      celltype = ct,
      n_NC = n_NC,
      n_AD = n_AD
    ) %>%
    relocate(celltype, n_NC, n_AD, feature)
  res
}
celltypes_to_run <- c("Excit", "Oligo", "Astro", "Inhib", "Endo", "Mic", "OPC", "All")
res_list <- lapply(celltypes_to_run, function(ct) run_de_one_celltype(seu, ct, "HERV_family"))
res_list <- Filter(Negate(is.null), res_list)
res_subfamily <- bind_rows(res_list)
res_list <- lapply(celltypes_to_run, function(ct) run_de_one_celltype(seu, ct, "HERV_group"))
res_list <- Filter(Negate(is.null), res_list)
res_group <- bind_rows(res_list)
res_list <- lapply(celltypes_to_run, function(ct) run_de_one_celltype(seu, ct, "HERV", TRUE, keep_top_percent = 1, expr_pct_threshold = 0.01))
res_list <- Filter(Negate(is.null), res_list)
res_locus <- bind_rows(res_list)
write.csv(res_subfamily, file.path(res_dir, "hERV_family_celltype.csv"), row.names = F)
write.csv(res_group, file.path(res_dir, "hERV_group_celltype.csv"), row.names = F)
write.csv(res_locus, file.path(res_dir, "hERV_locus_celltype.csv"), row.names = F)
```

#### 总结hERV位点指标

```r
# 从gtf文件中提取位点染色体位置等信息
get_attr <- function(x, key) {
  pat <- paste0(key, ' "([^"]+)"')
  m <- regexec(pat, x, perl = TRUE)
  rr <- regmatches(x, m)
  out <- rep(NA_character_, length(x))
  ok <- lengths(rr) >= 2
  out[ok] <- vapply(rr[ok], `[`, character(1), 2)
  out
}
lines <- readLines("/public/home/wangtianhao/Desktop/STAR_ref/transcripts.gtf", warn = FALSE)
lines <- lines[!grepl("^\\s*#", lines)]
lines <- lines[nzchar(lines)]
gtf <- fread(
  text = lines,
  sep = "\t",
  header = FALSE,
  quote = "",
  fill = TRUE,
  data.table = TRUE
)
gtf <- gtf[, 1:9]
setnames(gtf, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))
gtf[, gene_id := get_attr(attribute, "gene_id")]
gtf[, locus := get_attr(attribute, "locus")]
gtf[, intModel := get_attr(attribute, "intModel")]
gtf[, repName := get_attr(attribute, "repName")]
gtf[, repFamily := get_attr(attribute, "repFamily")]
gtf[, geneRegion := get_attr(attribute, "geneRegion")]
gtf[, category_gtf := get_attr(attribute, "category")]
gene_dt <- unique(
  gtf[feature == "gene", .(
    locus_id = gene_id,
    chr = seqname,
    start = as.integer(start),
    end = as.integer(end),
    strand,
    locus,
    intModel,
    category_gtf
  )][!is.na(locus_id)],
  by = "locus_id"
)
gene_dt[, length_bp := end - start + 1L]
gene_dt[, coord := paste0(chr, ":", start, "-", end)]
# 从locus2family表中提取位点所属类别信息
family_dt <- fread("/public/home/GENE_proc/wth/metadata/hERV/hERV_locus2family.csv")
family_dt <- unique(family_dt, by = "locus_id")
setnames(family_dt, old = c("group", "class"), new = c("erv_group", "erv_class"))
anno_dt <- merge(gene_dt, family_dt, by = "locus_id", all = TRUE)
anno_dt[, category := fifelse(!is.na(category), category, category_gtf)]
anno_dt[, category_gtf := NULL]
# DE分析
DefaultAssay(seu) <- "HERV"
Idents(seu) <- "group"
de_res <- FindMarkers(
  object = seu,
  ident.1 = "AD",
  ident.2 = "NC",
  assay = "HERV",
  test.use = "MAST",
  min.pct = 0.001,
  logfc.threshold = 0
)
de_dt <- as.data.table(de_res, keep.rownames = "locus_id")[
  , .(locus_id, avg_log2FC, p_val_adj, pct.1, pct.2)
]
de_dt$AD_pct <- de_dt$pct.1
de_dt$NC_pct <- de_dt$pct.2
de_dt$locus_id <- gsub("-", "_", de_dt$locus_id)
de_dt <- de_dt %>% select(locus_id, avg_log2FC, p_val_adj, AD_pct, NC_pct)
# 表达细胞数/比例
cnt_mat <- GetAssayData(seu, assay = "HERV", layer = "counts")
data_mat <- GetAssayData(seu, assay = "HERV", layer = "data")
rownames(cnt_mat) <- gsub("-", "_", rownames(cnt_mat))
rownames(data_mat) <- gsub("-", "_", rownames(data_mat))
meta <- seu@meta.data[colnames(cnt_mat), c("group", "orig.ident", "dataset", "celltype"), drop = FALSE]
expr_bin <- cnt_mat > 0
ad_idx <- meta$group == "AD"
nc_idx <- meta$group == "NC"
expr_dt <- data.table(
  locus_id = rownames(cnt_mat),
  n_expr_cells_total = as.numeric(Matrix::rowSums(expr_bin)),
  expr_prop_total = as.numeric(Matrix::rowSums(expr_bin)) / ncol(cnt_mat),
  n_expr_cells_AD = as.numeric(Matrix::rowSums(expr_bin[, ad_idx, drop = FALSE])),
  expr_prop_AD = as.numeric(Matrix::rowSums(expr_bin[, ad_idx, drop = FALSE])) / sum(ad_idx),
  n_expr_cells_NC = as.numeric(Matrix::rowSums(expr_bin[, nc_idx, drop = FALSE])),
  expr_prop_NC = as.numeric(Matrix::rowSums(expr_bin[, nc_idx, drop = FALSE])) / sum(nc_idx)
)
expr_dt[, expr_unique_group := fifelse(
  n_expr_cells_AD > 0 & n_expr_cells_NC == 0, "AD_only",
  fifelse(n_expr_cells_AD == 0 & n_expr_cells_NC > 0, "NC_only",
          fifelse(n_expr_cells_AD > 0 & n_expr_cells_NC > 0, "both", "none"))
)]
# TPM
len_vec <- anno_dt$length_bp[match(rownames(cnt_mat), anno_dt$locus_id)]
valid_len <- !is.na(len_vec) & len_vec > 0
cnt_valid <- cnt_mat[valid_len, , drop = FALSE]
len_kb <- len_vec[valid_len] / 1000
rpk_mat <- Diagonal(x = 1 / len_kb) %*% cnt_valid
rpk_sum <- Matrix::colSums(rpk_mat)
rpk_sum[rpk_sum == 0] <- 1
tpm_mat <- rpk_mat %*% Diagonal(x = 1e6 / rpk_sum)
tpm_mat <- as.matrix(tpm_mat)
rownames(tpm_mat) <- rownames(cnt_mat)[valid_len]
colnames(tpm_mat) <- colnames(cnt_mat)
# 按orig.ident/dataset/celltype求平均
merge_many <- function(lst) {
  Reduce(function(a, b) merge(a, b, by = "locus_id", all = TRUE, sort = FALSE), lst)
}
mean_by_group <- function(mat, groups, prefix) {
  groups <- as.character(groups)
  keep <- !is.na(groups)
  mat <- mat[, keep, drop = FALSE]
  groups <- groups[keep]
  group_names <- unique(groups)
  res_list <- lapply(group_names, function(g) {
    idx <- which(groups == g)
    if (length(idx) == 1) {
      as.numeric(mat[, idx])
    } else {
      Matrix::rowMeans(mat[, idx, drop = FALSE])
    }
  })
  mean_mat <- do.call(cbind, res_list)
  if (is.null(dim(mean_mat))) {
    mean_mat <- matrix(mean_mat, ncol = 1)
  }
  rownames(mean_mat) <- rownames(mat)
  colnames(mean_mat) <- paste0(prefix, group_names)
  out <- data.frame(
    locus_id = rownames(mean_mat),
    mean_mat,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}
make_avg_block <- function(var_name) {
  grp <- meta[[var_name]]
  x1 <- mean_by_group(cnt_mat, grp, paste0("counts_", var_name, "__"))
  x2 <- mean_by_group(data_mat, grp, paste0("data_", var_name, "__"))
  x3 <- mean_by_group(tpm_mat, grp, paste0("TPM_", var_name, "__"))
  data.table::as.data.table(merge_many(list(x1, x2, x3)))
}
avg_orig_ident <- make_avg_block("orig.ident")
avg_dataset <- make_avg_block("dataset")
avg_celltype <- make_avg_block("celltype")
avg_all <- data.table(
  locus_id = rownames(cnt_mat),
  counts_all_cells = as.numeric(Matrix::rowMeans(cnt_mat)),
  data_all_cells = as.numeric(Matrix::rowMeans(data_mat))
)
avg_all_tpm <- data.table(
  locus_id = rownames(tpm_mat),
  TPM_all_cells = as.numeric(Matrix::rowMeans(tpm_mat))
)
avg_all <- merge(avg_all, avg_all_tpm, by = "locus_id", all = TRUE)
# 合并
summary_dt <- merge_many(list(data.table(locus_id = rownames(cnt_mat)), anno_dt, de_dt, expr_dt)) %>% 
  arrange(locus_id) %>% 
  select(c(
    "locus_id",
    "coord", "strand", "length_bp",  # 染色体位置
    "category", "erv_class", "subfamily_raw",,  # 所属家族
    "avg_log2FC", "p_val_adj", "AD_pct", "NC_pct",  # DE结果
    "n_expr_cells_total", "expr_prop_total",  # 有表达的细胞数和表达比例
    "expr_unique_group"  # 是否唯一表达
  ))
summary_dt <- merge_many(list(summary_dt, avg_all, avg_orig_ident, avg_dataset, avg_celltype))
write.csv(summary_dt, file = file.path(res_dir, "locus_summary.csv"), row.names = F)
```

### 细胞亚群中的hERV家族表达情况

#### 拆分各细胞类群

```r
get_subCelltype <- function(ct, pc.use=1:30, harmony=F, resolution=0.3){
  sub <- subset(seu, subset = celltype == ct)
  DefaultAssay(sub) <- "RNA"
  sub <- FindVariableFeatures(sub, nfeatures = 2000)
  sub <- ScaleData(
    sub,
    vars.to.regress = c("nCount_RNA", "percent_mito"),
    features = VariableFeatures(sub),
  )
  sub <- RunPCA(
    sub,
    npcs = 50,
    features = VariableFeatures(sub),
  )
  pdf(file.path(res_dir, paste0(ct, "_ElbowPlot.pdf")), width = 8, height = 8)
  print(ElbowPlot(sub, ndims = 50))
  dev.off()
  if(harmony){
    sub <- RunHarmony(
      object = sub,
      group.by.vars = "sample_uid",
      reduction.use = "pca",
      dims.use = pc.use,
      assay.use = "RNA",
    )
    sub <- RunUMAP(sub, reduction = "harmony", dims = pc.use)
    sub <- FindNeighbors(sub, reduction = "harmony", dims = pc.use)
  } else {
    sub <- RunUMAP(sub, dims = pc.use)
    sub <- FindNeighbors(sub, dims = pc.use)
  }
  sub <- FindClusters(sub, resolution = c(0.01,0.05,0.1,0.2,0.3,0.4,0.5))
  pdf(file.path(res_dir, paste0(ct, "_umap_sample_group.pdf")), width = 24, height = 8)
  print(DimPlot(sub, group.by = "orig.ident", raster = TRUE) | DimPlot(sub, group.by = "dataset", raster = TRUE) | DimPlot(sub, group.by = "group", raster = TRUE))
  dev.off()
  pdf(file.path(res_dir, paste0(ct, "_clustree.pdf")), width = 12, height = 12)
  print(clustree(sub@meta.data, prefix = "RNA_snn_res."))
  dev.off()
  sub$subcluster <- sub[[paste0("RNA_snn_res.", resolution), drop = TRUE]]
  Idents(sub) <- sub$subcluster
  pdf(file.path(res_dir, paste0(ct, "_umap_cluster.pdf")), width = 12, height = 12)
  print(DimPlot(sub, reduction = "umap", group.by = "subcluster", label = TRUE))
  dev.off()
  return(sub)
}
astro <- get_subCelltype("Astro", pc.use = 1:25, harmony = T, resolution = 0.1)
excit <- get_subCelltype("Excit", pc.use = 1:30, harmony = T, resolution = 0.3)
oligo <- get_subCelltype("Oligo", pc.use = 1:30, harmony = T, resolution = 0.5)
inhib <- get_subCelltype("Inhib", pc.use = 1:30, harmony = T, resolution = 0.2)
endo <- get_subCelltype("Endo", pc.use = 1:20, harmony = T, resolution = 0.01)
mic <- get_subCelltype("Mic", pc.use = 1:25, harmony = T, resolution = 0.3)
OPC <- get_subCelltype("OPC", pc.use = 1:30, harmony = T, resolution = 0.4)
saveRDS(
  astro, 
  file = "/public/home/GENE_proc/wth/my_data/astro.rds", 
  compress = "xz"
)
saveRDS(
  excit, 
  file = "/public/home/GENE_proc/wth/my_data/excit.rds", 
  compress = "xz"
)
saveRDS(
  oligo, 
  file = "/public/home/GENE_proc/wth/my_data/oligo.rds", 
  compress = "xz"
)
saveRDS(
  inhib, 
  file = "/public/home/GENE_proc/wth/my_data/inhib.rds", 
  compress = "xz"
)
saveRDS(
  endo, 
  file = "/public/home/GENE_proc/wth/my_data/endo.rds", 
  compress = "xz"
)
saveRDS(
  mic, 
  file = "/public/home/GENE_proc/wth/my_data/mic.rds", 
  compress = "xz"
)
saveRDS(
  OPC, 
  file = "/public/home/GENE_proc/wth/my_data/OPC.rds", 
  compress = "xz"
)
```

#### 细胞亚群中的家族表达情况

此段为之后分析的共用代码：

```r
# Seurat对象
seu <- readRDS("/public/home/GENE_proc/wth/my_data/all_final.rds")
seu_GSE157827 <- subset(seu, dataset == "GSE157827")
seu_GSE174367 <- subset(seu, dataset == "GSE174367")
astro <- readRDS("/public/home/GENE_proc/wth/my_data/astro.rds")
oligo <- readRDS("/public/home/GENE_proc/wth/my_data/oligo.rds")
OPC <- readRDS("/public/home/GENE_proc/wth/my_data/OPC.rds")
excit <- readRDS("/public/home/GENE_proc/wth/my_data/excit.rds")
inhib <- readRDS("/public/home/GENE_proc/wth/my_data/inhib.rds")
endo <- readRDS("/public/home/GENE_proc/wth/my_data/endo.rds")
mic <- readRDS("/public/home/GENE_proc/wth/my_data/mic.rds")
# 共用函数
my_print <- function(..., sep = ""){  # 输出时间和信息
  formatted_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  str <- paste(as.vector(list(...)), collapse = sep)
  message(paste0(formatted_time, "  ", str))
}
identity_seu <- function(seu, target_col = "celltype"){  # 判断该Seurat对象是哪种细胞类型
  celltype_name <- unique(as.character(seu[[target_col]][, 1]))
  celltype_name <- ifelse(length(celltype_name) == 1, celltype_name[1], "All")
  return(celltype_name)
}
```

```r
seu_list <- list()
seu_list$astro <- astro
seu_list$oligo <- oligo
seu_list$OPC <- opc
seu_list$excit <- excit
seu_list$inhib <- inhib
seu_list$endo <- endo
seu_list$mic <- mic
herv_subcluster_dotplot <- function(seu){
  ct <- unique(seu$celltype)
  if(length(ct) == 1) ct <- ct[1]
  else ct <- "ALL"
  Idents(seu) <- seu$subcluster
  DefaultAssay(seu) <- herv_assay
  DotPlot(
    seu,
    features = plot_subfam,
    group.by = "subcluster",
    assay = herv_assay
  ) + 
    RotatedAxis() + 
    ggtitle(ct) +
    theme(plot.title = element_text(hjust = 0.5))
}
herv_subcluster_pheatmap <- function(seu, is_plot = T){
  ct <- unique(seu$celltype)
  if(length(ct) == 1) ct <- ct[1]
  else ct <- "ALL"
  Idents(seu) <- seu$subcluster
  DefaultAssay(seu) <- herv_assay
  avg_mat <- AverageExpression(
    seu,
    assays = herv_assay,
    features = plot_subfam,
    group.by = "subcluster",
    slot = "data"
  )[[herv_assay]]
  avg_mat_z <- t(scale(t(avg_mat)))
  avg_mat_z[is.na(avg_mat_z)] <- 0
  if(!is_plot) return(avg_mat_z)
  pheatmap::pheatmap(
    avg_mat_z,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 10,
    fontsize_col = 10,
    main = ct,
    border = F,
  )$gtable
}
herv_assay <- "HERV_family"
plot_subfam <- c("HERV3-int", "HERVH-int", "HERVE-int", "HERVL-int", "HERVK-int", "HERVK3-int", "HERVK9-int", "HERVK11-int", "HERVK13-int", "HERVK14-int", "HERVK22-int")
dotplot_list <- lapply(seu_list, herv_subcluster_dotplot)
pheatmap_list <- lapply(seu_list, herv_subcluster_pheatmap)
pdf(file.path(res_dir, "herv_subcluster_family_dotplot.pdf"), width = 24, height = 24)
print(wrap_plots(dotplot_list, ncol = 3))
dev.off()
pdf(file.path(res_dir, "herv_subcluster_family_pheatmap.pdf"), width = 24, height = 24)
wrap_plots(pheatmap_list, ncol = 3)
dev.off()
herv_assay <- "HERV_subfamily_scTE"
plot_subfam <- c("HERV3-int", "HERVH-int", "HERVE-int", "HERVL-int", "HERVK-int", "HERVK3-int", "HERVK9-int", "HERVK11-int", "HERVK13-int", "HERVK14-int", "HERVK22-int")
dotplot_list <- lapply(seu_list, herv_subcluster_dotplot)
pheatmap_list <- lapply(seu_list, herv_subcluster_pheatmap)
pdf(file.path(res_dir, "herv_subcluster_family_scTE_dotplot.pdf"), width = 24, height = 24)
print(wrap_plots(dotplot_list, ncol = 3))
dev.off()
pdf(file.path(res_dir, "herv_subcluster_family_scTE_pheatmap.pdf"), width = 24, height = 24)
print(wrap_plots(pheatmap_list, ncol = 3))
dev.off()
herv_assay <- "HERV_subfamily_LTR_scTE"
plot_subfam <- c("HERV3", "HERVH", "HERVE", "HERVL", "HERVK", "HERVK3", "HERVK9", "HERVK11", "HERVK13", "HERVK14", "HERVK22")
dotplot_list <- lapply(seu_list, herv_subcluster_dotplot)
pheatmap_list <- lapply(seu_list, herv_subcluster_pheatmap)
pdf(file.path(res_dir, "herv_subcluster_family_LTR_scTE_dotplot.pdf"), width = 24, height = 24)
print(wrap_plots(dotplot_list, ncol = 3))
dev.off()
pdf(file.path(res_dir, "herv_subcluster_family_LTR_scTE_pheatmap.pdf"), width = 24, height = 24)
print(wrap_plots(pheatmap_list, ncol = 3))
dev.off()
herv_assay <- "HERV_group"
plot_subfam <- c("ERVL", "ERV1", "ERV3", "HERVE", "HERVF", "HERVH", "HERVI", "HERVK", "HERVL", "HERVW")
dotplot_list <- lapply(seu_list, herv_subcluster_dotplot)
pheatmap_list <- lapply(seu_list, herv_subcluster_pheatmap)
pdf(file.path(res_dir, "herv_subcluster_group_dotplot.pdf"), width = 24, height = 24)
print(wrap_plots(dotplot_list, ncol = 3))
dev.off()
pdf(file.path(res_dir, "herv_subcluster_group_pheatmap.pdf"), width = 24, height = 24)
print(wrap_plots(pheatmap_list, ncol = 3))
dev.off()
```

亚群的AD/NC分布和样本分布

```r
summarize_subcluster_dist <- function(seu, subcluster_col = "subcluster", group_col = "group", sample_col = "orig.ident", digits = 1) {
  ct <- unique(seu$celltype)
  if(length(ct) == 1) ct <- ct[1]
  else ct <- "ALL"
  meta <- seu@meta.data
  df <- meta %>%
    dplyr::select(
      subcluster = all_of(subcluster_col),
      group      = all_of(group_col),
      sample     = all_of(sample_col)
    ) %>%
    mutate(
      subcluster = as.character(subcluster),
      group      = as.character(group),
      sample     = as.character(sample)
    ) %>%
    filter(group %in% c("AD", "NC"))
  overall_group <- df %>%
    count(group, name = "n") %>%
    mutate(pct = n / sum(n))
  overall_ad_pct <- overall_group$pct[overall_group$group == "AD"]
  overall_nc_pct <- overall_group$pct[overall_group$group == "NC"]
  if (length(overall_ad_pct) == 0) overall_ad_pct <- 0
  if (length(overall_nc_pct) == 0) overall_nc_pct <- 0
  cell_n_df <- df %>%
    count(subcluster, name = "cell_n")
  group_stat <- df %>%
    count(subcluster, group, name = "n") %>%
    group_by(subcluster) %>%
    mutate(
      cluster_total = sum(n),
      pct = n / cluster_total
    ) %>%
    ungroup()
  group_df <- group_stat %>%
    select(subcluster, group, pct) %>%
    pivot_wider(
      names_from = group,
      values_from = pct,
      values_fill = 0
    )
  if (!"AD" %in% colnames(group_df)) group_df$AD <- 0
  if (!"NC" %in% colnames(group_df)) group_df$NC <- 0
  group_df <- group_df %>%
    mutate(
      AD_overall_pct = overall_ad_pct,
      NC_overall_pct = overall_nc_pct,
      AD_pct_diff = AD - AD_overall_pct,
      NC_pct_diff = NC - NC_overall_pct
    ) %>%
    mutate(
      AD_pct         = sprintf(paste0("%.", digits, "f%%"), AD * 100),
      NC_pct         = sprintf(paste0("%.", digits, "f%%"), NC * 100),
      AD_overall_pct = sprintf(paste0("%.", digits, "f%%"), AD_overall_pct * 100),
      NC_overall_pct = sprintf(paste0("%.", digits, "f%%"), NC_overall_pct * 100),
      AD_pct_diff    = sprintf(paste0("%+.", digits, "f%%"), AD_pct_diff * 100),
      NC_pct_diff    = sprintf(paste0("%+.", digits, "f%%"), NC_pct_diff * 100)
    ) %>%
    select(
      subcluster,
      AD_pct, NC_pct,
      AD_overall_pct, NC_overall_pct,
      AD_pct_diff, NC_pct_diff
    )
  top_sample_df <- df %>%
    count(subcluster, sample, name = "n") %>%
    group_by(subcluster) %>%
    mutate(pct = n / sum(n)) %>%
    arrange(subcluster, desc(n), desc(pct), sample) %>%
    summarise(
      top_samples = list(sprintf(paste0("%s(%.", digits, "f%%)"), sample, pct * 100)),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      top1_sample = ifelse(length(top_samples) >= 1, top_samples[[1]], NA),
      top2_sample = ifelse(length(top_samples) >= 2, top_samples[[2]], NA),
      top3_sample = ifelse(length(top_samples) >= 3, top_samples[[3]], NA)
    ) %>%
    ungroup() %>%
    select(subcluster, top1_sample, top2_sample, top3_sample)
  res <- cell_n_df %>%
    left_join(group_df, by = "subcluster") %>%
    left_join(top_sample_df, by = "subcluster") %>%
    arrange(subcluster)
  write.csv(res, file = file.path(res_dir, paste0(ct, "_summarize_subcluster.csv")), row.names = F)
}
lapply(seu_list, summarize_subcluster_dist)
```

#### 选定目标亚群进行差异表达分析

```r
subcluster_DE <- function(seu, target_clusters, logfc_threshold=0.1, min_pct=0.1){
  Idents(seu) <- seu$subcluster
  celltype_name <- identity_seu(seu)
  for (cl in target_clusters) {
    my_print(paste0("running ", celltype_name, "-", cl))
    other_ids <- setdiff(levels(Idents(seu)), cl)
    deg <- FindMarkers(
      seu,
      ident.1 = cl,
      ident.2 = other_ids,
      assay = "RNA",
      slot = "data",
      test.use = "MAST",
      # latent.vars = c("nCount_RNA", "percent_mito"),
      min.pct = min_pct,
      logfc.threshold = logfc_threshold
    )
    deg <- subset(deg, p_val_adj < 0.05)
    deg$gene <- rownames(deg)
    deg <- deg %>%
      select(gene, avg_log2FC, p_val_adj) %>%
      arrange(-avg_log2FC)
    deg_up <- subset(deg, avg_log2FC > logfc_threshold) %>%
      arrange(-avg_log2FC)
    write.csv(deg, file.path(res_dir, paste0(celltype_name, "_DEG_subcluster_", cl, ".csv")), row.names = FALSE)
    write.csv(deg_up, file.path(res_dir, paste0(celltype_name, "_DEG_subcluster_", cl, "_up.csv")), row.names = FALSE)
  }
}
subcluster_DE(astro, c(0, 1))
subcluster_DE(oligo, c(3, 5, 9))
subcluster_DE(excit, c(6, 10))  # 6下调
subcluster_DE(inhib, c(2))
subcluster_DE(OPC, c(6, 8))
subcluster_DE(endo, c(1, 4))
subcluster_DE(mic, c(2, 3, 6))  # 2AD多
```

### hERV-gene共表达网络

#### 构建网络

```r
# 此段代码不建议在R控制台中运行，可能出现unexpected symbol的bug，可以在Rstudio或用Rscript命令运行
get_module_size <- function(seu, wgcna_name) {
  df <- GetModules(seu, wgcna_name = wgcna_name)
  gene_col <- intersect(c("gene_name", "gene", "feature"), colnames(df))[1]
  mod_col  <- intersect(c("module", "color"), colnames(df))[1]
  out <- as.data.frame(table(df[[mod_col]]), stringsAsFactors = FALSE)
  colnames(out) <- c("module", "n_genes")
  return(out[order(out$n_genes, decreasing = TRUE), ])
}
construct_network <- function(seu, herv_assay, target_herv, mix=FALSE, filter_fraction=0.1, deep_split=4, min_module_size=20, merge_cut_height=0.05, soft_power=NA){
  celltype_name <- identity_seu(seu)
  plot_name <- paste0(celltype_name, "_", herv_assay, ifelse(mix, "_mix_", "_"), filter_fraction, "_", deep_split, "_", min_module_size, "_", merge_cut_height, "_gene_WGCNA")
  my_print(paste(
    paste0("=====", celltype_name, "====="),
    paste0("- herv assay：", herv_assay),
    paste0("- target herv：", paste(target_herv, collapse = ", ")),
    paste0("- mix herv：", mix),
    paste0("- filterFraction：", filter_fraction),
    paste0("- deepSplit：", deep_split),
    paste0("- minModuleSize：", min_module_size),
    paste0("- mergeCutHeight：", merge_cut_height),
    paste0("- res prefix：", plot_name),
    sep="\n"
  ))
  if(!mix){
    DefaultAssay(seu) <- "RNA"
    seu$AD_binary <- ifelse(seu$group == "AD", 1, 0)
    mat <- LayerData(seu, assay = herv_assay, layer = "data")[target_herv, , drop = FALSE]
    mat <- t(as.matrix(mat))
    colnames(mat) <- paste0("herv_", make.names(target_herv))
    seu <- AddMetaData(seu, metadata = as.data.frame(mat))
    my_print("start SetupForWGCNA")
    seu <- SetupForWGCNA(
      seu,
      gene_select = "fraction",
      fraction = filter_fraction,
      assay = "RNA",
      wgcna_name = "gene_net"
    )
    my_print("start MetacellsByGroups")
    seu <- tryCatch({  # 解决报错`error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': no applicable method for @ applied to an object of class "NULL"`
      MetacellsByGroups(
        seurat_obj = seu,
        group.by = c("celltype", "orig.ident"),
        ident.group = "celltype",
        assay = "RNA",
        layer = "data",
        reduction = "pca",
        dims = 1:30,
        k = 25,
        max_shared = 10,
        min_cells = 50,
        mode = "average",
        wgcna_name = "gene_net"
      )
    }, error = function(e) {
      my_print("change k/max_shared/min_cells...")
      MetacellsByGroups(
        seurat_obj = seu,
        group.by = c("celltype", "orig.ident"),
        ident.group = "celltype",
        assay = "RNA",
        layer = "data",
        reduction = "pca",
        dims = 1:30,
        k = 10,
        max_shared = 5,
        min_cells = 20,
        mode = "average",
        wgcna_name = "gene_net"
      )
    })
    my_print("start NormalizeMetacells")
    seu <- NormalizeMetacells(seu, wgcna_name = "gene_net")
    my_print("start SetDatExpr")
    seu <- SetDatExpr(
      seu,
      group_name = celltype_name,
      group.by = "celltype",
      use_metacells = TRUE,
      assay = "RNA",
      layer = "data",
      wgcna_name = "gene_net"
    )
    my_print("start TestSoftPowers")
    seu <- TestSoftPowers(
      seu,
      networkType = "signed",
      wgcna_name = "gene_net"
    )
    power_table <- GetPowerTable(seu)
    soft_power <- power_table$Power[which(power_table$SFT.R.sq >= 0.8)[1]]
    if (is.na(soft_power)) {
      soft_power <- power_table$Power[which.max(power_table$SFT.R.sq)]
    }
    # 如果soft_power为1，后面基因无法聚类，则ModuleTraitCorrelation会报错
    plot_list <- PlotSoftPowers(seu)
    pdf(file.path(res_dir, paste0(plot_name, "_soft_power.pdf")), width = 16, height = 16)
    p <- patchwork::wrap_plots(plot_list, ncol = 2)
    print(p)
    dev.off()
    my_print("start ConstructNetwork")
    seu <- ConstructNetwork(
      seu,
      soft_power = soft_power,
      networkType = "signed",
      overwrite_tom = TRUE,
      tom_name = paste0(plot_name, "_TOM"),
      deepSplit = deep_split,
      minModuleSize = min_module_size,
      mergeCutHeight = merge_cut_height,
      pamStage = FALSE,
      wgcna_name = "gene_net"
    )
    my_print("start ModuleEigengenes")
    seu <- ModuleEigengenes(
      seu,
      group.by.vars = "orig.ident",
      assay = "RNA",
      wgcna_name = "gene_net"
    )
    my_print("start ModuleTraitCorrelation")
    trait_cols <- c(
      paste0("herv_", make.names(target_herv)),
      "AD_binary"
    )
    seu <- ModuleTraitCorrelation(
      seu,
      traits = trait_cols,
      group.by = "celltype",
      wgcna_name = "gene_net"
    )
    mt <- GetModuleTraitCorrelation(seu)
    pdf(file.path(res_dir, paste0(plot_name, "_trait_corr.pdf")), width = 16, height = 16)
    p <- PlotModuleTraitCorrelation(seu, label = "fdr", label_symbol = "stars")
    print(p)
    dev.off()
    my_print("start ModuleConnectivity")
    seu <- ModuleConnectivity(
      seu,
      group.by = "celltype",
      group_name = celltype_name,
      assay = "RNA",
      wgcna_name = "gene_net"
    )
    my_print("summary...")
    write.csv(get_module_size(seu, "gene_net"), file=file.path(res_dir, paste0(plot_name, "_n_genes.csv")))
    my_print("saving net_obj")
    save(seu, file = file.path(res_dir, paste0(plot_name, ".RData")))
  } else {
    rna_counts <- LayerData(seu, assay = "RNA", layer = "counts")
    feature_rowsum <- Matrix::rowSums(rna_counts > 0)
    retained_f_low <- feature_rowsum > ncol(seu) * filter_fraction
    genes_keep <- names(feature_rowsum)[retained_f_low]
    seu[["RNA"]] <- subset(seu[["RNA"]], features = genes_keep)
    rna_counts <- LayerData(seu, assay = "RNA", layer = "counts")
    rna_data  <- LayerData(seu, assay = "RNA", layer = "data")
    herv_counts <- LayerData(seu, assay = herv_assay, layer = "counts")[target_herv, , drop = FALSE]
    herv_data <- LayerData(seu, assay = herv_assay, layer = "data")[target_herv, , drop = FALSE]
    rownames(herv_counts) <- paste0("HERV__", rownames(herv_counts))
    rownames(herv_data) <- paste0("HERV__", rownames(herv_data))
    mixed_counts <- rbind(rna_counts, herv_counts)
    mixed_data <- rbind(rna_data, herv_data)
    seu[["RNA_HERV"]] <- CreateAssayObject(counts = mixed_counts)
    LayerData(seu, assay = "RNA_HERV", layer = "data") <- mixed_data
    DefaultAssay(seu) <- "RNA_HERV"
    mixed_features <- rownames(seu[["RNA_HERV"]])
    my_print("start SetupForWGCNA")
    seu <- SetupForWGCNA(
      seu,
      gene_select = "custom",
      gene_list = mixed_features,
      assay = "RNA_HERV",
      wgcna_name = "mixed_net"
    )
    my_print("start MetacellsByGroups")
    seu <- tryCatch({
      MetacellsByGroups(
        seurat_obj = seu,
        group.by = c("celltype", "orig.ident"),
        ident.group = "celltype",
        assay = "RNA_HERV",
        layer = "data",
        reduction = "pca",
        dims = 1:30,
        k = 25,
        max_shared = 10,
        min_cells = 50,
        mode = "average",
        wgcna_name = "mixed_net"
      )
    }, error = function(e) {
      my_print("change k/max_shared/min_cells...")
      MetacellsByGroups(
        seurat_obj = seu,
        group.by = c("celltype", "orig.ident"),
        ident.group = "celltype",
        assay = "RNA_HERV",
        layer = "data",
        reduction = "pca",
        dims = 1:30,
        k = 10,
        max_shared = 5,
        min_cells = 20,
        mode = "average",
        wgcna_name = "mixed_net"
      )
    })
    my_print("start NormalizeMetacells")
    seu <- NormalizeMetacells(seu, wgcna_name = "mixed_net")
    my_print("start SetDatExpr")
    seu <- SetDatExpr(
      seu,
      group_name = celltype_name,
      group.by = "celltype",
      use_metacells = TRUE,
      assay = "RNA_HERV",
      layer = "data",
      wgcna_name = "mixed_net"
    )
    print(grep("^HERV__", rownames(seu)))
    my_print("start TestSoftPowers")
    seu <- TestSoftPowers(
      seu,
      networkType = "signed",
      wgcna_name = "mixed_net"
    )
    power_table <- GetPowerTable(seu)
    soft_power <- power_table$Power[which(power_table$SFT.R.sq >= 0.8)[1]]
    if (is.na(soft_power)) {
      soft_power <- power_table$Power[which.max(power_table$SFT.R.sq)]
    }
    plot_list <- PlotSoftPowers(seu)
    pdf(file.path(res_dir, paste0(plot_name, "_soft_power.pdf")), width = 16, height = 16)
    p <- patchwork::wrap_plots(plot_list, ncol = 2)
    print(p)
    dev.off()
    my_print("start ConstructNetwork")
    seu <- ConstructNetwork(
      seu,
      assay = "RNA_HERV",
      soft_power = soft_power,
      networkType = "signed",
      overwrite_tom = TRUE,
      tom_name = paste0(plot_name, "_TOM"),
      deepSplit = deep_split,
      minModuleSize = min_module_size,
      mergeCutHeight = merge_cut_height,
      pamStage = FALSE,
      wgcna_name = "mixed_net"
    )
    my_print("start ModuleEigengenes")
    seu <- ScaleData(
      seu,
      assay = "RNA_HERV",
      vars.to.regress = c("nCount_RNA", "percent_mito"),
      features = rownames(seu[["RNA_HERV"]])
    )
    seu <- ModuleEigengenes(
      seu,
      # group.by.vars = "orig.ident",
      assay = "RNA_HERV",
      wgcna_name = "mixed_net"
    )
    my_print("start ModuleConnectivity")
    seu <- ModuleConnectivity(
      seu,
      group.by = "celltype",
      group_name = celltype_name,
      assay = "RNA_HERV",
      wgcna_name = "mixed_net"
    )
    my_print("summary...")
    write.csv(get_module_size(seu, "mixed_net"), file=file.path(res_dir, paste0(plot_name, "_n_genes.csv")))
    modules_df <- GetModules(seu)
    herv_module_df <- modules_df %>%
      filter(grepl("^HERV--", gene_name)) %>%
      mutate(herv = sub("^HERV--", "", gene_name)) %>%
      select(herv, module)
    write.csv(herv_module_df, file=file.path(res_dir, paste0(plot_name, "_herv_module.csv")), row.names=F)
    my_print("saving net_obj")
    save(seu, file = file.path(res_dir, paste0(plot_name, ".RData")))
  }
}
herv_family_assay <- "HERV_family"
target_herv_family <- c("HERV3-int", "HERVH-int", "HERVE-int", "HERVL-int", "HERVK-int", "HERVK3-int", "HERVK9-int", "HERVK11-int", "HERVK13-int", "HERVK14-int", "HERVK22-int")
herv_group_assay <- "HERV_group"
target_herv_group <- c("HERVK", "HERVH", "HERVE", "ERV3", "HERVF", "HERVW", "HERVL", "ERV1", "ERVL", "HERVI", "HERVP", "HERVS")
construct_network(astro, herv_group_assay, target_herv_group, min_module_size=50, merge_cut_height=0.2)
construct_network(astro, herv_family_assay, target_herv_family, min_module_size=50, merge_cut_height=0.2)
construct_network(oligo, herv_group_assay, target_herv_group)
construct_network(oligo, herv_family_assay, target_herv_family)
construct_network(excit, herv_group_assay, target_herv_group, min_module_size=50, merge_cut_height=0.2)
construct_network(excit, herv_family_assay, target_herv_family, min_module_size=50, merge_cut_height=0.2)
construct_network(inhib, herv_group_assay, target_herv_group, min_module_size=30, merge_cut_height=0.1)
construct_network(inhib, herv_family_assay, target_herv_family, min_module_size=30, merge_cut_height=0.1)
construct_network(OPC, herv_group_assay, target_herv_group, min_module_size=30, merge_cut_height=0.1)
construct_network(OPC, herv_family_assay, target_herv_family, min_module_size=30, merge_cut_height=0.1)
construct_network(endo, herv_group_assay, target_herv_group, min_module_size=30, merge_cut_height=0.1)
construct_network(endo, herv_family_assay, target_herv_family, min_module_size=30, merge_cut_height=0.1)
construct_network(mic, herv_group_assay, target_herv_group, min_module_size=50, merge_cut_height=0.2)
construct_network(mic, herv_family_assay, target_herv_family, min_module_size=50, merge_cut_height=0.2)
```

```r
# 所有细胞类型整合一起画图
my_plot <- function(net_obj_name){
  my_print(paste0(net_obj_name))
  my_print("- reading RData")
  load(file.path(res_dir, net_obj_name))
  my_print("- drawing plot")
  PlotModuleTraitCorrelation(seu, label = "fdr", label_symbol = "stars", combine = FALSE)[[2]]
}
net_obj_name_list <- list(
  "Astro_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData",
  "Oligo_HERV_family_0.1_4_20_0.05_gene_WGCNA.RData",
  "OPC_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData",
  "Excit_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData",
  "Inhib_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData",
  "Endo_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData",
  "Mic_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData"
)
plot_list <- lapply(net_obj_name_list, my_plot)
pdf(file.path(res_dir, "HERV_family_subcluster_network_corr.pdf"), width = 32, height = 24)
print(wrap_plots(plot_list, ncol = 3))
dev.off()
net_obj_name_list <- list(
  "Astro_HERV_group_0.1_4_50_0.2_gene_WGCNA.RData",
  "Oligo_HERV_group_0.1_4_20_0.05_gene_WGCNA.RData",
  "OPC_HERV_group_0.1_4_30_0.1_gene_WGCNA.RData",
  "Excit_HERV_group_0.1_4_50_0.2_gene_WGCNA.RData",
  "Inhib_HERV_group_0.1_4_30_0.1_gene_WGCNA.RData",
  "Endo_HERV_group_0.1_4_30_0.1_gene_WGCNA.RData",
  "Mic_HERV_group_0.1_4_50_0.2_gene_WGCNA.RData"
)
plot_list <- lapply(net_obj_name_list, my_plot)
pdf(file.path(res_dir, "HERV_group_subcluster_network_corr.pdf"), width = 32, height = 24)
print(wrap_plots(plot_list, ncol = 3))
dev.off()
```

#### 计算每个细胞类型的各亚群的模块分数

```r
calc_module_score <- function(seu_obj_name, net_obj_name, modules_use = NA){
  my_print(paste0(seu_obj_name, " ", net_obj_name))
  my_print("- reading RData")
  load(file.path(res_dir, net_obj_name))
  seu_net <- seu
  seu <- readRDS(file.path("/public/home/GENE_proc/wth/my_data/", seu_obj_name))
  celltype_name <- identity_seu(seu)
  res_name <- paste0(celltype_name, "_subcluster_module_scores")
  my_print("- AddModuleScore...")
  # 取出模块基因
  modules_df <- GetModules(seu_net, wgcna_name = "gene_net")
  modules_df <- subset(modules_df, module != "grey")
  if(length(modules_use)==1) if(is.na(modules_use)) modules_use <- as.character(unique(modules_df$module))
  modules_df <- subset(modules_df, module %in% modules_use)
  module_list <- split(modules_df$gene_name, modules_df$module)
  module_list <- lapply(module_list, unique)
  module_list <- lapply(module_list, intersect, rownames(seu[["RNA"]]))
  # 给每个细胞打模块分数
  for (m in modules_use) {
    seu <- AddModuleScore(
      seu,
      features = list(module_list[[m]]),
      assay = "RNA",
      name = paste0(m, "_score_")
    )
  }
  score_cols <- setNames(paste0(modules_use, "_score_1"), modules_use)
  Idents(seu) <- seu[["subcluster"]][, 1]
  # 各亚群的模块平均分数
  score_avg <- seu@meta.data[, c("subcluster", unname(score_cols)), drop = FALSE]
  colnames(score_avg)[1] <- "subcluster"
  score_avg <- score_avg %>%
    group_by(subcluster) %>%
    summarise(across(everything(), mean))
  write.csv(score_avg, file.path(res_dir, paste0(res_name, ".csv")), row.names = FALSE)
  score_mat <- as.matrix(score_avg[, -1])
  rownames(score_mat) <- score_avg$subcluster
  colnames(score_mat) <- names(score_cols)
  pdf(file.path(res_dir, paste0(res_name, ".pdf")), width = 8, height = 8)
  p <- pheatmap(
    scale(score_mat),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 10,
    fontsize_col = 10,
    main = celltype_name,
    border = F
  )
  draw(p, heatmap_legend_side = "right", padding = unit(c(5, 5, 5, 5), "mm"))
  dev.off()
  return(as.ggplot(p))
}
# 注：使用hERV group和family的RData都一样，因为这个只影响hERV~gene的相关性，它们的模块都是一样的
astro <- calc_module_score("astro.rds", "Astro_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData")
oligo <- calc_module_score("oligo.rds", "Oligo_HERV_family_0.1_4_20_0.05_gene_WGCNA.RData")
excit <- calc_module_score("excit.rds", "Excit_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData")
inhib <- calc_module_score("inhib.rds", "Inhib_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData")
OPC <- calc_module_score("OPC.rds", "OPC_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData")
endo <- calc_module_score("endo.rds", "Endo_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData")
mic <- calc_module_score("mic.rds", "Mic_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData")
plot_list <- list(astro, oligo, OPC, excit, inhib, endo, mic)
pdf(file.path(res_dir, "subcluster_ModuleScore_pheatmap.pdf"), width = 24, height = 24)
print(wrap_plots(plot_list, ncol = 3))
dev.off()
```

#### 对每个模块进行GO分析

```r
run_go_one <- function(gene_vec, universe_genes = NULL, ont = "BP") {
  gene_vec <- unique(na.omit(gene_vec))
  my_print(paste0("  - n_gene:", length(gene_vec), " n_universe_gene:", length(universe_genes)))
  ego <- enrichGO(
    gene          = gene_vec,
    universe      = universe_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
  if (is.null(ego)){
    my_print("  - no res")
    return(NULL)
  }
  ego
}
mat_to_df <- function(mat, value_name) {
  df <- as.data.frame(as.table(as.matrix(mat)), stringsAsFactors = FALSE)
  colnames(df) <- c("trait", "module", value_name)
  df
}
run_go <- function(net_obj_name, target_modules = NA, ont = "BP", show_category = 15, n_hubs = 30) {
  my_print(net_obj_name)
  my_print("- reading RData")
  load(file.path(res_dir, net_obj_name))
  celltype_name <- identity_seu(seu)
  modules_df <- GetModules(seu, wgcna_name = "gene_net")
  # 背景基因设置为网络中的所有基因
  universe_genes <- unique(modules_df$gene_name)
  universe_genes <- universe_genes[!is.na(universe_genes)]
  modules_df <- subset(modules_df, module != "grey")
  if(length(target_modules)==1){
    if(is.na(target_modules)){
      celltype_name <- paste0(celltype_name, "_all")
      target_modules <- as.character(unique(modules_df$module))
    }
  }
  # 目标模块的基因
  module_gene_list <- lapply(target_modules, function(m) {
    unique(modules_df[modules_df$module == m, "gene_name", drop = TRUE])
  })
  names(module_gene_list) <- target_modules
  go_res_plot_list <- list()
  go_res_list <- list()
  for (m in target_modules) {
    my_print(paste0("- go module: ", m))
    ego <- run_go_one(
      gene_vec = module_gene_list[[m]],
      universe_genes = universe_genes,
      ont = ont
    )
    if (is.null(ego)) next
    go_df <- as.data.frame(ego)
    go_res_list[[m]] <- go_df
    res_name <- paste0(celltype_name, "_GO_", ont, "_", m)
    write.csv(
      go_df,
      file.path(res_dir, paste0(res_name, ".csv")),
      row.names = FALSE
    )
    go_res_plot_list[[m]] <- dotplot(ego, showCategory = show_category) +
      ggtitle(res_name)
  }
  pdf(
    file.path(res_dir, paste0(celltype_name, "_GO_", ont, ".pdf")), 
    width = 15, 
    height = ((length(target_modules)-1)%/%2+1)*8
  )
  print(wrap_plots(go_res_plot_list, ncol = 2))
  dev.off()
  # 取hub genes
  hub_df <- GetHubGenes(seu, n_hubs = n_hubs)
  module_hub_list <- lapply(target_modules, function(m) {
    unique(hub_df[hub_df$module == m, "gene_name", drop = TRUE])
  })
  names(module_hub_list) <- target_modules
  # 汇总模块相关性和GO的结果
  mt <- GetModuleTraitCorrelation(seu)
  cor_df <- mat_to_df(mt$cor$all_cells,  "cor")
  fdr_df <- mat_to_df(mt$fdr$all_cells,  "fdr")
  sig_tbl <- left_join(cor_df, fdr_df, by = c("module", "trait")) %>%
    arrange(module, desc(abs(cor)))
  top_go_tbl <- lapply(target_modules, function(m) {
    ego <- go_res_list[[m]]
    if (is.null(ego)) return(NULL)
    go_df <- as.data.frame(ego)
    go_df <- go_df[order(go_df$p.adjust), ]
    trait_df <- sig_tbl %>%
      filter(module == m) %>%
      arrange(desc(abs(cor)))
    data.frame(
      module = m,
      trait = paste(trait_df$trait, collapse = "; "),
      cor = paste(round(trait_df$cor, 3), collapse = "; "),
      n_gene = length(module_gene_list[[m]]),
      hub_gene = paste(module_hub_list[[m]], collapse = "; "),
      GO = paste(go_df$Description[1:30], collapse = "; "),
      GO_padj = paste(go_df$p.adjust[1:30], collapse = "; "),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  write.csv(top_go_tbl, file.path(res_dir, paste0(celltype_name, "_module_trait_GO_", ont, "_summary.csv")), row.names = FALSE)
}
run_go("Astro_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData", c("blue", "red"))
run_go("Oligo_HERV_family_0.1_4_20_0.05_gene_WGCNA.RData", c("turquoise", "yellow", "blue"))
run_go("Excit_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData", c("blue", "green", "magenta", "pink", "black", "red"))
run_go("Inhib_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData", c("brown"))
run_go("OPC_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData", c("turquoise", "yellow"))
run_go("Endo_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData", c("blue"))
run_go("Mic_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData", c("blue", "turquoise", "brown", "pink"))
run_go("Mic_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData")  # mic是否有炎症相关模块
run_go("Astro_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData")  # 有没有免疫相关模块
run_go("Oligo_HERV_family_0.1_4_20_0.05_gene_WGCNA.RData")
run_go("Excit_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData")
run_go("Inhib_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData")
run_go("OPC_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData")
run_go("Endo_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData")
```

#### 亚群的DEG与模块基因取交集的GO分析

```r
run_DEG_go <- function(net_obj_name, target_module_map, up_or_down = "up", ont = "BP", show_category = 15, is_save_go_res = F) {
  if(length(up_or_down)==1) up_or_down <- rep(up_or_down, length(target_module_map))
  my_print(net_obj_name)
  my_print("- reading RData")
  load(file.path(res_dir, net_obj_name))
  celltype_name <- identity_seu(seu)
  modules_df <- GetModules(seu, wgcna_name = "gene_net")
  universe_genes <- unique(modules_df$gene_name)
  universe_genes <- universe_genes[!is.na(universe_genes)]
  modules_df <- subset(modules_df, module != "grey")
  res_df <- list()
  res_plot <- list()
  for (i in 1:length(target_module_map)) {
    cl <- names(target_module_map)[i]
    deg_all <- read.csv(file.path(res_dir, paste0(celltype_name, "_DEG_subcluster_", cl, ".csv")))
    if(up_or_down[i]=="up"){
      deg <- filter(deg_all, avg_log2FC>0) %>% pull(gene)
    } else {
      deg <- filter(deg_all, avg_log2FC<0) %>% pull(gene)
    }
    target_module <- target_module_map[i]
    target_modules_gene <- modules_df %>% 
      filter(module == target_module) %>% 
      pull(gene_name)
    gene_use <- intersect(deg, target_modules_gene)
    ego <- run_go_one(gene_use, universe_genes, ont)
    if (is.null(ego)) return(NULL)
    res_plot[[paste0(cl, "_", target_module)]] <- dotplot(ego, showCategory = show_category) +
      ggtitle(paste0(celltype_name, " cluster", cl, " DEG - ", target_module, " gene GO_", ont))
    go_df <- as.data.frame(ego)
    go_df <- go_df[order(go_df$p.adjust), ]
    if(is_save_go_res){
      res_name <- paste0(celltype_name, "_DEG-GO_", ont, "_", target_module)
      write.csv(
        go_df,
        file.path(res_dir, paste0(res_name, ".csv")),
        row.names = FALSE
      )
    }
    res_df[[paste0(cl, "_", target_module)]] <- data.frame(
      cluster = cl,
      module = target_module,
      up_or_down = up_or_down[i],
      overlap_n = length(gene_use),
      overlap_ratio = length(gene_use)/length(target_modules_gene),
      n_module_gene = length(target_modules_gene),
      n_universe_gene = length(universe_genes),
      GO = paste(go_df$Description[1:30], collapse = "; "),
      GO_padj = paste(go_df$p.adjust[1:30], collapse = "; "),
      gene_use = paste(gene_use, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  res_name <- paste(celltype_name, "DEG_GO", ont, sep = "_")
  pdf(
    file.path(res_dir, paste0(res_name, ".pdf")), 
    width = 15, 
    height = ((length(target_module_map)-1)%/%2+1)*8
  )
  print(wrap_plots(res_plot, ncol = 2))
  dev.off()
  write.csv(bind_rows(res_df), file.path(res_dir, paste0(res_name, "_summary.csv")), row.names = FALSE)
}
run_DEG_go("Astro_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData", c(
  "0" = "blue",
  "1" = "red"
))
run_DEG_go("Oligo_HERV_family_0.1_4_20_0.05_gene_WGCNA.RData", c(
  "5" = "turquoise",
  "3" = "yellow",
  "9" = "blue"
))
run_DEG_go("Excit_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData", c(
  "6" = "black",
  "6" = "red",
  "6" = "pink",
  "6" = "green"
), up_or_down = c("up", "up", "down", "down"))
run_DEG_go("Mic_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData", c(
  "2" = "turquoise",
  "2" = "brown",
  "2" = "pink",
  "2" = "blue"
), up_or_down = c("up", "up", "up", "down"), is_save_go_res = T)
```

#### 相关性分析补充

```r
# 相当于提取了构建共表达网络的ModuleTraitCorrelation部分
herv_gene_module_corr <- function(seu, herv_assay, target_herv, is_plot = T){
  celltype_name <- identity_seu(seu)
  DefaultAssay(seu) <- "RNA"
  seu$AD_binary <- ifelse(seu$group == "AD", 1, 0)
  mat <- LayerData(seu, assay = herv_assay, layer = "data")[target_herv, , drop = FALSE]
  mat <- t(as.matrix(mat))
  colnames(mat) <- paste0("herv_", make.names(target_herv))
  seu <- AddMetaData(seu, metadata = as.data.frame(mat))
  trait_cols <- c(
    paste0("herv_", make.names(target_herv)),
    "AD_binary"
  )
  seu <- ModuleTraitCorrelation(
    seu,
    traits = trait_cols,
    group.by = "celltype",
    wgcna_name = "gene_net"
  )
  mt <- GetModuleTraitCorrelation(seu)
  # pdf(file.path(res_dir, paste(celltype_name, herv_assay, "trait_corr.pdf", sep = "_")), width = 16, height = 16)
  # p <- PlotModuleTraitCorrelation(seu, label = "fdr", label_symbol = "stars")
  # print(p)
  # dev.off()
  if(is_plot) return(PlotModuleTraitCorrelation(seu, label = "fdr", label_symbol = "stars", combine = FALSE)[[2]])
  else return(seu)
}
get_plot <- function(net_obj_name_list, target_herv_assay, target_herv){
  plot_list <- lapply(net_obj_name_list, function(net_obj_name){
    my_print(net_obj_name)
    my_print("- reading RData")
    load(file.path(res_dir, net_obj_name))
    my_print("- calc corr")
    herv_gene_module_corr(seu, target_herv_assay, target_herv, T)
  })
  pdf(file.path(res_dir, paste0(target_herv_assay, "_subcluster_network_corr.pdf")), width = 32, height = 24)
  print(wrap_plots(plot_list, ncol = 3))
  dev.off()
  return(plot_list)
}
target_herv_locus <- c("HERV3-1p34.3", "HML3-1q32.1b", "HML1-1q32.1", "HML4-16p13.3", "MER41-13q13.3", "MER61-21q21.1a", "MER61-21q21.1b", "HERVW-15q21.2", "HARLEQUIN-6p21.31")
target_herv_family <- c("HERV3-int", "HERVH-int", "HERVE-int", "HERVL-int", "HERVK-int", "HERVK3-int", "HERVK9-int", "HERVK11-int", "HERVK13-int", "HERVK14-int", "HERVK22-int")
target_herv_family_LTR <- c("HERV3", "HERVH", "HERVE", "HERVL", "HERVK", "HERVK3", "HERVK9", "HERVK11", "HERVK13", "HERVK14", "HERVK22")
net_obj_name_list <- list(
  "scTE_Astro_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData",
  "scTE_Oligo_HERV_family_0.1_4_20_0.05_gene_WGCNA.RData",
  "scTE_OPC_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData",
  "scTE_Excit_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData",
  "scTE_Inhib_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData",
  "scTE_Endo_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData",
  "scTE_Mic_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData"
)
get_plot(net_obj_name_list, "HERV", target_herv_locus)
get_plot(net_obj_name_list, "HERV_family", target_herv_family)
get_plot(net_obj_name_list, "HERV_subfamily_scTE", target_herv_family)
get_plot(net_obj_name_list, "HERV_subfamily_LTR_scTE", target_herv_family_LTR)
```

#### 位点级hERV与gene相关性

```r
summary_locus_WGCNA <- function(seu, fdr_sh = 0.2, is_get_cor_mat = F){
  ct <- identity_seu(seu)
  wgcna_name <- "gene_net"
  my_print("- calc herv~module cor")
  seu <- herv_gene_module_corr(seu, herv_locus_assay, target_herv_locus, F)
  my_print("- summary cor res")
  mt <- GetModuleTraitCorrelation(seu)
  cor_df <- mat_to_df(mt$cor$all_cells, "cor")
  fdr_df <- mat_to_df(mt$fdr$all_cells, "fdr")
  if(is_get_cor_mat) return(list(cor_df, fdr_df))
  pval_df <- mat_to_df(mt$pval$all_cells, "pval")
  mt_tbl <- cor_df %>%
    left_join(pval_df, by = c("trait", "module")) %>%
    left_join(fdr_df, by = c("trait", "module")) %>%
    mutate(celltype = ct) %>%
    select(celltype, trait, module, cor, pval, fdr)
  locus_map <- tibble(
    locus = rownames(GetAssayData(seu, assay = herv_locus_assay, layer = "data")),
    trait = paste0("herv_", make.names(locus))
  )
  pair_tbl <- mt_tbl %>%
    filter(grepl("^herv_", trait), module != "grey", fdr < fdr_sh) %>%
    inner_join(locus_map, by = "trait") %>%
    group_by(celltype, locus) %>%
    arrange(desc(abs(cor)), fdr) %>%
    slice_head(n = 2) %>%  # 取每个“细胞类型-位点”保留相关性绝对值最大的前2个module
    ungroup() %>%
    select(celltype, locus, trait, module, cor, pval, fdr)
  my_print("- get hub genes")
  hub_tbl <- GetHubGenes(seu, n_hubs = 50, wgcna_name = wgcna_name)
  hub_sel <- pair_tbl %>%
    select(celltype, locus, module, cor, fdr) %>%
    distinct() %>%
    inner_join(hub_tbl, by = "module") %>%
    arrange(celltype, locus, module)
  # 在对应celltype内、对应module里，找与该位点最相关的具体基因
  modules_tbl <- GetModules(seu, wgcna_name = wgcna_name) %>%
    select(gene_name, module) %>%
    distinct() %>%
    filter(module != "grey")
  gene_cor_list <- list()
  my_print("- calc herv~gene cor")
  for (i in seq_len(nrow(pair_tbl))) {
    ct <- pair_tbl$celltype[i]
    locus_i <- pair_tbl$locus[i]
    module_i <- pair_tbl$module[i]
    x <- as.numeric(GetAssayData(seu, assay = herv_locus_assay, layer = "data")[locus_i, ])
    genes_i <- modules_tbl$gene_name[modules_tbl$module == module_i]
    genes_i <- intersect(genes_i, rownames(GetAssayData(seu, assay = "RNA", layer = "data")))
    if (length(genes_i) == 0) next
    expr_i <- GetAssayData(seu, assay = "RNA", layer = "data")[genes_i, , drop = FALSE]
    cp <- WGCNA::corAndPvalue(
      x = t(as.matrix(expr_i)),
      y = x,
      method = "spearman"
    )
    tmp <- tibble(
      celltype = ct,
      locus = locus_i,
      module = module_i,
      gene = genes_i,
      rho = as.numeric(cp$cor),
      p = as.numeric(cp$p)
    ) %>%
      mutate(fdr = p.adjust(p, method = "fdr")) %>%
      arrange(desc(abs(rho)), fdr)
    gene_cor_list[[i]] <- tmp
  }
  gene_cor_tbl <- bind_rows(gene_cor_list)
  my_print("- summary...")
  top_gene_tbl <- gene_cor_tbl %>%  # 每个“位点-模块”里相关性最强的前50个基因
    filter(fdr < fdr_sh) %>%
    group_by(celltype, locus, module) %>%
    slice_max(order_by = abs(rho), n = 50, with_ties = FALSE) %>%
    ungroup()
  top_hub_tbl <- hub_sel %>%  # 每个“位点-模块”里前50个hub gene
    group_by(celltype, locus, module) %>%
    slice_head(n = 50) %>%
    ungroup()
  key_gene_tbl <- top_gene_tbl %>%  # 交集基因
    inner_join(
      top_hub_tbl %>% select(celltype, locus, module, gene_name),
      by = c("celltype", "locus", "module", "gene" = "gene_name")
    ) %>%
    arrange(celltype, locus, module, desc(abs(rho)), fdr)
  key_gene_tbl
}
net_obj_name_list <- list(
  "Astro_HERV_group_0.1_4_50_0.2_gene_WGCNA.RData",
  "Oligo_HERV_group_0.1_4_20_0.05_gene_WGCNA.RData",
  "OPC_HERV_group_0.1_4_30_0.1_gene_WGCNA.RData",
  "Excit_HERV_group_0.1_4_50_0.2_gene_WGCNA.RData",
  "Inhib_HERV_group_0.1_4_30_0.1_gene_WGCNA.RData",
  "Endo_HERV_group_0.1_4_30_0.1_gene_WGCNA.RData",
  "Mic_HERV_group_0.1_4_50_0.2_gene_WGCNA.RData"
)
herv_locus_assay <- "HERV"
target_herv_locus <- c("HERV3-1p34.3", "HML3-1q32.1b", "HML1-1q32.1", "HML4-16p13.3", "MER41-13q13.3", "MER61-21q21.1a", "MER61-21q21.1b", "HERVW-15q21.2", "HARLEQUIN-6p21.31")
res <- lapply(net_obj_name_list, function(net_obj_name){
  my_print(net_obj_name)
  my_print("- loading RData")
  load(file.path(res_dir, net_obj_name))
  summary_locus_WGCNA(seu)
}) %>% bind_rows()  
write.csv(res, file.path(res_dir, "locus_module_keygenes.csv"), row.names = FALSE)
```

### 线性拟合计算相关性

#### 家族级

```r
# 计算残差
residualize_vec <- function(y, meta) {
  fit <- lm(y ~ log10_nCount_RNA + percent_mito, data = meta)
  resid(fit)
}
# 线性拟合计算相关性，herv_features可以传多个，gene_features如果传多个会根据agg方法聚合成一个分数来计算
run_herv_geneset_coexpr_one_celltype <- function(seu_ct, herv_assay, herv_features, gene_features, agg = "mean", is_donor = FALSE, donor_var = "orig.ident") {
  ct <- identity_seu(seu_ct)
  level_name <- ifelse(is_donor, "donor", "cell")
  res_name <- paste(ct, herv_features[1], gene_features[1], agg, level_name, sep = " ")
  my_print(res_name, ":")
  herv_use <- intersect(rownames(seu_ct[[herv_assay]]), herv_features)
  gene_use <- intersect(rownames(seu_ct[["RNA"]]), gene_features)
  my_print("- hERV features used:", length(herv_use), " genes used:", length(gene_use))
  if(length(herv_use)==0 | length(gene_use)==0) return(NULL)
  my_print("- calculating residualize_vec")
  # 构建模型，设置协变量
  meta <- seu_ct@meta.data %>%
    mutate(log10_nCount_RNA = log10(nCount_RNA + 1))
  herv_mat <- GetAssayData(seu_ct, assay = herv_assay, layer = "data")[herv_use, , drop = FALSE]
  rna_mat  <- GetAssayData(seu_ct, assay = "RNA",  layer = "data")[gene_use, , drop = FALSE]
  herv_counts <- GetAssayData(seu_ct, assay = herv_assay, layer = "counts")[herv_use, , drop = FALSE]
  rna_counts  <- GetAssayData(seu_ct, assay = "RNA", layer = "counts")[gene_use, , drop = FALSE]
  herv_pct <- Matrix::rowMeans(herv_counts > 0)
  gene_set_pct <- mean(Matrix::colSums(rna_counts > 0) > 0)
  if (agg == "mean") {
    gene_score <- Matrix::colMeans(rna_mat)
  } else {
    gene_score <- Matrix::colSums(rna_mat)
  }
  gene_set_pct <- mean(  # 该组基因中任意一个表达>0即记为检测到
    Matrix::colSums(GetAssayData(seu_ct, assay = "RNA", layer = "counts")[gene_use, , drop = FALSE] > 0) > 0
  )
  # 残差化
  if (!is_donor) {
    my_print("- cell-level residualize")
    herv_resid <- t(apply(as.matrix(herv_mat), 1, residualize_vec, meta = meta))
    gene_resid <- residualize_vec(as.numeric(gene_score), meta = meta)
    n_used <- ncol(herv_resid)
  } else {
    my_print("- donor-level aggregation")
    donors <- as.character(meta[[donor_var]])
    donor_levels <- unique(donors)
    meta_donor <- meta %>%  # 协变量：对每个样本取平均
      mutate(donor = donors) %>%
      group_by(donor) %>%
      summarise(
        log10_nCount_RNA = mean(log10_nCount_RNA, na.rm = TRUE),
        percent_mito = mean(percent_mito, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      as.data.frame()
    rownames(meta_donor) <- meta_donor$donor
    herv_donor <- sapply(donor_levels, function(d) {  # hERV表达：每个样本内按细胞平均
      idx <- which(donors == d)
      Matrix::rowMeans(herv_mat[, idx, drop = FALSE])
    })
    if (is.vector(herv_donor)) {
      herv_donor <- matrix(herv_donor, nrow = length(herv_use), dimnames = list(herv_use, donor_levels))
    }
    colnames(herv_donor) <- donor_levels
    gene_donor <- sapply(donor_levels, function(d) {  # gene set表达：每个样本内按细胞平均
      idx <- which(donors == d)
      mean(gene_score[idx], na.rm = TRUE)
    })
    gene_donor <- as.numeric(gene_donor)
    names(gene_donor) <- donor_levels
    meta_donor <- meta_donor[colnames(herv_donor), , drop = FALSE]
    my_print("- donor-level residualize")
    herv_resid <- t(apply(as.matrix(herv_donor), 1, residualize_vec, meta = meta_donor))
    gene_resid <- residualize_vec(as.numeric(gene_donor), meta = meta_donor)
    n_used <- ncol(herv_resid)
  }
  my_print("- calculating corr")
  # 计算rho+p值
  res_list <- vector("list", length = nrow(herv_resid))
  for (i in seq_len(nrow(herv_resid))) {
    x <- as.numeric(herv_resid[i, ])
    y <- as.numeric(gene_resid)
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    if (length(x) < 10 || sd(x) == 0 || sd(y) == 0) {
      res_list[[i]] <- data.frame(
        feature = rownames(herv_resid)[i],
        gene_set = paste(gene_use, collapse = ";"),
        n_gene = length(gene_use),
        agg = agg,
        level = level_name,
        rho = NA_real_,
        p_val = NA_real_,
        n = length(x),
        herv_pct = herv_pct[rownames(herv_resid)[i]],
        gene_set_pct = gene_set_pct
      )
    } else {
      ct_res <- suppressWarnings(
        cor.test(x, y, method = "spearman", exact = FALSE)
      )
      res_list[[i]] <- data.frame(
        feature = rownames(herv_resid)[i],
        gene_set = paste(gene_use, collapse = ";"),
        n_gene = length(gene_use),
        agg = agg,
        level = level_name,
        rho = unname(ct_res$estimate),
        p_val = ct_res$p.value,
        n = length(x),
        herv_pct = herv_pct[rownames(herv_resid)[i]],
        gene_set_pct = gene_set_pct
      )
    }
  }
  bind_rows(res_list) %>%
    mutate(
      p_adj = p.adjust(p_val, method = "BH")
    ) %>%
    arrange(p_adj, p_val, desc(abs(rho)))
}
get_gene_herv_corr_res <- function(seu, herv_assay, target_herv, gene_panel, res_name = NA, agg = "mean", is_donor = FALSE, donor_var = "orig.ident", write_csv = T){
  if(is.na(res_name)) res_name <- identity_seu(seu)
  loop_tag <- names(gene_panel)
  if(is.null(loop_tag)) loop_tag <- 1:length(gene_panel)
  res <- lapply(loop_tag, function(sig) {
    row <- run_herv_geneset_coexpr_one_celltype(
      seu_ct = seu,
      herv_assay = herv_assay,
      herv_features = target_herv,
      gene_features = gene_panel[[sig]],
      agg = agg,
      is_donor = is_donor,
      donor_var = donor_var
    )
    if(is.null(row)) return(NULL)
    else return(row %>% mutate(signature = sig, .before = 1))
  }) %>% bind_rows()
  if(write_csv) write.csv(res, file = file.path(res_dir, paste0(res_name, "_geneset_", herv_assay, ifelse(is_donor, "_donor", ""), "_corr.csv")), row.names = F)
  else return(res)
}
```

**Mic的family-level**

```r
gene_panel_mic <- list(
  innate_immune = c("TYROBP", "CD74", "HLA-DRA"),
  inflammation_signal = c("TLR2", "MAP3K8", "EPSTI1"),
  stress_immune = c("HSP90AA1", "HSPH1", "DNAJB1"),
  homeostatic_microglia = c("P2RY12", "CX3CR1", "TMEM119"),
  IFN_STAT1 = c("STAT1", "ISG15", "IFIT3"),
  TLR8 = c("TLR8")
)
target_herv_family <- c("HERVK-int", "HERVK3-int", "HERVK9-int", "HERVK11-int", "HERVK13-int", "HERVK14-int", "HERVK22-int")
# Mic
get_gene_herv_corr_res(mic, "HERV_subfamily_scTE", target_herv_family, gene_panel_mic)
# 全体细胞
get_gene_herv_corr_res(seu, "HERV_subfamily_scTE", target_herv_family, gene_panel_mic)
# mic的cluster2
get_gene_herv_corr_res(subset(mic, subcluster=="2"), "HERV_subfamily_scTE", target_herv_family, gene_panel_mic, res_name = "mic2")
```

#### 位点级

**各细胞类型的locus-level**

```r
get_gene_herv_corr_res_plus <- function(seu, herv_assay, gene_panel_list, res_name = NA, agg = "mean", is_donor = FALSE, donor_var = "orig.ident", write_csv = T){
  if(is.na(res_name)) res_name <- identity_seu(seu)
  res <- lapply(names(gene_panel_list), function(herv){
    get_gene_herv_corr_res(
      seu,
      herv_assay = herv_assay,
      target_herv = herv,
      gene_panel = gene_panel_list[[herv]],
      res_name = res_name,
      agg = agg,
      is_donor = is_donor,
      donor_var = donor_var, 
      write_csv = F
    )
  }) %>% bind_rows()
  if(write_csv) write.csv(res, file = file.path(res_dir, paste0(res_name, "_geneset_", herv_assay, ifelse(is_donor, "_donor", ""), "_corr_plus.csv")), row.names = F)
  else return(res)
}
mic_target <- list(
  "HARLEQUIN-6p21.31" = list("IL4R", "STAT3", "HCK", "ADGRE2", "FGD2", "ALOX5", "ADAP2", "SH2B3"),
  "HERVW-15q21.2" = list("ITGAX", "LPAR5", "PLD1", "RASGEF1C", "FGD2", "ADAP2")
)
get_gene_herv_corr_res_plus(mic, "HERV", mic_target)
get_gene_herv_corr_res_plus(mic, "HERV", mic_target, is_donor = T)
oligo_target <- list(
  "HERVW-15q21.2" = list("GLDN", "FA2H", "TMEM63A", "DSCAML1", "PIK3C2B", "MICAL3"),
  "MER61-21q21.1a" = list("TMEM63A", "PIK3C2B", "MICAL3", "FNBP1", "GLDN", "FA2H"),
  "MER41-13q13.3" = list("GLDN", "TMEM63A", "FA2H", "ATP11A", "MICAL3", "PIK3C2B"),
  "HARLEQUIN-6p21.31" = list("PIK3C2B", "CORO2B", "SLC13A3", "TMEM63A", "GLDN")
)
get_gene_herv_corr_res_plus(oligo, "HERV", oligo_target)
get_gene_herv_corr_res_plus(oligo, "HERV", oligo_target, is_donor = T)
astro_target <- list(
  "HML3-1q32.1b" = list("STAT3", "PLSCR1", "MAP3K14", "OSMR", "COLEC12", "ABCC9", "TPST1"),
  "HML1-1q32.1" = list("STAT3", "PLSCR1", "MAP3K14", "COLEC12", "ABCC9", "TPST1", "OSMR"),
  "HML4-16p13.3" = list("NPAS3", "STAT3", "PLSCR1", "MAP3K14", "COLEC12")
)
get_gene_herv_corr_res_plus(astro, "HERV", astro_target)
get_gene_herv_corr_res_plus(astro, "HERV", astro_target, is_donor = T)
excit_target <- list(
  "HERV3-1p34.3" = list("THY1", "NCDN", "FAIM2", "SYNPO", "CAMK2B", "NRGN"),
  "HML4-16p13.3" = list("MTR", "TTC8", "MICU2", "TMEM241", "SPATA7"),
  "HML3-1q32.1b" = list("MTR", "TTC8", "MICU2", "SPATA7", "TMEM241"),
  "HML1-1q32.1" = list("GOLGA8B", "WSB1", "DYNC2LI1", "SCFD1", "MTR")
)
get_gene_herv_corr_res_plus(excit, "HERV", excit_target)
get_gene_herv_corr_res_plus(excit, "HERV", excit_target, is_donor = T)
inhib_target <- list(
  "HARLEQUIN-6p21.31" = list("COL4A1", "COL4A2", "AJAP1", "CRTAC1", "HAPLN1", "TOX2"),
  "HML4-16p13.3" = list("COL4A1", "COL4A2", "AJAP1", "CRTAC1", "HAPLN1", "TOX2"),
  "HML3-1q32.1b" = list("RASSF8", "CTXND1", "LIN7A", "TACR1", "GOLIM4", "AJAP1"),
  "HML1-1q32.1" = list("RASSF8", "LIN7A", "CRTAC1", "COL4A1", "AJAP1", "COL4A2")
)
get_gene_herv_corr_res_plus(inhib, "HERV", inhib_target)
get_gene_herv_corr_res_plus(inhib, "HERV", inhib_target, is_donor = T)
opc_target <- list(
  "HML3-1q32.1b" = list("TMEM63A", "PIK3C2B", "SLC13A3"),
  "HERVK11-2p12" = list("TMEM63A", "PIK3C2B", "SLC13A3")
)
get_gene_herv_corr_res_plus(opc, "HERV", opc_target)
get_gene_herv_corr_res_plus(opc, "HERV", opc_target, is_donor = T)
endo_target <- list(
  "HARLEQUIN-6p21.31" = list("STAT3", "PLSCR1", "MAP3K14"),
  "HML1-1q32.1" = list("COLEC12", "ABCC9", "TPST1")
)
get_gene_herv_corr_res_plus(endo, "HERV", endo_target)
get_gene_herv_corr_res_plus(endo, "HERV", endo_target, is_donor = T)
```

**Mic2的HERVK位点与TLR8+IFN_STAT1**

```r
summary_HERVK_locus <- function(seu_ct, herv_locus_assay = "HERV"){
  ct <- identity_seu(seu_ct)
  all_loci <- rownames(seu_ct[[herv_locus_assay]])
  hml_loci <- all_loci[grepl("^(HML|HERVK11)", all_loci)]
  my_print("- HML loci found:", length(hml_loci))
  herv_counts <- GetAssayData(seu_ct, assay = herv_locus_assay, layer = "counts")[hml_loci, , drop = FALSE]
  locus_pct <- Matrix::rowMeans(herv_counts > 0)
  locus <- names(locus_pct)
  data.frame(
    locus = locus,
    pct = as.numeric(locus_pct),
    family = case_when(
      grepl("^HML1", locus) ~ "HERVK14",
      grepl("^HML2", locus) ~ "HERVK",
      grepl("^HML3", locus) ~ "HERVK9",
      grepl("^HML4", locus) ~ "HERVK13",
      grepl("^HML5", locus) ~ "HERVK22",
      grepl("^HML6", locus) ~ "HERVK3",
      grepl("^HERVK11", locus) ~ "HERVK11",
      TRUE ~ NA_character_
    ),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(pct), locus)
}
mic2 <- subset(mic, subcluster=="2")
target_herv_locus <- summary_HERVK_locus(mic2) %>%
  filter(pct > 0.001) %>%
  pull(locus)
gene_panel <- list(
  IFN_STAT1 = c("STAT1", "ISG15", "IFIT3"),
  TLR8 = c("TLR8")
)
get_gene_herv_corr_res(mic2, "HERV", target_herv_locus, gene_panel, res_name = "mic2")
get_gene_herv_corr_res(mic2, "HERV", target_herv_locus, gene_panel, res_name = "mic2", is_donor = T)
```

#### 画图展示

**单个基因/基因集-单个hERV的相关性点图**

```r
plot_herv_gene_corr <- function(seu_ct, herv_assay, herv_feature, gene, res_name = NA, agg = "mean", is_donor = FALSE, donor_var = "orig.ident", is_plot = T) {
  if(is.na(res_name)) res_name <- identity_seu(seu_ct)
  my_print("calc corr: ", herv_feature[1], "-", gene[1])
  meta <- seu_ct@meta.data %>%
    mutate(log10_nCount_RNA = log10(nCount_RNA + 1))
  # hERV
  herv_mat <- GetAssayData(seu_ct, assay = herv_assay, layer = "data")
  if (!herv_feature %in% rownames(herv_mat)) stop("hERV feature not found in assay.")
  x_raw <- as.numeric(herv_mat[herv_feature, ])
  # gene/gene set
  rna_mat <- GetAssayData(seu_ct, assay = "RNA", layer = "data")
  if (length(gene) == 1 && gene %in% colnames(meta)) {
    y_raw <- meta[[gene]]
    gene_name <- gene
  } else {
    gene_use <- intersect(gene, rownames(rna_mat))
    if (length(gene_use) == 0) stop("gene not found in RNA assay.")
    y_mat <- rna_mat[gene_use, , drop = FALSE]
    y_raw <- if (agg == "mean") Matrix::colMeans(y_mat) else Matrix::colSums(y_mat)
    gene_name <- if (length(gene_use) == 1) gene_use else paste(gene_use, collapse = "+")
  }
  # cell-level/donor-level
  if (!is_donor) {
    x <- residualize_vec(x_raw, meta)
    y <- residualize_vec(as.numeric(y_raw), meta)
    df_plot <- data.frame(x = x, y = y)
  } else {
    donors <- as.character(meta[[donor_var]])
    meta_donor <- meta %>%
      mutate(donor = donors) %>%
      group_by(donor) %>%
      summarise(
        log10_nCount_RNA = mean(log10_nCount_RNA, na.rm = TRUE),
        percent_mito = mean(percent_mito, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      as.data.frame()
    rownames(meta_donor) <- meta_donor$donor
    donor_levels <- meta_donor$donor
    x_donor <- sapply(donor_levels, function(d) mean(x_raw[donors == d], na.rm = TRUE))
    y_donor <- sapply(donor_levels, function(d) mean(y_raw[donors == d], na.rm = TRUE))
    x <- residualize_vec(as.numeric(x_donor), meta_donor)
    y <- residualize_vec(as.numeric(y_donor), meta_donor)
    df_plot <- data.frame(x = x, y = y, donor = donor_levels)
  }
  ct_res <- suppressWarnings(cor.test(df_plot$x, df_plot$y, method = "spearman", exact = FALSE))
  my_print("- drawing")
  rho_txt <- round(unname(ct_res$estimate), 3)
  p_txt <- signif(ct_res$p.value, 3)
  if(!is_plot) return(list(df_plot, rho_txt, p_txt))
  ggplot(df_plot, aes(x = x, y = y)) +
    geom_point(size = 0.5, alpha = 0.35) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.7) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste0("rho=", rho_txt, "\np=", p_txt),
      hjust = -0.1, vjust = 1.2,
      size = 4.5
    ) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    labs(
      title = paste0(res_name, ":    ", gene_name, " ~ ", herv_feature),
      x = herv_feature,
      y = gene_name
    )
}
```

**同一细胞类型中多个基因/基因集（列表）与多个hERV（vector）**

```r
plot_herv_gene_corr_panel <- function(seu_ct, herv_assay, herv_feature, gene_list, res_name = NA, agg = c("mean", "sum"), is_donor = FALSE, donor_var = "orig.ident", ncol = length(herv_feature), is_plot = T) {
  plots <- list()
  for(gene in gene_list){
    for(herv in herv_feature){
      plots[[paste0(gene[1], "_", herv)]] <- plot_herv_gene_corr(seu_ct, herv_assay, herv, gene, res_name, agg, is_donor, donor_var)
    }
  }
  if(is_plot){
    pdf(
      file.path(res_dir, "herv_gene_corr_panel.pdf"), 
      width = 6*ncol, 
      height = ((length(plots)-1)%/%4+1)*6
    )
    print(wrap_plots(plots, ncol = ncol))
    dev.off()
  }
  else return(plots)
}
mic2 <- subset(mic, subcluster=="2")
gene_list <- list("TLR8")
top_tlr8_loci <- c(
  "HML2-3q12.3",
  "HML6-19q13.43b",
  "HML1-17q25.1",
  "HML2-22q11.23",
  "HML3-1q32.1b",
  "HML2-1q22"
)
plot_herv_gene_corr_panel(mic2, "HERV", top_tlr8_loci, gene_list, agg = "mean", is_donor = TRUE, ncol = 3)
```

**多细胞类型组合图**

```r
plot_herv_gene_corr_df <- function(seu_list, plot_specs, ncol = 4, res_name = "", title = "", is_plot = F){
  plots <- vector("list", nrow(plot_specs))
  for(i in seq_len(nrow(plot_specs))) {
    ct_i <- plot_specs$ct[i]
    herv_i <- plot_specs$herv[i]
    gene_i <- plot_specs$gene[i]
    p <- plot_herv_gene_corr(
      seu_ct = seu_list[[ct_i]],
      herv_assay = "HERV",
      herv_feature = herv_i,
      gene = gene_i,
      res_name = ifelse(ct_i=="Mic2", "Mic cluster2", NA),
      is_donor = TRUE,
    ) +
      geom_point(size = 1, alpha = 0.85) +
      labs(x = NULL, y = NULL) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(size = 13.5, face = "bold"),
        axis.title = element_blank(),
        axis.text = element_text(size = 9)
      )
    plots[[i]] <- p
  }
  res <- wrap_plots(plots, ncol = ncol)
  if(title != "") res <- res + plot_annotation(title = title)
  if(!is_plot) return(res)
  pdf(file.path(res_dir, paste0("herv_gene_donor_corr", res_name, ".pdf")), width = 24, height = ((length(plots)-1)%/%4+1)*6)
  print(res)
  dev.off()
}
seu_list <- list(
  Mic = mic,
  Oligo = oligo,
  Astro = astro,
  Excit = excit
)
plot_specs <- tibble::tribble(
  ~ct,     ~herv,                 ~gene,
  "Mic",   "HARLEQUIN-6p21.31",   "IL4R",
  "Mic",   "HARLEQUIN-6p21.31",   "ADGRE2",
  "Mic",   "HARLEQUIN-6p21.31",   "ALOX5",
  "Mic",   "HARLEQUIN-6p21.31",   "HCK",
  "Oligo", "HERVW-15q21.2",       "GLDN",
  "Oligo", "HERVW-15q21.2",       "PIK3C2B",
  "Oligo", "MER41-13q13.3",       "ATP11A",
  "Astro", "HML1-1q32.1",         "ABCC9",
  "Astro", "HML3-1q32.1b",        "ABCC9",
  "Astro", "HML3-1q32.1b",        "COLEC12",
  "Excit", "HML4-16p13.3",        "MTR",
  "Excit", "HML4-16p13.3",        "MICU2",
)
p_corr <- plot_herv_gene_corr_df(seu_list, plot_specs)
```

**气泡图**：

```r
plot_herv_gene_bubble <- function(cor_df, p_cutoff = 0.05) {
  df <- cor_df %>%
    filter(!is.na(rho), p_adj < p_cutoff) %>%
    group_by(signature) %>%
    arrange(abs(rho)) %>%
    mutate(
      sig_value = -log10(.data[["p_adj"]]),
      sig_value = ifelse(is.infinite(sig_value), max(sig_value[is.finite(sig_value)], na.rm = TRUE) + 1, sig_value)
    )
  df$signature <- factor(df$signature)
  df$subfamily <- factor(df$subfamily)
  ggplot(df, aes(x = signature, y = subfamily, size = abs(rho), color = rho)) +
    geom_point(aes(size = sig_value, color = rho), alpha = 0.9) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red") +
    theme_bw() +
    labs(
      x = NULL,
      y = NULL,
      color = "rho",
      size = "-log10(FDR)"
    )
}
cor_df <- read.csv(file.path(res_dir, paste0("mic2_geneset_hERV_locus_corr_donor.csv")))
pdf(file.path(res_dir, "herv_gene_corr_panel_dotplot.pdf"), width = 16, height = 16)
print(plot_herv_gene_bubble(cor_df, 0.15))
dev.off()
```

### 其它结果图

#### 质控+聚类图

```r
load("/public/home/GENE_proc/wth/my_data/all_final.RData"); seu=all; rm(all)
seu_before_qc <- readRDS("/public/home/GENE_proc/wth/my_data/all.rds")
dataset_cols <- c(
  "GSE157827" = "#F8766D",
  "GSE174367" = "#00BFC4"
)
plot_tag_fz <- 18
theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = 14, color = "black"),
      axis.text = element_text(size = 11, color = "black"),
      strip.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.tag = element_text(size = plot_tag_fz, face = "bold"),
    )
}
sample_order <- NULL
qc_cutoff_df <- bind_rows(
  tibble(
    dataset = "GSE157827",
    metric = c("nFeature_RNA", "nFeature_RNA", "nCount_RNA", "nCount_RNA", "percent_mito", "percent_ribo"),
    cutoff = c(200, 8000, 500, 40000, 14, 3),
    type = c("low", "high", "low", "high", "high", "high")
  ),
  tibble(
    dataset = "GSE174367",
    metric = c("nFeature_RNA", "nFeature_RNA", "nCount_RNA", "nCount_RNA", "percent_mito", "percent_ribo"),
    cutoff = c(300, 12000, 1000, 50000, 10, 2),
    type = c("low", "high", "low", "high", "high", "high")
  )
)
# 过滤前QC分布+阈值线
qc_metrics <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo")
meta_before <- seu_before_qc@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  select(cell, dataset, all_of(qc_metrics)) %>%
  mutate(dataset = factor(dataset, levels = names(dataset_cols)))
qc_long <- meta_before %>%
  pivot_longer(
    cols = all_of(qc_metrics),
    names_to = "metric",
    values_to = "value"
  )
qc_cutoff_plot <- qc_cutoff_df %>%
  mutate(
    dataset = factor(dataset, levels = names(dataset_cols)),
    metric = factor(metric, levels = qc_metrics),
    x = c(1,1,1,1,1,1, 2,2,2,2,2,2)
  ) %>%
  mutate(
    x0 = x - 0.32,
    x1 = x + 0.32
  )
pA <- ggplot(qc_long, aes(x = dataset, y = value, fill = dataset)) +
  geom_violin(scale = "width", trim = TRUE, color = "grey20", linewidth = 0.3) +
  geom_segment(
    data = qc_cutoff_plot,
    aes(x = x0, xend = x1, y = cutoff, yend = cutoff, linetype = type),
    inherit.aes = FALSE,
    color = "firebrick",
    linewidth = 0.6
  ) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = dataset_cols) +
  scale_linetype_manual(values = c("low" = "dashed", "high" = "solid")) +
  labs(x = NULL, y = NULL, linetype = "Cutoff") +
  theme_pub(base_size = 11) +
  theme(
    legend.position = "right",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(0.8, "lines")
  )
# 各样本QC前后细胞数对比
count_before <- seu_before_qc@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  count(dataset, orig.ident, name = "n_before")
count_after <- seu@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  count(dataset, orig.ident, name = "n_after")
count_df <- full_join(count_before, count_after, by = c("dataset", "orig.ident")) %>%
  mutate(
    n_before = ifelse(is.na(n_before), 0, n_before),
    n_after = ifelse(is.na(n_after), 0, n_after),
    retain_rate = ifelse(n_before == 0, NA, n_after / n_before)
  )
if (is.null(sample_order)) {
  sample_order <- count_df %>%
    distinct(dataset, orig.ident) %>%
    arrange(dataset, orig.ident) %>%
    pull(orig.ident)
}
count_df <- count_df %>%
  mutate(
    dataset = factor(dataset, levels = names(dataset_cols)),
    orig.ident = factor(orig.ident, levels = sample_order)
  ) %>%
  arrange(dataset, orig.ident)
label_pad <- 0.015 * max(count_df$n_before, na.rm = TRUE)
pB <- ggplot(count_df, aes(x = orig.ident)) +
  geom_col(aes(y = n_before), fill = "grey85", width = 0.9, color = "grey50", linewidth = 0.2) +
  geom_col(aes(y = n_after, fill = dataset), width = 0.65, color = "grey20", linewidth = 0.2) +
  geom_text(
    data = count_df %>% mutate(label_y = n_after + label_pad),
    aes(y = label_y, label = percent(retain_rate, accuracy = 1)),
    angle = 90, vjust = 0.5, hjust = 0, size = 3, color = "black"
  ) +
  facet_grid(~dataset, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = dataset_cols) +
  coord_cartesian(clip = "off") +
  labs(x = "Sample", y = "Number of cells") +
  theme_pub(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    strip.text = element_text(size = 15, face = "bold"),
    panel.spacing.x = unit(0.5, "lines")
  )
# 整合后UMAP按dataset着色
pC <- DimPlot(seu, group.by = "dataset", label = FALSE, raster = TRUE) +
  labs(title = NULL) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.text = element_text(size = 13),
    plot.tag = element_text(size = plot_tag_fz, face = "bold"),
    plot.tag.position = c(-0.03, 0.992),
  )
p_left <- (pA + labs(tag = "A")) / (pB + labs(tag = "B"))
p_final <- p_left | (pC + labs(tag = "C")) + 
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(20, 10, 10, 40)
  )
ggsave(
  file.path(res_dir, "qc.pdf"),
  p_final, width = 18, height = 12
)
```

#### hERV计数情况

```r
celltype_order <- c("Astro", "Endo", "Excit", "Inhib", "Mic", "Oligo", "OPC")
dataset_cols <- c(
  "GSE157827" = "#F8766D",
  "GSE174367" = "#00BFC4"
)
class_cols <- c(
  "ERVL" = "#0FA08A",
  "ERV1" = "#56B1C7",
  "ERVK" = "#F05A49"
)
top_n_subfamily <- 10
plot_tag_fz <- 20
plot_margin <- theme(
  plot.margin = margin(20, 10, 20, 20)  # 上右下左
)
theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.tag = element_text(size = plot_tag_fz, face = "bold"),
    )
}
meta <- seu@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    dataset = factor(dataset, levels = names(dataset_cols)),
    celltype = factor(celltype, levels = celltype_order)
  )
# 总体分布
df_A <- meta %>%
  select(dataset, nCount_HERV, nFeature_HERV, HERV_fraction) %>%
  pivot_longer(
    cols = c(nCount_HERV, nFeature_HERV, HERV_fraction),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("nCount_HERV", "nFeature_HERV", "HERV_fraction"),
      labels = c("nCount_HERV", "nFeature_HERV", "HERV_fraction")
    )
  )
p_A <- ggplot(df_A, aes(x = dataset, y = value, fill = dataset)) +
  geom_violin(scale = "width", trim = TRUE, color = "grey20", linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.3) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = dataset_cols) +
  labs(x = NULL, y = NULL, tag = "A") +
  theme_pub() +
  theme(
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 13), 
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12), 
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )
# 计算某assay中herv组成比例
make_comp_df <- function(seu, assay, top_n = NULL) {
  mat <- GetAssayData(seu, assay = assay, layer = "counts")
  df <- tibble(
    feature = rownames(mat),
    count = Matrix::rowSums(mat)
  ) %>%
    filter(count > 0) %>%
    arrange(desc(count))
  if(!is.null(top_n)) {
    top_feat <- head(df$feature, top_n)
    df <- df %>%
      mutate(feature = ifelse(feature %in% top_feat, feature, "Others")) %>%
      group_by(feature) %>%
      summarise(count = sum(count), .groups = "drop") %>%
      arrange(desc(count))
  }
  if("Others" %in% df$feature) {
    lev <- c(setdiff(df$feature, "Others"), "Others")
  } else {
    lev <- df$feature
  }
  df %>%
    mutate(
      assay = assay,
      feature = factor(feature, levels = lev),
      prop = count / sum(count),
      label = ifelse(prop >= 0.05, percent(prop, accuracy = 0.1), "")
    )
}
plot_comp_bar <- function(df, fill_cols, title, legend_title, show_y = TRUE) {
  p <- ggplot(df, aes(x = "All cells", y = prop, fill = feature)) +
    geom_col(width = 0.7, color = "grey20", linewidth = 0.3) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5), size = 4, color = "black"
    ) +
    scale_fill_manual(values = fill_cols) +
    scale_y_continuous(
      labels = label_percent(accuracy = 1),
      limits = c(0, 1), expand = c(0, 0)
    ) +
    labs(
      x = NULL,
      y = if(show_y) "Proportion of total HERV counts" else NULL,
      fill = legend_title
    ) +
    theme_pub(11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 12), 
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12), 
    )
  if(!show_y) {
    p <- p + theme(axis.title.y = element_blank())
  }
  p
}
df_C <- make_comp_df(seu, assay = "HERV_class_scTE", top_n = NULL)
df_D <- make_comp_df(seu, assay = "HERV_subfamily_scTE", top_n = top_n_subfamily)
subfamily_names <- setdiff(as.character(df_D$feature), "Others")
subfamily_cols <- setNames(scales::hue_pal()(length(subfamily_names)), subfamily_names)
if("Others" %in% as.character(df_D$feature)) {
  subfamily_cols <- c(subfamily_cols, "Others" = "#B3B3B3")
}
# HERV_class_scTE在全体细胞中的组成比例
p_C <- plot_comp_bar(
  df = df_C,
  fill_cols = class_cols,
  title = "HERV_class_scTE",
  legend_title = "Class",
  show_y = TRUE
) + labs(tag = "B")
# HERV_subfamily_scTE在全体细胞中的组成比例
p_D <- plot_comp_bar(
  df = df_D,
  fill_cols = subfamily_cols,
  title = "HERV_subfamily_scTE",
  legend_title = "Subfamily",
  show_y = FALSE
) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(
    legend.text = element_text(size = 8)
  )
# 拼图
p_final <- p_A / (p_C | p_D) +
  plot_annotation(
    theme = theme(
      plot.tag.position = c(0, 1),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(10, 5, 20, 15) 
    )
  )
ggsave(
  file.path(res_dir, "herv_counts.pdf"),
  p_final, width = 12, height = 12
)
```

#### 各细胞类型的hERV表达量和DE显著性

```r
my_heatmap <- function(seu, herv_assay, de_res_path){
  dataset <- identity_seu(seu, "dataset")
  mat <- GetAssayData(seu, assay = herv_assay, layer = "data")
  celltype_order <- c("Astro", "Oligo", "OPC", "Excit", "Inhib", "Mic", "Endo")
  group_order <- c("AD", "NC")
  meta <- seu@meta.data[, c("celltype", "group"), drop = FALSE]
  mat <- mat[, rownames(meta), drop = FALSE]
  meta$celltype <- factor(as.character(meta$celltype), levels = celltype_order)
  meta$group <- factor(as.character(meta$group), levels = group_order)
  mat <- mat[, rownames(meta), drop = FALSE]
  meta$cellgroup <- paste(meta$celltype, meta$group, sep = "_")
  col_order <- expand.grid(
    celltype = celltype_order,
    group = group_order,
    stringsAsFactors = FALSE
  )
  col_order$cellgroup <- paste(col_order$celltype, col_order$group, sep = "_")
  col_order <- col_order[order(
    match(col_order$celltype, celltype_order),
    match(col_order$group, group_order)
  ), , drop = FALSE]
  col_order <- col_order[col_order$cellgroup %in% meta$cellgroup, , drop = FALSE]
  avg_list <- lapply(col_order$cellgroup, function(x) {
    cells <- rownames(meta)[meta$cellgroup == x]
    Matrix::rowMeans(mat[, cells, drop = FALSE])
  })
  avg_mat <- do.call(cbind, avg_list)
  colnames(avg_mat) <- col_order$cellgroup
  rownames(avg_mat) <- rownames(mat)
  row_mean <- rowMeans(avg_mat)
  row_sd <- apply(avg_mat, 1, sd)
  min_mean <- 0.05
  top_n <- 40
  keep <- row_mean >= min_mean
  stat_df <- data.frame(
    family = rownames(avg_mat),
    mean_expr = row_mean,
    sd_expr = row_sd,
    stringsAsFactors = FALSE
  )
  stat_df <- stat_df[keep, ]
  stat_df <- stat_df[order(-stat_df$sd_expr, -stat_df$mean_expr), ]
  family_use <- head(stat_df$family, top_n)
  plot_mat <- avg_mat[family_use, , drop = FALSE]
  plot_mat_z <- t(scale(t(as.matrix(plot_mat))))
  plot_mat_z <- plot_mat_z[apply(plot_mat_z, 1, function(x) all(is.finite(x))), , drop = FALSE]
  peak_col <- max.col(plot_mat_z, ties.method = "first")
  peak_val <- apply(plot_mat_z, 1, max)
  ord <- order(peak_col, -peak_val)
  plot_mat_z <- plot_mat_z[ord, , drop = FALSE]
  family_anno <- read.table("/public/home/GENE_proc/wth/metadata/hERV/LTR2Class.tsv", header=T, sep="\t")
  family_anno$repName <- gsub("_", "-", family_anno$repName)
  mat_rowname <- rownames(plot_mat_z)
  if(!endsWith(mat_rowname[1], "-int")) mat_rowname <- paste0(mat_rowname, "-int")
  row_class <- family_anno$repClass[match(mat_rowname, family_anno$repName)]
  row_class[is.na(row_class)] <- "Unknown"
  class_cols <- c(
    "ERV1" = "#E6E6E6",
    "ERVK" = "#969696",
    "ERVL" = "#252525",
    "Unknown" = "#BDBDBD"
  )
  right_anno <- rowAnnotation(
    Class = row_class,
    col = list(Class = class_cols),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 0),
    simple_anno_size = unit(2, "mm")
  )
  celltype_cols <- setNames(hcl.colors(length(celltype_order), "Set 2"), celltype_order)
  group_cols <- c("NC" = "#2B6CB0", "AD" = "#B22222")
  col_order$celltype <- factor(col_order$celltype, levels = celltype_order)
  col_order$group <- factor(col_order$group, levels = group_order)
  top_anno <- HeatmapAnnotation(
    Celltype = col_order$celltype,
    Group = col_order$group,
    col = list(
      Celltype = celltype_cols,
      Group = group_cols
    ),
    annotation_legend_param = list(
      Celltype = list(
        at = celltype_order
      ),
      Group = list(
        at = group_order
      )
    ),
    annotation_name_gp = gpar(fontsize = 0),
    simple_anno_size = unit(2, "mm")
  )
  de_df <- read.csv(de_res_path)
  de_res <- de_df %>%
    filter(celltype %in% celltype_order) %>%
    select(celltype, feature, avg_log2FC, p_val_adj)
  de_res$feature <- gsub("_", "-", de_res$feature)
  de_res <- de_res[, c("celltype", "feature", "p_val_adj")]
  de_res <- de_res[de_res$celltype %in% celltype_order, ]
  de_res$star <- ""
  de_res$star[de_res$p_val_adj < 0.05] <- "*"
  de_res$star[de_res$p_val_adj < 0.01] <- "**"
  de_res$star[de_res$p_val_adj < 0.001] <- "***"
  sig_df <- pivot_wider(
    de_res,
    id_cols = feature,
    names_from = celltype,
    values_from = star
  )
  sig_mat_ct <- as.matrix(sig_df[, -1])
  rownames(sig_mat_ct) <- sig_df$feature
  sig_mat_ct[is.na(sig_mat_ct)] <- ""
  sig_mat_ct <- sig_mat_ct[rownames(plot_mat_z), , drop = FALSE]
  sig_mat <- matrix(
    "",
    nrow = nrow(plot_mat_z),
    ncol = ncol(plot_mat_z),
    dimnames = dimnames(plot_mat_z)
  )
  for (ct in celltype_order) {
    ad_col <- paste0(ct, "_AD")
    nc_col <- paste0(ct, "_NC")
    if (!(ct %in% colnames(sig_mat_ct))) next
    if (!(ad_col %in% colnames(plot_mat_z))) next
    if (!(nc_col %in% colnames(plot_mat_z))) next
    stars <- sig_mat_ct[rownames(plot_mat_z), ct]
    ad_vals <- plot_mat_z[, ad_col]
    nc_vals <- plot_mat_z[, nc_col]
    idx_ad <- stars != "" & ad_vals > nc_vals
    sig_mat[idx_ad, ad_col] <- stars[idx_ad]
    idx_nc <- stars != "" & nc_vals > ad_vals
    sig_mat[idx_nc, nc_col] <- stars[idx_nc]
  }
  ht <- Heatmap(
    plot_mat_z,
    name = "Z-score",
    col = colorRamp2(c(-2, 0, 2), c("#3B4CC0", "white", "#B40426")),
    top_annotation = top_anno,
    right_annotation = right_anno,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 10),
    border = FALSE,
    rect_gp = gpar(col = NA),
    # column_title = paste0("Average expression of hERV families"),
    # column_title_gp = gpar(fontsize = 16),
    cell_fun = function(j, i, x, y, w, h, fill) {
      lab <- sig_mat[i, j]
      if (lab != "") {
        grid.text(
          lab,
          x = x, y = y,
          gp = gpar(fontsize = 8, col = "black", fontface = "bold")
        )
      }
    }
  )
  pdf(file.path(res_dir, paste0(herv_assay, "_celltype_avgexp_heatmap.pdf")), width = 12, height = 8)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", padding = unit(c(5, 10, 5, 5), "mm"))
  dev.off()
}
my_heatmap(seu, "HERV_subfamily_scTE", "/public/home/GENE_proc/wth/res/hERV_family_scTE_celltype.csv")
```

#### DE结果图

```r
celltype_order <- c("All", "Astro", "Endo", "Excit", "Inhib", "Mic", "Oligo", "OPC")
celltype_order_noall <- c("Astro", "Endo", "Excit", "Inhib", "Mic", "Oligo", "OPC")
class_order <- c("ERV1", "ERVK", "ERVL", "HERV-total")
# 准备样本级family_DE数据
summary_family_donor <- function(seu, features, celltypes = c("Excit", "Inhib", "Mic", "Oligo"), is_plot = T){
  mat <- GetAssayData(seu, assay = "HERV_subfamily_scTE", layer = "data")[features, , drop = FALSE]
  meta <- seu@meta.data %>%
    transmute(
      cell = rownames(seu@meta.data),
      donor = .data[["orig.ident"]],
      group = factor(.data[["group"]], levels = c("NC", "AD")),
      celltype = factor(.data[["celltype"]], levels = celltype_order_noall)
    ) %>%
    filter(celltype %in% celltypes)
  df <- as.data.frame(as.matrix(t(mat))) %>%
    rownames_to_column("cell") %>%
    inner_join(meta, by = "cell") %>%
    pivot_longer(cols = all_of(features), names_to = "feature", values_to = "expr") %>%
    group_by(celltype, donor, group, feature) %>%
    summarise(
      expr = mean(expr, na.rm = TRUE),
      n_cell = n(),
      .groups = "drop"
    )
  stat_df <- df %>%
    group_by(celltype, feature) %>%
    summarise(
      n_NC = sum(group == "NC"),
      n_AD = sum(group == "AD"),
      p = tryCatch(
        wilcox.test(expr ~ group, exact = FALSE)$p.value,
        error = function(e) NA_real_
      ),
      delta = median(expr[group == "AD"], na.rm = TRUE) -
        median(expr[group == "NC"], na.rm = TRUE),
      y = max(expr, na.rm = TRUE) + 0.08,
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0(fmt_p(p), "\nΔ=", sprintf("%.3f", delta))
    )
  return(list(stat_df, df))
}
family_donor <- summary_family_donor(seu, features = c("HERVH-int", "HERV3-int", "HERVK-int", "HERVK11-int"))
save(family_donor, file = file.path(res_dir, "class_family_DE_donor.RData"))
# 准备细胞级class_DE数据
class_ct <- read.csv(file.path(res_dir, "hERV_class_scTE_celltype.csv")) %>%
  select(-matches("^X")) %>%
  mutate(
    celltype = factor(celltype, levels = celltype_order_noall),
    feature = factor(feature, levels = class_order)
  )
class_DE <- plot_class_violin(seu, class_ct, assay = "HERV_class_scTE", celltype_order = celltype_order_noall, class_order = class_order, is_plot = F)
save(class_DE, file = file.path(res_dir, "class_family_DE.RData"))
```

```r
celltype_order <- c("All", "Astro", "Endo", "Excit", "Inhib", "Mic", "Oligo", "OPC")
celltype_order_noall <- c("Astro", "Endo", "Excit", "Inhib", "Mic", "Oligo", "OPC")
class_order <- c("ERV1", "ERVK", "ERVL", "HERV-total")
# 加载数据
class_ct <- read.csv(file.path("C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/hERV_DE", "hERV_class_scTE_celltype.csv")) %>%
  select(-matches("^X")) %>%
  mutate(
    celltype = factor(celltype, levels = celltype_order_noall),
    feature = factor(feature, levels = class_order)
  )
class_all <- read.csv(file.path("C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/hERV_DE", "hERV_class_all.csv")) %>%
  rename(feature = X) %>%
  mutate(
    celltype = "All",
    celltype = factor(celltype, levels = celltype_order),
    feature = factor(feature, levels = class_order)
  )
family_de <- read.csv(file.path("C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/hERV_DE", "hERV_family_scTE_celltype.csv")) %>%
  mutate(
    celltype = factor(celltype, levels = celltype_order_noall),
    logFDR = pmin(-log10(p_val_adj + 1e-300), 30),
    abs_log2FC = abs(avg_log2FC)
  )
class_de <- bind_rows(class_all, class_ct %>% mutate(celltype = as.character(celltype))) %>%
  mutate(
    celltype = factor(celltype, levels = celltype_order),
    feature = factor(feature, levels = class_order),
    logFDR = pmin(-log10(p_val_adj + 1e-300), 30)
  )
load("C:\\Users\\17185\\Desktop\\fsdownload\\class_family_DE_donor.RData")
stat_df <- family_donor[[1]]
df <- family_donor[[2]]
load("C:\\Users\\17185\\Desktop\\fsdownload\\class_family_DE.RData")
df_expr <- class_DE[[1]]
anno_df <- class_DE[[2]]
p_to_star <- function(p){
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}
fmt_p <- function(p){
  dplyr::case_when(
    is.na(p) ~ "p=NA",
    p < 0.001 ~ paste0("p=", format(p, scientific = TRUE, digits = 2)),
    TRUE ~ paste0("p=", sprintf("%.3f", p))
  )
}
```

```r
plot_tag_fz <- 24
# class层级DE汇总bubble plot
class_de <- class_de %>%
  mutate(star = p_to_star(p_val_adj))
max_abs_class <- max(abs(class_de$avg_log2FC), na.rm = TRUE)
p_class_bubble <- ggplot(class_de, aes(x = celltype, y = feature)) +
  geom_point(aes(size = logFDR, color = avg_log2FC)) +
  geom_text(
    data = class_de %>% filter(star != ""),
    aes(label = star),
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +
  scale_size(range = c(2, 8), name = expression(-log[10](FDR))) +
  scale_color_gradient2(
    low = "#3B82F6", mid = "white", high = "#EF4444",
    midpoint = 0,
    limits = c(-max_abs_class, max_abs_class),
    name = "avg_log2FC"
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.6),
    axis.line.y.left = element_line(color = "black", linewidth = 0.6),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 16),
    plot.tag = element_text(size = plot_tag_fz, face = "bold"),
  )
# family层级DE汇总bubble plot
sel_family <- c("HERVH-int", "HERV3-int", "HERVK-int", "HERVK3-int", "HERVK9-int", "HERVK11-int", "HERVK13-int", "HERVK14-int", "HERVK22-int", "HERVIP10F-int", "ERVL-E-int", "ERVL-B4-int", "HERVKC4-int")
family_plot_df <- family_de %>%
  filter(feature %in% sel_family) %>%
  group_by(feature) %>%
  mutate(feature_score = max(abs(avg_log2FC), na.rm = TRUE) + 0.01 * sum(p_val_adj < 0.05)) %>%
  ungroup() %>%
  mutate(star = p_to_star(p_val_adj))
family_order <- family_plot_df %>%
  distinct(feature, feature_score) %>%
  arrange(feature_score) %>%
  pull(feature)
max_abs_family <- max(abs(family_plot_df$avg_log2FC), na.rm = TRUE)
p_family_bubble <- family_plot_df %>%
  mutate(feature = factor(feature, levels = family_order)) %>%
  ggplot(aes(x = celltype, y = feature)) +
  geom_point(aes(size = logFDR, color = avg_log2FC)) +
  geom_text(
    data = function(x) x %>% filter(star != ""),
    aes(label = star),
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +
  scale_size(range = c(1.8, 7), name = expression(-log[10](FDR))) +
  scale_color_gradient2(
    low = "#3B82F6", mid = "white", high = "#EF4444",
    midpoint = 0,
    limits = c(-max_abs_family, max_abs_family),
    name = "avg_log2FC"
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.6),
    axis.line.y.left = element_line(color = "black", linewidth = 0.6),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 13),
    plot.tag = element_text(size = plot_tag_fz, face = "bold"),
  )
# 效应值-显著性散点图
# - 每个点是某个cell type中的某个hERV family的DE结果
# - 蓝线是FDR=0.05，红线是avg_log2FC=±0.1
p_fc_sig <- family_de %>%
  ggplot(aes(x = avg_log2FC, y = logFDR)) +
  geom_point(alpha = 0.7, size = 1.8, color = "grey35") +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = 2, color = "#EF4444") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "#2563EB") +
  facet_wrap(~ celltype, ncol = 4, scales = "free_y") +
  labs(x = "avg_log2FC", y = expression(-log[10](FDR))) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.tag = element_text(size = plot_tag_fz, face = "bold"),
  )
# family-样本箱线图
stat_df <- df %>%
  group_by(celltype, feature) %>%
  summarise(
    p = tryCatch(wilcox.test(expr ~ group)$p.value, error = function(e) NA_real_),
    delta = mean(expr[group == "AD"], na.rm = TRUE) - mean(expr[group == "NC"], na.rm = TRUE),
    ymax = max(expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    p_txt = ifelse(is.na(p), "NA", sprintf("%.3f", p)),
    delta_txt = sprintf("%.3f", delta),
    label = paste0("p=", p_txt, "\n", "\u0394", "=", delta_txt),
    y = ymax * 1.08
  )
plot_family_donor_box <- ggplot(df, aes(x = group, y = expr, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.65) +
  geom_jitter(aes(color = group), width = 0.12, size = 1.8, alpha = 0.85) +
  geom_text(
    data = stat_df, aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE,
    size = 4.6, lineheight = 0.9, vjust = 0
  ) +
  facet_grid(feature ~ celltype, scales = "free_y") +
  scale_fill_manual(values = c("NC" = "#F8766D", "AD" = "#00BFC4")) +
  scale_color_manual(values = c("NC" = "#F8766D", "AD" = "#00BFC4")) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.22))
  ) +
  coord_cartesian(clip = "off") +
  labs(y = "Mean donor-level expression") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.tag = element_text(size = plot_tag_fz, face = "bold"),
    plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
  )
# class-细胞DE图
plot_class_violin <- ggplot(df_expr, aes(x = celltype, y = expr, fill = group)) +
  geom_violin(
    position = position_dodge(width = 0.85),
    trim = TRUE,
    scale = "width",
    linewidth = 0.25
  ) +
  facet_wrap(~ feature, ncol = 2, scales = "free_y") +
  geom_text(
    data = anno_df, aes(x = celltype, y = y, label = star),
    inherit.aes = FALSE, size = 5, vjust = 0
  ) +
  scale_fill_manual(values = c("NC" = "#F8766D", "AD" = "#00BFC4"), name = "Group") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = "Cell type", y = "Expression level") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    plot.tag = element_text(size = plot_tag_fz, face = "bold"),
  )
# 拼图
final_plot <- (p_class_bubble | p_family_bubble) / p_fc_sig / plot_family_donor_box +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 1.2, 2)) +
  theme(
    plot.tag.position = c(0, 1),
    plot.margin = margin(20, 10, 20, 20),
    plot.background = element_rect(fill = "white", colour = NA)
  )
ggsave(
  file.path(res_dir, "class_family_DE_summary.pdf"),
  final_plot, width = 16, height = 24, device = cairo_pdf
)
```

#### 细胞亚类UMAP图

```r
label_size <- 5
p_mic <- DimPlot(mic, reduction = "umap", group.by = "subcluster", label = TRUE, label.size = label_size, repel = TRUE, raster = TRUE) + labs(title = "Mic")
p_astro <- DimPlot(astro, reduction = "umap", group.by = "subcluster", label = TRUE, label.size = label_size, repel = TRUE, raster = TRUE) + labs(title = "Astro")
p_excit <- DimPlot(excit, reduction = "umap", group.by = "subcluster", label = TRUE, label.size = label_size, repel = TRUE, raster = TRUE) + labs(title = "Excit")
p_oligo <- DimPlot(oligo, reduction = "umap", group.by = "subcluster", label = TRUE, label.size = label_size, repel = TRUE, raster = TRUE) + labs(title = "Oligo")
final_plot <- ((p_mic | p_astro) / (p_excit | p_oligo)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.background = element_rect(fill = "white", colour = NA)
    )
  ) &
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(30, 30, 30, 30),
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = c(-0.08, 1.03)
  )
ggsave(
  file.path(res_dir, "subcluster_umap.pdf"),
  final_plot, width = 12, height = 12
)
p_mic <- DimPlot(mic, reduction = "umap", group.by = "group", label = F, raster = TRUE) + labs(title = "Mic")
p_astro <- DimPlot(astro, reduction = "umap", group.by = "group", label = F, raster = TRUE) + labs(title = "Astro")
p_excit <- DimPlot(excit, reduction = "umap", group.by = "group", label = F, raster = TRUE) + labs(title = "Excit")
p_oligo <- DimPlot(oligo, reduction = "umap", group.by = "group", label = F, raster = TRUE) + labs(title = "Oligo")
final_plot <- ((p_mic | p_astro) / (p_excit | p_oligo)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.background = element_rect(fill = "white", colour = NA)
    )
  ) &
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    plot.margin = margin(30, 30, 30, 30),
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = c(-0.08, 1.03)
  )
ggsave(
  file.path(res_dir, "subcluster_group_umap.pdf"),
  final_plot, width = 12, height = 12
)
```

#### 共表达网络组合图

```r
# herv在各亚群的表达量矩阵
herv_avg_mat_list = list()
herv_assay <- "HERV_subfamily_scTE"
plot_subfam <- c("HERV3-int", "HERVH-int", "HERVE-int", "HERVL-int", "HERVK-int", "HERVK3-int", "HERVK9-int", "HERVK11-int", "HERVK13-int", "HERVK14-int", "HERVK22-int")
herv_avg_mat_list[["astro"]] = herv_subcluster_pheatmap(astro, F)
herv_avg_mat_list[["oligo"]] = herv_subcluster_pheatmap(oligo, F)
herv_avg_mat_list[["opc"]] = herv_subcluster_pheatmap(opc, F)
herv_avg_mat_list[["excit"]] = herv_subcluster_pheatmap(excit, F)
herv_avg_mat_list[["inhib"]] = herv_subcluster_pheatmap(inhib, F)
herv_avg_mat_list[["endo"]] = herv_subcluster_pheatmap(endo, F)
herv_avg_mat_list[["mic"]] = herv_subcluster_pheatmap(mic, F)
# 各亚群模块分数
calc_module_score <- function(seu, net_obj_name, modules_use = NA){
  tmp <- seu
  load(file.path(res_dir, net_obj_name))
  seu_net <- seu
  seu <- tmp
  celltype_name <- identity_seu(seu)
  res_name <- paste0(celltype_name, "_subcluster_module_scores")
  modules_df <- GetModules(seu_net, wgcna_name = "gene_net")
  modules_df <- subset(modules_df, module != "grey")
  if(length(modules_use)==1) if(is.na(modules_use)) modules_use <- as.character(unique(modules_df$module))
  modules_df <- subset(modules_df, module %in% modules_use)
  module_list <- split(modules_df$gene_name, modules_df$module)
  module_list <- lapply(module_list, unique)
  module_list <- lapply(module_list, intersect, rownames(seu[["RNA"]]))
  for (m in modules_use) {
    seu <- AddModuleScore(
      seu,
      features = list(module_list[[m]]),
      assay = "RNA",
      name = paste0(m, "_score_")
    )
  }
  score_cols <- setNames(paste0(modules_use, "_score_1"), modules_use)
  Idents(seu) <- seu[["subcluster"]][, 1]
  score_avg <- seu@meta.data[, c("subcluster", unname(score_cols)), drop = FALSE]
  colnames(score_avg)[1] <- "subcluster"
  score_avg <- score_avg %>%
    group_by(subcluster) %>%
    summarise(across(everything(), mean))
  score_mat <- as.matrix(score_avg[, -1])
  rownames(score_mat) <- score_avg$subcluster
  colnames(score_mat) <- names(score_cols)
  return(score_mat)
}
module_avg_mat_list = list()
module_avg_mat_list$astro <- calc_module_score(astro, "Astro_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData")
module_avg_mat_list$oligo <- calc_module_score(oligo, "Oligo_HERV_family_0.1_4_20_0.05_gene_WGCNA.RData")
module_avg_mat_list$excit <- calc_module_score(excit, "Excit_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData")
module_avg_mat_list$inhib <- calc_module_score(inhib, "Inhib_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData")
module_avg_mat_list$opc <- calc_module_score(opc, "OPC_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData")
module_avg_mat_list$endo <- calc_module_score(endo, "Endo_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData")
module_avg_mat_list$mic <- calc_module_score(mic, "Mic_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData")
# 模块相关性
get_plot <- function(seu){
  target_herv <- c("HERV3-int", "HERVH-int", "HERVE-int", "HERVL-int", "HERVK-int", "HERVK3-int", "HERVK9-int", "HERVK11-int", "HERVK13-int", "HERVK14-int", "HERVK22-int")
  my_print("- calc corr")
  seu <- herv_gene_module_corr(seu, "HERV_subfamily_scTE", target_herv, F)
  mt <- GetModuleTraitCorrelation(seu)
  cor_df <- mat_to_df(mt$cor$all_cells, "cor")
  fdr_df <- mat_to_df(mt$fdr$all_cells, "fdr")
  return(list(cor_df, fdr_df))
}
net_obj_name_list <- list(
  "scTE_Astro_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData",
  "scTE_Oligo_HERV_family_0.1_4_20_0.05_gene_WGCNA.RData",
  "scTE_OPC_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData",
  "scTE_Excit_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData",
  "scTE_Inhib_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData",
  "scTE_Endo_HERV_family_0.1_4_30_0.1_gene_WGCNA.RData",
  "scTE_Mic_HERV_family_0.1_4_50_0.2_gene_WGCNA.RData"
)
corr_mat_list = list()
fdr_mat_list = list()
for(net_obj_name in net_obj_name_list){  # 接续前文summary_locus_WGCNA
  my_print(net_obj_name)
  my_print("- reading RData")
  load(file.path(res_dir, net_obj_name))
  ct <- tolower(identity_seu(seu))
  res <- get_plot(seu)
  corr_mat_list[[ct]] <- res[[1]]
  fdr_mat_list[[ct]] <- res[[2]]
}
wide2long <- function(df_list, values_from){
  lapply(df_list, function(df){
    m <- pivot_wider(df, names_from = "module", values_from = values_from) %>% 
      column_to_rownames("trait") %>%
      as.matrix(., drop = F)
    rownames(m) <- gsub("herv_", "", rownames(m))
    rownames(m) <- gsub("\\.", "-", rownames(m))
    return(m)
  })
}
corr_mat_list <- wide2long(corr_mat_list, "cor")
fdr_mat_list <- wide2long(fdr_mat_list, "fdr")
```

```r
fix_cluster_name <- function(x) {
  x <- as.character(x)
  gsub("g", "", x)
}
# hERV平均表达量热图
plot_herv_avg_heatmap <- function(mat, ct = "", herv_keep = NULL, cluster_order = NULL, highlight = NULL) {
  colnames(mat) <- fix_cluster_name(colnames(mat))
  if (!is.null(herv_keep)) {
    herv_keep <- intersect(herv_keep, rownames(mat))
    mat <- mat[herv_keep, , drop = FALSE]
  }
  if (!is.null(cluster_order)) {
    cluster_order <- intersect(cluster_order, colnames(mat))
    mat <- mat[, cluster_order, drop = FALSE]
  }
  ann_col <- NULL
  ann_colors <- NULL
  if (!is.null(highlight)) {
    ann_col <- data.frame(
      key_cluster = ifelse(colnames(mat) %in% highlight, "yes", "no"),
      row.names = colnames(mat)
    )
    ann_colors <- list(key_cluster = c(yes = "#d73027", no = "grey90"))
  }
  p <- pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 11,
    fontsize_col = 11,
    main = ct,
    border_color = NA,
    angle_col = 0,
    annotation_col = ann_col,
    annotation_colors = ann_colors
  )
  p$gtable
}
# 模块分数热图
plot_module_avg_heatmap <- function(mat, ct = "", module_keep = NULL, cluster_order = NULL, highlight = NULL) {
  mat <- t(as.matrix(mat))
  colnames(mat) <- fix_cluster_name(colnames(mat))
  if (!is.null(module_keep)) {
    module_keep <- intersect(module_keep, rownames(mat))
    mat <- mat[module_keep, , drop = FALSE]
  }
  if (!is.null(cluster_order)) {
    cluster_order <- intersect(cluster_order, colnames(mat))
    mat <- mat[, cluster_order, drop = FALSE]
  }
  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0
  ann_col <- NULL
  ann_colors <- NULL
  if (!is.null(highlight)) {
    ann_col <- data.frame(
      key_cluster = ifelse(colnames(mat_z) %in% highlight, "yes", "no"),
      row.names = colnames(mat_z)
    )
    ann_colors <- list(key_cluster = c(yes = "#d73027", no = "grey90"))
  }
  p <- pheatmap(
    mat_z,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 11,
    fontsize_col = 11,
    main = ct,
    border_color = NA,
    angle_col = 0,
    annotation_col = ann_col,
    annotation_colors = ann_colors
  )
  p$gtable
}
# hERV~module相关性热图
plot_corr_heatmap <- function(corr_mat, fdr_mat, ct = "", herv_keep = NULL, module_keep = NULL) {
  corr_mat <- as.matrix(corr_mat)
  fdr_mat <- as.matrix(fdr_mat)
  if (!is.null(herv_keep)) {
    herv_keep <- intersect(herv_keep, rownames(corr_mat))
    corr_mat <- corr_mat[herv_keep, , drop = FALSE]
    fdr_mat <- fdr_mat[herv_keep, , drop = FALSE]
  }
  if (!is.null(module_keep)) {
    module_keep <- intersect(module_keep, colnames(corr_mat))
    corr_mat <- corr_mat[, module_keep, drop = FALSE]
    fdr_mat <- fdr_mat[, module_keep, drop = FALSE]
  }
  lab <- matrix("", nrow = nrow(fdr_mat), ncol = ncol(fdr_mat), dimnames = dimnames(fdr_mat))
  lab[fdr_mat < 0.05] <- "*"
  lab[fdr_mat < 0.01] <- "**"
  lab[fdr_mat < 0.001] <- "***"
  p <- pheatmap(
    corr_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = lab,
    number_color = "black",
    fontsize_row = 11,
    fontsize_col = 11,
    fontsize_number = 9,
    main = ct,
    border_color = NA,
    angle_col = 0
  )
  p$gtable
}
# GO结果图
parse_ratio <- function(x) {
  sapply(strsplit(as.character(x), "/"), function(v) {
    as.numeric(v[1]) / as.numeric(v[2])
  })
}
plot_go_dot <- function(go_df_path_list, wide_or_long = c("wide", "long"), n_show = 6, cluster_keep = NA, module_name_list = NA, facet_by = "module") {
  # 如果是DEG_GO这种把p值写到一个单元格的就是wide模式，如果是模块的GO就是long模式
  if(wide_or_long == "wide"){
    if(length(go_df_path_list) != 1) return("length(go_df_path_list) != 1")
    go_df <- read.csv(go_df_path_list)
    if(length(cluster_keep) == 1) if(is.na(cluster_keep)) cluster_keep <- unique(go_df$cluster)
    go_plot <- go_df %>%
      filter(cluster %in% cluster_keep) %>%
      mutate(
        GO = str_split(GO, ";"),
        GO_padj = str_split(GO_padj, ";")
      ) %>%
      unnest(c(GO, GO_padj)) %>%
      mutate(
        GO = str_trim(GO),
        GO_padj = as.numeric(str_trim(GO_padj)),
        neglog10 = -log10(GO_padj),
        group_lab = paste0("g", cluster, " (", module, ")")
      ) %>%
      group_by(group_lab) %>%
      arrange(GO_padj, .by_group = TRUE) %>%
      slice_head(n = n_show) %>%
      ungroup() %>%
      mutate(
        GO = str_wrap(GO, width = 35),
        GO_show = reorder_within(GO, neglog10, group_lab)
      )
    p <- ggplot(go_plot, aes(x = neglog10, y = GO_show)) +
      geom_point(aes(color = GO_padj), size = 3) +
      facet_grid(rows = vars(group_lab), scales = "free_y", space = "free_y") +
      scale_y_reordered() +
      scale_color_gradient(
        low = "#d73027", high = "#4575b4",
        trans = "reverse",
        name = "FDR"
      ) +
      labs(x = "-log10(FDR)", y = NULL) +
      theme_bw(base_size = 10) +
      theme(
        strip.background = element_rect(fill = "grey90", color = "grey60"),
        strip.text.y.right = element_text(angle = 0, face = "bold"),
        strip.placement = "outside",
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.y = unit(0.4, "lines")
      )
  } else {
    if(length(module_name_list) == 1) if(is.na(module_name_list)){
      module_name_list <- lapply(go_df_path_list, function(path){
        strsplit(sub("\\.csv$", "", basename(path)), "_")[[1]][4]
      }) %>% unlist()
    }
    go_df <- lapply(1:length(go_df_path_list), function(i){
      read.csv(go_df_path_list[i]) %>%
        mutate(module = module_name_list[i])
    }) %>% bind_rows()
    go_df$group_lab <- factor(go_df[[facet_by]], levels = module_name_list)
    go_df$x_value <- parse_ratio(go_df$GeneRatio)
    x_lab <- "GeneRatio"
    go_plot <- go_df %>%
      group_by(group_lab) %>%
      arrange(p.adjust, .by_group = TRUE) %>%
      slice_head(n = n_show) %>%
      ungroup() %>%
      mutate(
        neglog10 = -log10(p.adjust),
        Description = str_wrap(Description, width = 35),
        GO_show = reorder_within(Description, x_value, group_lab)
      )
    p <- ggplot(go_plot, aes(x = x_value, y = GO_show)) +
      geom_point(aes(size = Count, color = p.adjust)) +
      facet_grid(rows = vars(group_lab), scales = "free_y", space = "free_y") +
      scale_y_reordered() +
      scale_color_gradient(
        low = "#d73027", high = "#4575b4",
        trans = "reverse",
        name = "p.adjust"
      ) +
      labs(x = x_lab, y = NULL, size = "Count") +
      theme_bw(base_size = 11) +
      theme(
        panel.background = element_rect(fill = "white", colour = "black"),  # 去除背景色
        legend.background = element_rect(fill = "white", colour = NA),
        legend.box.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        strip.background = element_rect(fill = "grey90", colour = "black"),  # 右侧注释
        strip.text.y.right = element_text(size = 12, face = "bold", angle = 0),  # 右侧注释字号
        strip.text.y = element_text(size = 2, face = "bold", angle = 0),
        axis.text.x = element_text(size = 11, colour = "black"),  # x/y轴字号
        axis.text.y = element_text(size = 11, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10),  # 图例字号
        legend.text = element_text(size = 9),
        legend.key = element_rect(fill = "white", colour = NA),  # 去除图例图标边框
        panel.grid.major.y = element_blank(),  # 去除网格线
        panel.grid.minor = element_blank(),
        plot.margin = margin(6, 10, 6, 6)  # 图边距
      )
  }
  return(p)
}
# 汇总
align_heatmap_grobs <- function(...) {
  grobs <- list(...)
  max_widths <- do.call(unit.pmax, lapply(grobs, function(g) g$widths))
  grobs <- lapply(grobs, function(g) {
    g$widths <- max_widths
    g
  })
  grobs
}
pad_plot <- function(p, l = 12, r = 6, t = 6, b = 6) {
  p + theme(
    plot.margin = margin(t, r, b, l),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )
}
wgcna_pheatmap_GO_plot <- function(celltype, go_df_path_list, herv_keep = NULL, cluster_order = NULL, highlight = NULL, module_keep = NULL, n_show = 6){
  g_herv <- plot_herv_avg_heatmap(
    herv_avg_mat_list[[celltype]],
    herv_keep = herv_keep,
    cluster_order = cluster_order,
  )
  g_module <- plot_module_avg_heatmap(
    module_avg_mat_list[[celltype]],
    module_keep = module_keep,
    cluster_order = cluster_order,
  )
  g_corr <- plot_corr_heatmap(
    corr_mat_list[[celltype]],
    fdr_mat_list[[celltype]],
    herv_keep = herv_keep,
    module_keep = module_keep
  )
  p_go <- plot_go_dot(
    go_df_path_list,
    wide_or_long = "longer",
    n_show = n_show
  )
  grobs_aligned <- align_heatmap_grobs(g_herv, g_module, g_corr)
  p_herv <- as.ggplot(grobs_aligned[[1]]) |> pad_plot(l = 14, r = 20, t = 6, b = 6)
  p_module <- as.ggplot(grobs_aligned[[2]]) |> pad_plot(l = 14, r = 20, t = 6, b = 6)
  p_corr <- as.ggplot(grobs_aligned[[3]]) |> pad_plot(l = 14, r = 20, t = 6, b = 6)
  p_left <- plot_grid(
    p_herv, p_module, p_corr,
    ncol = 1,
    labels = c("A", "B", "C"),
    label_size = 16,
    label_fontface = "bold",
    hjust = -0.3,
    vjust = 1.2,
    rel_heights = c(1, 1, 1)
  )
  p_right <- plot_grid(
    p_go,
    ncol = 1,
    labels = "D",
    label_size = 16,
    label_fontface = "bold",
    hjust = -0.3,
    vjust = 1.2
  )
  p_all <- plot_grid(
    p_left, p_right,
    ncol = 2,
    rel_widths = c(1, 1.25)
  ) + theme(
    plot.margin = margin(15, 0, 5, 15),
    plot.background = element_rect(fill = "white", colour = NA)
  )
  ggsave(
    file.path(res_dir, paste0("wgcna_", celltype, ".pdf")),
    p_all, width = 12, height = 12
  )
  return(p_all)
}
```

```r
target_herv <- c("AD_binary", "HERV3-int", "HERVH-int", "HERVL-int", "HERVK-int", "HERVK3-int", "HERVK9-int", "HERVK11-int", "HERVK13-int", "HERVK14-int", "HERVK22-int")
wgcna_pheatmap_GO_plot(
  celltype = "mic",
  herv_keep = target_herv,
  module_keep = c("blue", "turquoise", "brown", "pink"),
  cluster_order = 0:10,
  highlight = c("2"),
  go_df_path_list = c(
    "C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/hERV_gene_WGCNA/mic/DEG_GO/Mic_DEG-GO_BP_blue.csv",
    "C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/hERV_gene_WGCNA/mic/DEG_GO/Mic_DEG-GO_BP_pink.csv"
  ),
  n_show = 10
)
wgcna_pheatmap_GO_plot(
  celltype = "oligo",
  herv_keep = target_herv,
  module_keep = c("blue", "turquoise", "yellow", "green", "brown"),
  cluster_order = 0:9,
  highlight = c("3", "5", "9"),
  go_df_path_list = c(
    "C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/hERV_gene_WGCNA/oligo/GO/Oligo_GO_BP_blue.csv",
    "C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/hERV_gene_WGCNA/oligo/GO/Oligo_GO_BP_yellow.csv",
    "C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/hERV_gene_WGCNA/oligo/GO/Oligo_GO_BP_turquoise.csv"
  )
)
wgcna_pheatmap_GO_plot(
  celltype = "astro",
  herv_keep =target_herv,
  module_keep = c("blue", "red"),
  cluster_order = 0:6,
  highlight = c("0", "1"),
  go_df_path_list = c(
    "C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/hERV_gene_WGCNA/astro/GO/Astro_GO_BP_blue.csv",
    "C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/hERV_gene_WGCNA/astro/GO/Astro_GO_BP_red.csv"
  ),
  n_show = 9
)
wgcna_pheatmap_GO_plot(
  celltype = "excit",
  herv_keep = target_herv,
  module_keep = c("blue", "green", "magenta", "pink", "black", "red"),
  cluster_order = 0:15,
  highlight = c("6"),
  go_df_path_list = c(
    "C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/hERV_gene_WGCNA/excit/GO/Excit_GO_BP_black.csv",
    "C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/hERV_gene_WGCNA/excit/GO/Excit_GO_BP_red.csv",
    "C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/hERV_gene_WGCNA/excit/GO/Excit_GO_BP_pink.csv",
    "C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/hERV_gene_WGCNA/excit/GO/Excit_GO_BP_green.csv"
  ),
  n_show = 5
)
```

#### 共表达网络总述图

```r
df <- read_xlsx("C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/hERV_gene_WGCNA/result.xlsx")
func_axis_levels <- c(
  "glial_neural_support",
  "neural_synaptic_function",
  "translation_energy_metabolism",
  "membrane_signal_remodeling",
  "cytoskeleton_adhesion",
  "antiviral_innate_immune",
  "protein_folding_stress",
  "vascular_circulation_program",
  "DNA_replication_RNA_processing"
)
func_axis_labels <- c(
  glial_neural_support = "Glial / neural\nsupport",
  neural_synaptic_function = "Neural / synaptic\nfunction",
  translation_energy_metabolism = "Translation /\nenergy metabolism",
  membrane_signal_remodeling = "Membrane / signal\nremodeling",
  cytoskeleton_adhesion = "Cytoskeleton /\nadhesion",
  antiviral_innate_immune = "Antiviral /\ninnate immune",
  inflammation_TLR_cytokine = "Inflammation /\nTLR-cytokine",
  protein_folding_stress = "Protein folding /\nstress",
  vascular_circulation_program = "Vascular /\ncirculation",
  DNA_replication_RNA_processing = "DNA replication /\nRNA processing"
)
celltype_levels <- c("Mic", "Oligo", "Astro", "Excit", "OPC", "Inhib", "Endo")
df$celltype <- factor(df$celltype, levels = rev(celltype_levels))
df$func_axis <- factor(df$func_axis, levels = func_axis_levels)
df$main_label <- gsub("/", "\n", df$herv_label)
# 纵轴顺序
celltype_pos <- c(
  Mic   = 7.4,
  Oligo = 6.0,
  Astro = 5.0,
  Excit = 4.0,
  OPC   = 2.6,
  Inhib = 1.6,
  Endo  = 0.6
)
celltype_order_topdown <- c("Mic", "Oligo", "Astro", "Excit", "OPC", "Inhib", "Endo")
func_pos <- setNames(seq_along(func_axis_levels), func_axis_levels)
df <- df %>%
  mutate(
    x = unname(func_pos[func_axis]),
    y = unname(celltype_pos[as.character(celltype)]),
    cluster_label = ifelse(is.na(cluster_label) | cluster_label %in% c("", "NA"), "", as.character(cluster_label))
  )
df_strong <- df %>% filter(evidence_level == "strong")
df_cluster <- df %>% filter(cluster_label != "")
p <- ggplot(df, aes(x = x, y = y)) +
  geom_tile(
    aes(fill = plot_fill, alpha = plot_alpha),
    width = 0.96, height = 0.96,
    color = "grey80", linewidth = 0.5
  ) +
  geom_tile(
    data = df_strong,
    fill = NA, color = "black", linewidth = 0.9,
    width = 0.96, height = 0.96
  ) +
  geom_hline(
    yintercept = c(6.7, 3.3),
    linewidth = 0.9,
    color = "grey45"
  ) +
  geom_text(
    aes(label = main_label),
    nudge_y = -0.08,
    size = 3, lineheight = 0.88, color = "black"
  ) +
  geom_text(
    data = df_cluster,
    aes(x = x, y = y, label = cluster_label),
    inherit.aes = FALSE,
    nudge_y = 0.28,
    size = 3.6, color = "white", fontface = "bold"
  ) +
  annotate(
    "text",
    x = max(df$x) + 1.3, y = celltype_pos["Mic"],
    label = "Key\nfinding",
    angle = 0, fontface = "bold", size = 5
  ) +
  annotate(
    "text",
    x = max(df$x) + 1.3, y = mean(celltype_pos[c("Oligo", "Astro", "Excit")]),
    label = "Main\nfindings",
    angle = 0, fontface = "bold", size = 5
  ) +
  annotate(
    "text",
    x = max(df$x) + 1.3, y = mean(celltype_pos[c("OPC", "Inhib", "Endo")]),
    label = "Supportive\nfindings",
    angle = 0, fontface = "bold", size = 5
  ) +
  scale_fill_gradient2(
    low = "#3b4cc0", mid = "white", high = "#d73027",
    midpoint = 0,
    limits = c(-3, 3), 
    breaks = c(-3, -2, -1, 0, 1, 2, 3),
    name = "Association"
  ) +
  scale_alpha_identity() +
  scale_x_continuous(
    breaks = seq_along(func_axis_levels),
    labels = unname(func_axis_labels[func_axis_levels]),
    expand = expansion(add = c(0.15, 1.15))
  ) +
  scale_y_continuous(
    breaks = unname(celltype_pos[celltype_order_topdown]),
    labels = celltype_order_topdown,
    limits = c(0, 10.4),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.margin = margin(-80, 10, 15, 20)
  )
ggsave(file.path(res_dir, "wgcna_summary.pdf"), p, width = 10, height = 8)
```

#### 位点级相关性

```r
final_loci <- c("HERV3-1p34.3", "HML3-1q32.1b", "HML1-1q32.1", "HML4-16p13.3", "MER41-13q13.3", "MER61-21q21.1a", "MER61-21q21.1b", "HERVW-15q21.2", "HARLEQUIN-6p21.31")
celltype_order <- c("Mic", "Oligo", "Astro", "Excit", "Inhib", "OPC", "Endo")
# 接续前文的一些变量
de_df <- read.csv(file.path(res_dir, "hERV_locus_celltype.csv"))  # DE结果
corr_mat_list = list()  # 位点wgcna结果
fdr_mat_list = list()
for(net_obj_name in net_obj_name_list){  # 接续前文summary_locus_WGCNA部分
  my_print(net_obj_name)
  my_print("- loading RData")
  load(file.path(res_dir, net_obj_name))
  ct <- identity_seu(seu)
  res <- summary_locus_WGCNA(seu, is_get_cor_mat = T)
  corr_mat_list[[ct]] <- res[[1]]
  fdr_mat_list[[ct]] <- res[[2]]
}
wide2long <- function(df_list, values_from){
  lapply(df_list, function(df){
    m <- pivot_wider(df, names_from = "module", values_from = values_from) %>% 
      column_to_rownames("trait") %>%
      as.matrix(., drop = F)
    rownames(m) <- gsub("herv_", "", rownames(m))
    rownames(m) <- gsub("\\.(?=.*\\.)", "-", rownames(m), perl = TRUE)  # 将所有的.转成-，但不包括最后一个.
    return(m)
  })
}
corr_mat_list <- wide2long(corr_mat_list, "cor")
fdr_mat_list <- wide2long(fdr_mat_list, "fdr")
p_corr <- plot_herv_gene_corr_df(seu_list, plot_specs, is_plot = F)  # 4*3相关性线图
# 画图参数
base_size_main <- 12
title_size_main <- 14
axis_text_size_main <- 13
strip_text_size_main <- 12
legend_title_size_main <- 12
legend_text_size_main <- 10
tag_size_main <- 18
theme_main <- theme_bw(base_size = base_size_main) +
  theme(
    plot.title = element_text(size = title_size_main, face = "bold"),
    axis.title = element_blank(),
    axis.text = element_text(size = axis_text_size_main, color = "black"),
    strip.text = element_text(size = strip_text_size_main),
    legend.title = element_text(size = legend_title_size_main),
    legend.text = element_text(size = legend_text_size_main),
    plot.tag = element_text(size = tag_size_main, face = "bold"),
    plot.tag.position = c(0, 1)
  )
```

```r
# 位点级DE
de_plot_df <- de_df %>%
  filter(celltype %in% celltype_order, feature %in% final_loci) %>%
  mutate(
    celltype = factor(celltype, levels = rev(celltype_order)),
    feature = factor(feature, levels = final_loci),
    pct_max = pmax(pct.1, pct.2),
    sig_lab = case_when(
      p_val_adj < 0.001 ~ "***",
      p_val_adj < 0.01 ~ "**",
      p_val_adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )
p_de <- ggplot(de_plot_df, aes(x = feature, y = celltype)) +
  geom_point(
    aes(size = pct_max, fill = avg_log2FC),
    shape = 21, color = "black", stroke = 0.3
  ) +
  geom_text(aes(label = sig_lab), size = 3) +
  scale_fill_gradient2(low = "#3b4cc0", mid = "white", high = "#b40426", midpoint = 0) +
  scale_size(range = c(2.5, 10)) +
  theme_main +
  theme(
    axis.text.x = element_text(size = axis_text_size_main, angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(fill = "avg_log2FC", size = "max pct")
# wgcna
ct_keep <- c("Mic", "Oligo", "Astro", "Excit")
herv_keep <- c("HARLEQUIN-6p21.31", "HERVW-15q21.2", "MER41-13q13.3", "HML4-16p13.3", "HML1-1q32.1", "HML3-1q32.1b")
module_keep_list <- list(
  Mic = c("green", "turquoise"),
  Oligo = c("black", "turquoise"),
  Astro = c("green"),
  Excit = c("green")
)
df_list <- list()
for (ct in ct_keep) {
  corr_mat <- as.matrix(corr_mat_list[[ct]])
  fdr_mat <- as.matrix(fdr_mat_list[[ct]])
  herv_use <- intersect(herv_keep, rownames(corr_mat))
  module_use <- intersect(module_keep_list[[ct]], colnames(corr_mat))
  tmp <- expand.grid(
    locus = herv_use,
    module = module_use,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      cor = mapply(function(h, m) corr_mat[h, m], locus, module),
      fdr = mapply(function(h, m) fdr_mat[h, m], locus, module),
      label = case_when(
        fdr < 0.001 ~ "***",
        fdr < 0.01 ~ "**",
        fdr < 0.05 ~ "*",
        TRUE ~ ""
      ),
      celltype = ct,
      x_id = paste(ct, module, sep = "__")
    )
  df_list[[ct]] <- tmp
}
plot_df <- bind_rows(df_list)
x_levels <- unlist(lapply(ct_keep, function(ct) {
  paste(ct, module_keep_list[[ct]], sep = "__")
}))
x_levels <- x_levels[x_levels %in% unique(plot_df$x_id)]
plot_df <- plot_df %>%
  mutate(
    x_id = factor(x_id, levels = x_levels),
    locus = factor(locus, levels = rev(herv_keep))
  )
lim <- max(abs(plot_df$cor), na.rm = TRUE)
p_wgcna <- ggplot(plot_df, aes(x = x_id, y = locus, fill = cor)) +
  geom_tile(color = NA) +
  geom_text(aes(label = label), size = 3.5) +
  facet_grid(. ~ celltype, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = function(x) sub("^.*__", "", x)) +
  scale_fill_gradient2(
    low = "#3b4cc0",
    mid = "white",
    high = "#b40426",
    midpoint = 0,
    limits = c(-lim, lim)
  ) +
  theme_main +
  theme(
    axis.text.x = element_text(size = axis_text_size_main, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = axis_text_size_main),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "black"),
    panel.spacing.x = unit(0.35, "lines")
  ) +
  labs(fill = "cor")
# 拼图
p_de <- p_de + labs(tag = "A")
p_wgcna <- p_wgcna + labs(tag = "B")
p_corr_panel <- wrap_elements(full = p_corr) + labs(tag = "C")
final_plot <- (p_de | p_wgcna) / p_corr_panel +
  plot_layout(heights = c(1, 2.4)) &
  theme(
    plot.tag = element_text(size = tag_size_main, face = "bold"),
    plot.tag.position = c(0, 1),
    plot.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(fill = "white", colour = NA)
  )
ggsave(
  file.path(res_dir, "locus_DE_WGCNA_geneCorr.pdf"),
  final_plot, width = 16, height = 16
)
```

#### 基因组浏览器和GWAS分析

```sh
# 可以使用之前下载的gtf文件，也可使用一个更简洁的版本
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz
# gunzip gencode.v49.basic.annotation.gtf.gz
# grep -v '^#' gencode.v49.basic.annotation.gtf | \
# awk -F '\t' '$9 ~ /gene_type "protein_coding"/' \
# > gencode.v49.protein_coding_only.gtf
module load miniconda3/base
conda activate pygenometracks
pyGenomeTracks \
  --tracks /public/home/GENE_proc/wth/my_data/tracks.ini \
  --BED /public/home/GENE_proc/wth/my_data/plot_windows.bed \
  --width 18 \
  --trackLabelFraction 0.18 \
  --fontSize 10 \
  -o /public/home/GENE_proc/wth/res/locus_plot.pdf
```

```r
# 拼pdf
pdf_files <- list(
  harlequin = "/public/home/GENE_proc/wth/res/locus_plot_chr6-34948360-35054307.pdf",
  hervw = "/public/home/GENE_proc/wth/res/locus_plot_chr15-51307365-51410808.pdf",
  mer41 = "/public/home/GENE_proc/wth/res/locus_plot_chr13-39030844-39136958.pdf",
  hml4 = "/public/home/GENE_proc/wth/res/locus_plot_chr16-2610261-2720537.pdf",
  hml1_hml3 = "/public/home/GENE_proc/wth/res/locus_plot_chr1-205807563-205946160.pdf"
)
pdf_titles <- list(
  harlequin = "HARLEQUIN-6p21.31",
  hervw = "HERVW-15q21.2",
  mer41 = "MER41-13q13.3",
  hml4 = "HML4-16p13.3",
  hml1_hml3 = "HML1-1q32.1 & HML3-1q32.1b"
)
read_pdf_panel <- function(pdf_file, cut_bottom = 0, page = 1, density = 300) {
  img <- image_read_pdf(path = pdf_file, pages = page, density = density)
  info <- image_info(img)
  w <- info$width
  h <- info$height
  left <- 0
  top <- 0
  bottom <- round(h * cut_bottom)
  new_w <- w 
  new_h <- h - bottom
  img <- image_crop(img, geometry = geometry_area(width = new_w, height = new_h, x_off = left, y_off = top))
  ggdraw() + draw_image(img, x = 0, y = 0, width = 1, height = 1)
}
crop_list <- list(
  harlequin = 0.2,
  hervw = 0.2,
  mer41 = 0.2,
  hml4 = 0.2,
  hml1_hml3 = 0.2
)
p_browser_list <- lapply(names(pdf_files), function(nm) {
  read_pdf_panel(
    pdf_file = pdf_files[[nm]],
    cut_bottom = crop_list[[nm]],
  )
})
names(p_browser_list) <- names(pdf_files)
p_browser <- p_browser_list$harlequin / p_browser_list$hervw / p_browser_list$mer41 / p_browser_list$hml4 / p_browser_list$hml1_hml3
ggsave(
  file.path(res_dir, "p_browser.pdf"),
  p_browser, width = 16, height = 16
)
# 用pdf工具裁剪一下（调整各图间隔，最后导出图片）
png_files <- list(
  harlequin = "C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/herv_locus_genome_GWAS/genome_res/locus_plot_chr6-34848360-35154307.png",
  hervw = "C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/herv_locus_genome_GWAS/genome_res/locus_plot_chr15-51307365-51410808.png",
  mer41 = "C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/herv_locus_genome_GWAS/genome_res/locus_plot_chr13-38930844-39236958.png",
  hml4 = "C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/herv_locus_genome_GWAS/genome_res/locus_plot_chr16-2610261-2720537.png",
  hml1_hml3 = "C:/Users/17185/Desktop/hERV_AD_snRNA-seq_analyze/res/herv_locus_genome_GWAS/genome_res/locus_plot_chr1-205807563-205946160.png"
)
p_png_list <- lapply(names(png_files), function(nm) {
  ggdraw() + draw_image(png_files[[nm]], x = 0, y = 0, width = 1, height = 1)
})
names(p_png_list) <- names(png_files)
p_png_list$hml1_hml3 <- p_png_list$hml1_hml3 +
  theme(
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = c(0.05, 1),
    plot.margin = margin(20, 0, 0, -20),
    plot.background = element_rect(fill = "white", colour = NA)
  )
p_png_list$harlequin <- p_png_list$harlequin +
  theme(
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = c(0.05, 1),
    plot.margin = margin(0, 0, 0, -20),
    plot.background = element_rect(fill = "white", colour = NA)
  )
p_png_list$mer41 <- p_png_list$mer41 +
  theme(
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = c(0.05, 1),
    plot.margin = margin(0, 0, 0, -20),
    plot.background = element_rect(fill = "white", colour = NA)
  )
p_png_list$hml4 <- p_png_list$hml4 +
  theme(
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = c(0.05, 1),
    plot.margin = margin(0, 0, 0, -20),
    plot.background = element_rect(fill = "white", colour = NA)
  )
p_png_list$hervw <- p_png_list$hervw +
  theme(
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = c(0.05, 1),
    plot.margin = margin(0, 0, 0, -20),
    plot.background = element_rect(fill = "white", colour = NA)
  )
p_png <- p_png_list$hml1_hml3 / p_png_list$harlequin / p_png_list$mer41 / p_png_list$hervw / p_png_list$hml4 +
  plot_annotation(tag_levels = "A")
ggsave(
  file.path(res_dir, "genome_pos.pdf"),
  p_png, width = 12, height = 16
)
```

```r
# 接续前文plot_herv_gene_corr_df
plot_pairs_main <- tribble(
  ~ct, ~herv, ~gene,
  "Mic",   "HARLEQUIN-6p21.31", "ANKS1A",
  "Oligo", "MER41-13q13.3",     "NHLRC3",
  "Excit", "HML4-16p13.3",      "KCTD5"
)
celltype_order <- c("Mic", "Oligo", "Astro", "Excit")
gene_info <- rbind(plot_specs, plot_pairs_main) %>%
  mutate(celltype = factor(ct, levels = celltype_order)) %>%
  arrange(celltype)
# 细胞类型内DE分析
de_gene_df <- lapply(seu_list, function(seu){
  seu_ct <- identity_seu(seu)
  my_print("Running DE for: ", seu_ct)
  Idents(seu) <- "group"
  genes_use <- gene_info %>%
    filter(ct == seu_ct) %>%
    pull(gene) %>%
    unique()
  FindMarkers(
    seu,
    ident.1 = "AD",
    ident.2 = "NC",
    assay = "RNA",
    test.use = "MAST",
    features = genes_use,
    min.pct = 0,
    logfc.threshold = 0,
  ) %>%
    rownames_to_column("gene") %>%
    mutate(celltype = seu_ct)
}) %>% bind_rows()
write.csv(de_gene_df, file.path(res_dir, "supply_gene_DE.csv"), row.names = FALSE)
# 分析哪些gene符合预测（若某herv位点上调，herv~gene正相关，gene也上调，就符合预测）
gene_de <- read.csv(file.path(res_dir, "supply_gene_DE.csv"), stringsAsFactors = FALSE)
herv_de <- read.csv(file.path(res_dir, "hERV_locus_celltype.csv"), stringsAsFactors = FALSE) %>%
  filter(feature %in% c(
    "HARLEQUIN-6p21.31",
    "HERVW-15q21.2",
    "MER41-13q13.3",
    "HML1-1q32.1",
    "HML3-1q32.1b",
    "HML4-16p13.3"
  )) %>%
  select(celltype, feature, herv_logFC = avg_log2FC)
pair_df <- tribble(  # 已经确认的herv~gene相关性
  ~celltype, ~gene, ~herv, ~corr_sign,
  "Mic",   "ADGRE2",  "HARLEQUIN-6p21.31",  1,
  "Mic",   "IL4R",    "HARLEQUIN-6p21.31",  1,
  "Mic",   "ALOX5",   "HARLEQUIN-6p21.31",  1,
  "Mic",   "HCK",     "HARLEQUIN-6p21.31",  1,
  "Mic",   "ANKS1A",  "HARLEQUIN-6p21.31",  1,
  "Oligo", "GLDN",    "HERVW-15q21.2",      1,
  "Oligo", "NHLRC3",  "HERVW-15q21.2",      1,
  "Oligo", "ATP11A",  "MER41-13q13.3",     -1,
  "Oligo", "PIK3C2B", "HERVW-15q21.2",      1,
  "Astro", "ABCC9",   "HML1-1q32.1",        1,
  "Astro", "ABCC9",   "HML3-1q32.1b",       1,
  "Astro", "COLEC12", "HML3-1q32.1b",       1,
  "Excit", "MTR",     "HML4-16p13.3",       1,
  "Excit", "MICU2",   "HML4-16p13.3",       1,
  "Excit", "KCTD5",   "HML4-16p13.3",      -1
)
plot_df <- pair_df %>%
  left_join(
    herv_de,
    by = c("celltype", "herv" = "feature")
  ) %>%
  left_join(
    gene_de %>% select(celltype, gene, gene_logFC = avg_log2FC, p_val_adj, pct.1, pct.2),
    by = c("celltype", "gene")
  ) %>%
  mutate(
    expected_sign = sign(herv_logFC) * corr_sign,
    observed_sign = sign(gene_logFC),
    concordance = case_when(
      expected_sign == observed_sign ~ "concordant",
      TRUE ~ "discordant"
    )
  )
# 计算每个基因在对应细胞类型中的表达比例
genes_use <- unique(plot_df$gene)
rna_counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
pct_expr_df <- lapply(genes_use, function(g) {
  expr_vec <- as.numeric(rna_counts[g, ] > 0)
  tibble(
    celltype = seu$celltype,
    expr = expr_vec
  ) %>%
    group_by(celltype) %>%
    summarise(
      pct_expr = mean(expr),
      .groups = "drop"
    ) %>%
    mutate(gene = g)
}) %>%
  bind_rows()
plot_df <- plot_df %>%
  left_join(pct_expr_df, by = c("celltype", "gene"))
write.csv(plot_df, file.path(res_dir, "supply_plot_df.csv"), row.names = FALSE)
```

```r
# 接续前文plot_herv_gene_corr_df
p_corr <- plot_herv_gene_corr_df(seu_list, plot_pairs_main, ncol = 2, is_plot = F)
# 基因DE的bubble plot
plot_df <- read.csv(file.path(res_dir, "supply_plot_df.csv"))
bubble_df <- plot_df %>%
  group_by(celltype, gene) %>%
  summarise(
    gene_logFC = first(gene_logFC),
    p_val_adj = first(p_val_adj),
    pct_expr = first(pct_expr),
    concordance = if(any(concordance == "discordant")) "discordant" else "concordant",
    .groups = "drop"
  ) %>%
  mutate(
    star = p_to_star(p_val_adj)
  ) %>%
  mutate(
    celltype = factor(celltype, levels = celltype_order),
    gene = factor(gene, levels = gene)
  )
bubble_df <- bubble_df %>% filter(gene %in% plot_specs$gene)
max_abs_gene <- max(abs(bubble_df$gene_logFC), na.rm = TRUE)
p_bubble <- ggplot(bubble_df, aes(x = celltype, y = gene)) +
  geom_point(
    aes(size = pct_expr, fill = gene_logFC),
    shape = 21, color = "black", stroke = 0.45
  ) +
  geom_point(
    data = bubble_df %>% filter(concordance == "discordant"),
    size = 14, shape = 22, fill = NA, color = "#DC2626", stroke = 1.15,
    show.legend = FALSE
  ) +
  geom_text(
    data = bubble_df %>% filter(star != ""),
    aes(label = star),
    size = 5, color = "black"
  ) +
  scale_size(range = c(4, 13)) +
  scale_fill_gradient2(
    low = "#3B82F6", mid = "white", high = "#EF4444",
    midpoint = 0,
    limits = c(-max_abs_gene, max_abs_gene),
  ) +
  labs(x = NULL, y = NULL,) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.6),
    axis.line.y.left = element_line(color = "black", linewidth = 0.6),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    plot.margin = margin(8, 8, 8, 8)
  )
# 拼图
# p_corr_panel <- wrap_elements(full = p_corr) + labs(tag = "A")
# p_bubble <- p_bubble + labs(tag = "B")
# p_final <- (p_corr_panel | p_bubble) &
#   theme(
#     plot.tag = element_text(size = 22, face = "bold"),
#     plot.tag.position = c(0, 1.08),
#     plot.margin = margin(40, 10, 20, 20),
#     plot.background = element_rect(fill = "white", colour = NA)
#   )
# ggsave(
#   file.path(res_dir, "genome_gene2.pdf"),
#   p_final, width = 18, height = 8
# )
ggsave(
  file.path(res_dir, "genome_gene2.pdf"),
  p_bubble + theme(plot.margin = margin(20, 10, 20, 20)), 
  width = 8, height = 6
)
```

#### 关键基因补充分析

```r
seu_list <- list(
  Astro = astro, 
  Oligo = oligo, 
  OPC = opc, 
  Excit = excit, 
  Inhib = inhib, 
  Endo = endo, 
  Mic = mic
)
target_herv <- c("HERVW-15q21.2", "HML1-1q32.1", "HML4-16p13.3")
target_gene <- c("GLDN", "ABCC9", "MTR")
celltype_order <- c("Astro", "Oligo", "OPC", "Excit", "Inhib", "Mic", "Endo")
herv_de_res <- lapply(seu_list, function(seu){
  ct <- identity_seu(seu)
  Idents(seu) <- "group"
  FindMarkers(
    seu,
    ident.1 = "AD", ident.2 = "NC",
    assay = "HERV",
    features = target_herv,
    test.use = "MAST",
    logfc.threshold = 0, min.pct = 0
  ) %>%
  rownames_to_column("feature") %>%
  mutate(celltype = ct)
}) %>% bind_rows()
herv_bubble_df <- herv_de_res %>%
  mutate(
    pct_max = pmax(pct.1, pct.2),
    star = p_to_star(p_val_adj),
  )
write.csv(herv_bubble_df, file.path(res_dir, "supply2_herv_DE.csv"), row.names = FALSE)
gene_de_res <- lapply(seu_list, function(seu){
  ct <- identity_seu(seu)
  Idents(seu) <- "group"
  FindMarkers(
    seu,
    ident.1 = "AD", ident.2 = "NC",
    assay = "RNA",
    features = target_gene,
    test.use = "MAST",
    logfc.threshold = 0, min.pct = 0
  ) %>%
  rownames_to_column("feature") %>%
  mutate(celltype = ct)
}) %>% bind_rows()
rna_counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
pct_expr_df <- lapply(target_gene, function(g) {
  expr_vec <- as.numeric(rna_counts[g, ] > 0)
  tibble(
    celltype = seu$celltype,
    expr = expr_vec
  ) %>%
    group_by(celltype) %>%
    summarise(
      pct_expr = mean(expr),
      .groups = "drop"
    ) %>%
    mutate(feature = g)
}) %>% bind_rows()
gene_bubble_df <- gene_de_res %>%
  left_join(pct_expr_df, by = c("celltype", "feature")) %>%
  mutate(star = p_to_star(p_val_adj))
write.csv(gene_bubble_df, file.path(res_dir, "supply2_gene_DE.csv"), row.names = FALSE)
target_pair <- list(
  "HERVW-15q21.2" = list("GLDN"),
  "HML1-1q32.1" = list("ABCC9"),
  "HML4-16p13.3" = list("MTR")
)
corr_res <- lapply(seu_list, function(seu){
  ct <- identity_seu(seu)
  get_gene_herv_corr_res_plus(seu, "HERV", target_pair, write_csv = F) %>%
    mutate(celltype = ct)
}) %>% bind_rows()
write.csv(corr_res, file.path(res_dir, "supply2_corr_res.csv"), row.names = FALSE)
corr_res_donor <- lapply(seu_list, function(seu){
  ct <- identity_seu(seu)
  get_gene_herv_corr_res_plus(seu, "HERV", target_pair, write_csv = F, is_donor = T) %>%
    mutate(celltype = ct)
}) %>% bind_rows()
write.csv(corr_res_donor, file.path(res_dir, "supply2_corr_res_donor.csv"), row.names = FALSE)
```

```r
tag_size_main <- 20
theme_bubble <- theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.6),
    axis.line.y.left   = element_line(color = "black", linewidth = 0.6),
    axis.line.x.top    = element_blank(),
    axis.line.y.right  = element_blank(),
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    plot.margin = margin(8, 8, 8, 8)
  )
# herv的DE结果
herv_bubble_df <- read.csv(file.path(res_dir, "supply2_herv_DE.csv")) %>%
  mutate(
    celltype = factor(celltype, levels = celltype_order),
    feature = factor(feature, levels = target_herv),
  )
max_abs_herv <- max(abs(herv_bubble_df$avg_log2FC), na.rm = TRUE)
p_herv_bubble <- ggplot(herv_bubble_df, aes(x = celltype, y = feature)) +
  geom_point(
    aes(size = pct_max, fill = avg_log2FC),
    shape = 21, color = "black", stroke = 0.45
  ) +
  geom_text(aes(label = star), size = 4.5) +
  scale_size(range = c(4, 12), name = "max pct") +
  scale_fill_gradient2(
    low = "#3B82F6", mid = "white", high = "#EF4444",
    midpoint = 0,
    limits = c(-max_abs_herv, max_abs_herv),
    name = "herv_log2FC"
  ) +
  labs(x = NULL, y = NULL) +
  theme_bubble
# gene的DE结果
gene_bubble_df <- read.csv(file.path(res_dir, "supply2_gene_DE.csv")) %>%
  mutate(
    celltype = factor(celltype, levels = celltype_order),
    feature = factor(feature, levels = target_gene),
  )
max_abs_gene <- max(abs(gene_bubble_df$avg_log2FC), na.rm = TRUE)
p_gene_bubble <- ggplot(gene_bubble_df, aes(x = celltype, y = feature)) +
  geom_point(
    aes(size = pct_expr, fill = avg_log2FC),
    shape = 21, color = "black", stroke = 0.45
  ) +
  geom_text(aes(label = star), size = 4.5) +
  scale_size(range = c(4, 12), name = "pct_expr") +
  scale_fill_gradient2(
    low = "#3B82F6", mid = "white", high = "#EF4444",
    midpoint = 0,
    limits = c(-max_abs_gene, max_abs_gene),
    name = "gene_log2FC"
  ) +
  labs(x = NULL, y = NULL) +
  theme_bubble
# 相关性图3*2
plot_specs <- tibble::tribble(
  ~ct,     ~herv,           ~gene,
  "Excit", "HERVW-15q21.2", "GLDN",
  "Mic",   "HERVW-15q21.2", "GLDN",
  "Astro", "HERVW-15q21.2", "GLDN",
  "Oligo", "HML4-16p13.3",  "MTR",
  "Astro", "HML4-16p13.3",  "MTR",
  "Inhib", "HML4-16p13.3",  "MTR",
)
p_corr <- plot_herv_gene_corr_df(seu_list, plot_specs, ncol = 3, is_plot = F)
# 拼图
p_gene_bubble <- p_gene_bubble + labs(tag = "A")
p_herv_bubble <- p_herv_bubble + labs(tag = "B")
p_corr_panel <- wrap_elements(full = p_corr) + labs(tag = "C")
final_plot <- (p_gene_bubble | p_herv_bubble) / p_corr_panel +
  plot_layout(heights = c(1, 2)) &
  theme(
    plot.tag = element_text(size = tag_size_main, face = "bold"),
    plot.tag.position = c(0, 1),
    plot.margin = margin(15, 15, 15, 15),
    plot.background = element_rect(fill = "white", colour = NA)
  )
ggsave(
  file.path(res_dir, "supply_gene_herv_DE_corr.pdf"),
  final_plot, width = 12, height = 12
)
```

### 结果表

#### 超家族-亚家族对应关系

```r
# 每个亚家族都属于哪个超家族，亚家族间以、分隔
ltr_class <- read.table("/public/home/GENE_proc/wth/metadata/hERV/LTR2Class.tsv", header=T, sep="\t") %>%
  mutate(repName = gsub("_", "-", repName))
features_use <- rownames(seu[["HERV_subfamily_scTE"]])
superfamily_subfamily_table <- ltr_class %>%
  filter(repName %in% features_use) %>%
  distinct(repClass, repName) %>%
  arrange(repClass, repName) %>%
  group_by(repClass) %>%
  summarise(
    subfamily = paste(repName, collapse = "、"),
    .groups = "drop"
  ) %>%
  rename(
    superfamily = repClass
  ) %>%
  mutate(superfamily = factor(superfamily, levels = c("ERV1", "ERVL", "ERVK"))) %>%
  arrange(superfamily)
write.csv(superfamily_subfamily_table, file.path(res_dir, "class_subfamily.csv"), row.names = FALSE)
```
