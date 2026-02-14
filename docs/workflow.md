flowchart TD
    classDef data fill:#2d2d2d,stroke:#fff,stroke-width:2px,color:#fff;
    classDef proc fill:#000,stroke:#fff,stroke-width:1px,color:#fff;
    classDef rstep fill:#1a1a1a,stroke:#ccc,stroke-width:1px,color:#eee;

    subgraph DATA ["INPUT DATA"]
        A[("Raw FASTQ")]:::data
        B[("Reference Transcriptome")]:::data
    end

    subgraph QC ["PRE-PROCESSING"]
        direction TB
        A --> C(FastQC):::proc
        C --> D(Trimmomatic):::proc
        D --> E(FastQC Post):::proc
    end

    subgraph QUANT ["QUANTIFICATION"]
        direction TB
        B --> F(Salmon Index):::proc
        F --> G(Salmon Quant):::proc
        D --> G
        G --> H[("quant.sf")]:::data
    end

    subgraph R_ANALYSIS ["ANALYSIS R"]
        direction TB
        H --> I(tximport):::rstep
        I --> J(DESeq2 Object):::rstep
        
        J --> K{PCA}:::rstep
        
        J --> L(Differential Gene Expression):::rstep
        L --> M[("DEGs Tables")]:::data
        M --> N(GSEA / KEGG / GO):::rstep
    end
