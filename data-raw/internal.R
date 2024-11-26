library(vroom)
library(dplyr)
library(fgsea)

PAM50genes <- vroom("data-raw/PAM50_symbol_entrez_2.csv",
                    col_types = cols(.default = col_character()))

PAM50genes <- PAM50genes %>%
  dplyr::select(entrez_id = ENTREZID, hugo_symbol = SYMBOL, ensembl_id = Gene.stable.ID,
         ensembl_id_version = Gene.stable.ID, hgnc_symbol_gene_id) %>%
  dplyr::mutate(ensembl_id = gsub("\\.(\\d)+", "", ensembl_id)) %>%
  as.data.frame()

hallmark_signatures <- gmtPathways("data-raw/h_rtk.all.v7.2.symbols.gmt")

breast_signatures <- gmtPathways("data-raw/c2.cgp.breast.v6.1.symbols.gmt")
signature_genes <- unique(c(unlist(breast_signatures, use.names = FALSE), unlist(hallmark_signatures, use.names = FALSE)))

bccb_signatures <- list(hallmark = hallmark_signatures, c2cgp_breast = breast_signatures)
genes_ESR1=c("ESR1","GATA3","FOXA1","KMT2C","KMT2D")
genes_AKT=c("EIF4EBP1","AKT1","AKT2","AKT3","AKT1S1","DEPDC5","DEPTOR","INPP4B","MAPKAP1","MLST8",
            "MTOR","NPRL2","NPRL3","PDK1","PIK3CA","PIK3CB","PIK3R1","PIK3R2","PIK3R3","PPP2R1A","PTEN",
            "RHEB","RICTOR","RPTOR","RPS6","RPS6KB1","STK11","TSC1","TSC2","PRR5","FLCN","PIK3C2B")
genes_MAPK=c("ABL1","SOS1","GRB2","PTPN11","KRAS","HRAS","NRAS","RIT1","ARAF","BRAF","RAF1","RAC1","MAP2K1",
             "MAP2K2","MAPK1","NF1","RASA1","CBL","ERRFI1","ERF","IRS2","CBLB","CBLC","IRS1","SOS2","SHC1",
             "SHC2","SHC3","SHC4","RASGRP1","RASGRP2","RASGRP3","RASGRP4","RAPGEF1","RAPGEF2","RASGRF1",
             "RASGRF2","FNTA","FNTB","RCE1","ICMT","MRAS","PLXNB1","MAPK3","ARHGAP35","RASA2","RASA3",
             "RASAL1","RASAL2","RASAL3","SPRED1","SPRED2","SPRED3","DAB2IP","SHOC2","PPP1CA","SCRIB",
             "PIN1","KSR1","KSR2","PEBP1","PEA15")
genes_RTK=c("EGFR","ERBB2","ERBB3","ERBB4","PDGFRA","PDGFRB","MET","FGFR1","FGFR2","FGFR3","FGFR4","FLT3",
            "ALK","RET","ROS1","KIT","IGF1R","NTRK1","NTRK2","NTRK3","JAK2","INSR","INSRR","FGF3")
genes_Prol=c("CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN2C","CCND1","CCND2","CCND3","CCNE1","CCNE2","CDK2",
             "CDK4","CDK6","RB1","E2F1","E2F3","AURKA","FAT1")
genes_Inmune=c("PDCD1", "CD274", "CTLA4")
genes_Onco <- vroom("data-raw/CDK46_Seth_gene_list_16177.ordered_3.plus_genes_from_Daniel.txt", delim = "\t") %>%
  pull("Hugo_Symbol")

bccb_genes <- unique(c(PAM50genes$hugo_symbol, genes_ESR1, genes_AKT, genes_MAPK,
                       genes_Prol, genes_RTK, genes_Inmune, genes_Onco))


genes <- vroom("data-raw/gene_dictionary.csv") %>%
  pull(final_name)
genes <- sample(grep("ENSG", genes, value = TRUE, invert = TRUE), 10000)

usethis::use_data(PAM50genes, bccb_genes, bccb_signatures, signature_genes, internal = TRUE, overwrite = TRUE)
