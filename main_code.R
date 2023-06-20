library(dplyr)
library(biomaRt) 
library(forcats)
library(patchwork)
library(ggplot2)
library(data.table)
library(tictoc)
library(fastman)
library(ggcorrplot)
library(corrplot)
library(matrixStats)
library(mvtnorm)
library (parallel)
library(CompQuadForm)
library(UpSetR)
library(stringr)
library(MTAR)
library(plyr)
#Const 
ALPHA <- 5e-8
SMALL <- 2.6e-324 # Not use .Machine$double.xmin.
set.seed(5)
# T2D  replaced with DIAMANTE-EUR. The latter has more data
# All pazoki data is the same with the previous data (GGT, ALT, GPL)
# Ghodsian_NAFLD data is the same with prev NAFLD 
# Mbatchou AST is chosen

## LOADING DATA ##
# ===== T2D =====
raw_t2d <-  data.table::fread("data/DIAMANTE-EUR.sumstat.txt", header = TRUE)
raw_t2d$`Fixed-effects_p-value` <- as.numeric(raw_t2d$`Fixed-effects_p-value`)
raw_t2d <- dplyr::rename(raw_t2d, CHR='chromosome(b37)', BP='position(b37)', 
                         SNP=chrposID,  EAF=effect_allele_frequency,
                         A1=effect_allele, A2=other_allele, 
                         BETA='Fixed-effects_beta', SE='Fixed-effects_SE',
                         P='Fixed-effects_p-value')
raw_t2d[, c("SNP", "A1", "A2") := 
          .(mapply(function(str) substr(str, 4, nchar(str)), SNP),
            mapply(toupper, A1),
            mapply(toupper, A2)
          )]

# ===== Fasting Glucagon ===== 
raw_gluc <-  data.table::fread("data/Stinson_FastingGlucagon_hg19.tsv", header = TRUE)
raw_gluc <- dplyr::rename(raw_gluc, CHR=chrom, BP=pos,
                         A1=Allele1, A2=Allele2, BETA=Effect, SE=StdErr,
                         P=`P-value`, EAF=Freq1)
raw_gluc[, c("A1", "A2") := 
          .(mapply(toupper, A1),
            mapply(toupper, A2)
          )]

# ===== BMI ===== 
raw_bmi <-  data.table::fread("data/Loh_BodyMassIndex_hg19.tsv", header = TRUE)
raw_bmi <- dplyr::rename(raw_bmi, CHR=chr, BP=pos,
                         A1=Allele1, A2=Allele2, BETA=Effect, SE=StdErr,
                         P=`P-value`, EAF=Freq1)
raw_bmi[, c("P", "A1", "A2") := 
           .(mapply(as.numeric, P),
             mapply(toupper, A1),
             mapply(toupper, A2)
           )]

# ===== ALT ===== (alanine transaminase ALT GCST90013405)
raw_alt <-  data.table::fread("data/Pazoki_ALT_hg19.tsv", header = TRUE)
raw_alt <- dplyr::rename(raw_alt, CHR=CHROM, BP=POS,
                         A1=EA, A2=OA, BETA=ES)

# ===== ALP ===== (Alkaline Phosphatat)
raw_apl <-  data.table::fread("data/Pazoki_APL_hg19.tsv", header = TRUE)
raw_apl <- dplyr::rename(raw_apl, CHR=CHROM, BP=POS,
                         A1=EA, A2=OA, BETA=ES)

# ===== GGT ===== gamma-glutamyl transferase (GGT GCST90013407) 
raw_ggt <-  data.table::fread("data/Pazoki_GGT_hg19.tsv", header = TRUE)
raw_ggt <- dplyr::rename(raw_ggt, CHR=CHROM, BP=POS,
                         A1=EA, A2=OA, BETA=ES)

# ===== NAFLD ===== (Non-alcoholic fatty liver disease)
raw_nafld <-  data.table::fread("data/Ghodsian_NAFLD_hg19.tsv", header = TRUE)
raw_nafld <- dplyr::rename(raw_nafld, CHR=CHROM, BP=POS,A1=EA, A2=OA, BETA=ES)


### ===== AST ===== (aspartate aminotransferase)
# raw_ast <-  data.table::fread("data/GCST90013996.tsv", header = TRUE)
# raw_ast <- dplyr::rename(raw_ast, CHR=chromosome, BP=base_pair_location, P=p_value,
#                          A1=effect_allele, A2=other_allele, BETA=beta,
#                          SE=standard_error)

# ===== Mbatchou AST ===== (aspartate aminotransferase)
mbc_raw_ast <-  data.table::fread("data/Mbatchou_AST_hg19.tsv", header = TRUE)
mbc_raw_ast <- dplyr::rename(mbc_raw_ast, CHR=CHROM, BP=POS,
                             A1=EA, A2=OA, BETA=ES)

# # ===== Mahajan AST ===== (aspartate aminotransferase)
# mhj_raw_ast <-  data.table::fread("data/Mahajan_AST_hg19.tsv", header = TRUE)
# mhj_raw_ast <- dplyr::rename(mhj_raw_ast, CHR=CHROM, BP=POS,
#                              A1=EA, A2=OA, BETA=ES)


## DATA TIDYING ##
## Sanity check and data tidying
gwas_sanity <- function(gwas_data, data_name) {
  # Set key for faster lookup
  gwas_data |>
    data.table::setkey(CHR, BP, A1, A2)
  
  # Filtering sex chromosome 
  if (!is.numeric(gwas_data$CHR)){
    gwas_data <- gwas_data[CHR %in% c(1:22)]
    gwas_data$CHR <- as.numeric(gwas_data$CHR)
  }
  else{
    gwas_data <- gwas_data[CHR<=22]
  }
  
  # Check duplicate SNPs
  if (nrow(gwas_data[duplicated(gwas_data, by = key(gwas_data))]) > 0) {
    gwas_data <- gwas_data[!duplicated(gwas_data, by = key(gwas_data))]
  }
  
  # Replace zero P-value with SMALL
  gwas_data[P==0, P := SMALL]
  
  # Retrieve SNP coloumn
  if (!("SNP" %in% colnames(gwas_data))) {
    gwas_data$SNP <- paste(gwas_data$CHR,gwas_data$BP,sep=":")
  }
  
  # Retrieve Z-score (B/SE)
  if (!(paste("Z", data_name, sep="_") %in% colnames(gwas_data))) {
    gwas_data$Z <- (gwas_data$BETA / gwas_data$SE)
    gwas_data <- dplyr::rename(gwas_data, !!sym(paste("Z", data_name, sep="_")):=Z)
  }
  
  return(gwas_data)
}

tic()
raw_t2d <- gwas_sanity(raw_t2d, "t2d") 
raw_nafld <- gwas_sanity(raw_nafld, "nafld")
raw_ggt <- gwas_sanity(raw_ggt, "ggt")
raw_alt <- gwas_sanity(raw_alt, "alt")
raw_apl <- gwas_sanity(raw_apl, "apl")
raw_gluc <- gwas_sanity(raw_gluc, "gluc")
raw_bmi <- gwas_sanity(raw_bmi, "bmi")

mbc_raw_ast <- gwas_sanity(mbc_raw_ast, "ast")
# mhj_raw_ast <- gwas_sanity(mhj_raw_ast, "ast")
toc()

## CREATING PLOTS ## 
## Number of Significant (5e-08) SNPs per traits ##
total_signf_snps <- data.table::data.table(gwas_data=c("GGT",
                                                       "ALT",
                                                       "ALP",
                                                       "NAFLD",
                                                       "T2D",
                                                       "AST",
                                                       "GLUC",
                                                       "BMI"
), 
snp_count=c(
  raw_ggt[P<5e-08, length(BP)],
  raw_alt[P<5e-08, length(BP)],
  raw_apl[P<5e-08, length(BP)],
  raw_nafld[P<5e-08, length(BP)],
  raw_t2d[P<5e-08, length(BP)],
  mbc_raw_ast[P<5e-08, length(BP)],
  raw_gluc[P<5e-08, length(BP)],
  raw_bmi[P<5e-08, length(BP)]
)
)
ggplot(total_signf_snps, aes(x=reorder(gwas_data, -snp_count), y=snp_count)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  geom_text(aes(label=snp_count), position=position_dodge(width=0.9), vjust=-0.25) + scale_y_continuous(limits=c(0,27500000)) +
  theme_minimal() +
  labs(title="Significant SNPs for Each GWAS Summary Statistics",
       y = "SNPs Count", x = "GWAS Data")

## CREATING PLOTS ## 
## Number of SNPs per traits ##
tic()
total_snps <- data.table::data.table(gwas_data=c("GGT",
                                                 "ALT",
                                                 "ALP",
                                                 "NAFLD",
                                                 "T2D",
                                                 "AST",
                                                 "GLUC",
                                                 "BMI"
), 
snp_count=c(
  raw_ggt[, length(BP)],
  raw_alt[, length(BP)],
  raw_apl[, length(BP)],
  raw_nafld[, length(BP)],
  raw_t2d[, length(BP)],
  mbc_raw_ast[, length(BP)],
  raw_gluc[, length(BP)],
  raw_bmi[, length(BP)]
)
)
ggplot(total_snps, aes(x=reorder(gwas_data, -snp_count), y=snp_count)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  geom_text(aes(label=snp_count), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_minimal() +
  labs(title="Total SNPs for Each GWAS Summary Statistics", 
       y = "SNP Count", x = "GWAS Data")


# Get Signnificant SNPs within a window
get_top_snps <- function(msnp, annotationWinMb=0.5, annotatePval=5e-08){
  if (!any(msnp$P<annotatePval)) {
    return(msnp[0,])
  } else {
    msnp <- msnp[P<annotatePval]
    print(ncol(msnp))
    m0=NULL; sunc=unique(msnp$CHR);
    for (i in sunc) { # loop over CHR
      m2=msnp[msnp$CHR==i,]; f=order(-log10(m2$P),decreasing=TRUE); 
      m2=m2[f,]; m1=m2[1,,drop=FALSE];  
      if (nrow(m2)>1) { 
        for (j in 2:nrow(m2)) {  
          f=abs(m1$BP-m2$BP[j]); 
          if (min(f)>(annotationWinMb*10^(6))) { 
            m1=rbind(m1,m2[j,,drop=FALSE]); 
          } 
        } 
      }
      m0=rbind(m0,m1); # identify top SNP within annotationWinMb window for each CHR
    }
    msnp=m0 # store information of only the top SNP within window for each CHR
    rm(m0,m1,m2)
    return(msnp)
  }
}

## CREATING PLOTS ## 
## Number of Windowed (500KB) Significant (5e-08) SNPs per traits ##
total_signf_snps <- data.table::data.table(gwas_data=c("GGT",
                                                       "ALT",
                                                       "APL",
                                                       "NAFLD",
                                                       "T2D",
                                                       "AST",
                                                       "GLUC",
                                                       "BMI"
), 
snp_count=c(
  get_top_snps(raw_ggt)[, length(BP)],
  get_top_snps(raw_alt)[, length(BP)],
  get_top_snps(raw_apl)[, length(BP)],
  get_top_snps(raw_nafld)[, length(BP)],
  get_top_snps(raw_t2d)[, length(BP)],
  get_top_snps(mbc_raw_ast)[, length(BP)],
  get_top_snps(raw_gluc)[, length(BP)],
  get_top_snps(raw_bmi)[, length(BP)]
)
)
ggplot(total_signf_snps, aes(x=reorder(gwas_data, -snp_count), y=snp_count)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  geom_text(aes(label=snp_count), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_minimal() +
  coord_cartesian(ylim = c(0,1800)) +
  labs(title="Significant Loci for Each GWAS Data",
       subtitle="Loci within 500KB distance considered as one locus",
       y = "Loci Count", x = "GWAS Data")
toc()


## CRATE PLOT ##
## Create manhattan plot for each trait ## 
# Print function
print_manhattan_plot <- function(DT, title, lower=0, upper=200, cex.text=0.4, WinMb=3){
  if (upper < -log10(5e-08)) {
    upper = 15
  } else if (upper == -log10(SMALL)) {
    upper = 350
  } else {
    upper = upper + 20
  }
  fastman::fastman(DT, maxP=upper, annotateHighlight=TRUE, annotationWinMb=WinMb,
                   main=title, suggestiveline=-log10(5e-06),
                   genomewideline=-log10(5e-08), annotatePval=-log10(5e-08), 
                   colAbovePval=TRUE, ylim=c(lower,upper), cex.axis=1, 
                   cex.text=cex.text)
}

cat(min(raw_ggt$P) ,"\n",
    min(raw_alt$P) ,"\n",
    min(raw_apl$P) ,"\n",
    min(raw_nafld$P) ,"\n",
    min(raw_t2d$P) ,"\n",
    min(mbc_raw_ast$P) ,"\n",
    min(raw_gluc$P) ,"\n",
    min(raw_bmi$P) ,"\n")

raw_ggt[P == 0, P := SMALL]
raw_alt[P == 0, P := SMALL]
raw_apl[P == 0, P := SMALL]
raw_nafld[P == 0, P := SMALL]
raw_t2d[P == 0, P := SMALL]
mbc_raw_ast[P == 0, P := SMALL]
raw_bmi[P == 0, P := SMALL]

tic()
print_manhattan_plot(raw_ggt, "Manhattan Plot for GGT GWAS Data", upper=-log10(min(raw_ggt$P)), cex.text=0.9, WinMb=66)
print_manhattan_plot(raw_alt, "Manhattan Plot for ALT GWAS Data", upper=-log10(min(raw_alt$P)), cex.text=0.9, WinMb=66)
print_manhattan_plot(raw_apl, "Manhattan Plot for ALP GWAS Data", upper=-log10(min(raw_apl$P)), cex.text=0.9, WinMb=66)
print_manhattan_plot(raw_nafld, "Manhattan Plot for NAFLD GWAS Data", upper=-log10(min(raw_nafld$P)), cex.text=0.9)
print_manhattan_plot(raw_t2d, "Manhattan Plot for T2D GWAS Data", upper=-log10(min(raw_t2d$P)), cex.text=0.9, WinMb=54)
print_manhattan_plot(mbc_raw_ast, "Manhattan Plot for AST GWAS Data", upper=-log10(min(mbc_raw_ast$P)), cex.text=0.9, WinMb=60)
print_manhattan_plot(raw_gluc, "Manhattan Plot for Fasting Glucagon GWAS Data", upper=-log10(min(raw_gluc$P)), cex.text=0.9)
print_manhattan_plot(raw_bmi, "Manhattan Plot for BMI GWAS Data", upper=-log10(min(raw_bmi$P)), cex.text=0.9, WinMb=66)
toc()


## CREATE JOINED Z MATRIX ##
create_z <- function(df1, df2) {
  result <- merge(df1, df2, all = TRUE)
  
  # Flipping the RHS
  df2$"A1_flip" <- df2$A2
  df2$"A2_flip" <- df2$A1
  
  df2 |>
    data.table::setkey(CHR, BP, A1_flip, A2_flip)
  
  flip <- (result[df2, on = .(CHR==CHR, BP == BP, A1 == A1_flip, A2 == A2_flip), nomatch=NULL])
  
  y <- colnames(df2)[length(df2)-2]
  flip[[y]] = flip[[paste0("i.",y)]]
  flip2 <- data.table::setkey(data.table::copy(flip), CHR, BP, i.A1, i.A2)
  
  result <- result[!flip]
  result <- result[!flip2]
  sel = c(1:length(result))
  result <- data.table::rbindlist(list(result, flip[, ..sel]))
  # Clearning memory
  rm (df2, flip, y, flip2, sel)
  return(result)
}

#Check all of the combination of pairwise data 
# tic()
# ggt_t2d <- create_z(raw_ggt, raw_t2d)
# ggt_t2d_count <- ggt_t2d[!is.na(Z_ggt) & !is.na(Z_t2d), length(CHR)]
# alt_t2d <- create_z(raw_alt, raw_t2d)
# alt_t2d_count <- alt_t2d[!is.na(Z_alt) & !is.na(Z_t2d), length(CHR)]
# apl_t2d <- create_z(raw_apl, raw_t2d)
# apl_t2d_count <- apl_t2d[!is.na(Z_apl) & !is.na(Z_t2d), length(CHR)]
# nafld_t2d <- create_z(raw_nafld, raw_t2d)
# nafld_t2d_count <- nafld_t2d[!is.na(Z_nafld) & !is.na(Z_t2d), length(CHR)]
# 
# mbc_ast_t2d <- create_z(mbc_raw_ast, raw_t2d)
# mbc_ast_t2d_count <- mbc_ast_t2d[!is.na(Z_ast) & !is.na(Z_t2d), length(CHR)]
# mbc_ast_ggt <- create_z(mbc_raw_ast, raw_ggt)
# mbc_ast_ggt_count <- mbc_ast_ggt[!is.na(Z_ast) & !is.na(Z_ggt), length(CHR)]
# mbc_ast_alt <- create_z(mbc_raw_ast, raw_alt)
# mbc_ast_alt_count <- mbc_ast_alt[!is.na(Z_ast) & !is.na(Z_alt), length(CHR)]
# mbc_ast_apl <- create_z(mbc_raw_ast, raw_apl)
# mbc_ast_apl_count <- mbc_ast_apl[!is.na(Z_ast) & !is.na(Z_apl), length(CHR)]
# mbc_ast_nafld <- create_z(mbc_raw_ast, raw_nafld)
# mbc_ast_nafld_count <- mbc_ast_nafld[!is.na(Z_ast) & !is.na(Z_nafld), length(CHR)]
# 
# mhj_ast_t2d <- create_z(mhj_raw_ast, raw_t2d)
# mhj_ast_t2d_count <- mhj_ast_t2d[!is.na(Z_ast) & !is.na(Z_t2d), length(CHR)]
# mhj_ast_ggt <- create_z(mhj_raw_ast, raw_ggt)
# mhj_ast_ggt_count <- mhj_ast_ggt[!is.na(Z_ast) & !is.na(Z_ggt), length(CHR)]
# mhj_ast_alt <- create_z(mhj_raw_ast, raw_alt)
# mhj_ast_alt_count <- mhj_ast_alt[!is.na(Z_ast) & !is.na(Z_alt), length(CHR)]
# mhj_ast_apl <- create_z(mhj_raw_ast, raw_apl)
# mhj_ast_apl_count <- mhj_ast_apl[!is.na(Z_ast) & !is.na(Z_apl), length(CHR)]
# mhj_ast_nafld <- create_z(mhj_raw_ast, raw_nafld)
# mhj_ast_nafld_count <- mhj_ast_nafld[!is.na(Z_ast) & !is.na(Z_nafld), length(CHR)]
# 
# 
# 
# pairwise_Z = data.table(
#   pairwise_name = c("ggt_t2d","alt_t2d","apl_t2d","nafld_t2d",
#                     "mbc_ast_t2d","mbc_ast_ggt","mbc_ast_alt","mbc_ast_apl",
#                     "mbc_ast_nafld","mhj_ast_t2d","mhj_ast_ggt","mhj_ast_alt",
#                     "mhj_ast_apl","mhj_ast_nafld"),
#   count_non_na = c(ggt_t2d_count,alt_t2d_count,apl_t2d_count,
#                    nafld_t2d_count,mbc_ast_t2d_count,mbc_ast_ggt_count,
#                    mbc_ast_alt_count,mbc_ast_apl_count,mbc_ast_nafld_count,
#                    mhj_ast_t2d_count,mhj_ast_ggt_count,mhj_ast_alt_count,
#                    mhj_ast_apl_count,mhj_ast_nafld_count )
# )
# print(pairwise_Z)
# toc()

# --- Choosing the Mbatchou AST ---
# Create the joined sum stat data
tic()
ggt_t2d <- create_z(raw_ggt[, .(CHR, BP, A1, A2, Z_ggt)], raw_t2d[, .(CHR, BP, A1, A2, Z_t2d)] )
setkey(ggt_t2d, CHR, BP, A1, A2)

ggt_t2d_nalfd <- create_z(ggt_t2d, raw_nafld[, .(CHR, BP, A1, A2, Z_nafld)])
setkey(ggt_t2d_nalfd, CHR, BP, A1, A2)
rm(ggt_t2d)

ggt_t2d_nalfd_ast <- create_z(ggt_t2d_nalfd, mbc_raw_ast[, .(CHR, BP, A1, A2, Z_ast)])
setkey(ggt_t2d_nalfd_ast, CHR, BP, A1, A2)
rm(ggt_t2d_nalfd)

ggt_t2d_nalfd_ast_alt <- create_z(ggt_t2d_nalfd_ast, raw_alt[, .(CHR, BP, A1, A2, Z_alt)])
setkey(ggt_t2d_nalfd_ast_alt, CHR, BP, A1, A2)
rm(ggt_t2d_nalfd_ast)

ggt_t2d_nalfd_ast_alt_apl <- create_z(ggt_t2d_nalfd_ast_alt, raw_apl[, .(CHR, BP, A1, A2, Z_apl)])
setkey(ggt_t2d_nalfd_ast_alt_apl, CHR, BP, A1, A2)
rm(ggt_t2d_nalfd_ast_alt)

ggt_t2d_nalfd_ast_alt_apl_gluc <- create_z(ggt_t2d_nalfd_ast_alt_apl, raw_gluc[, .(CHR, BP, A1, A2, Z_gluc)])
setkey(ggt_t2d_nalfd_ast_alt_apl_gluc, CHR, BP, A1, A2)

ggt_t2d_nalfd_ast_alt_apl_gluc_bmi <- create_z(ggt_t2d_nalfd_ast_alt_apl_gluc, raw_bmi[, .(CHR, BP, A1, A2, Z_bmi)])
setkey(ggt_t2d_nalfd_ast_alt_apl_gluc_bmi, CHR, BP, A1, A2)
rm(ggt_t2d_nalfd_ast_alt_apl_gluc)
toc()


## EMPIRICAL CORRELATION MATRIX ##                               
tic()
ggt_t2d_nalfd_ast_alt_apl_gluc_bmi <- 
  dplyr::relocate(ggt_t2d_nalfd_ast_alt_apl_gluc_bmi, Z_t2d, .before=Z_ggt)        
sel_range <- c(5 : length(ggt_t2d_nalfd_ast_alt_apl_gluc_bmi))
z_mat <- data.matrix(ggt_t2d_nalfd_ast_alt_apl_gluc_bmi[, ..sel_range])
z_mat_cut_05 <- data.matrix((ggt_t2d_nalfd_ast_alt_apl_gluc_bmi[abs(Z_t2d)<0.5 & abs(Z_ggt)<0.5 & abs(Z_nafld)<0.5 & abs(Z_ast)<0.5 & abs(Z_alt)<0.5 & abs(Z_apl)<0.5 & abs(Z_gluc)<0.5 & abs(Z_bmi)<0.5])[, ..sel_range])
z_mat_cut_1 <- data.matrix((ggt_t2d_nalfd_ast_alt_apl_gluc_bmi[abs(Z_t2d)<1 & abs(Z_ggt)<1 & abs(Z_nafld)<1 & abs(Z_ast)<1 & abs(Z_alt)<1 & abs(Z_apl)<1 & abs(Z_gluc)<1 & abs(Z_bmi)<1])[, ..sel_range])
z_mat_cut_2 <- data.matrix((ggt_t2d_nalfd_ast_alt_apl_gluc_bmi[abs(Z_t2d)<2 & abs(Z_ggt)<2 & abs(Z_nafld)<2 & abs(Z_ast)<2 & abs(Z_alt)<2 & abs(Z_apl)<2 & abs(Z_gluc)<2 & abs(Z_bmi)<2])[, ..sel_range])

# Using the naive correlation matrix
z_cor <- cor(z_mat, use = "complete.obs")
z_cor_cut_05 <- cor(z_mat_cut_05, use = "complete.obs")
z_cor_cut_1 <- cor(z_mat_cut_1, use = "complete.obs")
z_cor_cut_2 <- cor(z_mat_cut_2, use = "complete.obs")
# Plot correlation
corrplot(z_cor,method = 'color',type = 'upper', addCoef.col = 'black', 
         col = COL2('PRGn'), title = "Correlation Matrix",
         number.digits = 4,  mar=c(0,0,1,0)#Fix the title position
         )
corrplot(z_cor_cut_05,method = 'color',type = 'upper', addCoef.col = 'black', 
         col = COL2('PRGn'), title = "|Z|<0.5 Correlation Matrix",
         number.digits = 3,  mar=c(0,0,1,0)#Fix the title position
)
corrplot(z_cor_cut_1,method = 'color',type = 'upper', addCoef.col = 'black', 
         col = COL2('PRGn'), title = "|Z|<1 Correlation Matrix",
         number.digits = 3,  mar=c(0,0,1,0)#Fix the title position
)
corrplot(z_cor_cut_2,method = 'color',type = 'upper', addCoef.col = 'black', 
         col = COL2('PRGn'), title = "|Z|<2 Correlation Matrix",
         number.digits = 3,  mar=c(0,0,1,0)#Fix the title position
)
# cut z_mat[, c(1,2,3,5,8)]
toc()

tic()
omnibus_test <- function(Z, Sigma) {
  # Q follows chisq distribution
  Q <- rowSums(Z %*% solve(Sigma) * Z)
  # # Return the significant index
  # tm = qchisq(5e-8,dim(z_cor)[1], lower=FALSE)
  # return(which(Q>tm))
  p_omni <- pchisq(Q, dim(Sigma)[1], lower=FALSE)
  rm(Q)
  return(p_omni)
}

sumZ_test <- function(Z, Sigma) {
  sum_of_z = (matrixStats::rowSums2(Z)^2)/sum(Sigma)
  sumZ_P <- pchisq(sum_of_z, 1, lower.tail = FALSE)
  rm(sum_of_z)
  return(sumZ_P)
}


sumZ2_test <- function(Z, Sigma) {
  eig_vals = eigen(Sigma,sym=TRUE,only.val=TRUE)$val
  sum_of_z_2 = matrixStats::rowSums2(Z^2)
  sumZ2_P <- CompQuadForm::liu(sum_of_z_2, eig_vals)
  return(sumZ2_P)
}

pc_test <- function(Z, Sigma){
  eig_decomp <- eigen(Sigma)
  eig_vals <- eig_decomp$values
  eig_vecs <- eig_decomp$vectors
  
  pc_P <- pchisq((Z%*%eig_vecs[,1])^2/eig_vals[1], 1, lower.tail = FALSE)
  rm(eig_decomp, eig_vals, eig_vecs)
  return(c(pc_P))
}



minP_test <- function(Z, Sigma) {
  n_traits = ncol(Sigma)
  bound = matrixStats::rowMaxs(abs(Z))
  
  CLUSTERS <- makeCluster(parallel::detectCores()-1)
  minP <- parSapply(CLUSTERS, bound, function(i){
    if (is.na(i)) {
      return(NA)
    } else{
      return(1 - mvtnorm::pmvnorm(lower = rep(-i, n_traits),
                              upper = rep(i,n_traits),
                              sigma = Sigma,
                              keepAttr = FALSE)
      )
    }
  })
  stopCluster(CLUSTERS)
  return(minP)
}

# Parallel of MUSAT with socket
par_amatz <- function(Z, Sigma) {
  result <- c()
  row_num <- nrow(Z) 
  for (i in seq(1, row_num, by = 7000)) {
    # Define the end index for the current chunk
    end_idx <- min(i + 7000 - 1, row_num)
    tic()
    print(paste("Processing Row", i, "to", end_idx))
    
    CLUSTERS <- makeCluster(parallel::detectCores()-1)
    at_P <- parSapply(CLUSTERS, i:end_idx, function(row, Sigma){
      # if (any(is.na(Z[row,]))) {
      #   # cat(row, "NA")
      #   return(NA)
      # } else{
      #   # cat(row, "calculate MTAR")
      #   return(MTAR::emats(Z[row,], Sigma)$p.value)
      # }
      return(MTAR::amatz(Z[row,], Sigma)$p.value['A'])
    }, Sigma=Sigma)
    stopCluster(CLUSTERS)
    result <- append(result, at_P)
    toc()
  }
  return(result)
}

# Parallel of AT with socket
par_emats <- function(Z, Sigma) {
  CLUSTERS <- makeCluster(parallel::detectCores()-4)
  at_P <- parSapply(CLUSTERS, 1:nrow(Z), function(row, Sigma){
    if (any(is.na(Z[row,]))) {
      # cat(row, "NA")
      return(NA)
    } else{
      # cat(row, "calculate MTAR")
      return(MTAR::emats(Z[row,], Sigma)$p.value)
    }
  }, Sigma=Sigma)
  stopCluster(CLUSTERS)
  return(at_P)
}

# Parallel of AT with socket
par_emats3 <- function(Z, Sigma) {
  CLUSTERS <- makeCluster(parallel::detectCores()-4)
  at_P <- parSapply(CLUSTERS, 1:nrow(Z), function(row, Sigma){
    return(MTAR::emats(Z[row,], Sigma)$p.value)
  }, Sigma=Sigma)
  stopCluster(CLUSTERS)
  return(at_P)
}

par_emats2 <- function(Z, Sigma) {
  result <- c()
  row_num <- nrow(Z) 
  for (i in seq(1, row_num, by = 7000)) {
    # Define the end index for the current chunk
    end_idx <- min(i + 7000 - 1, row_num)
    tic()
    print(paste("Processing Row", i, "to", end_idx))
    
    CLUSTERS <- makeCluster(parallel::detectCores()-1)
    at_P <- parSapply(CLUSTERS, i:end_idx, function(row, Sigma){
      # if (any(is.na(Z[row,]))) {
      #   # cat(row, "NA")
      #   return(NA)
      # } else{
      #   # cat(row, "calculate MTAR")
      #   return(MTAR::emats(Z[row,], Sigma)$p.value)
      # }
      return(MTAR::emats(Z[row,], Sigma)$p.value)
    }, Sigma=Sigma)
    stopCluster(CLUSTERS)
    result <- append(result, at_P)
    toc()
  }
  return(result)
}

for (row in c(1:17001)) {
  result <- c()
  if (any(is.na(z_mat[row,]))) {
    # cat(row, "NA")
    cat("NA ")
  } else{
    cat("Z vector", z_mat[row,])
    cat(MTAR::emats(z_mat[row,], z_cor)$p.value)
  }
}


non_NA_z_mat <-as.matrix(na.omit(data.table::data.table(z_mat)))
non_NA_DT <- ggt_t2d_nalfd_ast_alt_apl_gluc_bmi[(!is.na(Z_ggt) 
                                                 & !is.na(Z_t2d) & !is.na(Z_nafld) 
                                                 & !is.na(Z_ast) & !is.na(Z_alt) 
                                                 & !is.na(Z_apl) & !is.na(Z_gluc) 
                                                 & !is.na(Z_bmi))]

tic()
P_at_1m <- par_emats2(non_NA_z_mat[1:1000000,], z_cor)
toc()
P_at_1m_DT <- get_SNP_DT_ref(P_at_1m, non_NA_DT[, .(CHR, BP, A1, A2)][1:1000000,])
P_at_1m_DT[P==0, P := SMALL]
print_manhattan_plot(P_at_1m_DT, "Manhattan Plot for Joint GWAS Data (AT)", upper=-log10(min(P_at_1m_DT[!is.na(P)]$P)), cex.text=0.9, WinMb=66)
# est. 807600

ggplot(melt(as.data.table(non_NA_z_mat[P_at_1m_DT$P==0,])), aes(x=variable, y=value)) + geom_boxplot()
tic()
P_at_11m <- par_emats2(non_NA_z_mat[1000000:1100000,], z_cor)
P_at_11m <- P_at_11m[2:length(P_at_11m)]
toc()
tic()
P_at_12m <- par_emats2(non_NA_z_mat[1100001:2000000,], z_cor)
toc()
tic()
P_at_3m <- par_emats2(non_NA_z_mat[2000001:3000000,], z_cor)
toc()



tic()
# Return the P-value
# All rows are kept to get the SNPs id
P_omni <- omnibus_test(z_mat, z_cor)
##P_at <- par_emats(z_mat, z_cor) # P_at <- par_emats(z_mat[9000:10000,], z_cor)x
 P_sumZ <- sumZ_test(z_mat, z_cor)
P_sumZ2 <- sumZ2_test(z_mat, z_cor)
P_pc <- pc_test(z_mat, z_cor)
minP <-minP_test(z_mat, z_cor)
minP2 <-minP_test(z_mat, z_cor)
toc()

# Get the SNP ref for multi trait summ stat for further analysis
get_SNP_DT_ref <- function(summ_P_stat, SNP_ref){
  DT <- data.table(SNP=paste(SNP_ref$CHR, SNP_ref$BP, sep = ":"), CHR=SNP_ref$CHR, BP=SNP_ref$BP, P=summ_P_stat)
  DT |>
    setkey(SNP)
  return(DT)
}

P_omni_DT <- get_SNP_DT_ref(P_omni, ggt_t2d_nalfd_ast_alt_apl_gluc_bmi)
P_omni_DT[is.infinite(P) | P==0, P := SMALL]

minP_DT <- get_SNP_DT_ref(minP, ggt_t2d_nalfd_ast_alt_apl_gluc_bmi)
minP_DT[is.infinite(P) | P==0, P := SMALL]

P_sumZ_DT <- get_SNP_DT_ref(P_sumZ, ggt_t2d_nalfd_ast_alt_apl_gluc_bmi)
P_sumZ_DT[is.infinite(P) | P==0, P := SMALL]

P_sumZ2_DT <- get_SNP_DT_ref(P_sumZ2, ggt_t2d_nalfd_ast_alt_apl_gluc_bmi)
P_sumZ2_DT[is.infinite(P) | P==0, P := SMALL]

P_pc_DT <- get_SNP_DT_ref(P_pc, ggt_t2d_nalfd_ast_alt_apl_gluc_bmi)
P_pc_DT[is.infinite(P) | P==0, P := SMALL]

cat(min(P_omni_DT[!is.na(P), P]) ,"\n",
    min(minP_DT[!is.na(P), P]) ,"\n",
    min(P_sumZ_DT[!is.na(P), P]) ,"\n",
    min(P_sumZ2_DT[!is.na(P), P]) ,"\n",
    min(P_pc_DT[!is.na(P), P]) ,"\n")

rm(P_omni)
rm(minP)
rm(P_sumZ)
rm(P_sumZ2)
rm(P_pc)

## PLOTTING THE MANHATTAN PLOT ##
print_manhattan_plot(P_omni_DT, "Manhattan Plot for Joint GWAS Data (OT)", upper=-log10(min(P_omni_DT[!is.na(P)]$P)), cex.text=0.9, WinMb=66)
print_manhattan_plot(minP_DT[-log10(P)<300], "Manhattan Plot for Joint GWAS Data (minP)", upper=-log10(min(minP_DT[-log10(P)<300][!is.na(P)]$P)), cex.text=0.9, WinMb=66)
print_manhattan_plot(P_sumZ_DT, "Manhattan Plot for Joint GWAS Data (SZ)", upper=-log10(min(P_sumZ_DT[!is.na(P)]$P)), cex.text=0.9, WinMb=66)
print_manhattan_plot(P_sumZ2_DT, "Manhattan Plot for Joint GWAS Data (SZ2)", upper=-log10(min(P_sumZ2_DT[!is.na(P)]$P)), cex.text=0.9, WinMb=66)
print_manhattan_plot(P_pc_DT, "Manhattan Plot for Joint GWAS Data (ET)", upper=-log10(min(P_pc_DT[!is.na(P)]$P)), cex.text=0.9, WinMb=66)




# Create the significant loci DT (500KB)  
P_omni_signf_loci <- get_top_snps(P_omni_DT[!is.na(P)])
minP_signf_loci <- get_top_snps(minP_DT[!is.na(P)])
SZ_signf_loci <- get_top_snps(P_sumZ_DT[!is.na(P)])
SZ2_signf_loci <- get_top_snps(P_sumZ2_DT[!is.na(P)])
ET_signf_loci <- get_top_snps(P_pc_DT[!is.na(P)])

# Plotting the barchart of loci distribution 
join_signf_loci <- data.table::data.table(gwas_data=c("OT",
                                                 # "minP",
                                                 "SZ",
                                                 "SZ2",
                                                 "ET"
), 
snp_count=c(
  P_omni_signf_loci[, length(BP)],
  # minP_signf_loci[, length(BP)],
  SZ_signf_loci[, length(BP)],
  SZ2_signf_loci[, length(BP)],
  ET_signf_loci[, length(BP)]
)
)
ggplot(join_signf_loci2, aes(x=reorder(gwas_data, -snp_count), y=snp_count)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  geom_text(aes(label=snp_count), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_minimal() +
  coord_cartesian(ylim = c(0,1800)) +
  labs(title="Significant Loci for Each Joint Test",
       subtitle="Loci within 500KB distance considered as one locus",
       y = "Loci Count", x = "Joint Test Method") 


# Plotting the venn
upset(fromList(list(OT=P_omni_signf_loci[,SNP], 
                    # minP=minP_signf_loci[,SNP],
                    SZ=SZ_signf_loci[,SNP],
                    SZ2=SZ2_signf_loci[,SNP],
                    ET=ET_signf_loci[,SNP])),
      order.by="freq",
      mainbar.y.label="Loci intersection",
      sets.x.label="Significant Loci per Method")


setkey(P_omni_signf_loci, CHR, BP)
setkey(minP_signf_loci, CHR, BP)
setkey(SZ_signf_loci, CHR, BP)
setkey(SZ2_signf_loci, CHR, BP)
setkey(ET_signf_loci, CHR, BP)


## Find Interesting Loci ## 
## Loci which detected from joint analyses but not in all original GWAS 

# Non significant SNPs
nss_ggt <- raw_ggt[P>=5e-08, SNP]
nss_alt <- raw_alt[P>=5e-08, SNP]
nss_nafld <- raw_nafld[P>=5e-08, SNP]
nss_t2d <- raw_t2d[P>=5e-08, SNP]
nss_ast <- mbc_raw_ast[P>=5e-08, SNP]
nss_gluc <- raw_gluc[P>=5e-08, SNP]
nss_bmi <- raw_bmi[P>=5e-08, SNP]
nss_apl <- raw_apl[P>=5e-08, SNP]

# non significant snips in all gwas data 
comb_nss <- Reduce(intersect, list(nss_alt, nss_ggt, nss_apl, nss_nafld, 
                                   nss_t2d, nss_ast, nss_gluc, nss_bmi))

# raw_ggt[SNP==comb_nss[1000]]; raw_alt[SNP==comb_nss[1000]]; raw_nafld[SNP==comb_nss[1000]];
# raw_t2d[SNP==comb_nss[1000]]; mbc_raw_ast[SNP==comb_nss[1000]]; raw_ggt[SNP==comb_nss[1000]];

P_omni_signf_loci$SNP <- paste(P_omni_signf_loci$CHR,P_omni_signf_loci$BP,sep=":")
minP_signf_loci$SNP <- paste(minP_signf_loci$CHR,minP_signf_loci$BP,sep=":")
SZ_signf_loci$SNP <- paste(SZ_signf_loci$CHR,SZ_signf_loci$BP,sep=":")
SZ2_signf_loci$SNP <- paste(SZ2_signf_loci$CHR,SZ2_signf_loci$BP,sep=":")
ET_signf_loci$SNP <- paste(ET_signf_loci$CHR,ET_signf_loci$BP,sep=":")

interest_omni <- Reduce(intersect, list(P_omni_signf_loci$SNP, comb_nss))
interest_minP <- Reduce(intersect, list(minP_signf_loci$SNP, comb_nss))
interest_SZ <- Reduce(intersect, list(SZ_signf_loci$SNP, comb_nss))
interest_SZ2 <- Reduce(intersect, list(SZ2_signf_loci$SNP, comb_nss))
interest_ET <- Reduce(intersect, list(ET_signf_loci$SNP, comb_nss))

interest_all <- Reduce(intersect, list(interest_omni, interest_minP, 
                                       interest_SZ, interest_SZ2, interest_ET))

new_loci <- data.table(OT=length(interest_omni), minP=length(interest_minP),
                       SZ=length(interest_SZ), SZ2=length(interest_SZ2),
                       ET=length(interest_ET))




# rm(nss_ggt, nss_alt, nss_nafld, nss_t2d, nss_ast, nss_gluc, nss_bmi, nss_apl, comb_nss)



z_ref <- copy(ggt_t2d_nalfd_ast_alt_apl_gluc_bmi)
z_ref$SNP <-  paste(z_ref$CHR,z_ref$BP,sep=":")
## Choosing loci of interest to be shown
# Chose from glucagon
sorted_interest_gluc <- merge(raw_gluc[SNP %in% interest_omni, c("SNP", "A1", "A2", "P")][order(P)][1:10], 
                              P_omni_DT[, c("SNP", "P")], by="SNP")
sorted_interest_gluc <- merge(sorted_interest_gluc, 
                              z_ref[, c("SNP", "Z_t2d", "Z_ggt", "Z_nafld", "Z_ast",
                                        "Z_alt", "Z_apl","Z_gluc", "Z_bmi")],
                              by="SNP")
# Chose from t2d
sorted_interest_t2d <- merge(raw_t2d[SNP %in% interest_omni, c("SNP", "A1", "A2", "P")][order(P)][1:10], 
                              P_omni_DT[, c("SNP", "P")], by="SNP")
sorted_interest_t2d <- merge(sorted_interest_t2d, 
                              z_ref[, c("SNP", "Z_t2d", "Z_ggt", "Z_nafld", "Z_ast",
                                        "Z_alt", "Z_apl","Z_gluc", "Z_bmi")],
                              by="SNP")
colnames(sorted_interest_t2d)[4:5] <- c("P_t2d", "P_OT")

# Plot the Z stat contribution
temp_plot <- sorted_interest_t2d[, c("SNP","Z_t2d", "Z_ggt", "Z_nafld", "Z_ast",
                                      "Z_alt", "Z_apl","Z_gluc", "Z_bmi")]
names(temp_plot)[2:9] <- c("t2d", "ggt", "nafld", "ast", "alt", "alp", "gluc", "bmi")
temp_plot <- melt(temp_plot, 
                  id.vars = "SNP", variable.name = "data", values.name = "z_score")[order(-rank(SNP), value)]  
temp_plot[data == "apl", data := "alp"] 
temp_plot[, data := toupper(data)] 

ggplot(temp_plot[SNP=="2:17626707"], aes(x=SNP, y=value, fill=reorder(data, -value))) + 
  geom_bar(position="dodge", stat="identity", alpha=0.75)  + 
  geom_text(aes(label=round(value,2)), position=position_dodge(width=0.9)) +
  theme_minimal() +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.01)) +
  labs(title="Z-score Plot of Interested Locus",
       subtitle="Locus 2:17626707",
       y = "Z Score", x = "Loci",
       fill = "")


names(sorted_interest_gluc)[4:5] <- c("P_gluc", "P_OT")




# Chose from t2d
sorted_interest_t2d <- merge(raw_t2d[SNP %in% interest_omni, c("SNP", "A1", "A2", "P")][order(-P)][1:10], 
                              P_omni_DT[, c("SNP", "P")], by="SNP")
names(sorted_interest_t2d)[4:5] <- c("P_t2d", "P_OT")
sorted_interest_t2d <- merge(sorted_interest_t2d, 
                              z_ref[, c("SNP", "Z_t2d", "Z_ggt", "Z_nafld", "Z_ast",
                                        "Z_alt", "Z_apl","Z_gluc", "Z_bmi")],
                              by="SNP")
# Plot the Z stat contribution
temp_plot <- sorted_interest_t2d[, c("SNP","Z_t2d", "Z_ggt", "Z_nafld", "Z_ast",
                                      "Z_alt", "Z_apl","Z_gluc", "Z_bmi")]
names(temp_plot)[2:9] <- c("t2d", "ggt", "nafld", "ast", "alt", "apl", "gluc", "bmi")
temp_plot <- melt(temp_plot, 
                  id.vars = "SNP", variable.name = "data", values.name = "z_score")[order(-rank(SNP), value)]  

ggplot(temp_plot[SNP=="7:95140031"], aes(x=SNP, y=value, fill=reorder(data, -value))) + 
  geom_bar(position="dodge", stat="identity", alpha=0.75)  + 
  geom_text(aes(label=round(value,2)), position=position_dodge(width=0.9)) +
  theme_minimal() +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.01)) +
  labs(title="Z-score Plot of Interested Loc",
       subtitle="Loci 7:95140031 within 500 KB",
       y = "Z Score", x = "Loci",
       fill = "")

## Get Gene
snp_mart <- biomaRt::useEnsembl(biomart="ENSEMBL_MART_SNP", 
                                host="grch37.ensembl.org", 
                                dataset="hsapiens_snp")
gene_ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_mart <- biomaRt::useMart('ensembl', dataset='hsapiens_gene_ensembl')
get_gene <- function(val) {
  chr <- str_split_1(val, ":")[1]
  bp <- str_split_1(val, ":")[2]
  gene <- getBM(
    attributes = c(
                   "external_gene_name"),
    filters =  c("biotype", "chromosome_name","start","end"), 
    values =  list("protein_coding", chr,bp,bp),
    mart = gene_ensembl,
    verbose = TRUE)
  return(gene[[1]])
}

sorted_interest_gene_gluc <- sapply(sorted_interest_gluc[, SNP], get_gene)
sorted_interest_gene_gluc <- plyr::ldply(sorted_interest_gene_gluc, data.frame)
names(sorted_interest_gene_gluc) <- c("SNP", "Gene")

sorted_interest_gene_t2d<- sapply(sorted_interest_t2d[, SNP], get_gene)
sorted_interest_gene_t2d <- plyr::ldply(sorted_interest_gene_t2d, data.frame)
names(sorted_interest_gene_t2d) <- c("SNP", "Gene")

# Joining gene information 
interest_gene_result <- merge(x = sorted_interest_t2d, 
                              y = sorted_interest_gene_t2d, by = "SNP")
fwrite(interest_gene_result, "interest_gene_result.csv")


# Locus Zoom Plot 
locus_zoom_plot <- function(df, locus, data_name, sel_chr, minBP, maxBP) {

  df_sel_region <- dplyr::filter(df,CHR==sel_chr, between(BP, minBP, maxBP))
  snp_mart <- biomaRt::useEnsembl(biomart="ENSEMBL_MART_SNP", 
                                  host="grch37.ensembl.org", 
                                  dataset="hsapiens_snp")
  gene_ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  gene_mart <- biomaRt::useMart('ensembl', dataset='hsapiens_gene_ensembl')
  # Range Gene 
  out_gene <- getBM(
    attributes = c("ensembl_gene_id", 
                   "external_gene_name",
                   "gene_biotype",
                   "chromosome_name", 
                   "start_position", 
                   "end_position"),
    filters =  c("chromosome_name","start","end"), 
    values =  list(sel_chr,minBP,maxBP),
    mart = gene_ensembl)
  
  plot.range <- c(minBP,maxBP)
  out_gene <- out_gene %>% mutate(gene_biotype_fac = fct_relevel(as.factor(gene_biotype), 
                                                                 "protein_coding"), external_gene_name = fct_reorder2(external_gene_name, 
                                                                                                                      start_position, gene_biotype_fac, .desc = TRUE))
  p2 <- ggplot(data = out_gene) + 
    geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, colour = gene_biotype_fac, group = gene_biotype_fac)) +
    coord_flip() + ylab("") +
    ylim(plot.range) + 
    geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name, colour = gene_biotype_fac), fontface = 2, alpha = I(0.7), hjust = "right", size= 2.5) + 
    labs(title = "", subtitle = paste0("Genes"), 
         # caption = paste0("Data source: ", gene_ensembl@host, " + Data set: ", gene_ensembl@dataset), 
         color = "Gene Biotype") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), 
          strip.text.y = element_text(angle = 0),
          legend.position="bottom", 
          panel.grid.major.y = element_blank()) + 
    expand_limits(y=c(-1, 1)) +
    scale_color_manual(values = c("black", metafolio::gg_color_hue(nlevels(out_gene$gene_biotype_fac)-1)))
  
  p1 <- ggplot(data = df_sel_region) + 
    geom_point(aes(BP, -log10(P)), shape = 1)  +
    geom_point(data=df[SNP==locus], 
               aes(BP, -log10(P)), 
               color='red',
               size=3) +
    labs(title = paste("Locuszoomplot of", data_name), subtitle = paste("Chromosome", sel_chr, "from", format(minBP, big.mark = "'"), "to", format(maxBP, big.mark = "'"), "bp"))
  print(p1)
  
  p1 <- p1 + geom_hline(yintercept=-log10(5e-06), linetype="dashed", color = "blue") + geom_hline(yintercept=-log10(5e-08), linetype="dashed", color = "red")
  
  p1b <- p1 + xlab("") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlim(plot.range)
  show(p1b + p2 + plot_layout(ncol = 1, heights = c(6, 6)))
}

# Locus Zoom Plot 
locus_zoom_plot_wo_gene <- function(df, locus, data_name, sel_chr, minBP, maxBP) {
  
  df_sel_region <- dplyr::filter(df,CHR==sel_chr, between(BP, minBP, maxBP))
  
  plot.range <- c(minBP,maxBP)
  
  p1 <- ggplot(data = df_sel_region) + 
    geom_point(aes(BP, -log10(P)), shape = 1)  +
    geom_point(data=df[SNP==locus], 
               aes(BP, -log10(P)), 
               color='red',
               size=3) +
    labs(title = paste("Locuszoomplot of", data_name), subtitle = paste("Chromosome", sel_chr, "from", format(minBP, big.mark = "'"), "to", format(maxBP, big.mark = "'"), "bp")) + 
    geom_hline(yintercept=-log10(5e-06), linetype="dashed", color = "blue") + geom_hline(yintercept=-log10(5e-08), linetype="dashed", color = "red")
  
  print(p1)
}

## Gene of interest in comment
## window 500000
# 1:1092367 
# 1:7727854 

## Plot gene of interest
locus_zoom_plot(P_omni_DT, "2:7727854", "OT", 2, 7227854, 8227854)
locus_zoom_plot(P_omni_DT, "7:95140031", "OT", 7, 94640031, 95640031)


locus_zoom_plot(P_omni_DT, "2:17626707", "OT", 2, 17126707, 18126707)
locus_zoom_plot_wo_gene(raw_gluc, "2:17626707", "Glucagon", 2, 17126707, 18126707)
locus_zoom_plot_wo_gene(raw_bmi, "2:17626707", "BMI", 2, 17126707, 18126707)


locus_zoom_plot(P_omni_DT, "19:7244884", "OT", 19, 6744884, 7744884)
locus_zoom_plot_wo_gene(raw_gluc, "19:7244884", "T2D", 19, 6744884, 7744884)
locus_zoom_plot_wo_gene(raw_bmi, "19:7244884", "BMI", 19, 6744884, 7744884)
locus_zoom_plot_wo_gene(mbc_raw_ast, "19:7244884", "AST", 19, 6744884, 7744884)

 ### POWER TEST 
z_cor_eig_vects <- eigen(z_cor,symmetric=TRUE)$vectors
delta <- c(6,0,0,0,0,0,0,0)
signal <- delta %*%z_cor_eig_vects
omnibus_test(signal, z_cor)

