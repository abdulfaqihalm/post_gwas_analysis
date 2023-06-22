library(tidyverse)
library(forcats)
library(patchwork)
library(data.table)
library(tictoc)
library(fastman)
library(ggcorrplot)
library(corrplot)
library(matrixStats)
library(mvtnorm)
library(parallel)
library(UpSetR)
library(MTAR)

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


# ===== ALT ===== (alanine transaminase ALT GCST90013405)
raw_alt <-  data.table::fread("data/Pazoki_ALT_hg19.tsv", header = TRUE)
raw_alt <- dplyr::rename(raw_alt, CHR=CHROM, BP=POS,
                         A1=EA, A2=OA, BETA=ES)

# ===== ALP ===== (Alkaline Phosphatat)
raw_apl <-  data.table::fread("data/Pazoki_APL_hg19.tsv", header = TRUE)
raw_apl <- dplyr::rename(raw_apl, CHR=CHROM, BP=POS,
                         A1=EA, A2=OA, BETA=ES)


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
raw_t2d <- gwas_sanity(raw_t2d, "t2d") 
raw_alt <- gwas_sanity(raw_alt, "alt")
raw_apl <- gwas_sanity(raw_apl, "apl")


## CREATING PLOTS ## 
## Number of Significant (5e-08) SNPs per traits ##
total_signf_snps <- data.table::data.table(gwas_data=c("T2D",
                                                       "ALT",
                                                       "ALP"
), 
snp_count=c(
  raw_t2d[P<5e-08, length(BP)],
  raw_alt[P<5e-08, length(BP)],
  raw_apl[P<5e-08, length(BP)]
)
)


ggplot(total_signf_snps, aes(x=reorder(gwas_data, -snp_count), y=snp_count)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  geom_text(aes(label=snp_count), position=position_dodge(width=0.9), vjust=-0.25) + 
  theme_minimal() +
  labs(title="Significant SNPs for Each GWAS Summary Statistics",
       y = "SNPs Count", x = "GWAS Data")

## CREATING PLOTS ## 
## Number of SNPs per traits ##
total_snps <- data.table::data.table(gwas_data=c("T2D",
                                                 "ALT",
                                                 "ALP"
), 
snp_count=c(
  raw_t2d[, length(BP)],
  raw_alt[, length(BP)],
  raw_apl[, length(BP)]
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
total_signf_snps <- data.table::data.table(gwas_data=c("T2D",
                                                       "ALT",
                                                       "APL"
), 
snp_count=c(
  get_top_snps(raw_t2d)[, length(BP)],
  get_top_snps(raw_alt)[, length(BP)],
  get_top_snps(raw_apl)[, length(BP)]
)
)

ggplot(total_signf_snps, aes(x=reorder(gwas_data, -snp_count), y=snp_count)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  geom_text(aes(label=snp_count), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_minimal() +
  labs(title="Significant Loci for Each GWAS Data",
       subtitle="Loci within 500KB distance considered as one locus",
       y = "Loci Count", x = "GWAS Data")

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

cat(min(raw_t2d$P) ,"\n",
    min(raw_alt$P) ,"\n",
    min(raw_apl$P) ,"\n")

raw_t2d[P == 0, P := SMALL]
raw_alt[P == 0, P := SMALL]
raw_apl[P == 0, P := SMALL]

print_manhattan_plot(raw_t2d, "Manhattan Plot for GGT GWAS Data", upper=-log10(min(raw_t2d$P)), cex.text=0.9, WinMb=66)
print_manhattan_plot(raw_alt, "Manhattan Plot for ALT GWAS Data", upper=-log10(min(raw_alt$P)), cex.text=0.9, WinMb=66)
print_manhattan_plot(raw_apl, "Manhattan Plot for ALP GWAS Data", upper=-log10(min(raw_apl$P)), cex.text=0.9, WinMb=66)

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

# Create the joined sum stat data
tic()
alt_t2d <- create_z(raw_alt[, .(CHR, BP, A1, A2, Z_alt)], raw_t2d[, .(CHR, BP, A1, A2, Z_t2d)] )
setkey(alt_t2d, CHR, BP, A1, A2)

alt_t2d_apl <- create_z(alt_t2d, raw_apl[, .(CHR, BP, A1, A2, Z_apl)])
setkey(alt_t2d_apl, CHR, BP, A1, A2)
rm(alt_t2d)
toc()


## EMPIRICAL CORRELATION MATRIX ##                               
tic()
sel_range <- c(5 : length(alt_t2d_apl))
z_mat <- data.matrix(alt_t2d_apl[, ..sel_range])
z_mat_cut_05 <- data.matrix((alt_t2d_apl[abs(Z_t2d)<0.5 & abs(Z_alt)<0.5 & abs(Z_apl)<0.5])[, ..sel_range])
z_mat_cut_1 <- data.matrix((alt_t2d_apl[abs(Z_t2d)<1 & abs(Z_alt)<1 & abs(Z_apl)<1])[, ..sel_range])
z_mat_cut_2 <- data.matrix((alt_t2d_apl[abs(Z_t2d)<2 & abs(Z_alt)<2 & abs(Z_apl)<2])[, ..sel_range])

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

# for (row in c(1:17001)) {
#   result <- c()
#   if (any(is.na(z_mat[row,]))) {
#     # cat(row, "NA")
#     cat("NA ")
#   } else{
#     cat("Z vector", z_mat[row,])
#     cat(MTAR::emats(z_mat[row,], z_cor)$p.value)
#   }
# }


non_NA_z_mat <-as.matrix(na.omit(data.table::data.table(z_mat)))
non_NA_DT <- alt_t2d_apl[(!is.na(Z_t2d) 
                          & !is.na(Z_alt) 
                          & !is.na(Z_apl))]

tic()
# Return the P-value
# All rows are kept to get the SNPs id
P_omni <- omnibus_test(z_mat, z_cor)
 P_sumZ <- sumZ_test(z_mat, z_cor)
P_sumZ2 <- sumZ2_test(z_mat, z_cor)
P_pc <- pc_test(z_mat, z_cor)
toc()

# Get the SNP ref for multi trait summ stat for further analysis
get_SNP_DT_ref <- function(summ_P_stat, SNP_ref){
  DT <- data.table(SNP=paste(SNP_ref$CHR, SNP_ref$BP, sep = ":"), CHR=SNP_ref$CHR, BP=SNP_ref$BP, P=summ_P_stat)
  DT |>
    setkey(SNP)
  return(DT)
}

P_omni_DT <- get_SNP_DT_ref(P_omni, alt_t2d_apl)
P_omni_DT[is.infinite(P) | P==0, P := SMALL]

P_sumZ_DT <- get_SNP_DT_ref(P_sumZ, alt_t2d_apl)
P_sumZ_DT[is.infinite(P) | P==0, P := SMALL]

P_sumZ2_DT <- get_SNP_DT_ref(P_sumZ2, alt_t2d_apl)
P_sumZ2_DT[is.infinite(P) | P==0, P := SMALL]

P_pc_DT <- get_SNP_DT_ref(P_pc, alt_t2d_apl)
P_pc_DT[is.infinite(P) | P==0, P := SMALL]

cat(min(P_omni_DT[!is.na(P), P]) ,"\n",
    min(P_sumZ_DT[!is.na(P), P]) ,"\n",
    min(P_sumZ2_DT[!is.na(P), P]) ,"\n",
    min(P_pc_DT[!is.na(P), P]) ,"\n")

rm(P_omni)
rm(P_sumZ)
rm(P_sumZ2)
rm(P_pc)

## PLOTTING THE MANHATTAN PLOT ##
print_manhattan_plot(P_omni_DT, "Manhattan Plot for Joint GWAS Data (OT)", upper=-log10(min(P_omni_DT[!is.na(P)]$P)), cex.text=0.9, WinMb=66)
print_manhattan_plot(P_sumZ_DT, "Manhattan Plot for Joint GWAS Data (SZ)", upper=-log10(min(P_sumZ_DT[!is.na(P)]$P)), cex.text=0.9, WinMb=66)
print_manhattan_plot(P_sumZ2_DT, "Manhattan Plot for Joint GWAS Data (SZ2)", upper=-log10(min(P_sumZ2_DT[!is.na(P)]$P)), cex.text=0.9, WinMb=66)
print_manhattan_plot(P_pc_DT, "Manhattan Plot for Joint GWAS Data (ET)", upper=-log10(min(P_pc_DT[!is.na(P)]$P)), cex.text=0.9, WinMb=66)


# Create the significant loci DT (500KB)  
P_omni_signf_loci <- get_top_snps(P_omni_DT[!is.na(P)])
SZ_signf_loci <- get_top_snps(P_sumZ_DT[!is.na(P)])
SZ2_signf_loci <- get_top_snps(P_sumZ2_DT[!is.na(P)])
ET_signf_loci <- get_top_snps(P_pc_DT[!is.na(P)])

# Plotting the barchart of loci distribution 
join_signf_loci <- data.table::data.table(gwas_data=c("OT","SZ","SZ2","ET"), 
snp_count=c(
  P_omni_signf_loci[, length(BP)],
  SZ_signf_loci[, length(BP)],
  SZ2_signf_loci[, length(BP)],
  ET_signf_loci[, length(BP)]
)
)

ggplot(join_signf_loci, aes(x=reorder(gwas_data, -snp_count), y=snp_count)) + 
  geom_bar(stat="identity", fill="steelblue") + 
  geom_text(aes(label=snp_count), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_minimal() +
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
nss_alt <- raw_alt[P>=5e-08, SNP]
nss_apl <- raw_apl[P>=5e-08, SNP]

# non significant snips in all gwas data 
comb_nss <- Reduce(intersect, list(nss_alt, nss_apl))

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