#!/usr/bin/env Rscript
#######################################################################################
############ This sciprt is a easy pipeline
############ To perform two sample mendelian randomization using TwoSampleMR R package
############ Make sure all required R packages are downloaded.
#######################################################################################
library(TwoSampleMR)
library(ieugwasr)
library(tidyverse)
library(MRPRESSO)
library(data.table)
library(glue)
library(argparser)


###### Get arguments from command line #############
p <- arg_parser("Perform mendelian randomization analysis using TwoSampleMR")
p <- add_argument(p, "--exp", help = "exposure gwas summary")
p <- add_argument(p, "--otc", help = "outcome gwas summary")
p <- add_argument(p, "--nexp", help = "sample size of exposure")
p <- add_argument(p, "--notc", help = "sample size of outcome")
p <- add_argument(p, "--prefix", help = "output prefix for results", default = "TMR_results")
p <- add_argument(p, "--expname", help = "exposure name", default = "exposure")
p <- add_argument(p, "--otcname", help = "outcome name", default = "outcome")
p <- add_argument(p, "--proxy", help = "use proxy snps", flag = TRUE)
p <- add_argument(p, "--online", help = "search proxy online", flag =TRUE)
p <- add_argument(p, "--pop", help = "specify LD population <default = EUR>", default = "EUR")
p <- add_argument(p, "--rsq", help = "specify LD r-square <default = 0.8>", default = 0.8)
p <- add_argument(p, "--bidirect", help = "option to perform bidirectional MR", flag = TRUE)
p <- add_argument(p, "--bfile", help = "plink format file to use LD reference panel", default = "")
p <- add_argument(p, "--plink", help = "path to plink", default = "")
p <- add_argument(p, "--piv", help = "pvalue threshold for IV", default = 5e-8)
p <- add_argument(p, "--rpiv", help = "pvalue threshold for reverse IV", default = 5e-8)
p <- add_argument(p, "--esnp", help = "column name of rsid in exposure gwas", default = "variant_id")
p <- add_argument(p, "--osnp", help = "column name of rsid in outcome gwas", default = "variant_id")
p <- add_argument(p, "--ebeta", help = "column name of beta in exposure gwas", default = "beta")
p <- add_argument(p, "--obeta", help = "column name of beta in outcome gwas", default = "beta")
p <- add_argument(p, "--ese", help = "column name of se in exposure gwas", default = "standard_error")
p <- add_argument(p, "--ose", help = "column name of se in outcome gwas", default = "standard_error")
p <- add_argument(p, "--echr", help = "column name of chr in exposure gwas", default = "chromosome")
p <- add_argument(p, "--ochr", help = "column name of chr in outcome gwas", default = "chromosome")
p <- add_argument(p, "--epos", help = "column name of snp position in exposure gwas", default = "base_pair_location")
p <- add_argument(p, "--opos", help = "column name of snp position in outcome gwas", default = "base_pair_location")
p <- add_argument(p, "--eea", help = "column name of effect allele in exposure gwas", default = "effect_allele")
p <- add_argument(p, "--oea", help = "column name of effect allele in outcome gwas", default = "effect_allele")
p <- add_argument(p, "--eoa", help = "column name of other allele in exposure gwas", default = "other_allele")
p <- add_argument(p, "--ooa", help = "column name of other allele in outcome gwas", default = "other_allele")
p <- add_argument(p, "--eeaf", help ="column name of effect allele frequency in exposure gwas", default = "effect_allele_frequency")
p <- add_argument(p, "--oeaf", help ="column name of effect allele frequency in outcome gwas", default = "effect_allele_frequency")
p <- add_argument(p, "--epval", help = "column name of pvalue in exposure gwas", default = "p_value")
p <- add_argument(p, "--opval", help = "column name of pvalue in outcome gwas", default = "p_value")
argv <- parse_args(p)
out.tsv.path <- paste(argv$prefix, ".tsv", sep = "")
out.tsv.nullpath <- paste(argv$prefix, ".null.tsv", sep = "")
out.reverse.tsv.path <- paste(argv$prefix, "_reverse.tsv", sep = "")
out.reverse.tsv.nullpath <- paste(argv$prefix, "_reverse.null.tsv", sep = "")
out.df <- data.frame()
out.reverse.df <- data.frame()
if (!argv$bidirect & (file.exists(out.tsv.path) | file.exists(out.tsv.nullpath))) {
    cat("File exits! Skip.\n")
    q(save = "no")
} else if (argv$bidirect & (file.exists(out.tsv.path) | file.exists(out.tsv.nullpath)) & (file.exists(out.reverse.tsv.path) | file.exists(out.reverse.tsv.nullpath))) {
    cat("File exits! Skip.\n")
    q(save = "no")
}
######## function for get proxy #########
get_proxy <- function(exp_dat, otc_dat, r2 = 0.8, pop = "EUR", online = TRUE, bin_plink, ref_panel, outfile = "proxy") {
            exp.snp = exp_dat$SNP
            otc.snp = otc_dat$SNP
            otc_dat = otc_dat %>% mutate(
                target_snp.outcome = NA,
                proxy_snp.outcome = NA,
                target_a1.outcome = NA,
                target_a2.outcome = NA,
                proxy_a1.outcome = NA,
                proxy_a2.outcome = NA,
                proxy.outcome = NA
                )
    miss.expsnps = setdiff(exp.snp, intersect(exp.snp, otc.snp))
    remain.otcsnps = setdiff(otc.snp, intersect(exp.snp, otc.snp))
    for (miss.snp in miss.expsnps) {
        if (online) {
            library(LDlinkR)
            proxies.df = LDproxy(snp = miss.snp, pop = pop, token = "6fb632e022ef")
            top_proxies.df = proxies.df %>% filter(R2 >= r2, Distance > 0, RS_Number %in% remain.otcsnps) %>% arrange(-R2) %>% head(1)
        } else {
            outfilepath = shQuote(outfile)
            system(glue("{bin_plink} -bfile {ref_panel} --r2 --ld-snp {miss.snp} --ld-window 1000000 --ld-window-kb 1000 --ld-window-r2 {r2} --out {outfilepath}_{miss.snp}_proxy"))
            proxies.df = read_table(glue("{outfile}_{miss.snp}_proxy.ld"))
            top_proxies.df = proxies.df %>% filter(SNP_A != SNP_B, SNP_B %in% remain.otcsnps) %>% arrange(-R2) %>% head(1) %>% rename(RS_Number = SNP_B)
        }
        if (nrow(top_proxies.df) > 0) {
            cat("Found proxies of ", miss.snp,"\n", sep = "")
            rsid = top_proxies.df$RS_Number
            effect_allele.exposure = (exp_dat %>% filter(SNP == miss.snp) %>% head(1))$effect_allele.exposure
            other_allele.exposure = (exp_dat %>% filter(SNP == miss.snp) %>% head(1))$other_allele.exposure
            if (online) {
                proxy.allele1 = (str_replace_all((top_proxies.df %>% head(1))$Alleles, "[()]", "") %>% str_split("/"))[[1]][1]
                proxy.allele2 = (str_replace_all((top_proxies.df %>% head(1))$Alleles, "[()]", "") %>% str_split("/"))[[1]][2]
                target.allele1 = (((top_proxies.df %>% head(1))$Correlated_Alleles  %>% str_split(","))[[1]][1] %>% str_split("="))[[1]][1]
                target.allele2 = (((top_proxies.df %>% head(1))$Correlated_Alleles  %>% str_split(","))[[1]][1] %>% str_split("="))[[1]][2]
            } else {
                allele_info <- system(glue('{bin_plink} -bfile {ref_panel} --ld {miss.snp} {rsid} | grep \'In phase\' | awk \'{{print substr($NF,1,1)"\\n"substr($NF,4,1)"\\n"substr($NF,2,1)"\\n"substr($NF,5,1)}}\''), intern = TRUE)
                target.allele1 = allele_info[1]
                target.allele2 = allele_info[2]
                proxy.allele1 = allele_info[3]
                proxy.allele2 = allele_info[4]
            }
            if (!identical(sort(c(target.allele1, target.allele2)), sort(c(effect_allele.exposure, other_allele.exposure)))){
                cat("Alleles are not consistent for ", miss.snp, "\n", sep = "")
                next
            }
            otc_dat[otc_dat$SNP == rsid,]$target_snp.outcome = miss.snp
            otc_dat[otc_dat$SNP == rsid,]$proxy_snp.outcome = rsid
            otc_dat[otc_dat$SNP == rsid,]$target_a1.outcome = ifelse(otc_dat[otc_dat$SNP == rsid,]$effect_allele.outcome == proxy.allele1, target.allele1, target.allele2)
            otc_dat[otc_dat$SNP == rsid,]$target_a2.outcome = ifelse(otc_dat[otc_dat$SNP == rsid,]$other_allele.outcome == proxy.allele2, target.allele2, target.allele1)
            otc_dat[otc_dat$SNP == rsid,]$proxy_a1.outcome = ifelse(otc_dat[otc_dat$SNP == rsid,]$effect_allele.outcome == proxy.allele1, proxy.allele1, proxy.allele2)
            otc_dat[otc_dat$SNP == rsid,]$proxy_a2.outcome = ifelse(otc_dat[otc_dat$SNP == rsid,]$other_allele.outcome == proxy.allele2, proxy.allele2, proxy.allele1)
            otc_dat[otc_dat$SNP == rsid,]$proxy.outcome = TRUE
            otc_dat[otc_dat$SNP == rsid,]$effect_allele.outcome = ifelse(otc_dat[otc_dat$SNP == rsid,]$effect_allele.outcome == proxy.allele1, target.allele1, target.allele2)
            otc_dat[otc_dat$SNP == rsid,]$other_allele.outcome = ifelse(otc_dat[otc_dat$SNP == rsid,]$other_allele.outcome == proxy.allele2, target.allele2, target.allele1)
            otc_dat[otc_dat$SNP == rsid,]$SNP = miss.snp
            file.remove(glue("{outfilepath}_{miss.snp}_proxy.ld"))
            file.remove(glue("{outfilepath}_{miss.snp}_proxy.log"))
            # system(glue("rm {outfilepath}_{miss.snp}_proxy.ld {outfilepath}_{miss.snp}_proxy.log"))
        } else {
            cat("Not found proxies of ", miss.snp, "\n", sep = "")
        }
    }
    return(otc_dat %>% filter(SNP %in% exp.snp))
}
########## perform MR #####################
performMR <- function(exp_dat, otc_dat, outprefix, expname = "exposure", otcname = "outcome", piv = 5e-8) {
    if (argv$proxy) {
        cat("Searching proxy snps...\n")
        otc_dat <- get_proxy(
            exp_dat,
            otc_dat,
            r2 = as.numeric(argv$rsq),
            online = argv$online,
            pop = argv$pop,
            bin_plink = argv$plink,
            ref_panel = argv$bfile,
            outfile = outprefix
            )
    }
    otc_dat <- otc_dat %>% mutate(outcome = argv$otcname, id.outcome = argv$otcname, originalname.outcome = argv$otcname)
    cat("Harmonising data......\n")
    exp_otc_dat <- harmonise_data(exp_dat, otc_dat)
    if (sum(exp_otc_dat$mr_keep) >= 4) {
            exp_otc_mrPresso <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure",
                                            OUTLIERtest = TRUE, DISTORTIONtest = FALSE, data = exp_otc_dat %>% filter(mr_keep == TRUE), SignifThreshold = 0.05)
            outlier.df <- exp_otc_mrPresso$`MR-PRESSO results`$`Outlier Test`
        if (!is.null(outlier.df) && nrow(outlier.df) > 0) {
            write.table(outlier.df, paste(argv$prefix, "_outlierPresso.tsv", sep = ""))
            row.names(outlier.df) <- (exp_otc_dat %>% filter(mr_keep == TRUE))$SNP
            exp_otc_dat <- exp_otc_dat %>% filter(SNP %in% row.names(outlier.df %>% filter(Pvalue > 0.05)))
        } else {
            cat("No outliers found or detected.\n")
        }
    }
    if (!is.null(exp_otc_dat) && nrow(exp_otc_dat) > 0 && sum(exp_otc_dat$mr_keep) >= 1) {
        exp_otc_dat <- exp_otc_dat %>% mutate(samplesize.exposure = as.numeric(argv$nexp), samplesize.outcome = as.numeric(argv$notc))
        cat("Performing MR...\n")
        exp_otc_res <- mr(exp_otc_dat, method_list = c("mr_wald_ratio", "mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_median", "mr_weighted_mode"))
        exp_otc_res <- generate_odds_ratios(exp_otc_res)
        exp_otc_res <- exp_otc_res %>% mutate(beta_95_CI = paste(round(b, 3), " [", round(lo_ci, 3), "-", round(up_ci, 3), "]", sep = ""),
                                            OR_95_CI = paste(round(or, 3), " [", round(or_lci95, 3), "-", round(or_uci95, 3), "]", sep = ""))
        exp_otc_res_single <- mr_singlesnp(exp_otc_dat)
        p1 <- mr_scatter_plot(exp_otc_res, exp_otc_dat)
        p2 <- mr_forest_plot(exp_otc_res_single)
        p3 <- mr_funnel_plot(exp_otc_res_single)
        res_loo <- mr_leaveoneout(exp_otc_dat)
        p4 <- mr_leaveoneout_plot(res_loo)
        ggsave(p1[[1]], file = glue("{outprefix}_scatterplot.pdf"), width = 7, height = 7)
        ggsave(p2[[1]], file = glue("{outprefix}_forestplot.pdf"), width = 7, height = 7)
        ggsave(p3[[1]], file = glue("{outprefix}_funnelplot.pdf"), width = 7, height = 7)
        ggsave(p4[[1]], file = glue("{outprefix}_leaveoneoutplot.pdf"), width = 7, height = 7)
        write_tsv(exp_otc_dat, glue("{outprefix}_harmonise.tsv"))
        write_tsv(exp_otc_res_single, glue("{outprefix}_single.tsv"))
        write_tsv(res_loo, glue("{outprefix}_leaveoneout.tsv"))
        out.df <- data.frame(
            exposure = expname,
            outcome = otcname,
            ivw_or_wald.b = (exp_otc_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")))$b,
            nsnp = (exp_otc_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")))$nsnp,
            ivw_or_wald.se = (exp_otc_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")))$se,
            ivw_or_wald.b.95ci = (exp_otc_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")))$beta_95_CI,
            ivw_or_wald.or.95ci = (exp_otc_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")))$OR_95_CI,
            ivw_or_wald.pval = (exp_otc_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")))$pval
            )
        cat("Finished MR\n")
        ####### sensitive test ##############
        cat("Performing sensitive MR...\n")
        # 1. Heterogeneity test
        if (sum(exp_otc_dat$mr_keep) >= 2) {
            exp_otc_mrHet <- mr_heterogeneity(exp_otc_dat, method_list = "mr_ivw")
            q.het <- exp_otc_mrHet$Q_pval
            I2 <- (exp_otc_mrHet  %>% mutate(I2 = (Q-Q_df)/Q*100) %>% filter(method == "Inverse variance weighted"))$I2
            out.df$q.het <- q.het
            out.df$I2 <- max(0, I2)
        } else {
            out.df$q.het <- NA
            out.df$I2 <- NA
        }
        # 2. Horizontal pleiotropy test
        if (sum(exp_otc_dat$mr_keep) >= 3) {
            out.df$eggReg.b <- (exp_otc_res %>% filter(method == "MR Egger"))$b
            out.df$eggReg.se <- (exp_otc_res %>% filter(method == "MR Egger"))$se
            out.df$eggReg.b.95ci <- (exp_otc_res %>% filter(method == "MR Egger"))$beta_95_CI
            out.df$eggReg.or.95ci <- (exp_otc_res %>% filter(method == "MR Egger"))$OR_95_CI
            out.df$eggReg.pval <- (exp_otc_res %>% filter(method == "MR Egger"))$pval
            out.df$weightMed.b <- (exp_otc_res %>% filter(method == "Weighted median"))$b
            out.df$weightMed.se <- (exp_otc_res %>% filter(method == "Weighted median"))$se
            out.df$weightMed.b.95ci <- (exp_otc_res %>% filter(method == "Weighted median"))$beta_95_CI
            out.df$weightMed.or.95ci <- (exp_otc_res %>% filter(method == "Weighted median"))$OR_95_CI
            out.df$weightMed.pval <- (exp_otc_res %>% filter(method == "Weighted median"))$pval
            out.df$weightMode.b <- (exp_otc_res %>% filter(method == "Weighted mode"))$b
            out.df$weightMode.se <- (exp_otc_res %>% filter(method == "Weighted mode"))$se
            out.df$weightMode.b.95ci <- (exp_otc_res %>% filter(method == "Weighted mode"))$beta_95_CI
            out.df$weightMode.or.95ci <- (exp_otc_res %>% filter(method == "Weighted mode"))$OR_95_CI
            out.df$weightMode.pval <- (exp_otc_res %>% filter(method == "Weighted mode"))$pval
            out.df$simpleM.b <- (exp_otc_res %>% filter(method == "Simple median"))$b
            out.df$simpleM.se <- (exp_otc_res %>% filter(method == "Simple median"))$se
            out.df$simpleM.b.95ci <- (exp_otc_res %>% filter(method == "Simple median"))$beta_95_CI
            out.df$simpleM.or.95ci <- (exp_otc_res %>% filter(method == "Simple median"))$OR_95_CI
            out.df$simpleM.pval <- (exp_otc_res %>% filter(method == "Simple median"))$pval
            exp_otc_mrEggerIntercept <- mr_pleiotropy_test(exp_otc_dat)
            out.df$p.eggerIntercept <- exp_otc_mrEggerIntercept$pval
            out.df$eggerIntercept <- exp_otc_mrEggerIntercept$egger_intercept
            out.df$se.eggerIntercept <- exp_otc_mrEggerIntercept$se
        } else {
            out.df$p.eggerIntercept <- NA
            out.df$eggerIntercept <- NA
            out.df$se.eggerIntercept <- NA
            out.df$eggReg.b <- NA
            out.df$eggReg.se <- NA
            out.df$eggReg.b.95ci <- NA
            out.df$eggReg.or.95ci <- NA
            out.df$eggReg.pval <- NA
            out.df$weightMed.b <- NA
            out.df$weightMed.se <- NA
            out.df$weightMed.b.95ci <- NA
            out.df$weightMed.or.95ci <- NA
            out.df$weightMed.pval <- NA
            out.df$weightMode.b <- NA
            out.df$weightMode.se <- NA
            out.df$weightMode.b.95ci <- NA
            out.df$weightMode.or.95ci <- NA
            out.df$weightMode.pval <- NA
            out.df$simpleM.b <- NA
            out.df$simpleM.se <- NA
            out.df$simpleM.b.95ci <- NA
            out.df$simpleM.or.95ci <- NA
            out.df$simpleM.pval <- NA
        }
        # 3. MR-presso, global pleiotropy test and outlier test
        out.df$presso.b <- NA
        out.df$presso.se <- NA
        out.df$presso.b.95ci <- NA
        out.df$presso.pval <- NA
        out.df$presso.or <- NA
        out.df$presso.or.95ci <- NA
        if (sum(exp_otc_dat$mr_keep) >= 4) {
            out.df$p.globalPleiotropy <- exp_otc_mrPresso$`MR-PRESSO results`$`Global Test`$Pvalue
            out.df$RSSobs <- exp_otc_mrPresso$`MR-PRESSO results`$`Global Test`$RSSobs
            index <- ifelse(is.na(exp_otc_mrPresso$`Main MR results`[2, 3]), 1, 2)
            out.df$presso.b <- exp_otc_mrPresso$`Main MR results`$`Causal Estimate`[index]
            out.df$presso.se <- exp_otc_mrPresso$`Main MR results`$Sd[index]
            out.df$presso.pval <- exp_otc_mrPresso$`Main MR results`$`P-value`[index]
            out.df$presso.or <- exp(exp_otc_mrPresso$`Main MR results`$`Causal Estimate`[index])
            lci <- exp_otc_mrPresso$`Main MR results`$`Causal Estimate`[index] - 1.96 * exp_otc_mrPresso$`Main MR results`$Sd[index]
            uci <- exp_otc_mrPresso$`Main MR results`$`Causal Estimate`[index] + 1.96 * exp_otc_mrPresso$`Main MR results`$Sd[index]
            out.df$presso.b.95ci <- paste("[", lci, "-", uci, "]", sep = "")
            or.lci <- exp(lci)
            or.uci <- exp(uci)
            out.df$presso.or.95ci <- paste("[", or.lci, "-", or.uci, "]", sep = "")
            # 4. Same direction for all method
            same_direction <- all(c(exp_otc_res$b, out.df$presso.b) > 0) | all(c(exp_otc_res$b, out.df$presso.b) < 0)
        } else {
            out.df$p.globalPleiotropy <- NA
            out.df$RSSobs <- NA
            same_direction <- all(exp_otc_res$b > 0) | all(exp_otc_res$b < 0)
        }
        out.df$same_direction <- same_direction


        # 5. Reverse causation test
        exp_otc_Steiger <- directionality_test(exp_otc_dat)
        (exp_otc_Steiger$snp_r2.exposure < exp_otc_Steiger$snp_r2.outcome) & exp_otc_Steiger$steiger_pval < 0.05
        out.df$snp_r2.exposure <- exp_otc_Steiger$snp_r2.exposure
        out.df$snp_r2.outcome <- exp_otc_Steiger$snp_r2.outcome
        out.df$steiger_pval <- exp_otc_Steiger$steiger_pval
        out.df$reverse_caution <- (exp_otc_Steiger$snp_r2.exposure < exp_otc_Steiger$snp_r2.outcome) & exp_otc_Steiger$steiger_pval < 0.05
        if (out.df$nsnp == 1) {
            out.df$pass_alltests <- ifelse((out.df$ivw_or_wald.pval < 0.05) & (!out.df$reverse_caution), TRUE, FALSE)
        } else if (out.df$nsnp == 2) {
            out.df$pass_alltests <- ifelse((out.df$ivw_or_wald.pval < 0.05) & !(out.df$I2 > 50 & out.df$q.het < 0.05) & (!out.df$reverse_caution), TRUE, FALSE)
        } else if (out.df$nsnp == 3) {
            out.df$pass_alltests <- ifelse((out.df$ivw_or_wald.pval < 0.05) & !(out.df$I2 > 50 & out.df$q.het < 0.05) & (out.df$p.eggerIntercept > 0.05) &(out.df$same_direction) & (!out.df$reverse_caution), TRUE, FALSE)
        } else {
            out.df$pass_alltests <- ifelse((out.df$ivw_or_wald.pval < 0.05) & !(out.df$I2 > 50 & out.df$q.het < 0.05) & (out.df$p.eggerIntercept > 0.05) &(out.df$same_direction) & (!out.df$reverse_caution) & (out.df$p.globalPleiotropy > 0.05), TRUE, FALSE)
        }
        cat("Finished sensitive test\n")
    } else {
        out.df = data.frame()
    }
   return(out.df)
}


########## Main ###################

# read data
exp_gwas <- fread(argv$exp)
otc_gwas <- fread(argv$otc)
exp_gwas.filter <- exp_gwas %>% filter(get(argv$epval) < argv$piv & str_detect(get(argv$esnp), "^rs\\d+"))

# format data
if (!is.null(exp_gwas.filter) && nrow(exp_gwas.filter) > 0) {
    exp_dat <- format_data(
        exp_gwas.filter,
        snp_col = argv$esnp,
        beta_col = argv$ebeta,
        se_col = argv$ese,
        chr_col = argv$echr,
        pos_col = argv$epos,
        effect_allele_col = argv$eea,
        other_allele_col = argv$eoa,
        eaf_col = argv$eeaf,
        pval_col = argv$epval,
        type = "exposure")
    otc_dat <- format_data(
        otc_gwas,
        snp_col = argv$osnp,
        beta_col = argv$obeta,
        se_col = argv$ose,
        chr_col = argv$ochr,
        pos_col = argv$opos,
        effect_allele_col = argv$oea,
        other_allele_col = argv$ooa,
        eaf_col = argv$oeaf,
        pval_col = argv$opval,
        type = "outcome"
    )
    # exp_dat <- exp_dat %>% filter(pval.exposure < argv$piv)
    rsid_pval_exp <- exp_dat %>% select(SNP, pval.exposure) %>% rename(rsid = SNP, pval = pval.exposure) %>% filter(!is.na(rsid) & nzchar(rsid))
    if (nrow(rsid_pval_exp) > 0){
        cat("Clumping...\n")
        rsid_pval_exp.clump <- tryCatch({ieugwasr::ld_clump_local(rsid_pval_exp,
                                            clump_kb = 10000,
                                            clump_r2 = 0.001,
                                            clump_p = 1,
                                            plink_bin = argv$plink,
                                            bfile = argv$bfile)},
                                        error = function(e){
                                                print(e)
                                                return(data.frame())
                                                })
        if (nrow(rsid_pval_exp.clump) > 0){
                exp_dat <- exp_dat %>% filter(SNP %in% rsid_pval_exp.clump$rsid)
                cat("Finished clump\n")
                cat("Perform MR...\n")
                out.df <- performMR(exp_dat, otc_dat, argv$prefix, expname = argv$expname, otcname = argv$otcname, piv = as.numeric(argv$piv))
        } else {
                out.df <- data.frame()
        }
    } else{
        writeLines("Not found IVs.",out.tsv.nullpath)
    }

} else {
    writeLines("nrow(exp_gwas.filter) = 0", out.tsv.nullpath)
}

# reverse MR
# if (nrow(out.df) > 0) {
#     out.df$reverse_MR_sig <- NA
#     out.df$reverse_MR_allpass <- NA
# } else {
#     writeLines("out.df=0. Not found IVs.",out.tsv.nullpath)
# }

otc_gwas.filter <- otc_gwas %>% filter(get(argv$opval) < argv$rpiv, str_detect(get(argv$osnp), "^rs\\d+"))
if (argv$bidirect) {
    if (!is.null(otc_gwas.filter) && nrow(otc_gwas.filter) > 0) {
        exp_rev_dat <- format_data(
        otc_gwas.filter,
        snp_col = argv$osnp,
        beta_col = argv$obeta,
        se_col = argv$ose,
        chr_col = argv$ochr,
        pos_col = argv$opos,
        effect_allele_col = argv$oea,
        other_allele_col = argv$ooa,
        eaf_col = argv$oeaf,
        pval_col = argv$opval,
        type = "exposure")
        otc_rev_dat <- format_data(
            exp_gwas,
            snp_col = argv$esnp,
            beta_col = argv$ebeta,
            se_col = argv$ese,
            chr_col = argv$echr,
            pos_col = argv$epos,
            effect_allele_col = argv$eea,
            other_allele_col = argv$eoa,
            eaf_col = argv$eeaf,
            pval_col = argv$epval,
            type = "outcome")
        # exp_rev_dat <- exp_rev_dat %>% filter(pval.exposure < argv$rpiv)
        rev_rsid_pval_exp <- exp_rev_dat %>% select(SNP, pval.exposure) %>% rename(rsid = SNP, pval = pval.exposure) %>% filter(!is.na(rsid) & nzchar(rsid))
        if (nrow(rev_rsid_pval_exp) > 0) {
            cat("Clumping...\n")
            rev_rsid_pval_exp.clump <- tryCatch({ieugwasr::ld_clump_local(rev_rsid_pval_exp,
                                                    clump_kb = 10000,
                                                    clump_r2 = 0.001,
                                                    clump_p = 1,
                                                    plink_bin = argv$plink,
                                                    bfile = argv$bfile)},
                                                error = function(e) {
                                                        print("No sig clump")
                                                        return(data.frame())})
            if (nrow(rev_rsid_pval_exp.clump) > 0){
                exp_rev_dat <- exp_rev_dat %>% filter(SNP %in% rev_rsid_pval_exp.clump$rsid)
                cat("Finished clump\n")
                cat("Perform reverse MR...\n")
                out.reverse.df <- performMR(exp_rev_dat, otc_rev_dat, paste(argv$prefix, "_reverse", sep = ""), expname = argv$otcname, otcname = argv$expname, piv = as.numeric(argv$rpiv))
                } else {
                        out.reverse.df = data.frame()
                }
        } else {
            writeLines("Not found IVs.",out.reverse.tsv.nullpath)
        }


   if (nrow(out.reverse.df) > 0 && "ivw_or_wald.pval" %in% names(out.reverse.df)) {
        out.reverse.df$reverse_MR_sig <- NA
        out.reverse.df$reverse_MR_allpass <- NA
        if (nrow(out.df) > 0) {
        out.reverse.df$reverse_MR_sig <- ifelse(out.df$ivw_or_wald.pval < 0.05, TRUE, FALSE)
        out.reverse.df$reverse_MR_allpass <- ifelse(out.df$pass_alltests, TRUE, FALSE)
        out.df$reverse_MR_sig <- ifelse(out.reverse.df$ivw_or_wald.pval < 0.05, TRUE, FALSE)
        out.df$reverse_MR_allpass <- ifelse(out.reverse.df$pass_alltests, TRUE, FALSE)
        }
    } else {
        writeLines("nrow(out.reverse.df) =0. Not found IVs.",out.reverse.tsv.nullpath)
        if (nrow(out.df) > 0) {
        out.df$reverse_MR_sig <- NA
        out.df$reverse_MR_allpass <- NA
        }
    }

    }
    else {
        if (nrow(out.df) > 0){
        out.df$reverse_MR_sig <- NA
        out.df$reverse_MR_allpass <- NA
        }
        writeLines("nrow(otc_gwas.filter) = 0. Not found IVs.",out.reverse.tsv.nullpath)
    }
}

# write to file
if (nrow(out.df) > 0) {
    write_tsv(out.df, out.tsv.path)
    cat("Successfully write results to ", out.tsv.path, "\n", sep = "")
} else {
    writeLines("nrow(out.df) = 0. Not found IVs.", out.tsv.nullpath)
}
if (argv$bidirect) {
    if (nrow(out.reverse.df) > 0) {
        write_tsv(out.reverse.df, out.reverse.tsv.path)
        cat("Successfully write reverse results to ", out.reverse.tsv.path, "\n", sep = "")
    } else {
        writeLines("nrow(out.reverse.df) = 0. Not found IVs.", out.reverse.tsv.nullpath)
    }
}
