#!/usr/bin/env Rscript
########################################################################
# Two-sample Mendelian Randomization pipeline using the TwoSampleMR package
# (bug-fixed & cleaned version)
#
# Make sure all required R packages are installed:
#   TwoSampleMR, ieugwasr, tidyverse, MRPRESSO, data.table, glue, argparser
#   (LDlinkR is only needed when --proxy + --online are both used)
########################################################################

suppressPackageStartupMessages({
    library(TwoSampleMR)
    library(ieugwasr)
    library(tidyverse)
    library(MRPRESSO)
    library(data.table)
    library(glue)
    library(argparser)
})

# =====================================================================
# 1. Command line arguments
# =====================================================================
p <- arg_parser("Perform mendelian randomization analysis using TwoSampleMR")
p <- add_argument(p, "--exp",      help = "exposure gwas summary")
p <- add_argument(p, "--otc",      help = "outcome gwas summary")
p <- add_argument(p, "--nexp",     help = "sample size of exposure",         type = "numeric")
p <- add_argument(p, "--notc",     help = "sample size of outcome",          type = "numeric")
p <- add_argument(p, "--prefix",   help = "output prefix for results",       default = "TMR_results")
p <- add_argument(p, "--expname",  help = "exposure name",                   default = "exposure")
p <- add_argument(p, "--otcname",  help = "outcome name",                    default = "outcome")
p <- add_argument(p, "--proxy",    help = "use proxy snps",                  flag    = TRUE)
p <- add_argument(p, "--online",   help = "search proxy online",             flag    = TRUE)
p <- add_argument(p, "--pop",      help = "LD population <default = EUR>",   default = "EUR")
p <- add_argument(p, "--rsq",      help = "LD r-square <default = 0.8>",     default = 0.8,  type = "numeric")
p <- add_argument(p, "--bidirect", help = "perform bidirectional MR",        flag    = TRUE)
p <- add_argument(p, "--bfile",    help = "plink LD reference panel prefix", default = "")
p <- add_argument(p, "--plink",    help = "path to plink",                   default = "")
p <- add_argument(p, "--piv",      help = "pvalue threshold for IV",         default = 5e-8, type = "numeric")
p <- add_argument(p, "--rpiv",     help = "pvalue threshold for reverse IV", default = 5e-8, type = "numeric")
# exposure GWAS column names
p <- add_argument(p, "--esnp",  help = "column name of rsid in exposure gwas",                    default = "variant_id")
p <- add_argument(p, "--ebeta", help = "column name of beta in exposure gwas",                    default = "beta")
p <- add_argument(p, "--ese",   help = "column name of se in exposure gwas",                      default = "standard_error")
p <- add_argument(p, "--echr",  help = "column name of chr in exposure gwas",                     default = "chromosome")
p <- add_argument(p, "--epos",  help = "column name of snp position in exposure gwas",            default = "base_pair_location")
p <- add_argument(p, "--eea",   help = "column name of effect allele in exposure gwas",           default = "effect_allele")
p <- add_argument(p, "--eoa",   help = "column name of other allele in exposure gwas",            default = "other_allele")
p <- add_argument(p, "--eeaf",  help = "column name of effect allele frequency in exposure gwas", default = "effect_allele_frequency")
p <- add_argument(p, "--epval", help = "column name of pvalue in exposure gwas",                  default = "p_value")
# outcome GWAS column names
p <- add_argument(p, "--osnp",  help = "column name of rsid in outcome gwas",                    default = "variant_id")
p <- add_argument(p, "--obeta", help = "column name of beta in outcome gwas",                    default = "beta")
p <- add_argument(p, "--ose",   help = "column name of se in outcome gwas",                      default = "standard_error")
p <- add_argument(p, "--ochr",  help = "column name of chr in outcome gwas",                     default = "chromosome")
p <- add_argument(p, "--opos",  help = "column name of snp position in outcome gwas",            default = "base_pair_location")
p <- add_argument(p, "--oea",   help = "column name of effect allele in outcome gwas",           default = "effect_allele")
p <- add_argument(p, "--ooa",   help = "column name of other allele in outcome gwas",            default = "other_allele")
p <- add_argument(p, "--oeaf",  help = "column name of effect allele frequency in outcome gwas", default = "effect_allele_frequency")
p <- add_argument(p, "--opval", help = "column name of pvalue in outcome gwas",                  default = "p_value")

argv <- parse_args(p)

# Belt-and-suspenders: force numeric for any threshold/sample-size arg.
# (argparser usually infers from `default`, but command-line strings like
#  "5e-8" can leak through and silently become lexicographic comparisons.)
argv$piv  <- as.numeric(argv$piv)
argv$rpiv <- as.numeric(argv$rpiv)
argv$rsq  <- as.numeric(argv$rsq)
argv$nexp <- as.numeric(argv$nexp)
argv$notc <- as.numeric(argv$notc)

# =====================================================================
# 2. Output paths & skip-if-already-done
# =====================================================================
out.tsv.path             <- paste0(argv$prefix, ".tsv")
out.tsv.nullpath         <- paste0(argv$prefix, ".null.tsv")
out.reverse.tsv.path     <- paste0(argv$prefix, "_reverse.tsv")
out.reverse.tsv.nullpath <- paste0(argv$prefix, "_reverse.null.tsv")

fwd_done <- file.exists(out.tsv.path)         || file.exists(out.tsv.nullpath)
rev_done <- file.exists(out.reverse.tsv.path) || file.exists(out.reverse.tsv.nullpath)

if ((!argv$bidirect && fwd_done) || (argv$bidirect && fwd_done && rev_done)) {
    cat("File exists! Skip.\n")
    q(save = "no")
}

# =====================================================================
# 3. Helpers
# =====================================================================

# 3.1 Wrap format_data() so we don't repeat the ten column names five times.
#     `type` controls whether the resulting frame is treated as exposure/outcome
#     by TwoSampleMR — orthogonal to which GWAS the columns came from.
fmt_exp_side <- function(gwas, type) {
    format_data(data.frame(gwas),
        snp_col           = argv$esnp,  beta_col          = argv$ebeta,
        se_col            = argv$ese,   chr_col           = argv$echr,
        pos_col           = argv$epos,  effect_allele_col = argv$eea,
        other_allele_col  = argv$eoa,   eaf_col           = argv$eeaf,
        pval_col          = argv$epval, type              = type)
}
fmt_otc_side <- function(gwas, type) {
    format_data(data.frame(gwas),
        snp_col           = argv$osnp,  beta_col          = argv$obeta,
        se_col            = argv$ose,   chr_col           = argv$ochr,
        pos_col           = argv$opos,  effect_allele_col = argv$oea,
        other_allele_col  = argv$ooa,   eaf_col           = argv$oeaf,
        pval_col          = argv$opval, type              = type)
}

# 3.2 LD-clump an exposure data frame; returns NULL if nothing remains.
clump_exposure <- function(exp_dat) {
    rsid_pval <- exp_dat %>%
        select(SNP, pval.exposure) %>%
        rename(rsid = SNP, pval = pval.exposure) %>%
        filter(!is.na(rsid) & nzchar(rsid))
    if (nrow(rsid_pval) == 0) return(NULL)

    cat("Clumping...\n")
    clumped <- tryCatch(
        ieugwasr::ld_clump_local(rsid_pval,
            clump_kb = 10000, clump_r2 = 0.001, clump_p = 1,
            plink_bin = argv$plink, bfile = argv$bfile),
        error = function(e) { print(e); data.frame() }
    )
    if (nrow(clumped) == 0) return(NULL)
    cat("Finished clump\n")
    exp_dat %>% filter(SNP %in% clumped$rsid)
}

# 3.3 Push one MR method's stats into the result row. Returns NA when the
#     method has no row (e.g. Egger isn't run when nsnp < 3).
add_method_cols <- function(out.df, mr_res, prefix, method_name) {
    sub <- mr_res %>% filter(method == method_name)
    has <- nrow(sub) > 0
    out.df[[paste0(prefix, ".b")]]       <- if (has) sub$b          else NA
    out.df[[paste0(prefix, ".se")]]      <- if (has) sub$se         else NA
    out.df[[paste0(prefix, ".b.95ci")]]  <- if (has) sub$beta_95_CI else NA
    out.df[[paste0(prefix, ".or.95ci")]] <- if (has) sub$OR_95_CI   else NA
    out.df[[paste0(prefix, ".pval")]]    <- if (has) sub$pval       else NA
    out.df
}

# =====================================================================
# 4. Proxy SNP search
# =====================================================================
get_proxy <- function(exp_dat, otc_dat, r2 = 0.8, pop = "EUR", online = TRUE,
                      bin_plink, ref_panel, outfile = "proxy") {
    exp.snp <- exp_dat$SNP
    otc.snp <- otc_dat$SNP
    otc_dat <- otc_dat %>% mutate(
        target_snp.outcome = NA, proxy_snp.outcome = NA,
        target_a1.outcome  = NA, target_a2.outcome = NA,
        proxy_a1.outcome   = NA, proxy_a2.outcome  = NA,
        proxy.outcome      = NA
    )
    miss.expsnps   <- setdiff(exp.snp, intersect(exp.snp, otc.snp))
    remain.otcsnps <- setdiff(otc.snp, intersect(exp.snp, otc.snp))

    for (miss.snp in miss.expsnps) {
        if (online) {
            library(LDlinkR)
            proxies.df     <- LDproxy(snp = miss.snp, pop = pop, token = "6fb632e022ef")
            top_proxies.df <- proxies.df %>%
                filter(R2 >= r2, Distance > 0, RS_Number %in% remain.otcsnps) %>%
                arrange(-R2) %>% head(1)
        } else {
            # shQuote is for the shell command; bare `outfile` for file.remove later.
            outfilepath_sh <- shQuote(outfile)
            system(glue(
                "{bin_plink} -bfile {ref_panel} --r2 --ld-snp {miss.snp} ",
                "--ld-window 1000000 --ld-window-kb 1000 --ld-window-r2 {r2} ",
                "--out {outfilepath_sh}_{miss.snp}_proxy"
            ))
            proxies.df     <- read_table(glue("{outfile}_{miss.snp}_proxy.ld"))
            top_proxies.df <- proxies.df %>%
                filter(SNP_A != SNP_B, SNP_B %in% remain.otcsnps) %>%
                arrange(-R2) %>% head(1) %>% rename(RS_Number = SNP_B)
        }

        if (nrow(top_proxies.df) == 0) {
            cat("Not found proxies of ", miss.snp, "\n", sep = "")
            next
        }
        cat("Found proxies of ", miss.snp, "\n", sep = "")

        rsid   <- top_proxies.df$RS_Number
        ea_exp <- (exp_dat %>% filter(SNP == miss.snp) %>% head(1))$effect_allele.exposure
        oa_exp <- (exp_dat %>% filter(SNP == miss.snp) %>% head(1))$other_allele.exposure

        # Pull target/proxy allele pairs from whichever source we used.
        if (online) {
            alleles_str   <- str_replace_all(top_proxies.df$Alleles[1], "[()]", "")
            proxy_pair    <- str_split(alleles_str, "/")[[1]]
            proxy.allele1 <- proxy_pair[1]
            proxy.allele2 <- proxy_pair[2]
            corr_first    <- str_split(top_proxies.df$Correlated_Alleles[1], ",")[[1]][1]
            target_pair   <- str_split(corr_first, "=")[[1]]
            target.allele1 <- target_pair[1]
            target.allele2 <- target_pair[2]
        } else {
            allele_info <- system(glue(
                "{bin_plink} -bfile {ref_panel} --ld {miss.snp} {rsid} ",
                "| grep 'In phase' ",
                "| awk '{{print substr($NF,1,1)\"\\n\"substr($NF,4,1)\"\\n\"",
                       "substr($NF,2,1)\"\\n\"substr($NF,5,1)}}'"
            ), intern = TRUE)
            target.allele1 <- allele_info[1]
            target.allele2 <- allele_info[2]
            proxy.allele1  <- allele_info[3]
            proxy.allele2  <- allele_info[4]
        }

        if (!identical(sort(c(target.allele1, target.allele2)),
                       sort(c(ea_exp,         oa_exp)))) {
            cat("Alleles are not consistent for ", miss.snp, "\n", sep = "")
            next
        }

        # Update the outcome row that carries the proxy SNP info.
        row_sel <- otc_dat$SNP == rsid
        ea_otc  <- otc_dat[row_sel, ]$effect_allele.outcome
        oa_otc  <- otc_dat[row_sel, ]$other_allele.outcome
        otc_dat[row_sel, "target_snp.outcome"]    <- miss.snp
        otc_dat[row_sel, "proxy_snp.outcome"]     <- rsid
        otc_dat[row_sel, "target_a1.outcome"]     <- ifelse(ea_otc == proxy.allele1, target.allele1, target.allele2)
        otc_dat[row_sel, "target_a2.outcome"]     <- ifelse(oa_otc == proxy.allele2, target.allele2, target.allele1)
        otc_dat[row_sel, "proxy_a1.outcome"]      <- ifelse(ea_otc == proxy.allele1, proxy.allele1,  proxy.allele2)
        otc_dat[row_sel, "proxy_a2.outcome"]      <- ifelse(oa_otc == proxy.allele2, proxy.allele2,  proxy.allele1)
        otc_dat[row_sel, "proxy.outcome"]         <- TRUE
        otc_dat[row_sel, "effect_allele.outcome"] <- ifelse(ea_otc == proxy.allele1, target.allele1, target.allele2)
        otc_dat[row_sel, "other_allele.outcome"]  <- ifelse(oa_otc == proxy.allele2, target.allele2, target.allele1)
        otc_dat[row_sel, "SNP"]                   <- miss.snp

        # Bug fix: file.remove takes a literal path, not a shell-quoted one.
        if (!online) {
            for (ext in c(".ld", ".log")) {
                f <- glue("{outfile}_{miss.snp}_proxy{ext}")
                if (file.exists(f)) file.remove(f)
            }
        }
    }
    otc_dat %>% filter(SNP %in% exp.snp)
}

# =====================================================================
# 5. Core MR routine
# =====================================================================
performMR <- function(exp_dat, otc_dat, outprefix,
                      expname = "exposure", otcname = "outcome", piv = 5e-8) {

    if (argv$proxy) {
        cat("Searching proxy snps...\n")
        otc_dat <- get_proxy(exp_dat, otc_dat,
            r2 = argv$rsq, online = argv$online, pop = argv$pop,
            bin_plink = argv$plink, ref_panel = argv$bfile, outfile = outprefix)
    }
    otc_dat <- otc_dat %>% mutate(
        outcome              = otcname,
        id.outcome           = otcname,
        originalname.outcome = otcname
    )

    cat("Harmonising data......\n")
    exp_otc_dat <- harmonise_data(exp_dat, otc_dat)

    # ---- MR-PRESSO outlier removal (only meaningful with >=4 IVs) ----
    exp_otc_mrPresso <- NULL
    if (sum(exp_otc_dat$mr_keep) >= 4) {
        exp_otc_mrPresso <- tryCatch(
            mr_presso(
                BetaOutcome     = "beta.outcome",  BetaExposure = "beta.exposure",
                SdOutcome       = "se.outcome",    SdExposure   = "se.exposure",
                OUTLIERtest     = TRUE,            DISTORTIONtest = FALSE,
                data            = exp_otc_dat %>% filter(mr_keep == TRUE),
                SignifThreshold = 0.05),
            error = function(e) {
                cat("MR-PRESSO failed: ", conditionMessage(e), "\n", sep = "")
                NULL
            }
        )
        if (!is.null(exp_otc_mrPresso)) {
            outlier.df <- exp_otc_mrPresso$`MR-PRESSO results`$`Outlier Test`
            if (!is.null(outlier.df) && nrow(outlier.df) > 0) {
                # Bug fix: assign SNP rownames BEFORE writing, write a real TSV,
                # and use this run's prefix so reverse MR doesn't clobber forward.
                row.names(outlier.df) <- (exp_otc_dat %>% filter(mr_keep == TRUE))$SNP
                write.table(outlier.df,
                            paste0(outprefix, "_outlierPresso.tsv"),
                            sep = "\t", quote = FALSE, col.names = NA)
                keep_snps  <- row.names(outlier.df %>% filter(Pvalue > 0.05))
                exp_otc_dat <- exp_otc_dat %>% filter(SNP %in% keep_snps)
            } else {
                cat("No outliers found or detected.\n")
            }
        }
    }

    # ---- Need at least one IV after outlier removal ----
    if (is.null(exp_otc_dat) || nrow(exp_otc_dat) == 0 || sum(exp_otc_dat$mr_keep) < 1) {
        return(data.frame())
    }

    exp_otc_dat <- exp_otc_dat %>% mutate(
        samplesize.exposure = argv$nexp,
        samplesize.outcome  = argv$notc
    )

    # ---- MR estimation ----
    cat("Performing MR...\n")
    exp_otc_res <- mr(exp_otc_dat,
        method_list = c("mr_wald_ratio", "mr_ivw", "mr_egger_regression",
                        "mr_weighted_median", "mr_simple_median", "mr_weighted_mode"))
    exp_otc_res <- generate_odds_ratios(exp_otc_res) %>% mutate(
        beta_95_CI = paste0(round(b,  3), " [", round(lo_ci,    3), "-", round(up_ci,    3), "]"),
        OR_95_CI   = paste0(round(or, 3), " [", round(or_lci95, 3), "-", round(or_uci95, 3), "]")
    )

    # ---- Single-SNP analyses + diagnostic plots ----
    exp_otc_res_single <- mr_singlesnp(exp_otc_dat)
    res_loo <- mr_leaveoneout(exp_otc_dat)
    ggsave(mr_scatter_plot(exp_otc_res, exp_otc_dat)[[1]],
           file = glue("{outprefix}_scatterplot.pdf"),     width = 7, height = 7)
    ggsave(mr_forest_plot(exp_otc_res_single)[[1]],
           file = glue("{outprefix}_forestplot.pdf"),      width = 7, height = 7)
    ggsave(mr_funnel_plot(exp_otc_res_single)[[1]],
           file = glue("{outprefix}_funnelplot.pdf"),      width = 7, height = 7)
    ggsave(mr_leaveoneout_plot(res_loo)[[1]],
           file = glue("{outprefix}_leaveoneoutplot.pdf"), width = 7, height = 7)
    write_tsv(exp_otc_dat,        glue("{outprefix}_harmonise.tsv"))
    write_tsv(exp_otc_res_single, glue("{outprefix}_single.tsv"))
    write_tsv(res_loo,            glue("{outprefix}_leaveoneout.tsv"))

    # ---- Primary estimate (Wald / IVW) ----
    primary <- exp_otc_res %>%
        filter(method %in% c("Wald ratio", "Inverse variance weighted"))
    out.df <- data.frame(
        exposure            = expname,
        outcome             = otcname,
        ivw_or_wald.b       = primary$b,
        nsnp                = primary$nsnp,
        ivw_or_wald.se      = primary$se,
        ivw_or_wald.b.95ci  = primary$beta_95_CI,
        ivw_or_wald.or.95ci = primary$OR_95_CI,
        ivw_or_wald.pval    = primary$pval
    )
    cat("Finished MR\n")
    cat("Performing sensitive MR...\n")

    # ---- Sensitivity 1: Heterogeneity (need >=2 IVs) ----
    if (sum(exp_otc_dat$mr_keep) >= 2) {
        het <- mr_heterogeneity(exp_otc_dat, method_list = "mr_ivw")
        I2  <- (het %>% mutate(I2 = (Q - Q_df) / Q * 100) %>%
                 filter(method == "Inverse variance weighted"))$I2
        out.df$q.het <- het$Q_pval
        out.df$I2    <- max(0, I2)
    } else {
        out.df$q.het <- NA
        out.df$I2    <- NA
    }

    # ---- Sensitivity 2: alternative methods + Egger intercept (need >=3) ----
    method_map <- list(
        eggReg     = "MR Egger",
        weightMed  = "Weighted median",
        weightMode = "Weighted mode",
        simpleM    = "Simple median"
    )
    for (nm in names(method_map)) {
        out.df <- add_method_cols(out.df, exp_otc_res, nm, method_map[[nm]])
    }
    if (sum(exp_otc_dat$mr_keep) >= 3) {
        eggInt <- mr_pleiotropy_test(exp_otc_dat)
        out.df$p.eggerIntercept  <- eggInt$pval
        out.df$eggerIntercept    <- eggInt$egger_intercept
        out.df$se.eggerIntercept <- eggInt$se
    } else {
        out.df$p.eggerIntercept  <- NA
        out.df$eggerIntercept    <- NA
        out.df$se.eggerIntercept <- NA
    }

    # ---- Sensitivity 3: MR-PRESSO summary + same-direction check ----
    out.df$presso.b       <- NA
    out.df$presso.se      <- NA
    out.df$presso.b.95ci  <- NA
    out.df$presso.pval    <- NA
    out.df$presso.or      <- NA
    out.df$presso.or.95ci <- NA
    if (!is.null(exp_otc_mrPresso)) {
        out.df$p.globalPleiotropy <- exp_otc_mrPresso$`MR-PRESSO results`$`Global Test`$Pvalue
        out.df$RSSobs             <- exp_otc_mrPresso$`MR-PRESSO results`$`Global Test`$RSSobs
        # Row 2 is the outlier-corrected estimate when it exists; row 1 otherwise.
        idx  <- ifelse(is.na(exp_otc_mrPresso$`Main MR results`[2, 3]), 1, 2)
        b_p  <- exp_otc_mrPresso$`Main MR results`$`Causal Estimate`[idx]
        se_p <- exp_otc_mrPresso$`Main MR results`$Sd[idx]
        out.df$presso.b       <- b_p
        out.df$presso.se      <- se_p
        out.df$presso.pval    <- exp_otc_mrPresso$`Main MR results`$`P-value`[idx]
        out.df$presso.or      <- exp(b_p)
        lci <- b_p - 1.96 * se_p
        uci <- b_p + 1.96 * se_p
        out.df$presso.b.95ci  <- paste0("[", lci,      "-", uci,      "]")
        out.df$presso.or.95ci <- paste0("[", exp(lci), "-", exp(uci), "]")
        same_direction <- all(c(exp_otc_res$b, b_p) > 0) | all(c(exp_otc_res$b, b_p) < 0)
    } else {
        out.df$p.globalPleiotropy <- NA
        out.df$RSSobs             <- NA
        same_direction <- all(exp_otc_res$b > 0) | all(exp_otc_res$b < 0)
    }
    out.df$same_direction <- same_direction

    # ---- Sensitivity 4: Steiger directionality / reverse-causation flag ----
    steig <- directionality_test(exp_otc_dat)
    out.df$snp_r2.exposure <- steig$snp_r2.exposure
    out.df$snp_r2.outcome  <- steig$snp_r2.outcome
    out.df$steiger_pval    <- steig$steiger_pval
    out.df$reverse_caution <- (steig$snp_r2.exposure < steig$snp_r2.outcome) &
                              (steig$steiger_pval < 0.05)

    # ---- Sensitivity 5: overall pass flag (cumulative on nsnp tier) ----
    nsnp <- out.df$nsnp
    pass <- (out.df$ivw_or_wald.pval < 0.05) & (!out.df$reverse_caution)
    if (nsnp >= 2) pass <- pass & !(out.df$I2 > 50 & out.df$q.het < 0.05)
    if (nsnp >= 3) pass <- pass & (out.df$p.eggerIntercept > 0.05) & out.df$same_direction
    if (nsnp >= 4) pass <- pass & (out.df$p.globalPleiotropy > 0.05)
    out.df$pass_alltests <- ifelse(pass, TRUE, FALSE)

    cat("Finished sensitive test\n")
    out.df
}

# =====================================================================
# 6. Main: forward MR
# =====================================================================
cat("Reading exposure GWAS...\n")
exp_gwas <- fread(argv$exp)
cat("Reading outcome GWAS...\n")
otc_gwas <- fread(argv$otc)

out.df <- data.frame()
exp_gwas.filter <- exp_gwas %>%
    filter(get(argv$epval) < argv$piv & str_detect(get(argv$esnp), "^rs\\d+"))

if (nrow(exp_gwas.filter) > 0) {
    exp_dat <- fmt_exp_side(exp_gwas.filter, "exposure")
    otc_dat <- fmt_otc_side(otc_gwas,        "outcome")
    exp_dat_clumped <- clump_exposure(exp_dat)
    if (is.null(exp_dat_clumped)) {
        writeLines("Not found IVs.", out.tsv.nullpath)
    } else {
        cat("Perform MR...\n")
        out.df <- performMR(exp_dat_clumped, otc_dat, argv$prefix,
                            expname = argv$expname, otcname = argv$otcname,
                            piv = argv$piv)
    }
} else {
    writeLines("nrow(exp_gwas.filter) = 0", out.tsv.nullpath)
}

# =====================================================================
# 7. Main: reverse MR (optional)
# =====================================================================
out.reverse.df <- data.frame()

if (argv$bidirect) {
    otc_gwas.filter <- otc_gwas %>%
        filter(get(argv$opval) < argv$rpiv, str_detect(get(argv$osnp), "^rs\\d+"))

    if (nrow(otc_gwas.filter) > 0) {
        # Reverse direction: outcome GWAS becomes the exposure, and vice versa.
        exp_rev_dat <- fmt_otc_side(otc_gwas.filter, "exposure")
        otc_rev_dat <- fmt_exp_side(exp_gwas,        "outcome")
        exp_rev_dat_clumped <- clump_exposure(exp_rev_dat)
        if (is.null(exp_rev_dat_clumped)) {
            writeLines("Not found IVs.", out.reverse.tsv.nullpath)
        } else {
            cat("Perform reverse MR...\n")
            out.reverse.df <- performMR(
                exp_rev_dat_clumped, otc_rev_dat,
                paste0(argv$prefix, "_reverse"),
                expname = argv$otcname, otcname = argv$expname,
                piv = argv$rpiv
            )
        }

        # Cross-link forward and reverse: each side carries the OTHER side's
        # significance + overall-pass status under `reverse_MR_*` columns.
        if (nrow(out.reverse.df) > 0 && "ivw_or_wald.pval" %in% names(out.reverse.df)) {
            out.reverse.df$reverse_MR_sig     <- NA
            out.reverse.df$reverse_MR_allpass <- NA
            if (nrow(out.df) > 0) {
                out.reverse.df$reverse_MR_sig     <- ifelse(out.df$ivw_or_wald.pval < 0.05, TRUE, FALSE)
                out.reverse.df$reverse_MR_allpass <- ifelse(out.df$pass_alltests,           TRUE, FALSE)
                out.df$reverse_MR_sig             <- ifelse(out.reverse.df$ivw_or_wald.pval < 0.05, TRUE, FALSE)
                out.df$reverse_MR_allpass         <- ifelse(out.reverse.df$pass_alltests,           TRUE, FALSE)
            }
        } else {
            writeLines("nrow(out.reverse.df) = 0. Not found IVs.", out.reverse.tsv.nullpath)
            if (nrow(out.df) > 0) {
                out.df$reverse_MR_sig     <- NA
                out.df$reverse_MR_allpass <- NA
            }
        }
    } else {
        if (nrow(out.df) > 0) {
            out.df$reverse_MR_sig     <- NA
            out.df$reverse_MR_allpass <- NA
        }
        writeLines("nrow(otc_gwas.filter) = 0. Not found IVs.", out.reverse.tsv.nullpath)
    }
}

# =====================================================================
# 8. Write outputs
# =====================================================================
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
