#!/usr/bin/env Rscript
# Cross-checks are computational source/integration checks, not independent scientific verification.
options(stringsAsFactors = FALSE, digits = 17)
root <- normalizePath(getwd(), mustWork = TRUE)
out <- file.path(root, "outputs/t19_t24_20261017/t20")
rd <- function(p) read.delim(p, colClasses = "character", check.names = FALSE, quote = "", na.strings = "NA")
wr <- function(x,p) write.table(x,p,sep="\t",row.names=FALSE,quote=FALSE,na="NA")
checks <- list()
check <- function(name, pass, detail="") {
  checks[[length(checks)+1L]] <<- data.frame(check=name,pass=isTRUE(pass),detail=as.character(detail))
}
equal <- function(a,b,tol=2e-13) isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol,check.attributes=FALSE))
panel <- rd(file.path(out,"frozen_panel.tsv"))
result <- rd(file.path(out,"lookup_results.tsv"))
planned <- rd(file.path(out,"frozen_tests.tsv"))
freeze <- rd(file.path(out,"freeze_record.tsv"))
traits <- rd(file.path(root,"config/traits.tsv"))
manifest <- rd(file.path(root,"raw_data/manifests/finngen_R13_manifest.tsv"))
check("twelve_distinct_variants",nrow(panel)==12L && !anyDuplicated(panel$variant_id))
check("twentyfour_planned_results",nrow(planned)==24L && nrow(result)==24L)
check("exact_planned_test_family",setequal(paste(result$endpoint,result$variant_id),paste(planned$endpoint,planned$variant_id)))
check("freeze_before_queries",file.info(file.path(out,"freeze_record.tsv"))$mtime < file.info(file.path(out,"GRAVES_OPHT.raw_position_rows.tsv"))$mtime)
check("general_qc_contract_unchanged",any(grepl('maf_min_general: 0.005',readLines(file.path(root,"config/analysis.yml")),fixed=TRUE)))
extraction <- list()
for (ep in c("GRAVES_OPHT","E4_GRAVES_OPHT_STRICT")) {
  for (kind in c("raw","harmonized")) {
    input <- if (kind=="raw") file.path(root,"raw_data/finngen_r13/summary_stats",paste0("finngen_R13_",ep,".gz")) else file.path(root,"outcome_data/harmonized_gwas",paste0(ep,".harmonized.tsv.gz"))
    output <- file.path(out,paste0(ep,".",kind,"_sequential_panel.tsv"))
    cmd <- paste("gzip -cd",shQuote(input),"| awk -F '\t' -v",shQuote(paste0("mode=",kind)),"-f",shQuote(file.path(out,"select_panel.awk")),shQuote(file.path(out,"frozen_panel.tsv")),"- >",shQuote(output))
    start <- format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC")
    message("Sequential cross-check: ",ep," ",kind," ",start)
    status <- system2("bash",c("-o","pipefail","-c",shQuote(cmd)))
    extraction[[length(extraction)+1L]] <- data.frame(started_utc=start,finished_utc=format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC"),command=cmd,exit_status=status)
    wr(do.call(rbind,extraction),file.path(out,"sequential_extraction_commands.tsv"))
    check(paste(ep,kind,"sequential_extraction_exit_zero",sep="/"),status==0L)
    if (status != 0L) stop("Sequential extraction failed")
  }
  raw <- rd(file.path(out,paste0(ep,".raw_sequential_panel.tsv")))
  raw$variant_id <- paste0("chr",sub("^chr","",raw[["#chrom"]]),"_",raw$pos,"_",raw$ref,"_",raw$alt)
  harm <- rd(file.path(out,paste0(ep,".harmonized_sequential_panel.tsv")))
  tabix <- rd(file.path(out,paste0(ep,".raw_position_rows.tsv")))
  tabix$variant_id <- paste0("chr",sub("^chr","",tabix[["#chrom"]]),"_",tabix$pos,"_",tabix$ref,"_",tabix$alt)
  check(paste(ep,"sequential_raw_exactly_twelve_unique",sep="/"),nrow(raw)==12L && !anyDuplicated(raw$variant_id))
  check(paste(ep,"sequential_harmonized_exactly_twelve_unique",sep="/"),nrow(harm)==12L && !anyDuplicated(harm$variant_id))
  check(paste(ep,"all_panel_ids_present_raw_harmonized",sep="/"),setequal(raw$variant_id,panel$variant_id) && setequal(harm$variant_id,panel$variant_id))
  visible <- readLines(file.path(root,"outputs/t19_t24_20261017/sources",paste0(ep,"_visible.txt")))
  check(paste(ep,"official_population_noncase_definition",sep="/"),any(visible=="Controls for this endpoint are individuals that are not cases."))
  for (id in panel$variant_id) {
    key <- paste(ep,id,sep="/")
    rr <- result[result$endpoint==ep & result$variant_id==id,,drop=FALSE]
    x <- raw[raw$variant_id==id,,drop=FALSE]; y <- harm[harm$variant_id==id,,drop=FALSE]; t <- tabix[tabix$variant_id==id,,drop=FALSE]
    check(paste0(key,"/source_row_unique"),nrow(rr)==1L && nrow(x)==1L && nrow(y)==1L && nrow(t)==1L)
    if (any(c(nrow(rr),nrow(x),nrow(y),nrow(t))!=1L)) next
    check(paste0(key,"/tabix_equals_sequential_all_raw_fields"),identical(unname(as.list(x)),unname(as.list(t))))
    check(paste0(key,"/exact_position_ref_alt_build"),rr$genome_build=="GRCh38" && rr$chrom==x[["#chrom"]] && rr$pos==x$pos && rr$ref==x$ref && rr$alt==x$alt && rr$ref==y$ref && rr$alt==y$alt)
    check(paste0(key,"/unmodified_source_strings"),rr$beta_source_string==x$beta && rr$se_source_string==x$sebeta && rr$source_p_string==x$pval && rr$eaf_source_string==x$af_alt)
    check(paste0(key,"/raw_vs_harmonized_numeric_fields"),equal(c(x$beta,x$sebeta,x$pval,x$af_alt),c(y$beta,y$se,y$p,y$eaf)))
    check(paste0(key,"/ALT_effect_orientation"),rr$effect_allele==x$alt && rr$other_allele==x$ref && y$effect_allele==x$alt && y$other_allele==x$ref)
    b <- as.numeric(x$beta); s <- as.numeric(x$sebeta); p <- as.numeric(x$pval); f <- as.numeric(x$af_alt)
    pass <- is.finite(b) && is.finite(s) && s>0 && is.finite(p) && p>0 && p<=1 && is.finite(f) && f>0 && f<1 && min(f,1-f)>=0.005
    check(paste0(key,"/general_qc"),pass && rr$qc_status=="PASS_GENERAL_QC")
    # Fixed normal quantile provides a separate numerical path from production qnorm().
    independent_ci <- exp(c(b,b-1.959963984540054*s,b+1.959963984540054*s))
    check(paste0(key,"/OR_normal95CI_from_raw"),equal(independent_ci,c(rr$OR_ALT,rr$OR_lower95_normal_approx,rr$OR_upper95_normal_approx)))
    check(paste0(key,"/OR_CI_order_and_direction"),as.numeric(rr$OR_lower95_normal_approx)<as.numeric(rr$OR_ALT) && as.numeric(rr$OR_ALT)<as.numeric(rr$OR_upper95_normal_approx) && sign(as.numeric(rr$OR_ALT)-1)==sign(b))
    check(paste0(key,"/optional_Bonferroni_all24"),equal(rr$p_bonferroni24_optional_descriptive,min(1,p*24)) && as.logical(rr$source_p_below_0_05_over_24)==(p<0.05/24))
    meta <- traits[traits$trait_id==ep,,drop=FALSE]
    check(paste0(key,"/counts_match_config_and_harmonized"),rr$n_cases==meta$n_cases && rr$n_controls==meta$n_controls && rr$n_cases==y$n_cases && rr$n_controls==y$n_controls)
    for (tr in c("graves","aiht")) {
      d <- as.numeric(rr[[paste0(tr,"_discovery_beta_ALT")]])
      expected <- if (is.na(d)) "NOT_ASSESSED_DISCOVERY_CONTEXT_NOT_IN_PANEL_SOURCE" else if (b==0 || d==0) "ZERO_BETA" else if ((b>0)==(d>0)) "same" else "opposite"
      check(paste0(key,"/",tr,"_discovery_direction"),rr[[paste0("direction_vs_",tr,"_discovery")]]==expected)
    }
  }
}
report <- do.call(rbind,checks)
wr(report,file.path(out,"validation_checks.tsv"))
cat(sum(report$pass),"/",nrow(report)," computational/source/integration checks passed\n",sep="")
if (!all(report$pass)) { print(report[!report$pass,]); quit(status=1L) }
