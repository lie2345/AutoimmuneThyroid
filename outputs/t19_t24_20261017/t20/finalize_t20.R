#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
root <- normalizePath(getwd(),mustWork=TRUE)
out <- file.path(root,"outputs/t19_t24_20261017/t20")
rd <- function(p) read.delim(p,colClasses="character",check.names=FALSE,quote="",na.strings="NA")
wr <- function(x,p) write.table(x,p,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
rel <- function(p) sub(paste0("^",root,"/"),"",p)
sha <- function(p) strsplit(system2("sha256sum",shQuote(p),stdout=TRUE),"[[:space:]]+")[[1]][[1]]
freeze <- rd(file.path(out,"freeze_record.tsv"))
stopifnot(sha(file.path(out,"frozen_panel.tsv"))==freeze$panel_sha256,
  sha(file.path(out,"frozen_tests.tsv"))==freeze$tests_sha256)
lead_path <- file.path(root,"outputs/t19_t24_20261017/t19/tables/hla_regional_lead_variants.tsv")
status_path <- file.path(root,"outputs/t19_t24_20261017/t19/validation/status.tsv")
lead <- rd(lead_path); status <- rd(status_path)
expected <- c(E4_GRAVES_STRICT="chr6_32703430_T_G",E4_HYTHY_AI_STRICT="chr6_32659200_T_C")
stopifnot(nrow(lead)==2L,setequal(lead$trait_id,names(expected)),
  all(lead$variant_id==unname(expected[lead$trait_id])),
  status$value[status$dimension=="implementation_checks"]=="42/42_PASS",
  status$value[status$dimension=="alternative_checks"]=="33/33_PASS")
panel <- rd(file.path(out,"frozen_panel.tsv"))
lead <- lead[match(names(expected),lead$trait_id),]
stopifnot(identical(lead$variant_id,panel$variant_id[11:12]))
confirmation <- data.frame(confirmed_utc=format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC"),
  trait_id=lead$trait_id,variant_id=lead$variant_id,effect_allele=lead$effect_allele,other_allele=lead$other_allele,
  beta_discovery_ALT=lead$beta,source_p_discovery=lead$source_p,
  confirmation_status="T19_MARGINAL_LEAD_ID_AND_NUMERICAL_SOURCE_VALIDATION_CONFIRMED_BY_PRIMARY",
  lead_source=rel(lead_path),lead_source_sha256=sha(lead_path),status_source=rel(status_path),status_source_sha256=sha(status_path),
  frozen_panel_modified=FALSE,inference_strength="DESCRIPTIVE_ONLY_NOT_INDEPENDENT_SCIENTIFIC_VALIDATION")
wr(confirmation,file.path(out,"hla_selection_confirmation.tsv"))
checks <- rd(file.path(out,"validation_checks.tsv"))
stopifnot(nrow(checks)==329L,all(checks$pass=="TRUE"))
wr(data.frame(dimension=c("task_progress","source_verification","inference_strength","frozen_panel","planned_lookups","missing_or_excluded","computational_source_integration_checks","manifest_count_checks","HLA_selection_confirmation","independent_scientific_verification"),
  value=c("COMPLETED_BOUNDED_T20_LOOKUP","OFFICIAL_ENDPOINT_METADATA_AND_LOCAL_TARGET_ROWS_VERIFIED","DESCRIPTIVE_SAME_FINNGEN_DISCOVERY_SELECTED_PANEL","12_VARIANTS_UNCHANGED_FROM_PRELOOKUP_FREEZE","24_OF_24_EXACT_MATCH_GENERAL_QC_PASS","0","329/329_PASS","2/2_ENDPOINT_COUNTS_AGREE","T19_BOTH_LEADS_CONFIRMED_SEPARATELY_WITHOUT_PANEL_CHANGE","SEPARATE_PRIMARY_AGENT_GATE")),
  file.path(out,"status.tsv"))
source_manifest_path <- file.path(out,"source_checksum_manifest.tsv")
source_manifest <- rd(source_manifest_path)
for(p in c(lead_path,status_path)) if(!rel(p) %in% source_manifest$path) source_manifest <- rbind(source_manifest,
  data.frame(path=rel(p),bytes=as.character(file.info(p)$size),sha256=sha(p)))
wr(source_manifest,source_manifest_path)
files <- list.files(out,full.names=TRUE,recursive=TRUE)
files <- files[!file.info(files)$isdir & basename(files)!="artifact_checksum_manifest.tsv"]
wr(data.frame(path=rel(files),bytes=file.info(files)$size,sha256=vapply(files,sha,character(1))),
  file.path(out,"artifact_checksum_manifest.tsv"))
cat("T20 finalized. Frozen panel/tests unchanged; both T19 HLA lead IDs confirmed; 329/329 checks retained.\n")
