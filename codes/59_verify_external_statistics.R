# Independently reproduce Python numerical outputs using base R distributions.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
out <- args[[1]]
stopifnot(!file.exists(file.path(out,'independent_statistics_validation.tsv')))
x <- read.delim(file.path(out,'priority_variant_validation.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(x)==20, !anyDuplicated(paste(x$signal_unit_id,x$trait)))
y <- x[x$availability=='aligned_exact_target',]
alpha <- .05/20; margin <- log(1.1)
b <- y$external_beta; s <- y$external_se; f <- y$discovery_beta; fs <- y$discovery_se
pow <- function(e,se) pnorm(-qnorm(1-alpha/2)-abs(e)/se)+pnorm(abs(e)/se-qnorm(1-alpha/2))
d <- b-f; ds <- sqrt(s*s+fs*fs); q <- d*d/(ds*ds)
ep <- pmax(pnorm(-(b+margin)/s),pnorm((b-margin)/s))
expected <- list(normal_p=2*pnorm(-abs(b/s)),source_p_bonferroni=pmin(1,20*y$external_source_p),
 difference_beta=d,difference_se=ds,heterogeneity_q=q,heterogeneity_p=pchisq(q,df=1,lower.tail=FALSE),
 heterogeneity_p_bonferroni=pmin(1,20*pchisq(q,df=1,lower.tail=FALSE)),i2_descriptive=ifelse(q==0,0,pmax(0,(q-1)/q)),
 external_or=exp(b),external_or_lower95=exp(b-qnorm(.975)*s),external_or_upper95=exp(b+qnorm(.975)*s),
 discovery_or=exp(f),discovery_or_lower95=exp(f-qnorm(.975)*fs),discovery_or_upper95=exp(f+qnorm(.975)*fs),
 difference_ci_lower95=d-qnorm(.975)*ds,difference_ci_upper95=d+qnorm(.975)*ds,
 equivalence_p=ep,equivalence_p_bonferroni=pmin(1,ep*20),
 equivalence_beta_lower995=b-qnorm(1-alpha)*s,equivalence_beta_upper995=b+qnorm(1-alpha)*s,
 precision_power_or1_10=pow(margin,s),precision_power_half_discovery=pow(.5*f,s),precision_power_full_discovery=pow(f,s))
checks <- do.call(rbind,lapply(names(expected),function(n){
 delta <- abs(y[[n]]-expected[[n]])
 data.frame(check=n,passed=all(delta<=1e-10*pmax(1,abs(expected[[n]]))),max_absolute_error=max(delta))
}))
desired <- ifelse(b*f>0,'same',ifelse(b*f<0,'opposite','zero'))
catg <- ifelse(y$external_source_p<=alpha & desired=='same','same_direction_association_support',
 ifelse(y$external_source_p<=alpha & desired=='opposite','opposite_direction_association_support',
 ifelse(ep<=alpha,'no_primary_support_small_effect_equivalence',
 ifelse(pow(margin,s)<.8,'inconclusive_low_precision_OR1_10_scenario','inconclusive'))))
checks <- rbind(checks,data.frame(check=c('classifications','equivalence_flags','unavailable_not_called_failure','mde_power80','AIHT_n_missing','denominator_kept_20'),
 passed=c(all(y$classification==catg),all(as.logical(y$equivalence_flag)==(ep<=alpha)),
 all(x$classification[x$availability!='aligned_exact_target']=='untestable_missing_aligned_target'),
 all(abs(pow(log(y$minimum_detectable_or80),s)-.8)<1e-10),all(is.na(x$external_n_total[x$trait=='AIHT'])),nrow(x)==20),max_absolute_error=NA_real_))
write.table(checks,file.path(out,'independent_statistics_validation.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
print(checks)
stopifnot(all(checks$passed))
