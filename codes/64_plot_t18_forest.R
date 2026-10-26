# Display source effects only; no new association threshold or pooled estimate.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
out <- args[[1]]
d <- read.delim(file.path(out,'forest_data.tsv'),stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(all(d$available %in% c('True','False')))
d$available <- d$available == 'True'
stopifnot(nrow(d)==40, sum(d$available)==30, !anyDuplicated(paste(d$signal_unit_id,d$trait,d$cohort)))
lims <- range(c(d$lower95,d$upper95),na.rm=TRUE)
lims <- exp(log(lims)+c(-.10,.10))
ticks <- c(.5,.6,.7,.8,1,1.2,1.5,2,2.5,3)
ticks <- ticks[ticks>=lims[1] & ticks<=lims[2]]
draw <- function(trait) {
  x <- d[d$trait==trait,]
  ids <- sort(unique(x$signal_unit_id)); n <- length(ids)
  par(mar=c(4.5,14,4,14),oma=c(3,0,0,0),family='sans',las=1,xpd=NA)
  plot(1,1,type='n',log='x',xlim=lims,ylim=c(.5,n+1),axes=FALSE,xlab='',ylab='',xaxs='i',yaxs='i')
  segments(1,.55,1,n+.40,col='#888888',lty=2,lwd=1,xpd=FALSE)
  axis(1,at=ticks,labels=ticks,cex.axis=.9)
  mtext('Odds ratio (log scale); effect allele shown at left',side=1,line=2.6,cex=.85)
  title(main=paste(trait,'- FinnGen and BBJ representative variants'),line=2.3,cex.main=1.05)
  mtext('Source P threshold: 0.0025 (20 planned tests); discovery classes are not externally revalidated',side=3,line=.7,cex=.72)
  usr <- par('usr'); left <- 10^(usr[1]-.035*diff(usr[1:2])); right <- 10^(usr[2]+.035*diff(usr[1:2]))
  text(left,n+.75,'Representative / nearest gene; EA / OA',adj=1,cex=.73,font=2)
  text(right,n+.75,'BBJ source P / variant-level result',adj=0,cex=.73,font=2)
  for (i in seq_along(ids)) {
    r <- x[x$signal_unit_id==ids[i],]; y <- n-i+1
    g <- r[r$cohort=='FinnGen',]; b <- r[r$cohort=='BBJ',]
    gene <- ifelse(is.na(g$nearest_genes),'unannotated',g$nearest_genes)
    text(left,y+.13,paste(ids[i],gene),adj=1,cex=.85)
    text(left,y-.18,paste0(g$effect_allele,' / ',g$other_allele,' | ',g$discovery_class),adj=1,cex=.68,col='#555555')
    for (cohort in c('FinnGen','BBJ')) {
      z <- r[r$cohort==cohort,]; yy <- y+ifelse(cohort=='FinnGen',.14,-.14)
      col <- ifelse(cohort=='FinnGen','#376DA0','#222222'); pch <- ifelse(cohort=='FinnGen',15,16)
      if(z$available) {
        stopifnot(z$lower95>=lims[1],z$upper95<=lims[2])
        segments(z$lower95,yy,z$upper95,yy,col=col,lwd=1.5)
        segments(c(z$lower95,z$upper95),yy-.035,c(z$lower95,z$upper95),yy+.035,col=col)
        points(z$odds_ratio,yy,pch=pch,col=col,cex=.8)
      }
    }
    ptext <- ifelse(b$available,formatC(b$source_p,format='g',digits=3),'NA')
    tag <- if(b$external_classification=='same_direction_association_support') 'Same-direction support' else if(!b$available) 'Untestable (BBJ missing)' else 'Inconclusive; low precision*'
    text(right,y+.13,paste0('P = ',ptext),adj=0,cex=.79)
    text(right,y-.18,tag,adj=0,cex=.73)
  }
  legend('topright',legend=c('FinnGen','BBJ'),pch=c(15,16),col=c('#376DA0','#222222'),bty='n',horiz=TRUE,cex=.75,inset=c(0,-.02))
  mtext('Intervals: ordinary Wald 95% CI, not SPA-calibrated; source P determines association support, not CI exclusion of 1.',side=1,outer=TRUE,line=.4,cex=.72)
  mtext('*Fixed OR=1.10 scenario; all 10 evaluable tests have <80% precision-scenario power, including the 3 supported variants.',side=1,outer=TRUE,line=1.6,cex=.70)
}
for(trait in c('Graves','AIHT')) {
  paths <- file.path(out,paste0('forest_',trait,c('.svg','.png')))
  stopifnot(!any(file.exists(paths)))
  svg(paths[1],width=13,height=8,family='sans'); draw(trait); dev.off()
  png(paths[2],width=2600,height=1600,res=200,type='cairo'); draw(trait); dev.off()
}
writeLines(c('40 plot rows; 30 observed cohort effects; 10 explicit BBJ missing rows.',
             'Identical x limits across phenotypes; every interval fits; no pooled effects.'),file.path(out,'forest_checks.txt'))
