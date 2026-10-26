#!/usr/bin/env Rscript
# R-only rendering of the approved circular Figure 1B and Figure 2 revision.
# Inputs are frozen output extracts. No estimation or statistical tests occur.
suppressPackageStartupMessages({library(grid);library(jsonlite)})
args <- commandArgs(trailingOnly=FALSE)
script_arg <- sub("^--file=", "", args[grepl("^--file=",args)])
BASE <- normalizePath(file.path(dirname(script_arg),".."),mustWork=TRUE)
ROOT <- normalizePath(file.path(BASE,"../../../.."),mustWork=TRUE)
FIG <- file.path(BASE,"figures"); DATA <- file.path(BASE,"source_data"); REV <- file.path(BASE,"reviews")
invisible(lapply(c(FIG,DATA,REV),dir.create,recursive=TRUE,showWarnings=FALSE))
MM <- 72/25.4; W <- 180*MM; H <- 245*MM
INK <- "#18212B"; MUTED <- "#505965"; RULE <- "#D9DEE4"; FONT <- "DejaVu Sans"
class_order <- c("shared_strong_H4","distinct_signals_H3","shared_suggestive_H4","unresolved_low_purity","unresolved_insufficient_overlap")
styles <- list(
 shared_strong_H4=list(label="Strong H4",color="#0072B2",shape="circle",dash="solid"),
 distinct_signals_H3=list(label="Strong H3",color="#D55E00",shape="diamond",dash="solid"),
 shared_suggestive_H4=list(label="Suggestive H4",color="#0072B2",shape="open_circle",dash="solid"),
 unresolved_low_purity=list(label="Low purity",color="#777777",shape="cross",dash="dashed"),
 unresolved_insufficient_overlap=list(label="Low overlap",color="#777777",shape="triangle",dash="dotted"))
checks <- list()
check <- function(name,value,details=NULL) {
 checks[[length(checks)+1L]] <<- list(check=name,passed=isTRUE(value),details=details)
 if(!isTRUE(value)) stop(name,": ",paste(details,collapse="; "))
}
read_frozen <- function(name) read.delim(file.path(DATA,name),colClasses="character",check.names=FALSE,quote="",na.strings=NULL)
f1 <- read_frozen("figure1_all42_comparisons.tsv")
regions <- read_frozen("figure1_region_counts.tsv")
endpoints <- read_frozen("figure1_endpoints.tsv")
f2 <- read_frozen("figure2_all24_comparisons.tsv")
nodes <- read_frozen("figure2_nodes.tsv")
frozen_paths <- file.path(DATA,c("figure1_all42_comparisons.tsv","figure1_region_counts.tsv","figure1_endpoints.tsv","figure2_all24_comparisons.tsv","figure2_nodes.tsv"))
original_paths <- c(file.path(ROOT,"outcome_data/t06/gwas_gwas_susie_coloc.tsv"),file.path(ROOT,"outcome_data/t06/susie_signal_catalog.tsv"),file.path(ROOT,"outcome_data/t06/gwas_gwas_region_architecture.tsv"),file.path(ROOT,"config/traits.tsv"))
hash_files <- function(paths) setNames(vapply(paths,function(p)strsplit(system2("sha256sum",shQuote(p),stdout=TRUE)," ")[[1]][1],character(1)),paths)
hash_before <- hash_files(c(frozen_paths,original_paths))
check("all42 unique pair identities",nrow(f1)==42L&&!anyDuplicated(f1[c("pair_id","graves_cs","aiht_cs")]))
check("all16 region pairs",identical(regions$pair_id,sprintf("RP%02d",1:16)))
check("five class totals",identical(as.integer(table(factor(f1$classification,levels=class_order))),c(8L,11L,1L,17L,5L)))
check("15 unique Graves /16 AIHT regions",length(unique(f1$graves_region))==15L&&length(unique(f1$aiht_region))==16L)
check("chrX included",any(f1$pair_id=="RP16"&f1$chromosome=="X"))
check("all24 exact mixed comparisons",nrow(f2)==24L&&!anyDuplicated(f2[c("pair_id","graves_cs","aiht_cs")]))
check("27 exact nodes",nrow(nodes)==27L&&!anyDuplicated(nodes[c("pair_id","node_key")]))
check("mixed class totals",identical(as.integer(table(factor(f2$classification,levels=class_order))),c(6L,9L,0L,7L,2L)))
check("RP06 five unresolved",sum(f2$pair_id=="RP06"&grepl("^unresolved",f2$classification))==5L)
check("15 strong edges use local Graves CS1",all(f2$graves_cs[f2$classification%in%class_order[1:2]]=="1"))
raw_pairs<-read.delim(original_paths[1],colClasses="character",check.names=FALSE,quote="",na.strings=NULL)
raw_catalog<-read.delim(original_paths[2],colClasses="character",check.names=FALSE,quote="",na.strings=NULL)
for(i in seq_len(nrow(f1))) {
 r<-f1[i,];q<-raw_pairs[raw_pairs$pair_id==r$pair_id&raw_pairs$graves_cs==r$graves_cs&raw_pairs$aiht_cs==r$aiht_cs,]
 fields<-c("graves_region","aiht_region","classification","PP.H3","PP.H4","sensitivity_PP.H3","sensitivity_PP.H4")
 check(paste("frozen pair matches original",r$pair_id,r$graves_cs,r$aiht_cs),nrow(q)==1&&all(q[fields]==r[fields]))
}
for(i in seq_len(nrow(nodes))) {
 r<-nodes[i,];q<-raw_catalog[raw_catalog$trait_id==r$trait_id&raw_catalog$region==r$region&raw_catalog$cs==r$cs,]
 fields<-c("lead_variant","cs_size","cs_min_r2","low_purity")
 check(paste("frozen full node matches original",r$pair_id,r$side,r$cs),nrow(q)==1&&all(q[fields]==r[fields]))
}

esc <- function(s) {s<-gsub("&","&amp;",as.character(s),fixed=TRUE);s<-gsub("<","&lt;",s,fixed=TRUE);s<-gsub(">","&gt;",s,fixed=TRUE);gsub('"',"&quot;",s,fixed=TRUE)}
new_drawing <- function(slug,title,desc) {
 d<-new.env(parent=emptyenv());d$slug<-slug;d$ops<-list();d$bounds<-list()
 d$svg<-c(sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="245mm" viewBox="0 0 %.5f %.5f" role="img" aria-labelledby="%s-title %s-desc">',W,H,slug,slug),sprintf('<title id="%s-title">%s</title>',slug,esc(title)),sprintf('<desc id="%s-desc">%s</desc>',slug,esc(desc)))
 d
}
add <- function(d,svg,op) {d$svg<-c(d$svg,svg);d$ops[[length(d$ops)+1L]]<-op;invisible(NULL)}
rect <- function(d,x,y,w,h,fill="#FFFFFF",stroke=RULE,lw=.5,rx=0) add(d,sprintf('<rect x="%.4f" y="%.4f" width="%.4f" height="%.4f" rx="%.2f" fill="%s" stroke="%s" stroke-width="%.2f"/>',x,y,w,h,rx,fill,stroke,lw),list(type="rect",x=x,y=y,w=w,h=h,fill=fill,stroke=stroke,lw=lw,rx=rx))
line <- function(d,x1,y1,x2,y2,color=RULE,lw=.5,dash="solid") {
 ds<-switch(dash,dashed=' stroke-dasharray="4,2"',dotted=' stroke-dasharray="1,2"',"")
 add(d,sprintf('<line x1="%.4f" y1="%.4f" x2="%.4f" y2="%.4f" stroke="%s" stroke-width="%.2f"%s/>',x1,y1,x2,y2,color,lw,ds),list(type="polyline",x=c(x1,x2),y=c(y1,y2),color=color,lw=lw,dash=dash))
}
text <- function(d,x,y,label,size=8,color=INK,bold=FALSE,anchor="start") {
 width<-as.numeric(systemfonts::string_width(as.character(label),family=FONT,weight=if(bold)"bold" else "normal",size=size,res=72))
 left<-switch(anchor,start=x,middle=x-width/2,end=x-width)
 d$bounds[[length(d$bounds)+1L]]<-c(left=left,top=y-size,right=left+width,bottom=y+size*.25)
 add(d,sprintf('<text x="%.4f" y="%.4f" font-family="DejaVu Sans,Arial,Helvetica,sans-serif" font-size="%.2f" font-weight="%d" fill="%s" text-anchor="%s">%s</text>',x,y,size,if(bold)700L else 400L,color,anchor,esc(label)),list(type="text",x=x,y=y,label=label,size=size,color=color,bold=bold,anchor=anchor))
}
circle <- function(d,x,y,r,fill="none",stroke=RULE,lw=.6) add(d,sprintf('<circle cx="%.4f" cy="%.4f" r="%.4f" fill="%s" stroke="%s" stroke-width="%.2f"/>',x,y,r,fill,stroke,lw),list(type="circle",x=x,y=y,r=r,fill=fill,stroke=stroke,lw=lw))
polygon <- function(d,x,y,fill,stroke="none",lw=.8) add(d,sprintf('<polygon points="%s" fill="%s" stroke="%s" stroke-width="%.2f"/>',paste(sprintf("%.4f,%.4f",x,y),collapse=" "),fill,stroke,lw),list(type="polygon",x=x,y=y,fill=fill,stroke=stroke,lw=lw))
symbol <- function(d,x,y,shape,color,r=3) {
 if(shape%in%c("circle","open_circle")) circle(d,x,y,r,if(shape=="circle")color else "#FFFFFF",color,.9)
 else if(shape=="diamond")polygon(d,x+c(0,r*1.2,0,-r*1.2),y+c(-r*1.2,0,r*1.2,0),color)
 else if(shape=="square")rect(d,x-r,y-r,2*r,2*r,"#FFFFFF",color,.9)
 else if(shape=="triangle")polygon(d,x+c(0,r*1.15,-r*1.15),y+c(-r*1.2,r,r),"#FFFFFF",color,.9)
 else if(shape=="cross"){line(d,x-r,y-r,x+r,y+r,color,.9);line(d,x-r,y+r,x+r,y-r,color,.9)}
}
arrow <- function(d,x,y1,y2){line(d,x,y1,x,y2-4,MUTED,.9);polygon(d,x+c(-2.8,2.8,0),c(y2-4,y2-4,y2),MUTED)}
panel <- function(d,letter,title,y){text(d,16,y,letter,11,bold=TRUE);text(d,32,y,title,9.5,bold=TRUE)}
bezier_xy <- function(x,y,t=seq(0,1,length.out=101)) {
 weights<-rbind((1-t)^3,3*(1-t)^2*t,3*(1-t)*t^2,t^3)
 list(x=as.numeric(x%*%weights),y=as.numeric(y%*%weights))
}
bezier <- function(d,x,y,color,lw=1,dash="solid") {
 ds<-switch(dash,dashed=' stroke-dasharray="4,2"',dotted=' stroke-dasharray="1,2"',"")
 p<-bezier_xy(x,y)
 add(d,sprintf('<path d="M %.4f %.4f C %.4f %.4f %.4f %.4f %.4f %.4f" fill="none" stroke="%s" stroke-width="%.2f"%s/>',x[1],y[1],x[2],y[2],x[3],y[3],x[4],y[4],color,lw,ds),list(type="polyline",x=p$x,y=p$y,color=color,lw=lw,dash=dash))
}
render_grid <- function(d) {
 grid.newpage();pushViewport(viewport(xscale=c(0,W),yscale=c(H,0)))
 col<-function(z)if(z=="none")NA_character_ else z
 # R's lwd uses 1/96 inch, so multiply physical point width by 96/72.
 for(op in d$ops) {
  if(op$type=="rect") {
   if(op$rx>0)grid.roundrect(x=unit(op$x+op$w/2,"native"),y=unit(op$y+op$h/2,"native"),width=unit(op$w,"points"),height=unit(op$h,"points"),r=unit(op$rx,"points"),gp=gpar(fill=col(op$fill),col=col(op$stroke),lwd=op$lw*4/3))
   else grid.rect(x=unit(op$x+op$w/2,"native"),y=unit(op$y+op$h/2,"native"),width=unit(op$w,"points"),height=unit(op$h,"points"),gp=gpar(fill=col(op$fill),col=col(op$stroke),lwd=op$lw*4/3))
  }
  else if(op$type=="circle")grid.circle(x=unit(op$x,"native"),y=unit(op$y,"native"),r=unit(op$r,"points"),gp=gpar(fill=col(op$fill),col=col(op$stroke),lwd=op$lw*4/3))
  else if(op$type=="polygon")grid.polygon(x=unit(op$x,"native"),y=unit(op$y,"native"),gp=gpar(fill=col(op$fill),col=col(op$stroke),lwd=op$lw*4/3))
  else if(op$type=="polyline")grid.lines(x=unit(op$x,"native"),y=unit(op$y,"native"),gp=gpar(col=op$color,lwd=op$lw*4/3,lty=switch(op$dash,dashed="42",dotted="12","solid"),lineend="butt"))
  else if(op$type=="text")grid.text(op$label,x=unit(op$x,"native"),y=unit(op$y-op$size*.32,"native"),just=c(switch(op$anchor,start="left",middle="centre",end="right"),"centre"),gp=gpar(fontfamily=FONT,fontface=if(op$bold)"bold" else "plain",fontsize=op$size,col=op$color))
 }
 popViewport()
}
png_metadata <- function(path) {
 bytes<-readBin(path,"raw",n=file.info(path)$size)
 uint32<-function(i)readBin(bytes[i:(i+3L)],"integer",n=1L,size=4L,endian="big")
 p<-grepRaw(charToRaw("pHYs"),bytes,fixed=TRUE,all=TRUE)[1]
 list(width_px=uint32(17L),height_px=uint32(21L),dpi_x=uint32(p+4L)*.0254,dpi_y=uint32(p+8L)*.0254)
}
finish <- function(d) {
 d$svg<-c(d$svg,"</svg>");svg<-paste(d$svg,collapse="\n")
 writeLines(c('<?xml version="1.0" encoding="UTF-8"?>',svg),file.path(FIG,paste0(d$slug,".svg")),useBytes=TRUE)
 writeLines(paste0('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>',d$slug,' editable R-authored source</title><style>html,body{margin:0;background:white}svg{display:block;max-width:100%;height:auto}</style><body>',svg,'</body></html>'),file.path(FIG,paste0(d$slug,"_source.html")),useBytes=TRUE)
 grDevices::cairo_pdf(file.path(FIG,paste0(d$slug,".pdf")),width=180/25.4,height=245/25.4,family=FONT,onefile=TRUE);render_grid(d);dev.off()
 pngs<-list()
 for(dpi in c(600,180)) {
  path<-file.path(FIG,paste0(d$slug,if(dpi==180)"_preview" else "",".png"))
  ragg::agg_png(path,width=ceiling(180/25.4*dpi),height=ceiling(245/25.4*dpi),units="px",res=dpi,background="white");render_grid(d);dev.off()
  pngs[[as.character(dpi)]]<-png_metadata(path)
  info<-pngs[[as.character(dpi)]]
  check(paste(d$slug,"PNG dimensions",dpi),info$width_px==ceiling(180/25.4*dpi)&&info$height_px==ceiling(245/25.4*dpi))
  check(paste(d$slug,"PNG resolution metadata",dpi),abs(info$dpi_x-dpi)<.02&&abs(info$dpi_y-dpi)<.02)
 }
 font_output<-system2("pdffonts",shQuote(file.path(FIG,paste0(d$slug,".pdf"))),stdout=TRUE)
 check(paste(d$slug,"PDF fonts embedded"),length(font_output)>2L&&all(grepl("yes[[:space:]]+yes[[:space:]]+yes",font_output[-c(1,2)])))
 selfcheck<-system2("python3",c(shQuote("tools/diagram-design/self_check.py"),shQuote(file.path(FIG,paste0(d$slug,"_source.html")))),stdout=TRUE,stderr=TRUE)
 check(paste(d$slug,"editable HTML self_check"),is.null(attr(selfcheck,"status")),selfcheck)
 bounds<-do.call(rbind,d$bounds)
 check(paste(d$slug,"text within physical canvas"),all(bounds[,"left"]>=0&bounds[,"right"]<=W&bounds[,"top"]>=0&bounds[,"bottom"]<=H))
 list(figure=d$slug,width_mm=180,height_mm=245,minimum_text_pt=min(vapply(Filter(function(x)x$type=="text",d$ops),`[[`,numeric(1),"size")),svg_text_count=nrow(bounds),png_exports=pngs,pdffonts_output=font_output,html_selfcheck_output=selfcheck)
}
write_map <- function(name,x)write.table(x,file.path(DATA,name),sep="\t",quote=FALSE,row.names=FALSE,na="NA",fileEncoding="UTF-8")

# Figure 1A preserves the existing layout and content. Only canvas height grows.
d1<-new_drawing("figure1","Scope and outcomes of signal-pair comparisons","Existing endpoint scope and all42 tested component pairs displayed in16 equal organizational sectors. Sector angle and radius have no genomic-distance interpretation.")
rect(d1,0,0,W,H,"#FFFFFF","none")
panel(d1,"A","Scope of the existing FinnGen R13 comparison",23)
for(i in 1:2) {
 x<-c(16,274)[i];r<-endpoints[i,];title<-c("Graves’ disease","Autoimmune hypothyroidism")[i]
 rect(d1,x,40,220,72,"#F8FAFC",RULE,.6,3)
 text(d1,x+10,57,title,9,bold=TRUE);text(d1,x+10,73,r$trait_id,8)
 text(d1,x+10,91,paste0(format(as.integer(r$n_cases),big.mark=",")," cases / ",format(as.integer(r$n_controls),big.mark=",")," controls"),8.5)
 arrow(d1,x+110,112,136)
}
rect(d1,16,136,478,61,"#F3F7FA","#A9BBC9",.5,3)
text(d1,W/2,153,"16 overlapping fine-mapping region pairs",9,bold=TRUE,anchor="middle")
text(d1,W/2,169,"15 unique Graves regions and 16 unique AIHT regions",8.5,anchor="middle")
text(d1,W/2,185,"42 pairwise comparisons of detected association components",8.5,anchor="middle")
text(d1,16,214,"Finnish-enriched European ancestry · GRCh38 · summary statistics",8,MUTED)
panel(d1,"B","All tested signal pairs in 16 equal region-pair sectors",239)
cx<-W/2;cy<-450;outer<-143;inner<-66;label_r<-164
circle(d1,cx,cy,outer,"none",RULE,.6);circle(d1,cx,cy,inner,"none",RULE,.6)
mark_map<-list()
for(i in 1:16) {
 pid<-sprintf("RP%02d",i);angle<-pi/2-(i-.5)*2*pi/16
 boundary<-pi/2-(i-1)*2*pi/16
 line(d1,cx+inner*cos(boundary),cy-inner*sin(boundary),cx+outer*cos(boundary),cy-outer*sin(boundary),RULE,.6)
 rr<-regions[regions$pair_id==pid,];pair_rows<-f1[f1$pair_id==pid,]
 pair_rows<-pair_rows[order(as.integer(pair_rows$graves_cs),as.integer(pair_rows$aiht_cs)),]
 check(paste("radial order is numeric Graves CS then AIHT CS",pid),!is.unsorted(as.integer(pair_rows$graves_cs))&&all(vapply(split(as.integer(pair_rows$aiht_cs),pair_rows$graves_cs),function(z)!is.unsorted(z),logical(1))))
 for(j in seq_len(nrow(pair_rows))) {
  rad<-79+(j-1)*8.4;x<-cx+rad*cos(angle);y<-cy-rad*sin(angle);s<-styles[[pair_rows$classification[j]]]
  symbol(d1,x,y,s$shape,s$color,2.6)
  mark_map[[length(mark_map)+1L]]<-data.frame(pair_rows[j,],circular_sector=i,sector_angle_degrees=angle*180/pi,radial_slot=j,plot_x_pt=x,plot_y_pt=y,radius_pt=rad,shape=s$shape,color=s$color,stringsAsFactors=FALSE)
 }
 lx<-cx+label_r*cos(angle);ly<-cy-label_r*sin(angle)
 extra<-if(pid%in%c("RP08","RP09"))"*" else ""
 text(d1,lx,ly,paste0(pid,extra," · chr",rr$chromosome),8,bold=TRUE,anchor="middle")
 text(d1,lx,ly+11,paste0("n=",nrow(pair_rows)),8,MUTED,anchor="middle")
}
text(d1,cx,434,"42 comparisons",10,bold=TRUE,anchor="middle")
text(d1,cx,450,"16 region pairs",8.5,anchor="middle")
text(d1,cx,465,"8 H4 · 11 H3",8.5,anchor="middle")
text(d1,cx,480,"1 suggestive · 22 unresolved",8,anchor="middle")
line(d1,16,639,494,639)
for(i in 1:5) {s<-styles[[class_order[i]]];x<-c(22,116,210,315,402)[i];symbol(d1,x,650,s$shape,s$color,2.6);text(d1,x+7,653,s$label,8)}
text(d1,16,676,"Equal sectors; RP order is organizational, not genomic distance. One mark = one tested CS pair.",8,MUTED)
text(d1,16,690,"* RP08/RP09 share one Graves region. The chrX region pair is retained.",8,MUTED)
mark_map<-do.call(rbind,mark_map)
check("all42 marks mapped one-to-one",nrow(mark_map)==42L&&identical(paste(mark_map$pair_id,mark_map$graves_cs,mark_map$aiht_cs),paste(f1$pair_id,f1$graves_cs,f1$aiht_cs)))
write_map("figure1_circular_marks.tsv",mark_map)
m1<-finish(d1)

# Figure 2: six small chord plots; every source node is unique in its panel.
d2<-new_drawing("figure2","Tested component correspondence within six known risk regions","Six chord diagrams retain every one of24 tested component pairs and27 source nodes. Curves have equal width, no arrows and no causal, effect-size or genomic-distance meaning. Full node descriptions are alongside.")
rect(d2,0,0,W,H,"#FFFFFF","none")
pids<-c("RP01","RP02","RP06","RP07","RP10","RP12")
region_names<-c("PTPN22","RTKN2–ARID5B","CTLA4–ICOS","AC007563.2–DIRC3","ZAP70–TMEM131–COX5B","LPP")
node_map<-list();edge_map<-list();circle_radius<-33;gcx<-63
for(i in seq_along(pids)) {
 pid<-pids[i];top<-14+(i-1)*104;gcy<-top+51
 pn<-nodes[nodes$pair_id==pid,];pn<-pn[order(match(pn$side,c("graves","aiht")),as.integer(pn$cs)),]
 gs<-pn[pn$side=="graves",];as<-pn[pn$side=="aiht",]
 ga<-if(nrow(gs)==1)pi else c(150,210)*pi/180
 aa<-if(nrow(as)==4)c(72,24,-24,-72)*pi/180 else if(nrow(as)==3)c(60,0,-60)*pi/180 else c(40,-40)*pi/180
 pn$node_code<-paste0(ifelse(pn$side=="graves","G","A"),pn$cs)
 pn$angle_radian<-c(ga,aa);pn$plot_x_pt<-gcx+circle_radius*cos(pn$angle_radian);pn$plot_y_pt<-gcy-circle_radius*sin(pn$angle_radian)
 pn$plot_shape<-ifelse(pn$side=="graves","circle","square")
 node_map[[i]]<-pn
 circle(d2,gcx,gcy,circle_radius,"none","#DCE1E7",.6)
 pe<-f2[f2$pair_id==pid,];pe<-pe[order(!grepl("^unresolved",pe$classification),as.integer(pe$graves_cs),as.integer(pe$aiht_cs)),]
 for(j in seq_len(nrow(pe))) {
  e<-pe[j,];g<-pn[pn$side=="graves"&pn$cs==e$graves_cs,];a<-pn[pn$side=="aiht"&pn$cs==e$aiht_cs,]
  check(paste("one source node per edge",pid,e$graves_cs,e$aiht_cs),nrow(g)==1&&nrow(a)==1)
  x<-c(g$plot_x_pt,gcx-8,gcx+8,a$plot_x_pt);y<-c(g$plot_y_pt,gcy+(g$plot_y_pt-gcy)*.75,gcy+(a$plot_y_pt-gcy)*.75,a$plot_y_pt)
  s<-styles[[e$classification]];bezier(d2,x,y,s$color,1,s$dash)
  gt<-if(e$graves_cs=="2") .73 else .51;p<-bezier_xy(x,y,gt)
  symbol(d2,p$x,p$y,s$shape,s$color,2.1)
  edge_map[[length(edge_map)+1L]]<-data.frame(e,start_node=g$node_code,end_node=a$node_code,x0=x[1],y0=y[1],x1=x[2],y1=y[2],x2=x[3],y2=y[3],x3=x[4],y3=y[4],glyph_x=p$x,glyph_y=p$y,glyph_t=gt,edge_width_pt=1,edge_dash=s$dash,edge_shape=s$shape,edge_color=s$color,stringsAsFactors=FALSE)
 }
 for(j in seq_len(nrow(pn))) {
  n<-pn[j,];low<-n$low_purity=="TRUE";color<-if(low)"#9AA2AA" else "#26313B"
  symbol(d2,n$plot_x_pt,n$plot_y_pt,n$plot_shape,color,3.5)
  lx<-gcx+(circle_radius+11)*cos(n$angle_radian);ly<-gcy-(circle_radius+11)*sin(n$angle_radian)+2.5
  text(d2,lx,ly,n$node_code,8,color,bold=TRUE,anchor="middle")
 }
 text(d2,128,top+7,paste0(LETTERS[i],"  ",pid," · ",region_names[i]),9.5,bold=TRUE)
 ar<-regions[regions$pair_id==pid,]
 text(d2,128,top+20,paste0("Graves ",ar$graves_region),8,MUTED)
 text(d2,128,top+30,paste0("AIHT   ",ar$aiht_region),8,MUTED)
 for(j in seq_len(nrow(pn))) {
  n<-pn[j,];label<-paste0(n$node_code,"  ",n$lead_variant,"  n=",n$cs_size,"  min r²=",sprintf("%.3f",as.numeric(n$cs_min_r2))," ",if(n$low_purity=="TRUE")"low" else "pass")
  width<-as.numeric(systemfonts::string_width(label,family=FONT,size=8,res=72))
  check(paste("full descriptor fits",pid,n$node_code),128+width<499,paste(width,label))
  text(d2,128,top+44+(j-1)*10,label,8,if(n$low_purity=="TRUE")MUTED else INK)
 }
}
edge_map<-do.call(rbind,edge_map);node_map<-do.call(rbind,node_map)
check("all24 edges mapped one-to-one",nrow(edge_map)==24&&!anyDuplicated(edge_map[c("pair_id","graves_cs","aiht_cs")]))
check("all27 nodes mapped one-to-one",nrow(node_map)==27&&!anyDuplicated(node_map[c("pair_id","node_key")]))
check("all24 constant width",all(edge_map$edge_width_pt==1))
check("RP06 each unresolved edge retained",setequal(paste(edge_map$graves_cs[edge_map$pair_id=="RP06"&grepl("^unresolved",edge_map$classification)],edge_map$aiht_cs[edge_map$pair_id=="RP06"&grepl("^unresolved",edge_map$classification)]),c("1 2","2 1","2 2","2 3","2 4")))
line(d2,16,640,494,640)
for(i in seq_along(c(1,2,4,5))) {k<-class_order[c(1,2,4,5)[i]];s<-styles[[k]];x<-c(18,139,262,382)[i];line(d2,x,653,x+17,653,s$color,1,s$dash);symbol(d2,x+8.5,653,s$shape,s$color,2);text(d2,x+23,656,s$label,8)}
text(d2,16,676,"Neutral circles = Graves; neutral squares = AIHT. n = CS size; min r² and purity are from the source CS.",8,MUTED)
text(d2,16,690,"Each equal-width chord is a tested pair; length/crossing has no quantitative or causal meaning. GRCh38.",8,MUTED)
write_map("figure2_circular_edges.tsv",edge_map);write_map("figure2_circular_nodes.tsv",node_map)
m2<-finish(d2)

hash_after<-hash_files(c(frozen_paths,original_paths))
check("all frozen/original inputs unchanged",identical(hash_before,hash_after))
record<-list(status="PASS_R_RENDER_AND_SOURCE_CHECKS_PENDING_INDEPENDENT_REVIEW",generated_utc=format(Sys.time(),tz="UTC",usetz=TRUE),scope="R-only circular visual revision; no statistical estimation or tests rerun.",command="Rscript outputs/sci_writing_20261017/revisions/circular_v1/scripts/build_circular_figures.R",R_version=R.version.string,package_versions=lapply(c("ragg","systemfonts","jsonlite"),function(p)list(package=p,version=as.character(packageVersion(p)))),font_family=FONT,figures=list(m1,m2),source_hashes_sha256=as.list(hash_after),checks=checks,session_info=capture.output(sessionInfo()))
write_json(record,file.path(REV,"circular_build_validation.json"),pretty=TRUE,auto_unbox=TRUE,digits=15)
writeLines(c("# R circular-figure production record","","All plotting, layout, and SVG generation were implemented in R. Existing grid/cairo_pdf/ragg/systemfonts/jsonlite were used; no packages were installed. Original artifacts and frozen source extracts remain unchanged.","","Figure 1A retains its existing content and coordinates; the canvas expands to180×245mm. Figure 1B uses16 equal sectors, one mark per testedCS pair. Angular order and radius are layout only, not genome distance or effect size.","","Figure 2 includes24 equal-width tested-pair chords and27 unique component nodes. Graves nodes are neutral circles and AIHT nodes neutral squares. The source low-purity nodes remain gray. Full region definitions, full-region leads, CSsizes and purity/min-r2 descriptors are retained alongside. Edge colors/dashes/glyphs express only existing pair-level model/quality labels.","","Source extracts and plot coordinates are in source_data/*circular*.tsv. R executes only layout, display, source joins and validation, without changing inference. Visual/export and independent review results will be recorded separately."),file.path(REV,"circular_build_notes.md"))
cat(toJSON(list(figures=list(m1,m2),checks=length(checks),passed=sum(vapply(checks,`[[`,logical(1),"passed"))),auto_unbox=TRUE,pretty=TRUE),"\n")
