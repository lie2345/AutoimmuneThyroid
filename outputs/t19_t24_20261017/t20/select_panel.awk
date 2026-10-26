BEGIN { OFS="\t" }
NR==FNR { if (FNR>1) wanted[$2]=1; next }
FNR==1 { print; next }
{
  if (mode=="raw") {
    chrom=$1; sub(/^chr/,"",chrom)
    key="chr" chrom "_" $2 "_" $3 "_" $4
  } else { key=$1 }
  if (key in wanted) print
}
