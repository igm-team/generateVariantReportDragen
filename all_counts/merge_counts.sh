#!~/bin/bash


join -j2 denovo.counts hemizygous.counts > both
join -j2 -a 1 denovo.counts hemizygous.counts |awk 'NF==2{print $1" "$2" 0"}' > 1only
join -j2 -a 2 denovo.counts hemizygous.counts |awk 'NF==2{print $1" 0 "$2}' > 2only
cat both 1only 2only |sort -k1,1 > tmp

join -11 -22 tmp homozygous.counts > both
join -11 -22 -a 1 tmp homozygous.counts |awk 'NF==3{print $1" "$2" "$3" 0"}' > 1only
join -11 -22 -a 2 tmp homozygous.counts |awk 'NF==2{print $2" 0 0 "$1}' > 2only
cat both 1only 2only |sort -k1,1 > tmp

join -11 -22 tmp compound.heterozygote.counts > both
join -11 -22 -a 1 tmp compound.heterozygote.counts |awk 'NF==4{print $1" "$2" "$3" "$4" 0"}' > 1only
join -11 -22 -a 2 tmp compound.heterozygote.counts |awk 'NF==2{print $2" 0 0 0 "$1}' > 2only
cat both 1only 2only |sort -k1,1 > tmp

join -11 -22 tmp possible.compound.heterozygote.counts > both
join -11 -22 -a 1 tmp possible.compound.heterozygote.counts |awk 'NF==5{print $1" "$2" "$3" "$4" "$5" 0"}' > 1only
join -11 -22 -a 2 tmp possible.compound.heterozygote.counts |awk 'NF==2{print $2" 0 0 0 0 "$1}' > 2only
cat both 1only 2only |sort -k1,1 > tmp

cat <(echo "Gene All.dnm All.hem All.hom All.chet All.pchet") tmp > gene_counts.txt

