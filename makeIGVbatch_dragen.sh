#!/bin/bash
#Input Format
#Family.ID       Child   Mother  Father  Flag    Variant.ID Gene.Name

IN=$1

awk '{print $2"\n"$3"\n"$4}' <(sed '1d' $IN) |sort -u > $(dirname $IN)/$(basename ${IN%.*}).list

while read line; do  
    alignLoc="NULL"
    if [[ $alignLoc == "NULL" ]]
        then tmp=$(find \
        $(mysql --defaults-group-suffix=sequencedragen --defaults-file=/nfs/goldstein/software/Bioinformatics_Tools/generateVariantReport/.my.cnf -e "select AlignSeqFileLoc from dragen_qc_metrics inner join dragen_sample_metadata on dragen_sample_metadata.pseudo_prepid = dragen_qc_metrics.pseudo_prepid where dragen_sample_metadata.sample_name = '$line' ;" |head -n2|tail -n1|sed 's/_temp[0-9]\?/18/')/$line* \
        -name '*recal.bam')
    else
        tmp=$(find \
        $alignLoc \
        -name '*final.bam') #final
    fi

    if [ -z "${tmp//\/\//\/}" ]
        then 
            tmp=$(find $alignLoc -name '*.bam')
            if [ -z "${tmp//\/\//\/}" ]
            then >&2 echo "ERROR: Did not find BAM for $line"
            continue
            fi
    elif [ $(du -sh $tmp|wc -l) -gt 1 ]
        then >&2 echo "ERROR: Found more than one BAM for $line"
             >&2 echo "$tmp"
             continue
             fi
    
    if [ $alignLoc == "NULL" ];then echo $tmp
    elif [[ -f $tmp ]];then echo $tmp
    else echo "$(dirname $tmp)/$(basename $(readlink $tmp))";fi
    

done<$(dirname $IN)/$(basename ${IN%.*}).list|sort -u > $(dirname $IN)/$(basename ${IN%.*}).bamloc

mkdir $(dirname $IN)/BAMS
while read i;
    do
       PRO_LOC=$(grep -w $(echo $i|awk '{print $2}') $(dirname $IN)/$(basename ${IN%.*}).bamloc|head -n1)
       if [ "$PRO_LOC" = "" ] ;then PRO_LOC="NA";fi
       DAD_LOC=$(grep -w $(echo $i|awk '{print $3}') $(dirname $IN)/$(basename ${IN%.*}).bamloc|head -n1)
       if [ "$DAD_LOC" = "" ] ;then DAD_LOC="NA";fi
       MOM_LOC=$(grep -w $(echo $i|awk '{print $4}') $(dirname $IN)/$(basename ${IN%.*}).bamloc|head -n1)
       if [ "$MOM_LOC" = "" ] ;then MOM_LOC="NA";fi
       VAR_CHR=$(echo $i|awk '{print $6}')
       VAR_LOC=$(echo $i|awk '{print $7}')
       /nfs/goldstein/software/jdk1.8.0_05/bin/java -Xmx24g -jar /nfs/goldstein/software/GATK/GATK-3.6.0-ArchivedVersion-g89b7209-patched/GenomeAnalysisTK.jar -R /scratch/HS_Build37/BWA_INDEX_hs37d5/hs37d5.fa -T HaplotypeCaller -L $VAR_CHR:$((VAR_LOC-100))-$((VAR_LOC+100)) -I $PRO_LOC -o $(dirname $IN)/BAMS/$(echo $i|awk    '{print $2"."$6"."$7".vcf.gz"}') -stand_call_conf 20 -stand_emit_conf 20 --emitRefConfidence GVCF -GQB 5 -GQB 15 -GQB 20 -GQB 60 --variant_index_type LINEAR --variant_index_parameter 128000 --dbsnp /nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/dbsnp_138.b37.vcf.gz -nct 1 --bamOutput $(dirname $IN)/BAMS/$(echo $i|awk   '{print $2"."$6"."$7".bam"}') -forceActive -disableOptimizations
       /nfs/goldstein/software/jdk1.8.0_05/bin/java -Xmx24g -jar /nfs/goldstein/software/GATK/GATK-3.6.0-ArchivedVersion-g89b7209-patched/GenomeAnalysisTK.jar -R /scratch/HS_Build37/BWA_INDEX_hs37d5/hs37d5.fa -T HaplotypeCaller -L $VAR_CHR:$((VAR_LOC-100))-$((VAR_LOC+100)) -I $DAD_LOC -o $(dirname $IN)/BAMS/$(echo $i|awk    '{print $3"."$6"."$7".vcf.gz"}') -stand_call_conf 20 -stand_emit_conf 20 --emitRefConfidence GVCF -GQB 5 -GQB 15 -GQB 20 -GQB 60 --variant_index_type LINEAR --variant_index_parameter 128000 --dbsnp /nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/dbsnp_138.b37.vcf.gz -nct 1  --bamOutput $(dirname $IN)/BAMS/$(echo $i|awk   '{print $3"."$6"."$7".bam"}') -forceActive -disableOptimizations
       if [ $(samtools view $(dirname $IN)/BAMS/$(echo $i|awk '{print $3"."$6"."$7".bam"}')|wc -l) -eq 0 ];then samtools view -b $DAD_LOC $VAR_CHR:$((VAR_LOC-100))-$((VAR_LOC+100)) > $(dirname $IN)/BAMS/$(echo $i|awk   '{print $3"."$6"."$7".bam"}');samtools index $(dirname $IN)/BAMS/$(echo $i|awk '{print $3"."$6"."$7".bam"}');fi
       /nfs/goldstein/software/jdk1.8.0_05/bin/java -Xmx24g -jar /nfs/goldstein/software/GATK/GATK-3.6.0-ArchivedVersion-g89b7209-patched/GenomeAnalysisTK.jar -R /scratch/HS_Build37/BWA_INDEX_hs37d5/hs37d5.fa -T HaplotypeCaller -L $VAR_CHR:$((VAR_LOC-100))-$((VAR_LOC+100)) -I $MOM_LOC -o $(dirname $IN)/BAMS/$(echo $i|awk    '{print $4"."$6"."$7".vcf.gz"}') -stand_call_conf 20 -stand_emit_conf 20 --emitRefConfidence GVCF -GQB 5 -GQB 15 -GQB 20 -GQB 60 --variant_index_type LINEAR --variant_index_parameter 128000 --dbsnp /nfs/goldstein/goldsteinlab/software/GATK_bundle_2.8_b37/dbsnp_138.b37.vcf.gz -nct 1 --bamOutput $(dirname $IN)/BAMS/$(echo $i|awk   '{print $4"."$6"."$7".bam"}') -forceActive -disableOptimizations
       if [ $(samtools view $(dirname $IN)/BAMS/$(echo $i|awk '{print $4"."$6"."$7".bam"}')|wc -l) -eq 0 ];then samtools view -b $MOM_LOC $VAR_CHR:$((VAR_LOC-100))-$((VAR_LOC+100)) > $(dirname $IN)/BAMS/$(echo $i|awk   '{print $4"."$6"."$7".bam"}');samtools index $(dirname $IN)/BAMS/$(echo $i|awk '{print $4"."$6"."$7".bam"}');fi
       #samtools index $(dirname $IN)/BAMS/$(echo $i|awk '{print $4"."$6"."$7".bam"}')
        done < <(sed 's/ /_/g' $IN|sed 's/-/\t/'|sed 's/-/\t/'|sed '1d'|sort -k6,6n) #> $(dirname $IN)/$(basename ${IN%.*}).info.txt


#sed 's/\/nfs/\\\\10.73.50.80/' $(dirname $IN)/$(basename ${IN%.*}).bamloc |sed 's/\//\\/g' > $(dirname $IN)/$(basename ${IN%.*}).bamwinloc
find $(pwds)/$(dirname $IN)/BAMS -name '*.bam'|sed 's/\.\///'|sed 's/home\/[^\/]*/\/10.73.50.80\/homes/'|sed 's/\/nfs/\\\\10.73.50.80/' |sed 's/homes\/svaprojects/homes/'| sed 's/\//\\/g'|sed 's/\\\\ALIGNMENT/\\ALIGNMENT/g' > $(dirname $IN)/$(basename ${IN%.*}).bamwinloc
#sed 's/\/nfs/\\\\igm-avere.igm.cumc.columbia.edu/' $(dirname $IN)/$(basename ${IN%.*}).bamloc |sed 's/homes\/svaprojects/homes/' | sed 's/\//\\/g' > $(dirname $IN)/$(basename ${IN%.*}).bamwinloc

while read i
    do PRO_LOC=$(grep -w $(echo $i|awk '{print $2"."$6"."$7}') $(dirname $IN)/$(basename ${IN%.*}).bamwinloc|head -n1)
       if [ "$PRO_LOC" = "" ] ;then PRO_LOC="NA";fi
       DAD_LOC=$(grep -w $(echo $i|awk '{print $3"."$6"."$7}') $(dirname $IN)/$(basename ${IN%.*}).bamwinloc|head -n1)
       if [ "$DAD_LOC" = "" ] ;then DAD_LOC="NA";fi
       MOM_LOC=$(grep -w $(echo $i|awk '{print $4"."$6"."$7}') $(dirname $IN)/$(basename ${IN%.*}).bamwinloc|head -n1)
       if [ "$MOM_LOC" = "" ] ;then MOM_LOC="NA";fi
       VAR_CHR=$(echo $i|awk '{print $6}')
       VAR_LOC=$(echo $i|awk '{print $7}')
       VAR_GENE=$(echo $i|awk '{print $9}')
       PRO=$(echo $i|awk '{print $2}')
       echo "$PRO_LOC $DAD_LOC $MOM_LOC $VAR_CHR $VAR_LOC $((VAR_LOC-40)) $((VAR_LOC+40)) $PRO $VAR_GENE";done < <(sed 's/ /_/g' $IN|sed 's/-/\t/'|sed 's/-/\t/'|sed '1d'|sort -k6,6n) > $(dirname $IN)/$(basename ${IN%.*}).info.txt

       dir=$(pwd|sed 's/\/nfs/\\\\10.73.50.80/'|sed 's/home\/[^\/]*/10.73.50.80\/homes/'|sed 's/homes\/goldsteinlab/goldsteinlab/'|sed 's/homes\/svaprojects/svaprojects\/ns3116/'|sed 's/\//\\\\/g')
       awk -v dir="$dir" '{print "#"$8" "$4":"$5"\n""new\ngenome 1kg_v37\nload " $1"\nload " $2"\nload "$3"\nsnapshotDirectory \\"dir"\\IGV\ngoto "$4":"$6"-"$7"\nsort position\ncolor read strand\nsnapshot "$9"."$8"."$4"-"$5".png\ncollapse\nsnapshot "$9"."$8"."$4"-"$5".collapsed.png\n"}' <(sort -k9,9 $(dirname $IN)/$(basename ${IN%.*}).info.txt) >$(dirname $IN)/$(basename ${IN%.*}).batch
