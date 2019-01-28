#!/nfs/goldstein/software/R-3.0.1/bin/Rscript
library(rtf)

writeDNM <- function(dnm,rtf){
    #Change different colname in list-var-geno
    colnames(dnm) <- gsub("OMIM.Disease","OMIM.Disease",colnames(dnm))
    colnames(dnm) <- gsub("MGI_desc","MGI.Essential",colnames(dnm))
    colnames(dnm) <- gsub("HGMD.Class","HGMD.Class",colnames(dnm))
    colnames(dnm) <- gsub("Het.Binomial.P..child.","Het.Binomial.P",colnames(dnm))
    for(i in 1:dim(dnm)[1]){
        setFontSize(rtf,12)
        gene = gsub("'","",dnm[i,]$Gene.Name)
        addHeader(rtf,gene,subtitle=dnm[i,]$Variant.ID,font.size=14)
        startParagraph(rtf)
        addText(rtf,paste0("This is a ",dnm[i,]$Denovo.Flag," ",dnm[i,]$GT," ",gsub("_"," ",tolower(dnm[i,]$Effect))," variant in ",gene,". "))
        if((dnm[i,]$Ctrl.AF == 0 | is.na(dnm[i,]$Ctrl.AF)) & (dnm[i,]$Evs.All.Maf == 0 | is.na(dnm[i,]$Evs.All.Maf)) & (dnm[i,]$ExAC.global.af == 0 | is.na(dnm[i,]$ExAC.global.af)) & (dnm[i,]$gnomAD.Exome.global_AF == 0 | is.na(dnm[i,]$gnomAD.Exome.global_AF))){
            addText(rtf,"This variant is absent from internal and external control samples. ")}
        else{addText(rtf,paste0("This variant has a control AF of ",dnm[i,]$Ctrl.AF," in IGM controls, ",dnm[i,]$Evs.All.Maf," in EVS, ",dnm[i,]$ExAC.global.af," in ExAC, and ",dnm[i,]$gnomAD.Exome.global_AF," in gnomAD. "))}
        if(dnm[i,]$Effect == "missense_variant"){addText(rtf,paste0("It is a ",gsub("_"," ",dnm[i,]$Polyphen.Humvar.Prediction)," missense variant with a PolyPhen2 score of ",dnm[i,]$Polyphen.Humvar.Score,". "))}
        adj=" ";if(!is.na(dnm[i,]$X0.05._anypopn_RVIS.tile.ExAC.) & as.vector(dnm[i,]$X0.05._anypopn_RVIS.tile.ExAC.) <= 25){adj="n in"}
        addText(rtf,paste0(gene," is a",adj,"tolerant gene with an RVIS score of ",signif(as.double(as.vector(dnm[i,]$X0.05._anypopn_RVIS.tile.ExAC.)),digits=2),". "))
        if(!is.na(dnm[i,]$LoF.pLI.ExAC.) & as.vector(dnm[i,]$LoF.pLI.ExAC.) >= .9){addText(rtf,paste0(gene," is LoF intolerant gene with a PLI score of ",signif(as.double(as.vector(dnm[i,]$LoF.pLI.ExAC.)),digits=2),". "))}
        if(is.na(dnm[i,]$Gerp.RS.Score)){adj=" not"}else{adj=" very strongly";if(as.vector(dnm[i,]$Gerp.RS.Score)<5){adj=" strongly"};if(as.vector(dnm[i,]$Gerp.RS.Score)<4){adj=""};if(as.vector(dnm[i,]$Gerp.RS.Score)<2){adj=" weakly"};if(as.vector(dnm[i,]$Gerp.RS.Score)<0){adj=" not"}}
        addText(rtf,paste0("This site is",adj," conserved with a GERP++ RS score of ",dnm[i,]$Gerp.RS.Score,". "))
        if(!is.na(dnm[i,]$OMIM.Disease)){
            adj="";if(!is.na(dnm[i,]$MGI.Essential) & dnm[i,]$MGI.Essential == 1){adj=", and an essential gene"}
           addText(rtf,paste0(gene," is an OMIM disease gene associated with ",gsub(" \\|",", and", dnm[i,]$OMIM.Disease),adj,". "))}
        else if(!is.na(dnm[i,]$MGI.Essential) & dnm[i,]$MGI.Essential == 1){addText(rtf,paste0(gene," is an essential gene. "))}
        if(!is.na(dnm[i,]$ClinVar.ClinSig)){
           addText(rtf,paste0("This variant is listed as ",tolower(dnm[i,]$ClinVar.ClinSig)," in ClinVar for ",dnm[i,]$ClinVar.Disease,". "))}
        if(!is.na(dnm[i,]$HGMD.Class)){
            addText(rtf,paste0("This variant is listed as ",dnm[i,]$HGMD.Class," in HGMD for ",dnm[i,]$HGMD.Disease,". "))}
        if(grepl("0", dnm[i,]$X.Transcripts)){
            addText(rtf,"This is a non-CCDS variant.  ")}
        if(!is.na(dnm[i,]$Exon) & eval(parse(text=as.vector(dnm[i,]$Exon))) == 1){
            addText(rtf,"This variant is in the last exon of the transcript.  ")}
        if(is.na(dnm[i,]$ExAC.Sample.Covered.10x) | dnm[i,]$ExAC.Sample.Covered.10x < 10000){
            addText(rtf,paste0("This variant is only covered in ",dnm[i,]$ExAC.Sample.Covered.10x," samples in ExAC.  "))}
        if(dnm[i,]$Percent.Alt.Read.Binomial.P < 0.05){
            addText(rtf,paste0("This variant appears to mosaic with a binomial p-value of ",dnm[i,]$Percent.Alt.Read.Binomial.P,".  "))}
        addText(rtf,paste0("We have cases with ",dnm[i,]$All.dnm," de novo, ",dnm[i,]$All.hem," newly hemizygous, ",dnm[i,]$All.hom," newly homozygous, and ",dnm[i,]$All.chet+dnm[i,]$All.pchet," compouned het, tier 1 variants. "))
        addText(rtf,"\n")
        endParagraph(rtf)
    }
}  

writeCHET <- function(chet,rtf){
    #colnames(chet) <- gsub("OMIM.Disease","OMIM.Disease",colnames(chet))
    #colnames(chet) <- gsub("MGI_desc","MGI.Essential",colnames(chet))
    #colnames(chet) <- gsub("HGMD.Class","HGMD.Class",colnames(chet))
    for(i in 1:dim(chet)[1]){
        setFontSize(rtf,12)
        gene = gsub("'","",chet[i,]$Gene.Name.1)
        addHeader(rtf,gene,paste0(subtitle=chet[i,]$Variant.ID.1," , ",subtitle=chet[i,]$Variant.ID.2),font.size=14)
        startParagraph(rtf)
        addText(rtf,paste0("These are ",chet[i,]$Comp.Het.Flag," variants in ",gene,". The first is a ",gsub("_"," ",tolower(chet[i,]$Effect.1))," and the second is a ",gsub("_"," ",tolower(chet[i,]$Effect.2)), " variant. "))
        adj=" ";if(!is.na(chet[i,]$X0.05._anypopn_RVIS.tile.ExAC..1) & as.vector(chet[i,]$X0.05._anypopn_RVIS.tile.ExAC..1) <= 25){adj="n in"}
        addText(rtf,paste0(gene," is a",adj,"tolerant gene with an RVIS score of ",signif(as.double(as.vector(chet[i,]$X0.05._anypopn_RVIS.tile.ExAC..1)),digits=2),". "))
        if(!is.na(chet[i,]$LoF.pLI.ExAC..1) & as.vector(chet[i,]$LoF.pLI.ExAC..1) >= .9){addText(rtf,paste0(gene," is LoF intolerant gene with a PLI score of ",signif(as.double(as.vector(chet[i,]$LoF.pLI.ExAC..1)),digits=2),". "))}
        if(!is.na(chet[i,]$OMIM.Disease.1)){
           adj="";if(!is.na(chet[i,]$MGI.Essential.1) & chet[i,]$MGI.Essential.1 == 1){adj=", and an essential gene"}
           addText(rtf,paste0(gene," is an OMIM disease gene associated with ",gsub(" \\|",", and", chet[i,]$OMIM.Disease.1),adj,". "))}
        else if(!is.na(chet[i,]$MGI.Essential.1) & chet[i,]$MGI.Essential.1 == 1){addText(rtf,paste0(gene," is an essential gene. "))}
        addText(rtf,paste0("We have cases with ",chet[i,]$All.dnm," de novo, ",chet[i,]$All.hem," newly hemizygous, ",chet[i,]$All.hom," newly homozygous, and ",chet[i,]$All.chet+dnm[i,]$All.pchet," compouned het, tier 1 variants. "))
        addText(rtf,"\n")
        
#First variant
        if((chet[i,]$Ctrl.AF.1 == 0 | is.na(chet[i,]$Ctrl.AF.1)) & (chet[i,]$Evs.All.Maf.1 == 0 | is.na(chet[i,]$Evs.All.Maf.1)) & (chet[i,]$ExAC.global.af.1 == 0 | is.na(chet[i,]$ExAC.global.af.1)) & (chet[i,]$gnomAD.Exome.global_AF.1 == 0 | is.na(chet[i,]$gnomAD.Exome.global_AF.1))){
            addText(rtf,"The first variant is absent from internal and external control samples. ")}
        else{addText(rtf,paste0("The first variant has a control AF of ",chet[i,]$Ctrl.AF.1," in IGM controls, ",chet[i,]$Evs.All.Maf.1," in EVS, ",chet[i,]$ExAC.global.af.1," in ExAC, "),chet[i,]$gnomAD.Exome.global_AF.1," in gnomAD. ")}
        if(chet[i,]$Effect.1 == "missense_variant"){addText(rtf,paste0("It is a ",gsub("_"," ",chet[i,]$Polyphen.Humvar.Prediction.1)," missense variant with a PolyPhen2 score of ",chet[i,]$Polyphen.Humvar.Score.1,". "))}
        if(is.na(chet[i,]$Gerp.RS.Score.1)){adj=" not"}else{adj=" very strongly";if(as.vector(chet[i,]$Gerp.RS.Score.1)<5){adj=" strongly"};if(as.vector(chet[i,]$Gerp.RS.Score.1)<4){adj=""};if(as.vector(chet[i,]$Gerp.RS.Score.1)<2){adj=" weakly"};if(as.vector(chet[i,]$Gerp.RS.Score.1)<0){adj=" not"}}
        addText(rtf,paste0("This site is",adj," conserved with a GERP++ RS score of ",chet[i,]$Gerp.RS.Score.1,". "))
        if(!is.na(chet[i,]$ClinVar.ClinSig.1)){
           addText(rtf,paste0("This variant is listed as ",tolower(chet[i,]$ClinVar.ClinSig.1)," in ClinVar for",chet[i,]$ClinVar.Disease.1,". "))}
        if(!is.na(chet[i,]$HGMD.Class.1)){
            addText(rtf,paste0("This variant is listed as ",chet[i,]$HGMD.Class.1," in HGMD for ",chet[i,]$HGMD.Disease.1,". "))}
        
        if(grepl("0", chet[i,]$X.Transcripts.1)){
            addText(rtf,"This is a non-CCDS variant.  ")}
        if(!is.na(chet[i,]$Exon.1) & eval(parse(text=as.vector(chet[i,]$Exon.1))) == 1){
            addText(rtf,"This variant is in the last exon of the transcript.  ")}
        if(is.na(chet[i,]$ExAC.Sample.Covered.10x.1) | chet[i,]$ExAC.Sample.Covered.10x.1 < 10000){
            addText(rtf,paste0("This variant is only covered in ",dnm[i,]$ExAC.Sample.Covered.10x," samples in ExAC.  "))}
        if(chet[i,]$Percent.Alt.Read.Binomial.P.1 < 0.05){
            addText(rtf,paste0("This variant appears to mosaic with a binomial p-value of ",chet[i,]$Percent.Alt.Read.Binomial.P.1,".  "))} 
        if(!is.null(chet[i,]$gnomAD.Exome.filter.1)){
            addText(rtf,paste0("This variant is marked as a ",chet[i,]$gnomAD.Exome.filter.1," in gnomAD. "))}
        #second variant
        addText(rtf,"\n")
        if((chet[i,]$Ctrl.AF.2 == 0 | is.na(chet[i,]$Ctrl.AF.2)) & (chet[i,]$Evs.All.Maf.2 == 0 | is.na(chet[i,]$Evs.All.Maf.2)) & (chet[i,]$ExAC.global.af.2 == 0 | is.na(chet[i,]$ExAC.global.af.2)) & (chet[i,]$gnomAD.Exome.global_AF.2 == 0 | is.na(chet[i,]$gnomAD.Exome.global_AF.2))){
            addText(rtf,"The second variant is absent from internal and external control samples. ")}
        else{addText(rtf,paste0("The second variant has a control AF of ",chet[i,]$Ctrl.AF.2," in IGM controls, ",chet[i,]$Evs.All.Maf.2," in EVS, ",chet[i,]$ExAC.global.af.2," in ExAC, ",chet[i,]$gnomAD.Exome.global_AF.2," in gnomAD. "))}
        if(chet[i,]$Effect.2 == "missense_variant"){addText(rtf,paste0("It is a ",gsub("_"," ",chet[i,]$Polyphen.Humvar.Prediction.2)," missense variant with a PolyPhen2 score of ",chet[i,]$Polyphen.Humvar.Score.2,". "))}
        if(is.na(chet[i,]$Gerp.RS.Score.2)){adj=" not"}else{adj=" very strongly";if(as.vector(chet[i,]$Gerp.RS.Score.2)<5){adj=" strongly"};if(as.vector(chet[i,]$Gerp.RS.Score.2)<4){adj=""};if(as.vector(chet[i,]$Gerp.RS.Score.2)<2){adj=" weakly"};if(as.vector(chet[i,]$Gerp.RS.Score.2)<0){adj=" not"}}
        addText(rtf,paste0("This site is",adj," conserved with a GERP++ RS score of ",chet[i,]$Gerp.RS.Score.2,". "))
        if(!is.na(chet[i,]$ClinVar.ClinSig.2)){
           addText(rtf,paste0("This variant is listed as ",tolower(chet[i,]$ClinVar.ClinSig.2)," in ClinVar for ",chet[i,]$ClinVar.Disease.2,". "))}
        if(!is.na(chet[i,]$HGMD.Class.2)){
            addText(rtf,paste0("This variant is listed as ",chet[i,]$HGMD.Class.2," in HGMD for ",chet[i,]$HGMD.Disease.2,". "))}
        
        if(grepl("0", chet[i,]$X.Transcripts.2)){
            addText(rtf,"This is a non-CCDS variant.  ")}
        if(!is.na(chet[i,]$Exon.2) & eval(parse(text=as.vector(chet[i,]$Exon.2))) == 1){
            addText(rtf,"This variant is in the last exon of the transcript.  ")}
        if(is.na(chet[i,]$ExAC.Sample.Covered.10x.2) | chet[i,]$ExAC.Sample.Covered.10x.2 < 10000){
            addText(rtf,paste0("This variant is only covered in ",dnm[i,]$ExAC.Sample.Covered.10x," samples in ExAC.  "))}
        if(chet[i,]$Percent.Alt.Read.Binomial.P.2 < 0.05){
            addText(rtf,paste0("This variant appears to mosaic with a binomial p-value of ",chet[i,]$Percent.Alt.Read.Binomial.P.2,".  "))} 
        if(!is.null(chet[i,]$gnomAD.Exome.filter.2)){
            addText(rtf,paste0("This variant is marked as a ",chet[i,]$gnomAD.Exome.filter.2," in gnomAD. "))}
        addText(rtf,"\n")
        addText(rtf,"\n")
        endParagraph(rtf)
    }
}   

writeSummary <- function(dnm,hom,hem,chet,tier2,dir){
    dir.create(file.path(dir,"Sample_Summaries"),showWarnings=F)
    allChilds <- c(as.vector(dnm$Sample.Name),as.vector(hom$Sample.Name),as.vector(hem$Sample.Name),as.vector(chet$Sample.Name.1),as.vector(tier2$Sample.Name))
    if(length(allChilds) == 0){return}
    out <- file.path(dir,"Sample_Summaries",paste0(allChilds[1],".doc"))
    rtf <- RTF(out,width=8.5,height=11,omi=c(1,1,1,1),font.size=18)
    addHeader(rtf,allChilds[1],font.size=18)
    addHeader(rtf,"Tier 1",font.size=18)
    if(dim(dnm)[1] > 0){writeDNM(dnm,rtf)}
    if(dim(hem)[1] > 0){writeDNM(hem,rtf)}
    if(dim(hom)[1] > 0){writeDNM(hom,rtf)}
    if(dim(chet)[1] > 0){writeCHET(chet,rtf)}
    addHeader(rtf,"Tier 2",font.size=18)
    if(dim(tier2)[1] > 0){writeDNM(tier2,rtf)}
    done(rtf)
    return
}

writeNonTrioSummary <- function(samp.kv,samp.kv5 = samp.kv[FALSE,],samp.pdnm = samp.kv[FALSE,],samp.prec = samp.kv[FALSE,],samp.pchet = samp.kv[FALSE,],samp.lofd = samp.kv[FALSE,],samp.CVExact = samp.kv[FALSE,],dir=dir){
    dir.create(file.path(dir,"Sample_Summaries"),showWarnings=F)
    allSamps <- c(as.vector(samp.kv$Sample.Name),as.vector(samp.kv5$Sample.Name),as.vector(samp.pdnm$Sample.Name),as.vector(samp.prec$Sample.Name),as.vector(samp.pchet$Sample.Name))
    if(length(allSamps) == 0){return}
    out <- file.path(dir,"Sample_Summaries",paste0(allSamps[1],".nonTrio.doc"))
    rtf <- RTF(out,width=8.5,height=11,omi=c(1,1,1,1),font.size=18)
    addHeader(rtf,allSamps[1],font.size=18)
    if(dim(samp.kv)[1] > 0){writeDNM(samp.kv,rtf)}
    if(dim(samp.pdnm)[1] > 0){writeDNM(samp.pdnm,rtf)}
    if(dim(samp.kv5)[1] > 0){writeDNM(samp.kv5,rtf)}
    if(dim(samp.prec)[1] > 0){writeDNM(samp.prec,rtf)}
    if(dim(samp.pchet)[1] > 0){writeDNM(samp.pchet,rtf)}
    if(dim(samp.lofd)[1] > 0){writeDNM(samp.lofd,rtf)}
    if(dim(samp.CVExact)[1] > 0){writeDNM(samp.CVExact,rtf)}
    done(rtf)
    return
}
