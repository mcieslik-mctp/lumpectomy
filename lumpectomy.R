lumpy.vcf2tbl <- function(vcf) {
    ## INFO
    ## SVTYPE    1      String  Type of structural variant                         
    ## SVLEN     .      Integer Difference in length between REF and ALT alleles   
    ## END       1      Integer End position of the variant described in this re...
    ## STRANDS   .      String  Strand orientation of the adjacency in BEDPE for...
    ## IMPRECISE 0      Flag    Imprecise structural variation                     
    ## CIPOS     2      Integer Confidence interval around POS for imprecise var...
    ## CIEND     2      Integer Confidence interval around END for imprecise var...
    ## CIPOS95   2      Integer Confidence interval (95%) around POS for impreci...
    ## CIEND95   2      Integer Confidence interval (95%) around END for impreci...
    ## MATEID    .      String  ID of mate breakends                               
    ## EVENT     1      String  ID of event associated to breakend                 
    ## SECONDARY 0      Flag    Secondary breakend in a multi-line variants        
    ## SU        .      Integer Number of pieces of evidence supporting the vari...
    ## PE        .      Integer Number of paired-end reads supporting the varian...
    ## SR        .      Integer Number of split reads supporting the variant acr...
    ## BD        .      Integer Amount of BED evidence supporting the variant ac...
    ## EV        .      String  Type of LUMPY evidence contributing to the varia...
    ## PRPOS     .      String  LUMPY probability curve of the POS breakend        
    ## PREND     .      String  LUMPY probability curve of the END breakend        
    ## GENO
    ## GT  1      String  Genotype                                                 
    ## SU  1      Integer Number of pieces of evidence supporting the variant      
    ## PE  1      Integer Number of paired-end reads supporting the variant        
    ## SR  1      Integer Number of split reads supporting the variant             
    ## BD  1      Integer Amount of BED evidence supporting the variant            
    ## GQ  1      Integer Genotype quality                                         
    ## SQ  1      Float   Phred-scaled probability that this site is variant (no...
    ## GL  G      Float   Genotype Likelihood, log10-scaled likelihoods of the d...
    ## DP  1      Integer Read depth                                               
    ## RO  1      Integer Reference allele observation count, with partial obser...
    ## AO  A      Integer Alternate allele observations, with partial observatio...
    ## QR  1      Integer Sum of quality of reference observations                 
    ## QA  A      Integer Sum of quality of alternate observations                 
    ## RS  1      Integer Reference allele split-read observation count, with pa...
    ## AS  A      Integer Alternate allele split-read observation count, with pa...
    ## ASC A      Integer Alternate allele clipped-read observation count, with ...
    ## RP  1      Integer Reference allele paired-end observation count, with pa...
    ## AP  A      Integer Alternate allele paired-end observation count, with pa...
    ## AB  A      Float   Allele balance, fraction of observations from alternat...
    var.tbl <- data.table(
    ID=rownames(info(vcf)),
    CHR=as.character(seqnames(rowRanges(vcf))),
    BEG=start(rowRanges(vcf)),
    FILTER=rowRanges(vcf)$FILTER,
    as.data.table(lapply(info(vcf)[,c("END", "SVTYPE", "IMPRECISE", "EVENT", "SECONDARY")], unlist))
   ,
    as.data.table(lapply(c("SVLEN"="SVLEN"), function(xx) {
        unlist(ifelse(lengths(info(vcf)[,xx])==0, NA_integer_, as.list(info(vcf)[,xx])))
    }))
   ,
    as.data.table(lapply(c("MATEID"="MATEID"), function(xx) {
        unlist(ifelse(lengths(info(vcf)[,xx])==0, NA_character_, as.list(info(vcf)[,xx])))
    }))
   ,
    as.data.table(lapply(c("GT"="GT", "DP"="DP", "SU"="SU", "PE"="PE", "SR"="SR", "GQ"="GQ",
                           "SQ"="SQ", "RO"="RO", "AO"="AO", "QR"="QR",
                           "QA"="QA", "RS"="RS", "AS"="AS", "ASC"="ASC", "RP"="RP",
                           "AP"="AP", "AB"="AB"
                           ),
                         function(xx) unlist(geno(vcf)[[xx]][,1])))
    )
}

main <- function() {
    option_list = list(
    )
    parser = optparse::OptionParser(
      "Rscript lumpectomy.R' [options] vcf_file bed_file out_vcf_file out_csv_file",
      description=c("Extract SV calls.\n"),
      epilogue=c(
        "Written by Marcin Cieslik (mcieslik@med.umich.edu) ",
        "Michigan Center for Translational Pathology (c) 2017\n"),
      option_list=option_list
    )
    opt = optparse::parse_args(parser, positional_arguments=TRUE)
    ##
    if (length(opt$args) < 4) {
        optparse::print_help(parser)
        write("Some of vcf_file bed_file out_vcf_file out_csv_file are missing.\n", stderr())
        quit("no", 1)
    }
    vcf.fn <- opt$args[1]
    bed.fn <- opt$args[2]
    out.vcf.fn <- opt$args[3]
    out.csv.fn <- opt$args[4]
    if (!("VariantAnnotation" %in% installed.packages())) {
        write("Please install the Bioconductor VariantAnnotation package\n", stderr())
        quit("no", 1)
    }
    if (!("data.table" %in% installed.packages())) {
        write("Please install the data.table package\n", stderr())
        quit("no", 1)
    }
    
    suppressMessages(library(VariantAnnotation))
    suppressMessages(library(data.table))

    vcf <- readVcf(vcf.fn, "hg38")
    rng <- with(fread(bed.fn), GRanges(V1, IRanges(V2, V3)))
    ##
    tbl <- lumpy.vcf2tbl(vcf)
    tbl.loc <- tbl[SVTYPE %in% c("DUP", "DEL", "INV"), .(chr.5=CHR, pos.5=BEG, chr.3=CHR, pos.3=END, spn=SR, enc=PE,
                                                         type=SVTYPE, filter=FILTER, dp=DP, gt=GT, id=ID)]
    tmp.1 <- tbl[SVTYPE=="BND"][grep("_1",ID)]
    tmp.2 <- tbl[SVTYPE=="BND"][grep("_2",ID)]
    setkey(tmp.1, EVENT)
    setkey(tmp.2, EVENT)
    tbl.bnd <- tmp.1[tmp.2,.(chr.5=CHR, pos.5=BEG, chr.3=i.CHR, pos.3=i.BEG, spn=SR, enc=PE, type=SVTYPE, filter=FILTER, dp=DP, gt=GT, id=EVENT)]
    tbl.jnc <- rbind(tbl.loc, tbl.bnd)
    bpt.rng <- with(tbl.jnc, {
        bpt.5 <- GRanges(chr.5, IRanges(pos.5, pos.5))
        bpt.3 <- GRanges(chr.3, IRanges(pos.3, pos.3))
        (bpt.5 %over% rng) | (bpt.3 %over% rng)
    })
    tbl.rng <- tbl.jnc[bpt.rng]
    sel.bnd <- (info(vcf)$SVTYPE == "BND") & (info(vcf)$EVENT %in% tbl.rng$id)
    sel.ddi <- (info(vcf)$SVTYPE != "BND") & (rownames(info(vcf)) %in% tbl.rng$id)
    sel.all <- sel.bnd | sel.ddi
    vcf.rng <- vcf[sel.all]
    writeVcf(vcf.rng, out.vcf.fn)
    fwrite(tbl.rng, out.csv.fn)
}

main()
