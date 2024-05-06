#!/usr/bin/env nextflow

// git 8.1
// Create channel of linkage groups
Channel
    .from( ('01'..'09') + ('10'..'19') + ('20'..'24') )
    .map{ "LG" + it }
    .into{ lg_ch1; lg_ch2; lg_ch3 }

// git 8.2
// Open phased genotype data
Channel
    .fromFilePairs("/nfs/data/haex1482/shared/demo/vcf/phylo2e_phased2.vcf.{gz,gz.tbi}")
    .set{ vcf_msmc }

// git 8.3
// Open unphased genotype data to extract depth information
Channel
    .fromFilePairs("/nfs/data/haex1482/shared/demo/vcf/phylo2e_snpsfilt.vcf.{gz,gz.tbi}")
    .set{ vcf_depth }

// git 8.4
// Gather depth per individual
process gather_depth {
    label 'L_20g2h_split_by_sample'
    publishDir "depth", mode: 'copy'

    input:
    set vcfID, file( vcf ) from vcf_depth

    output:
    file( "depth_by_sample.txt" ) into depth_ch

    script:
    """
    ml VCFtools/0.1.16-GCC-13.1.0

    zcat < ${vcf[0]} |
       vcf-query -l |
       grep "tor\\|tab\\|tig" \
       > pop.txt

    vcftools \
       --gzvcf ${vcf[0]} \
       --remove pop.txt \
       --depth \
       --stdout > depth_by_sample.txt
   """
}

// git 8.5
// Create channel out of sequencing depth table
depth_ch
    .splitCsv(header:true, sep:"\t")
    .map{ row -> [ id:row.INDV, sites:row.N_SITES, depth:row.MEAN_DEPTH] }
    .map{ [it.id, it.sites, it.depth] }
    .set{ depth_by_sample_ch }

// git 8.6
// Create channel from bam files and add sample id
Channel
    .fromPath("/nfs/data/haex1482/shared/demo/bam/*_dedup.bam")
    .map{ file ->
        def key = file.name.toString()
        return tuple(key - "_dedup.bam", file) }
    .set{ sample_bams }

// git 8.7
// Combine sample bams and sequencing depth
sample_bams
    .join( depth_by_sample_ch )
    .set{ sample_bam_and_depth }

// git 8.8
// Multiply the sample channel by the linkage groups
sample_bam_and_depth
    .combine( vcf_msmc )
    .combine( lg_ch1 )
    .set{ samples_msmc }

// git 8.9
// Split vcf by individual
process split_vcf_by_individual {
    label 'L_20g15m_split_by_vcf'

    input:
    set val( id ), file( bam ), val( sites ), val( depth ), val( vcf_id ), file( vcf ), val( lg ) from samples_msmc

    output:
    set val( id ), val( lg ), file( bam ), val( depth ), file( "phased_mac2.${id}.${lg}.vcf.gz" ) into ( sample_vcf, sample_vcf2 )

    script:
    """
    ml GATK/4.4.0.0-GCCcore-13.1.0-Java-17

    ref="/nfs/data/haex1482/shared/demo/ref/Hpue_genome_unmasked_01.fasta"

    gatk --java-options "-Xmx10G" \
       SelectVariants \
       -R \${ref} \
       -V ${vcf[0]} \
       -sn ${id} \
       -L ${lg}\
       -O phased_mac2.${id}.${lg}.vcf.gz
    """
}

// git 8.10
// Create coverage mask from original mapped sequences
process bam_caller {
    label 'L_36g47h_bam_caller'
    publishDir "masks", mode: 'copy' , pattern: "*.coverage_mask.bed.gz"
    // conda "$HOME/miniconda2/envs/py3"

    input:
    set val( id ), val( lg ), file( bam ), val( depth ), file( vcf ) from sample_vcf

    output:
    set val( id ), val( lg ), file( "*.bam_caller.vcf.gz" ), file( "*.coverage_mask.bed.gz" ) into coverage_by_sample_lg

    script:
    """
    ml BCFtools/1.18-GCC-13.1.0
    ml SAMtools/1.18-GCC-13.1.0
    ml OpenSSL/1.1

    ref="/nfs/data/haex1482/shared/demo/ref/Hpue_genome_unmasked_01.fasta"

    samtools index ${bam}

    bcftools mpileup -q 25 -Q 20 -C 50 -r ${lg} -f \${ref} ${bam} | \
       bcftools call -c -V indels | \
       bamHamletCaller.py ${depth} ${id}.${lg}.coverage_mask.bed.gz | \
       gzip -c > ${id}.${lg}.bam_caller.vcf.gz
   """
}

// git 8.11
// Create segsites file
process generate_segsites {
    label 'L_20g15m_msmc_generate_segsites'
    publishDir "segsites", mode: 'copy' , pattern: "*.segsites.vcf.gz"

    input:
    set val( id ), val( lg ), file( bam ), val( depth ), file( vcf ) from sample_vcf2

    output:
    set val( id ), val( lg ), file( "*.segsites.vcf.gz" ), file( "*.covered_sites.bed.txt.gz" ) into segsites_by_sample_lg

    script:
    """
    zcat ${vcf} | \
       python2 ~/apps/msmc-tools/vcfAllSiteParser.py ${id} ${id}.${lg}.covered_sites.bed.txt.gz | \
       gzip -c > ${id}.${lg}.segsites.vcf.gz
    """
}

// git 8.12
// Assign samples randomly across MSMC and cross coalescence runs
process msmc_sample_grouping {
    label "L_loc_msmc_grouping"
    publishDir "setup", mode: 'copy'

    output:
    file( "msmc_grouping_phylo2e-n3.txt" ) into msmc_grouping
    file( "msmc_grouping_phylo2e_cc-m2.txt" ) into cc_grouping

    script:
    """
    ml R/4.3.1-foss-2023a

    base="\$WORK/demo/msmc"

    Rscript --vanilla \$base/R/sample_assignment_msmc_mh.R \
           \$base/R/distribute_samples_msmc_and_cc.R \
           \$base/R/cross_cc.R \
           \$base/R/sample_info_phylo2e.txt \
           msmc
   """
}

// git 8.13
// Read grouping into a channel
msmc_grouping
    .splitCsv(header:true, sep:"\t")
    .map{ row -> [ run:row.msmc_run, spec:row.spec, geo:row.geo, group_nr:row.group_nr, group_size:row.group_size, samples:row.samples ] }
    .set{ msmc_runs }

// git 8.14
// Wait for bam_caller and generate_segsites to finish
coverage_by_sample_lg.collect().map{ "coverage done!" }.into{ coverage_done; coverage_cc }
segsites_by_sample_lg.collect().map{ "segsites done!" }.into{ segsites_done; segsites_cc }

// git 8.15
// Attach masks to MSMC group assignment
lg_ch2
    .combine( msmc_runs )
    .combine( coverage_done )
    .combine( segsites_done )
    .map{ [it[0], it[1].run, it[1]] }
    .set{ msmc_grouping_after_segsites }

// git 8.16
// Generate MSMC input files (4 or 3 inds per species)
process generate_multihetsep {
    label "L_120g40h_msmc_generate_multihetsep"
    publishDir "input/run_${run}", mode: 'copy' , pattern: "*.multihetsep.txt"

    input:
    /* content msmc_gr: val( msmc_run ), val( spec ), val( geo ), val( group_nr ), val( group_size ), val( samples ) */
    /* [LG20, [msmc_run:45, spec:uni, geo:pan, group_nr:4, group_size:3, samples:ind1, ind2, ind3], coverage done!, segsites done!]*/
    set val( lg ), val( run ), msmc_gr from msmc_grouping_after_segsites

    output:
    set val( run ), val( lg ), val( msmc_gr.spec ), val( msmc_gr.geo ), val( msmc_gr.group_size ), file( "msmc_run.*.multihetsep.txt" ) into msmc_input_lg

    script:
    """
    base="\$WORK/demo/msmc"

    covdir="\$base/masks/"
    smp=\$(echo ${msmc_gr.samples} | \
        sed "s|, |\\n--mask=\${covdir}|g; s|^|--mask=\${covdir}|g" | \
        sed "s/\$/.${lg}.coverage_mask.bed.gz/g" | \
        echo \$( cat ) )

    segdir="\$base/segsites/"
    seg=\$(echo ${msmc_gr.samples} | \
        sed "s|, |\\n\${segdir}|g; s|^|\${segdir}|g" | \
        sed "s/\$/.${lg}.segsites.vcf.gz/g" | \
        echo \$( cat ) )

    generate_multihetsep.py \
        \$smp \
        --mask=/nfs/data/haex1482/shared/demo/map/Hpue_${lg}.mapmask.bed.txt.gz \
        \$seg > msmc_run.${msmc_gr.run}.${msmc_gr.spec}.${msmc_gr.geo}.${lg}.multihetsep.txt
   """
}

// git 8.17
// Collect all linkage groups for each run
msmc_input_lg
    .groupTuple()
    .set { msmc_input }

// git 8.18
// Run MSMC analysis
process msmc_run {
    label "L_190g100h_msmc_run"
    publishDir "output", mode: 'copy' , pattern: "*.final.txt"
    publishDir "loops", mode: 'copy' , pattern: "*.loop.txt"

    input:
    set msmc_run, lg , spec, geo, group_size, file( hetsep ) from msmc_input

    output:
    file("*.msmc2.*.txt") into msmc_output

    script:
    """
    nhap=\$(echo \$(seq 0 \$((${group_size[0]}*2-1))) | sed 's/ /,/g')
    infiles=\$( echo ${hetsep} )

    msmc2-linux \
        -t 32 \
        -m 0.00255 \
        -p 1*2+25*1+1*2+1*3 \
        -o run${msmc_run}.${spec[0]}.${geo[0]}.msmc2 \
        -I \${nhap} \
        \${infiles}
   """
}

// git 8.19
// Generate MSMC cross coalescence input files (2 inds x 2 species)
cc_grouping
    .splitCsv(header:true, sep:"\t")
    .map{ row -> [ run:row.run_nr, geo:row.geo, spec_1:row.spec_1, spec_2:row.spec_2, contrast_nr:row.contrast_nr, samples_1:row.samples_1, samples_2:row.samples_2 ] }
    .set { cc_runs }

// git 8.20
// Attach masks to cross coalescence group assignment
lg_ch3
    .combine( cc_runs )
    .combine( coverage_cc )
    .combine( segsites_cc )
    .map{[it[0], it[1].run, it[1]]}
    .set{ cc_grouping_after_segsites }

// git 8.21
// Create multihetsep files (combination of all 4 individuals)
process generate_multihetsep_cc {
    label "L_105g30h_cc_generate_multihetsep"
    publishDir "input_cc/run_cc_${run}", mode: 'copy' , pattern: "*.multihetsep.txt"

    input:
    /* content cc_gr: val( run_nr ), val( geo ), val( spec_1 ), val( spec_2 ), val( contrast_nr ), val( samples_1 ), val( samples_2 ) */
    set val( lg ), val( run ), cc_gr from cc_grouping_after_segsites

    output:
    set val( cc_gr.run ), val( lg ), val( cc_gr.spec_1 ), val( cc_gr.spec_2 ), val( cc_gr.geo ), val( cc_gr.contrast_nr ), val( cc_gr.samples_1 ), val( cc_gr.samples_2 ), file( "cc_run.*.multihetsep.txt" ) into cc_input_lg

    script:
    """
    base="\$WORK/demo/msmc"
    
    covdir="\$base/masks/"
    smp1=\$(echo ${cc_gr.samples_1}  | \
        sed "s|, |\\n--mask=\${covdir}|g; s|^|--mask=\${covdir}|g" | \
        sed "s/\$/.${lg}.coverage_mask.bed.gz/g" | \
        echo \$( cat ) )
    smp2=\$(echo ${cc_gr.samples_2}  | \
        sed "s|, |\\n--mask=\${covdir}|g; s|^|--mask=\${covdir}|g" | \
        sed "s/\$/.${lg}.coverage_mask.bed.gz/g" | \
        echo \$( cat ) )

    segdir="\$base/segsites/"
    seg1=\$(echo ${cc_gr.samples_1}  | \
        sed "s|, |\\n\${segdir}|g; s|^|\${segdir}|g" | \
        sed "s/\$/.${lg}.segsites.vcf.gz/g" | \
        echo \$( cat ) )
    seg2=\$(echo ${cc_gr.samples_2}  | \
        sed "s|, |\\n\${segdir}|g; s|^|\${segdir}|g" | \
        sed "s/\$/.${lg}.segsites.vcf.gz/g" | \
        echo \$( cat ) )

    generate_multihetsep.py \
        \${smp1} \
        \${smp2} \
         --mask=/nfs/data/haex1482/shared/demo/map/Hpue_${lg}.mapmask.bed.txt.gz \
        \${seg1} \
        \${seg2} \
        > cc_run.${run}.${cc_gr.spec_1}-${cc_gr.spec_2}.${cc_gr.contrast_nr}.${cc_gr.geo}.${lg}.multihetsep.txt
   """
}

// git 8.22
// Collect all linkage groups for each run
cc_input_lg
    .groupTuple()
    .set{ cc_input }
    // cc_input.view()

// git 8.23
// Run cross coalescence analysis
process cc_run {
    label "L_190g10ht24_cc_run"
    publishDir "output_cc", mode: 'copy', pattern: "*.final.txt"
    tag "${cc_run}-${geo[0]}:${spec1[0]}/${spec2[0]}"

    input:
    set cc_run, lg, spec1, spec2, geo, contr_nr, samples_1, samples_2, file( hetsep ) from cc_input

    output:
    file("*.final.txt") into cc_output

    script:
    """
    infiles=\$( echo ${hetsep} )
    pop1=\$( echo "${samples_1}" | sed 's/\\[//g; s/, /,/g; s/\\]//g' )
    pop2=\$( echo "${samples_2}" | sed 's/\\[//g; s/, /,/g; s/\\]//g' )

    msmc2-linux \
        -t 32 \
        -m 0.00255 \
        -p 1*2+25*1+1*2+1*3 \
        -o cc_run.${cc_run}.${spec1[0]}.cc \
        -I 0,1,2,3 \
        \${infiles}

    msmc2-linux \
        -t 32 \
        -m 0.00255 \
        -p 1*2+25*1+1*2+1*3 \
        -o cc_run.${cc_run}.${spec2[0]}.cc \
        -I 4,5,6,7 \
        \${infiles}

    msmc2-linux \
        -t 32 \
        -m 0.00255 \
        -p 1*2+25*1+1*2+1*3 \
        -o cc_run.${cc_run}.cross.cc \
        -I 0,1,2,3,4,5,6,7 \
        \${infiles}

    combineCrossCoal.py \
        cc_run.${cc_run}.cross.msmc.final.txt \
        cc_run.${cc_run}.${spec1[0]}.msmc.final.txt \
        cc_run.${cc_run}.${spec2[0]}.msmc.final.txt \
    > cc_run.${cc_run}.combined.final.txt
    """
}
