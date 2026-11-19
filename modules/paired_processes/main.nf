/*
 * Every process that consumes per-sample files takes ONE tuple keyed by SAMPLE.
 * Never feed a process several independent queue channels (e.g. a sample tuple + bare
 * trimmed-read paths): Nextflow pairs those by arrival order, not by sample, which
 * silently mislabels outputs (see CHANGELOG 2026-November).
 */

process fastp_paired {
    container 'us-east1-docker.pkg.dev/compute-workspace/omics-docker-repo/rnaseq2'

    publishDir "${params.out_bucket}/${SAMPLE}/logs/", mode: 'copy', pattern: '*.json'

    input:
    tuple val(SAMPLE), path(R1), path(R2)

    output:
    path("${SAMPLE}_fastp.json"), emit: fastp_log
    tuple val(SAMPLE), path("synced_clt_${SAMPLE}_R1.fastq.gz"), path("synced_clt_${SAMPLE}_R2.fastq.gz"), emit: trimmed

    script:
    """
    fastp \
            --in1 ${R1} \
            --in2 ${R2} \
            --out1 synced_clt_${SAMPLE}_R1.fastq.gz \
            --out2 synced_clt_${SAMPLE}_R2.fastq.gz \
            --length_required 35 \
            --cut_front_window_size 1 \
            --cut_front_mean_quality 13 \
            --cut_front \
            --cut_tail_window_size 1 \
            --cut_tail_mean_quality 13 \
            --cut_tail \
            -w 2 -y -r \
            -j ${SAMPLE}_fastp.json \
            --correction
    """

    stub:
    """
    sleep \$(( RANDOM % 5 ))
    echo "reads_of=\$(zcat -f ${R1} | head -1)" | gzip > synced_clt_${SAMPLE}_R1.fastq.gz
    echo "reads_of=\$(zcat -f ${R2} | head -1)" | gzip > synced_clt_${SAMPLE}_R2.fastq.gz
    echo '{}' > ${SAMPLE}_fastp.json
    """
}


process fastqc_paired {
    container 'us-east1-docker.pkg.dev/compute-workspace/omics-docker-repo/rnaseq2'

    publishDir "${params.out_bucket}/${SAMPLE}/logs/", mode: 'copy'

    input:
    tuple val(SAMPLE), path(R1), path(R2), path(trim_R1), path(trim_R2)

    output:
    path("raw_fastqc_${SAMPLE}")
    path("clean_fastqc_${SAMPLE}")

    script:
    """
    mkdir raw_fastqc_${SAMPLE}
    fastqc \
                        -o raw_fastqc_${SAMPLE} \
                        -f fastq \
                        -q \
                        -t ${task.cpus} \
                        ${R2} \
                        ${R1} ;

    mkdir clean_fastqc_${SAMPLE}
    fastqc \
                        -o clean_fastqc_${SAMPLE} \
                        -f fastq \
                        -q \
                        -t ${task.cpus} \
                        ${trim_R2} \
                        ${trim_R1} ;
    """

    stub:
    """
    mkdir raw_fastqc_${SAMPLE} clean_fastqc_${SAMPLE}
    """
}

process STAR_paired {
    container 'us-east1-docker.pkg.dev/compute-workspace/omics-docker-repo/rnaseq2'

    publishDir "${params.out_bucket}/${SAMPLE}/logs/STAR/", mode: 'copy', pattern: '*.{out,log}'
    publishDir "${params.out_bucket}/${SAMPLE}/STAR/", mode: 'copy', pattern: 'synced_ctl*.bam'

    input:
    tuple val(SAMPLE), path(trim_R1), path(trim_R2)
    path STARREF

    output:
    path("*")
    tuple val(SAMPLE), path("synced_ctl_${SAMPLE}Aligned.toTranscriptome.out.bam"), emit: Tr_bam

    script:
    """
    # Guard: refuse to align reads that belong to a different sample
    if [ "${trim_R1}" != "synced_clt_${SAMPLE}_R1.fastq.gz" ] || [ "${trim_R2}" != "synced_clt_${SAMPLE}_R2.fastq.gz" ]; then
        echo "SAMPLE MISMATCH: task ${SAMPLE} received ${trim_R1} ${trim_R2}" >&2
        exit 1
    fi

    STAR \
                        --twopassMode Basic \
                        --readFilesCommand zcat \
                        --runThreadN ${task.cpus}-2 \
                        --runMode alignReads \
                        --genomeDir ${STARREF} \
                        --readFilesIn ${trim_R1} ${trim_R2} \
                        --outSAMtype BAM Unsorted \
                        --outSAMstrandField intronMotif \
                        --outFileNamePrefix synced_ctl_${SAMPLE} \
                        --outFilterIntronMotifs RemoveNoncanonical \
                        --outReadsUnmapped Fastx \
			            --quantMode TranscriptomeSAM ;
    """

    stub:
    """
    if [ "${trim_R1}" != "synced_clt_${SAMPLE}_R1.fastq.gz" ] || [ "${trim_R2}" != "synced_clt_${SAMPLE}_R2.fastq.gz" ]; then
        echo "SAMPLE MISMATCH: task ${SAMPLE} received ${trim_R1} ${trim_R2}" >&2
        exit 1
    fi
    sleep \$(( RANDOM % 5 ))
    zcat ${trim_R1} > synced_ctl_${SAMPLE}Aligned.toTranscriptome.out.bam
    cp synced_ctl_${SAMPLE}Aligned.toTranscriptome.out.bam synced_ctl_${SAMPLE}Aligned.out.bam
    touch synced_ctl_${SAMPLE}Log.final.out
    """
}

process salmon_paired {
    container 'us-east1-docker.pkg.dev/compute-workspace/omics-docker-repo/rnaseq2'

    publishDir "${params.out_bucket}/${SAMPLE}/logs", mode: 'copy'

    input:
    path salmon_index
    tuple val(SAMPLE), path(sub_R1), path(sub_R2)

    output:
    path("salmon_quant.log")
    tuple val(SAMPLE), env(rsemSTRAND), emit: salmon_out


    shell:
    '''
    salmon quant \
                        -i !{salmon_index} \
                        -l A \
                        --seqBias \
                        --gcBias \
                        --posBias \
                        -1 !{sub_R1} \
                        -2 !{sub_R2} \
                        -p !{task.cpus} \
                        -o salmon_!{SAMPLE}
    mv salmon_!{SAMPLE}/logs/salmon_quant.log salmon_quant.log
    # Define strandedness by salmon output
    salmonOUTSTRAND="$(grep Automatically salmon_quant.log | tr ' ' '\t' | cut -f 12)"

    if [ "${salmonOUTSTRAND}" == "ISF" ]
    then
        rsemSTRAND="forward"
    elif [ "${salmonOUTSTRAND}" == "ISR" ]
    then
        rsemSTRAND="reverse"
    elif [ "${salmonOUTSTRAND}" == "IU" ]
    then
        rsemSTRAND="none"
    fi

    '''

    stub:
    """
    sleep \$(( RANDOM % 5 ))
    echo "strand_from=\$(zcat ${sub_R1})" > salmon_quant.log
    rsemSTRAND="none"
    """
}

process rsem_expr_paired {
    container 'us-east1-docker.pkg.dev/compute-workspace/omics-docker-repo/rnaseq2'

    publishDir(path: {"${params.out_bucket}/${SAMPLE}/RSEM"}, mode: 'copy', pattern: '*.results')
    publishDir(path: {"${params.out_bucket}/${SAMPLE}/logs"}, mode: 'copy', pattern: '*.log')

    input:
    path(rsem_path)
    tuple val(SAMPLE), path(Tr_bam), val(rsemSTRAND)

    output:
    path("*")

    shell:
    '''
    # Guard: refuse to quantify a BAM that belongs to a different sample
    if [ "!{Tr_bam}" != "synced_ctl_!{SAMPLE}Aligned.toTranscriptome.out.bam" ]; then
        echo "SAMPLE MISMATCH: task !{SAMPLE} received !{Tr_bam}" >&2
        exit 1
    fi

    # Prepare STAR bam file
    samtools view -H !{Tr_bam} > file1 ;
    samtools view -@ !{task.cpus} !{Tr_bam} | awk '{printf "%s", \$0 " "; getline; print}' | sort -S 20G -T ./ | tr ' ' '\n' > file2 ;
    cat file1 file2 | samtools view -@ !{task.cpus} -bS - > rsem_!{SAMPLE}_Aligned.toTranscriptome.out.bam ;

    # RSEM index is supplied weirdely as path_to_index/index_prefix (*.grp, *.chrlist, *idx.fa, etc.)
    # Learned from nf-core rna-seq how to supply rsem index
    INDEX="$(find -L ./ -name "*.grp" | sed 's/\\.grp\$//')"

    # Finally initiate rsem
    rsem-calculate-expression \
                        -p !{task.cpus} \
                        --paired-end \
                        --strandedness !{rsemSTRAND} \
                        --no-bam-output \
                        --estimate-rspd \
                        --alignments \
                        rsem_!{SAMPLE}_Aligned.toTranscriptome.out.bam \
                        ${INDEX} \
                        !{SAMPLE} >& rsem_!{SAMPLE}.log

    sed -i '1i Quantifying expression level for !{SAMPLE} with strandedness set to !{rsemSTRAND}' rsem_!{SAMPLE}.log
    # Clean up for MultiQC
    rm file1
    rm file2
    '''

    stub:
    """
    if [ "${Tr_bam}" != "synced_ctl_${SAMPLE}Aligned.toTranscriptome.out.bam" ]; then
        echo "SAMPLE MISMATCH: task ${SAMPLE} received ${Tr_bam}" >&2
        exit 1
    fi
    cat ${Tr_bam} > ${SAMPLE}.genes.results
    echo "Quantifying expression level for ${SAMPLE} with strandedness set to ${rsemSTRAND}" > rsem_${SAMPLE}.log
    """
}
