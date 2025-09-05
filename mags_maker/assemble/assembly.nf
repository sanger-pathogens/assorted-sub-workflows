include { METASPADES                 } from './modules/spades.nf'
include { MEGAHIT                    } from './modules/megahit.nf'
include { REMOVE_SMALL_CONTIGS;
          FIX_MEGAHIT_CONTIG_NAMING;
          SORT_CONTIGS               } from './modules/helper_scripts.nf'
include { BWA_INDEX;
          BWA                        } from './modules/bwa.nf'
include { SAM_TO_FASTQ               } from './modules/samtools.nf'
include { QUAST                      } from './modules/quast.nf'

/*
#############################################################################################################################################################
#
# This script is meant to be a comprehensive solution for producing the best metagenomic assembly given paired end reads from one or more samples.
# Ideally it should take in fully QC'd reads. First, the reads are assembled with metaSPAdes3.14, then all the reads that did not map back to the
# contigs are re-assembled with MEGAHIT (which works better on lower coverage contigs. The resulting assemblies are combined, sorted, and short 
# contigs are removed. 
#
# Author of original pipeline: German Uritskiy. I do not have any authorship of the many programs this pipeline uses.
# For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.
#
# Modified by Yan Shao: 1) Default hybrid assembly (metaspades followed by megahit) assemblers ; 2) Added --fastspades option (use if default metaspades runs out of time);
# 3) Disabled plots; 4) Disabled assembly QC by QUAST 
#
# Modified by Sam Dougan into nextflow :)
##############################################################################################################################################################
*/

workflow METAWRAP_ASSEMBLE {
    take:
    reads_ch

    main:

    if (!params.metaspades && !params.megahit) {
        log.warn("have to select at least one of --metaspades or --megahit")
    }

    if (params.metaspades) {
        METASPADES(reads_ch)
        | REMOVE_SMALL_CONTIGS
        REMOVE_SMALL_CONTIGS.out.long_contigs.set { metaspades_scaffolds }

        if (params.megahit) {
            BWA_INDEX(metaspades_scaffolds)
            | set { indexed_scaffolds }

            reads_ch.join(indexed_scaffolds)
            | BWA
            | SAM_TO_FASTQ
            | set { final_reads_ch }
        }

    } else {
        final_reads_ch = reads_ch
    }

    if (params.megahit) {
        MEGAHIT(final_reads_ch)
        | FIX_MEGAHIT_CONTIG_NAMING
        | set { megahit_contigs }
    }

    if (params.metaspades && params.megahit) {
        metaspades_scaffolds.join(megahit_contigs)
        | map { meta, contig1, contig2 -> tuple(meta, [contig1, contig2])}
        | set { final_contigs }
        
    } else if (params.metaspades) {
        metaspades_scaffolds
        | set { final_contigs }

    } else if (params.megahit) {
        megahit_contigs
        | set { final_contigs }
    }

    SORT_CONTIGS(final_contigs)

    if (params.assembly_stats) {
        QUAST(SORT_CONTIGS.out.sorted_contigs)
    }

    emit:
    SORT_CONTIGS.out.sorted_contigs
}