import groovy.json.JsonOutput

include { METASPADES                 } from './modules/spades.nf'
include { MEGAHIT                    } from './modules/megahit.nf'
include { REMOVE_SMALL_CONTIGS;
          FIX_MEGAHIT_CONTIG_NAMING;
          SORT_CONTIGS               } from './modules/helper_scripts.nf'
include { BWA_INDEX;
          BWA                        } from './modules/bwa.nf'
include { MAPPED_READS_TO_FASTQ      } from './modules/samtools.nf'
include { QUAST                      } from './modules/quast.nf'

include { GET_LABEL                  } from  '../../pipeline_resource_optimisation/modules/get_label.nf'

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

        if (!params.model_x_table) {
            error("--metaspades needs --model_x_table: a CSV of model features with one row per sample, keyed by a sample_ID column")
        }

        // One emission per row, keyed by whatever feature columns the CSV happens to have.
        // Every column except sample_ID becomes one JSON object for that sample.
        channel
        .fromPath ( params.model_x_table, checkIfExists: true )
        .splitCsv ( header:true, sep:',', strip:true )
        .map { row ->
            def features = row.findAll { column, value -> column != 'sample_ID' }
                .collectEntries { column, value -> [column, value?.isNumber() ? value as BigDecimal : value] }
            [ row.sample_ID, JsonOutput.toJson(features) ]
        }
        .set { model_x_ch }

        // Pair each sample's reads with its own row of features, so the prediction can be
        // tagged with the sample it belongs to and merged back into that sample's meta
        reads_ch
        .map { meta, read_1, read_2 -> [ meta.ID, meta, read_1, read_2 ] }
        .join ( model_x_ch, failOnMismatch: true, failOnDuplicate: true )
        .set { reads_with_features }

        GET_LABEL( reads_with_features.map { id, meta, read_1, read_2, features -> [ meta, features ] } )

        GET_LABEL.out.label
        .map { meta, label -> [ meta.ID, meta + [ memory_label: "${label.trim()} MB" ] ] }
        .join ( reads_with_features.map { id, meta, read_1, read_2, features -> [ id, read_1, read_2 ] } )
        .map { id, meta, read_1, read_2 -> [ meta, read_1, read_2 ] }
        .set { reads_with_label }

        METASPADES(reads_with_label)
        | REMOVE_SMALL_CONTIGS
        REMOVE_SMALL_CONTIGS.out.long_contigs.set { metaspades_scaffolds }

        REMOVE_SMALL_CONTIGS.out.warning_log.subscribe { file ->
            if (file.exists()) {
                file.readLines().each { line ->
                    log.info "[REMOVE_SMALL_CONTIGS] $line"
                }
            }
        }

        if (params.megahit) {
            BWA_INDEX(metaspades_scaffolds)
            | set { indexed_scaffolds }

            // reads_with_label, not reads_ch: the scaffolds carry the labelled meta, so the
            // join key has to be the labelled meta too
            reads_with_label.join(indexed_scaffolds)
            | BWA
            | MAPPED_READS_TO_FASTQ
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
