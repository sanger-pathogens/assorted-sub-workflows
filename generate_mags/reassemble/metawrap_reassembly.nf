include { COMBINE_BINS; SPLIT_READS } from "../modules/reassemble/helper_scripts.nf"
include { BWA_INDEX; BWA            } from '../modules/assemble/bwa.nf'
include { SPADES_REASSEMBLE         } from '../modules/reassemble/spades.nf'
include { REMOVE_SMALL_CONTIGS; 
		  COLLECT_BINS;
		  RENAME_ORIGINAL;
		  CHOOSE_BEST_BIN			} from '../modules/reassemble/helper_scripts.nf'
include { CHECKM; SUMMARISE_CHECKM	} from '../modules/reassemble/checkm1.nf'
include { CHECKM2					} from '../modules/reassemble/checkm2.nf'

/*
##############################################################################################################################################################
#
# Improves a set of bins by aligning reads back to the bins and reassembling them.
#
# Author of pipeline: German Uritskiy. I do not clain any authorship of the many programs this pipeline uses.
# For questions, bugs, and suggestions, contact me at guritsk1@jhu.edu.
#
# Modified by Yan Shao: 1) default checkM 50/5 to medium-quality ; 2) accept fastq.gz; 3) output clean up (only keep final bins and stats)
##############################################################################################################################################################
*/

def filesFromDir(Path dirPath) {
    return dirPath.toFile()
        .listFiles()
        .findAll { it.isFile() }
        .collect { it.toPath() }
}

workflow METAWRAP_REASSEMBLY {
    take:
    un_trusted_contigs
    reads
	//checkm_reports

    main:
	un_trusted_contigs
	| map {meta, fastas, checkm_report -> [meta, fastas] }
	| set { to_reassemble_bins }

    COMBINE_BINS(to_reassemble_bins)
    | BWA_INDEX
    | set { index }

    reads
    | join(index)
    | BWA
    | set { sam }

    to_reassemble_bins
    | join(sam)
    | SPLIT_READS
	| transpose
	| map { meta, path -> 
		def basename = path.getName() // e.g., "bin.9.strict_2.fastq.gz"

		def parts = basename.split('\\.')  // [ "bin", "9", "strict_2", "fastq", "gz" ]

		def bin = parts[1] as Integer
		def level = parts[2].split('_')[0] //remove read info (1 or 2)
		def join_name = "${meta.ID}_${bin}_${level}" // to join reads

		def bin_info = [ bin: bin, level: level ]

		[join_name, meta, bin_info, path]
	}
	| groupTuple
	| map { join_name, meta, bin_info, path ->
		//drop the join name and take just one of the now duplicated maps
		def single_meta = meta.first()
		def single_bin = bin_info.first()
		def key = "${single_meta.ID}_${single_bin.bin}"

		def (read_1, read_2) = path.sort { it.name }
		
		[key, single_meta, single_bin, read_1, read_2]
	}
	| set { split_reads }

	to_reassemble_bins
	| map { meta, directory_path ->
		def file_list = filesFromDir(directory_path)
		[meta, file_list]
	}
	| transpose
	| set { original_bins }

	original_bins
	| RENAME_ORIGINAL
	| set { renamed_original_bins }

	original_bins
	| map { meta, bin_path -> 
		def basename = bin_path.getName() // e.g., "bin.9.fasta"

		def parts = basename.split('\\.')  // [ "bin", "9", "fasta" ]

		def bin = parts[1] as Integer

		def bin_info = [ bin: bin ]

		def key = "${meta.ID}_${bin_info.bin}"

		[key, meta, bin_info, bin_path]
	}
	| combine(split_reads, by: 0)
	| map { joinKey, metaOne, binInfoSmall, binFasta, metaTwo, fullBinInfo, read_1, read_2 -> [metaOne, fullBinInfo, binFasta, read_1, read_2] }
	| SPADES_REASSEMBLE
	| REMOVE_SMALL_CONTIGS
	| mix(renamed_original_bins)
	| groupTuple
	| COLLECT_BINS
	| set { reassembled_bins }
	
	if (params.checkm1) {
        CHECKM(reassembled_bins)
        | SUMMARISE_CHECKM
        | set { checkm }
    } else {
        CHECKM2(reassembled_bins)
        | set { checkm }
    }

	checkm
	| CHOOSE_BEST_BIN
	| set { final_bins }

	emit:
	final_bins
}