include { COMBINE_BINS; SPLIT_READS } from "../modules/reassemble/helper_scripts.nf"
include { BWA_INDEX; BWA            } from '../modules/assemble/bwa.nf'
include { SPADES_REASSEMBLE         } from '../modules/reassemble/spades.nf'

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

def filesFromDir(dirPath) {
	def basePath = dirPath instanceof Path ? dirPath : file(dirPath.toString())
    return basePath.list()
        .findAll { basePath.resolve(it).isFile() }
        .collect { basePath.resolve(it) }
}

workflow METAWRAP_REASSEMBLY {
    take:
    un_trusted_contigs
    reads

    main:
    COMBINE_BINS(un_trusted_contigs)
    | BWA_INDEX
    | set { index }

    reads
    | join(index)
    | BWA
    | set { sam }

    un_trusted_contigs
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
		
		[key, meta, single_bin, read_1, read_2]
	}
	| set { split_reads }

	un_trusted_contigs
	| map { meta, directory_path ->
		def file_list = filesFromDir(directory_path)
		[meta, file_list]
	}
	| transpose
	| set { original_bins }

	original_bins
	| map { meta, bin_path -> 
		def basename = bin_path.getName() // e.g., "bin.9.fasta"

		def parts = basename.split('\\.')  // [ "bin", "9", "fasta" ]

		def bin = parts[1] as Integer

		def bin_info = [ bin: bin ]

		def key = "${meta.ID}_${bin_info.bin}"

		[key, meta, bin_info, bin_path]
	}
	| join(split_reads)
	| map { joinKey, metaOne, binInfoSmall, binFasta, metaTwo, fullBinInfo, read_1, read_2 -> [metaOne, fullBinInfo, binFasta, read_1, read_2] }
	| SPADES_REASSEMBLE
	| mix(original_bins)
	| groupTuple
	| view



}
/*
# removing short contigs and placing reassemblies in the final folder
comm "Finalizing reassemblies"
mkdir ${out}/reassembled_bins
for i in $( ls ${out}/reassemblies/ ); do
	spades_folder=${out}/reassemblies/$i
	bin_name=${spades_folder##}
	
	#remove shortest contigs (probably artifacts...)
	if [ -s ${out}/reassemblies/${bin_name}/scaffolds.fasta ]; then
		${SOFT}/rm_short_contigs.py $len\
		 ${out}/reassemblies/${bin_name}/scaffolds.fasta\
		 > ${out}/reassemblies/${bin_name}/long_scaffolds.fasta

		if [ -s ${out}/reassemblies/${bin_name}/long_scaffolds.fasta ]; then
			echo "$bin_name was reassembled! Processing..."
			mv ${out}/reassemblies/${bin_name}/long_scaffolds.fasta\
			${out}/reassembled_bins/${bin_name}.fa
		else
			comm "$bin_name was reassembled, but did not yeild contigs $len bp. It is possible there were not enough reads."
		fi
	else
		comm "$bin_name was not successfully reassembled. It is possible there were not enough reads."
	fi
done


if [[ $(ls ${out}/reassembled_bins/ | wc -l) -lt 1 ]]; then
	error "None of the bins were successfully reassembled. ${out}/reassembled_bins/ is empty."
else
	comm "Looks like the reassemblies went well. Now to see if they made the bins better or worse..."
fi



if [ "$run_checkm" = true ]; then
	########################################################################################################
	########################             RUN CHECKM ON REASSEMBLED BINS             ########################
	########################################################################################################
	announcement "RUN CHECKM ON REASSEMBLED BINS"

	# determine --pplacer_threads count. It is either the max thread count or RAM/4, whichever is higher
	ram_max=$(($mem / 40))
	if (( $ram_max < $threads )); then
		p_threads=$ram_max
	else
		p_threads=$threads
	fi
	comm "There is $mem RAM and $threads threads available, and each pplacer thread uses ~40GB, so I will use $p_threads threads for pplacer"

	# copy over original bins
	for base in $( ls ${out}/original_bins/ | grep "\.fa$" ); do 
		i=${out}/original_bins/$base
		cp $i ${out}/reassembled_bins/${base%.*}.orig.fa
	done

	comm "Running CheckM on best bins (reassembled and original)"
	if [[ -d ${out}/reassembled_bins.checkm ]]; then rm -r ${out}/reassembled_bins.checkm; fi
	mkdir ${out}/tmp
	checkm lineage_wf -x fa ${out}/reassembled_bins ${out}/reassembled_bins.checkm -t $threads --tmpdir ${out}/tmp --pplacer_threads $p_threads
	if [[ ! -s ${out}/reassembled_bins.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi
	${SOFT}/summarize_checkm.py ${out}/reassembled_bins.checkm/storage/bin_stats_ext.tsv | (read -r; printf "%s\n" "$REPLY"; sort) > ${out}/reassembled_bins.stats
	if [[ $? -ne 0 ]]; then error "Cannot make checkm summary file. Exiting."; fi
	rm -r ${out}/tmp


	########################################################################################################
        ########################          FINDING THE BEST VERSION OF EACH BIN          ########################
	########################################################################################################
	announcement "FINDING THE BEST VERSION OF EACH BIN"

	if [ ! -d ${out}/reassembled_best_bins ]; then mkdir ${out}/reassembled_best_bins; fi
	for i in $(${SOFT}/choose_best_bin.py ${out}/reassembled_bins.stats $comp $cont); do 
		echo "Copying best bin: $i"
		cp ${out}/reassembled_bins/${i}.fa ${out}/reassembled_best_bins 
	done
	
	o=$(ls -l ${out}/reassembled_best_bins | grep orig | wc -l)
	s=$(ls -l ${out}/reassembled_best_bins | grep strict | wc -l)
	p=$(ls -l ${out}/reassembled_best_bins | grep permissive | wc -l)

	announcement "Reassembly results are in! $s bins were improved with 'strict' reassembly, $p bins were improved with 'permissive' reassembly, and $o bins were not improved by any reassembly, and thus will stay the same."
	
	if [[ $(ls ${out}/reassembled_best_bins | wc -l) -gt 0 ]]; then 
		comm "Seems that the reassembly went well. You will find the final, best, reassembled bins in ${out}/reassembled_bins, and all intermediate files in ${out}/work_files (which we recomend you delete to save space after you confirm that the pipeline worked)"
		mkdir ${out}/work_files
		mv ${out}/reassembled_bins ${out}/work_files/
		mv ${out}/reassembled_bins.checkm ${out}/work_files/
		mv ${out}/reassembled_bins.stats ${out}/work_files/
		mv ${out}/reads_for_reassembly ${out}/work_files/
		mv ${out}/nanopore_reads_for_reassembly ${out}/work_files/
		mv ${out}/binned_assembly ${out}/work_files/
		mv ${out}/reassemblies ${out}/work_files/
		rm -r ${out}/original_bins
		mv ${out}/reassembled_best_bins ${out}/reassembled_bins 
	else
		error "there are no good bins found in ${out}/reassembled_best_bins - something went wrong with choosing the best bins between the reassemblies."
	fi


	comm "Re-running CheckM on the best reasembled bins."
	if [[ -d ${out}/reassembled_bins.checkm ]]; then rm -r ${out}/reassembled_bins.checkm; fi
	mkdir ${out}/tmp
        checkm lineage_wf -x fa ${out}/reassembled_bins ${out}/reassembled_bins.checkm -t $threads --tmpdir ${out}/tmp --pplacer_threads $p_threads
        if [[ ! -s ${out}/reassembled_bins.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi
	rm -r ${out}/tmp
        comm "Finalizing CheckM stats..."
        ${SOFT}/summarize_checkm.py ${out}/reassembled_bins.checkm/storage/bin_stats_ext.tsv | (read -r; printf "%s\n" "$REPLY"; sort) > ${out}/reassembled_bins.stats
        if [[ $? -ne 0 ]]; then error "Cannot make checkm summary file. Exiting."; fi

#        comm "Making CheckM plot of ${out}/reassembled_bins bins"
#        checkm bin_qa_plot -x fa ${out}/reassembled_bins.checkm ${out}/reassembled_bins ${out}/reassembled_bins.plot
#        if [[ ! -s ${out}/reassembled_bins.plot/bin_qa_plot.png ]]; then warning "Something went wrong with making the CheckM plot. Exiting."; fi
#        mv ${out}/reassembled_bins.plot/bin_qa_plot.png ${out}/reassembled_bins.png
#        rm -r ${out}/reassembled_bins.plot
	
	comm "you will find the info on the final reassembled bins in ${out}/reassembled_bins.stats, and a figure summarizing it in ${out}/reassembled_bins.png"

#	comm "making reassembly N50, compleiton, and contamination summary plots."
#	head -n 1 ${out}/work_files/reassembled_bins.stats > ${out}/original_bins.stats
#	grep orig ${out}/work_files/reassembled_bins.stats >> ${out}/original_bins.stats
#	${SOFT}/plot_reassembly.py $out $comp $cont ${out}/reassembled_bins.stats ${out}/original_bins.stats
#	if [[ $? -ne 0 ]]; then error "Something went wrong with plotting the reassembly summary plots. Exiting..."; fi
fi

##### CLEAN UP #####
comm "cleaning up to save space..."
rm -r  ${out}/work_files
rm -r  ${out}/reassembled_bins.checkm

comm "you will find the final bins in ${out}/reassembled_bins"
########################################################################################################
########################    REASSEMBLY PIPELINE SUCCESSFULLY FINISHED!!!        ########################
########################################################################################################
announcement "BIN REASSEMBLY PIPELINE SUCCESSFULLY FINISHED!!!"
*/
