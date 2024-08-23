include { SEARCH_LIBRARY; KMER2READ_DIST; GENERATE_KMER_DISTRIBUTION } from '../modules/kraken2bracken/bracken'

workflow BRACKEN_BUILD {
    // Unfortunately, due to checks in bracken_build script, softlinking to the kraken2_db subdirs doesn't work.
    // Copying dir structure and hardlinking is not flexible enough (cannot hardlink between filesystems).
    // Recursive copy of all data will work but eats up space temporarily.
    // Run individual bracken scripts is the most flexible option, but less pretty (see https://github.com/jenniferlu717/Bracken#step-1-generate-the-bracken-database-file-databasexmerskmer_distrib-1)
    take:
    kraken2_db

    main:
    SEARCH_LIBRARY(kraken2_db)
    SEARCH_LIBRARY.out.kraken2_library.dump(tag: "kraken2_library").set { ch_kraken2_library }

    KMER2READ_DIST(ch_kraken2_library)
    KMER2READ_DIST.out.kmer2read_distr.dump(tag: "kmer2read_distr").set { ch_kmer2read_distr }

    GENERATE_KMER_DISTRIBUTION(ch_kmer2read_distr)
    GENERATE_KMER_DISTRIBUTION.out.kmer_distrib.dump(tag: "kmer_distrib").set { ch_kmer_distrib }

    emit:
    ch_kmer_distrib = GENERATE_KMER_DISTRIBUTION.out.kmer_distrib
}