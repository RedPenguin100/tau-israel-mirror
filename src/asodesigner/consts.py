import os
from pathlib import Path

PROJECT_PATH = Path('/tmp') / '.cache' / 'asodesigner'

TEST_PATH = PROJECT_PATH

DATA_PATH = PROJECT_PATH / 'data'

# tmp folder - for files that are dumped to disk during calculation
TMP_PATH = PROJECT_PATH / 'tmp'

# Experiments
EXPERIMENT_RESULTS = PROJECT_PATH / 'experiment_results'
CACHE_DIR = PROJECT_PATH / 'cache'

# Yeast
YEAST_DATA = DATA_PATH / 'yeast'
YEAST_FASTA_PATH = YEAST_DATA / 'GCF_000146045.2_R64_genomic.fna'
YEAST_GFF_PATH = YEAST_DATA / 'genomic.gff'
YEAST_GFF_DB_PATH = YEAST_DATA / 'yeast_gff.db'
YEAST_FIVE_PRIME_UTR = YEAST_DATA / 'SGD_all_ORFs_5prime_UTRs.fsa'
YEAST_THREE_PRIME_UTR = YEAST_DATA / 'SGD_all_ORFs_3prime_UTRs.fsa'
YEAST_README = YEAST_DATA / 'README.md'

# Human
HUMAN_DATA = DATA_PATH / 'human'
HUMAN_V34 = HUMAN_DATA / 'human_v34'
HUMAN_GTF_BASIC_GZ = HUMAN_V34 / 'gencode.v34.basic.annotation.gtf.gz'
HUMAN_GTF_BASIC = HUMAN_V34 / 'gencode.v34.basic.annotation.gtf'
HUMAN_GTF = HUMAN_GTF_BASIC
HUMAN_GTF_GZ = HUMAN_GTF_BASIC_GZ
HUMAN_GFF = HUMAN_GTF  # backwards-compatible alias
HUMAN_GFF_GZ = HUMAN_GTF_GZ
HUMAN_TRANSCRIPTS_FASTA_GZ = HUMAN_V34 / 'gencode.v34.transcripts.fa.gz'
HUMAN_TRANSCRIPTS_FASTA = HUMAN_V34 / 'gencode.v34.transcripts.fa'
HG38_CACHE_DIR = CACHE_DIR / 'genomes' / 'hg38'
HUMAN_GENOME_FASTA_GZ = HG38_CACHE_DIR / 'hg38.fa.gz'
HUMAN_GENOME_FASTA = HG38_CACHE_DIR / 'hg38.fa'
HUMAN_GENOME_FASTA = HUMAN_V34 / 'GRCh38.p13.genome.fa' ## TO EDIT
HUMAN_DB_PATH = HUMAN_V34 / 'dbs'
HUMAN_DB = HUMAN_DB_PATH / 'human_gff.db'
HUMAN_DB_BASIC = HUMAN_DB_PATH / 'human_gff_basic.db'
HUMAN_DB_BASIC_INTRONS = HUMAN_DB_PATH / 'human_gff_basic_introns.db'
HUMAN_DB_BASIC_INTRONS_GZ = HUMAN_DB_PATH / 'human_gff_basic_introns.db.gz'

# External
EXTERNAL_PATH = PROJECT_PATH / 'external'
RISEARCH_PATH = EXTERNAL_PATH / 'risearch'
RISEARCH1_PATH = RISEARCH_PATH / 'RIsearch1'
RISEARCH1_BINARY_PATH = RISEARCH1_PATH  / 'RIsearch'


if __name__ == "__main__":
	print(HUMAN_DB_BASIC_INTRONS)
	print(HUMAN_GENOME_FASTA)