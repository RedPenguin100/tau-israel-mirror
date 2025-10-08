import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass

from fuzzysearch import find_near_matches
from tqdm import tqdm

from asodesigner.consts import EXPERIMENT_RESULTS
from asodesigner.timer import Timer
from asodesigner.util import get_antisense


class LocusInfoOld:
    def __init__(self, seq=None):
        # Default empty fields
        self.exons = []
        self.introns = []
        self.exon_indices = []
        self.intron_indices = []
        self.stop_codons = []
        self.five_prime_utr = ""
        self.three_prime_utr = ""
        self.exon_concat = None
        self.full_mrna = None
        self.cds_start = None
        self.cds_end = None
        self.strand = None
        self.gene_type = None
        self.utr_indices = []

        # If a sequence is provided, create a simple gene with one exon
        if seq is not None:
            self.exons = [seq]
            self.exon_indices = [(0, len(seq) - 1)]
            self.cds_start = 0
            self.cds_end = len(seq) - 1
            self.strand = "+"
            self.gene_type = "unknown"
            self.exon_concat = seq
            self.full_mrna = seq

    def __repr__(self):
        print("LocusInfo:")
        for field, value in self.__dict__.items():
            print(f"  {field}: {value}")


def get_simplified_fasta_dict(fasta_dict):
    simplified_fasta_dict = dict()
    for locus_tag, locus_info in fasta_dict.items():
        simplified_fasta_dict[locus_tag] = str(locus_info.upper().seq)
    return simplified_fasta_dict


def validated_get_simplified_fasta_dict(fasta_dict, simplified_fasta_dict):
    if simplified_fasta_dict is None and fasta_dict is None:
        raise ValueError('Either simplified_fasta_dict or fasta_dict must be specified')

    if simplified_fasta_dict is None:
        return get_simplified_fasta_dict(fasta_dict)
    return simplified_fasta_dict


def process_watson_crick_differences(args):
    idx, l, aso_sense, locus_to_data = args
    matches_per_distance = [0, 0, 0, 0]

    for locus_tag, locus_info in locus_to_data.items():
        matches = find_near_matches(aso_sense, locus_info, max_insertions=0, max_deletions=0, max_l_dist=3)
        for match in matches:
            matches_per_distance[match.dist] += 1
            if match.dist == 0:
                print(locus_tag)

    # Return a tuple containing the starting index, current l, and match counts
    return (idx, l, matches_per_distance[0],
            matches_per_distance[1], matches_per_distance[2], matches_per_distance[3])


def validate_organism(organism: str):
    organisms = ['human', 'yeast']
    if organism not in organisms:
        raise ValueError(f'Organism={organism} must be in {organisms}')


def parallelize_function(function, tasks, max_threads=None):
    """
    :param function: To be parallelized
    :param tasks: to be submitted to function
    :param max_threads: pass None to use all cores
    :return: results of parallel operation
    """
    results = []
    with Timer() as t:
        with ProcessPoolExecutor(max_workers=max_threads) as executor:
            futures = [executor.submit(function, task) for task in tasks]

            for future in tqdm(as_completed(futures), total=len(futures), desc='Processing'):
                results.append(future.result())
    print(f"Parallel task done in: {t.elapsed_time}s")
    return results


class Task:
    def __init__(self, sense_start, sense_length, sense, simplified_fasta_dict, target_cache_filename):
        self.sense_start = sense_start
        self.sense_length = sense_length
        self.sense = sense
        self.simplified_fasta_dict = simplified_fasta_dict
        self.target_cache_filename = target_cache_filename
        # Settings. TODO: consider moving to separate class
        self.minimum_score = 900
        self.parsing_type = '2'
        self.binary_cutoff = -20

    def get_sense(self):
        return self.sense

    def get_antisense(self):
        return get_antisense(self.get_sense())


@dataclass
class ResultHybridization:
    sense_start: int
    sense_length: int
    total_hybridization_candidates: int
    total_hybridization_energy: int
    total_hybridization_max_sum: int
    total_hybridization_binary_sum: int

    @staticmethod
    def results_to_result_dict(results, experiment):
        results_dict = dict()
        for result in results:
            start = result.sense_start
            length = result.sense_length
            total_hybridization_candidates: int
            total_hybridization_energy: int
            total_hybridization_max_sum: int
            total_hybridization_binary_sum: int

            antisense = experiment.get_aso_antisense_by_index(idx=start, length=length)
            results_dict[antisense] = (
                result.total_hybridization_candidates, result.total_hybridization_energy,
                result.total_hybridization_max_sum, result.total_hybridization_binary_sum)

        return results_dict
