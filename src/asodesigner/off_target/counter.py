#!/usr/bin/env python3
import json
import argparse
from collections import Counter

def main():
    parser = argparse.ArgumentParser(description="Count k-values in genome_hits.json per sequence")
    parser.add_argument("-i", "--input", default="genome_hits.json", help="Input JSON file")
    parser.add_argument("-o", "--output", default="genome_hits_count.json", help="Output JSON file")
    args = parser.parse_args()

    with open(args.input) as f:
        genome_hits = json.load(f)

    result = {}
    for seq, chr_dict in genome_hits.items():
        counter = Counter()
        for entries in chr_dict.values():  # each chromosome
            for _, _, k in entries:
                counter[k] += 1
        # make list [count_0, count_1, count_2]
        result[seq] = [counter.get(i, 0) for i in range(3)]

    with open(args.output, "w") as out:
        json.dump(result, out, indent=2)

if __name__ == "__main__":
    main()
