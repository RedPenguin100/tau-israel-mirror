#!/usr/bin/env python3
import os
import json
import argparse

def collect_json_files(root_dir):
    return [
        os.path.join(dirpath, f)
        for dirpath, _, files in os.walk(root_dir)
        for f in files if f.endswith(".json")
    ]

def merge_json_files(json_files):
    merged = {}
    for file in json_files:
        with open(file) as f:
            data = json.load(f)

        for seq, chr_dict in data.items():
            if seq not in merged:
                merged[seq] = {}
            # each chr_dict has only one chromosome
            chr_name, entries = next(iter(chr_dict.items()))
            if chr_name not in merged[seq]:
                merged[seq][chr_name] = []
            merged[seq][chr_name].extend(entries)
    return merged

def main():
    parser = argparse.ArgumentParser(description="Merge genome hit JSON files into genome_hits.json")
    parser.add_argument("-d", "--dir", required=True, help="Directory containing JSON files")
    args = parser.parse_args()

    json_files = collect_json_files(args.dir)
    if not json_files:
        print("No JSON files found in the given directory.")
        return

    merged = merge_json_files(json_files)

    with open("genome_hits.json", "w") as out:
        json.dump(merged, out, indent=2)

if __name__ == "__main__":
    main()
