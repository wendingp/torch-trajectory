import argparse
import json
import logging
import os
import pickle
import random
import traceback
from pprint import pprint

from scipy.stats.mstats_basic import spearmanr
from tqdm import tqdm


def plot():
    pass


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--original_loc", default="./original_top_k.json")
    parser.add_argument("--transformed_loc", default="./transformed_top_k.json")
    args = parser.parse_args()

    with open(args.original_loc) as fin1:
        original_top_k = json.load(fin1)
    with open(args.transformed_loc) as fin2:
        transformed_top_k = json.load(fin2)

    coefficient, p_val = spearmanr(original_top_k, transformed_top_k)
    pprint(original_top_k)
    pprint(transformed_top_k)
    print(coefficient, p_val)


if __name__ == "__main__":
    main()
