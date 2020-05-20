import argparse
import json
import logging
import os
import pickle
import random
import traceback
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


def get_coordinates(traj):
    return [(p['x'], p['y']) for p in traj]


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--original_loc", default="./data/original_traj.json")
    parser.add_argument("--transformed_loc", default="./data/transformed_traj.json")
    args = parser.parse_args()

    with open(args.original_loc) as fin1:
        original_traj = json.load(fin1)[0]
    with open(args.transformed_loc) as fin2:
        transformed_traj = json.load(fin2)[0]

    pprint(original_traj)

    original_coords = get_coordinates(original_traj)
    transformed_coords = get_coordinates(transformed_traj)

    x_original, y_original = zip(*original_coords)
    x_trans, y_trans = zip(*transformed_coords)
    fig, axs = plt.subplots(2)
    axs[0].scatter(x_original, y_original)
    axs[1].scatter(x_trans, y_trans)
    plt.show()
    print(len(x_original), len(x_trans))
    print(set(transformed_coords) >= set(original_coords))
    print(set(transformed_coords) - set(original_coords))
    print(set(original_coords) - set(transformed_coords))

if __name__ == "__main__":
    main()
