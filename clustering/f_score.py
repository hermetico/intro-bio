import numpy as np
import argparse
import math



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True, type=str, default=None, help="File containing results")
    parser.add_argument("-g", "--gold", required=True, type=str, default=1.0, help="File containing gold standards")


    args = parser.parse_args()

    results = read_results(args.file)
    results = read_gold(args.file)