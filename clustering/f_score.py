import numpy as np
import argparse


def read_input( filename ):
    data = []
    cluster_data = []
    current_cluster = 0
    with open(filename) as f:
        for line in f:
            parts = line.rstrip().split('\t')
            if len(parts) == 2:
                if int(parts[1]) != current_cluster:
                    data.append(cluster_data)
                    current_cluster = int(parts[1])
                    cluster_data = []
                cluster_data.append(parts[0])

        if len(cluster_data) > 0:
            data.append(cluster_data)
    return data


def mapping_clusters( cla, clb):
    mapping = {}
    for i, data_a in enumerate(cla):

        commons = [0] * len(clb)
        for j, data_b in enumerate(clb):
            commons[j] = sum([ 1 for a in data_a if a in data_b])

        mapping[i] = np.array(commons).argmax()
    return mapping

def compute_f_score(cla, clb, mappings):

    for i, cl in enumerate(cla):
        mapp = clb[mappings[i]]
        TP = sum([1 for a in cl if a in mapp]) * 1.
        FP = sum([1 for a in cl if a not in mapp]) * 1.
        FN = sum([1 for a in mapp if a not in cl]) * 1.
        #TN = sum([1 for a in mapp if a in cl]) * 1.

    precission = TP / (TP + FP)
    recall = TP / (TP + FN)
    return  2 * precission * recall / (precission + recall)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True, type=str, default=None, help="File containing results")
    parser.add_argument("-g", "--gold", required=True, type=str, default=None, help="File containing gold standards")


    args = parser.parse_args()

    results = read_input(args.file)
    gold = read_input(args.gold)
    if len(results) != len(gold):
        print "There aren't the same number of clusters"
        exit(1)
    mappings = mapping_clusters(results, gold)
    f_score = compute_f_score(results, gold, mappings)
    print f_score