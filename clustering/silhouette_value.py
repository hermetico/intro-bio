import numpy as np
import argparse


def read_clusters( filename ):
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


def read_expressions( filename ):
    data = {}
    with open(filename) as f:
        for line in f:
            parts = line.rstrip().split('\t')
            if len(parts) > 1:
                data[parts[0]] = [long(val) for val in parts[1:]]
    return data


def map_data( clusters, data ):
    clusters_ids = []
    clusters_data = []
    ids = 0
    for cluster in clusters:
        current_cluster_ids = []
        for patient in cluster:
            clusters_data.append(np.array(data[patient]))
            current_cluster_ids.append(ids)
            ids += 1
        clusters_ids.append(np.array(current_cluster_ids))
    return clusters_ids, np.array(clusters_data)


def centers( c_ids, c_values ):
    '''Returns the centers of each cluster'''
    return np.array([ np.mean(c_values[cluster], axis=0) for cluster in c_ids])

def euclidean(a, b):
    """Calculate the euclidean distance for two points"""
    return np.linalg.norm(a - b)


def closes_cluster( c_centers ):
    distances = []
    for center in c_centers:
        distances.append( np.array([euclidean(center, b) for b in c_centers]))

    closest = []
    for i, distance in enumerate(distances):
        distance[i] = np.inf # put to infinity for the same one
        closest.append(np.argmin(distance))

    return closest


def compute_distances( values ):
    distances = [[ euclidean(gene_a, gene_b)
                   for gene_b in values]
                 for gene_a in values]
    return np.array(distances)


def compute_a_values( clusters, distances):
    results = []
    for cluster in clusters:
        for gene_idx in cluster:
            c_distances = np.array([ distances[gene_idx, other_gene_idx ]
                                for other_gene_idx in cluster
                                if gene_idx != other_gene_idx])
            results.append(c_distances.mean())
    return np.array(results)


def compute_b_values(clusters, distances, closest_clusters):
    results = []
    for i, cluster in enumerate(clusters):
        neighbor_cluster = clusters[closest_clusters[i]]
        for gene_idx in cluster:
            c_distances = np.array([distances[gene_idx, other_gene_idx]
                                    for other_gene_idx in neighbor_cluster])
            results.append(c_distances.mean())
    return np.array(results)


def compute_s_values( a_values, b_values):
    return np.array([ (b_values[i] - a_values[i]) / max((a_values[i], b_values[i] )) for i in xrange(a_values.shape[0])])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True, type=str, default=None, help="File containing results")
    parser.add_argument("-g", "--genes", required=True, type=str, default=None, help="File containing gene expressions standards")


    args = parser.parse_args()

    clusters = read_clusters(args.file)
    data = read_expressions(args.genes)
    ids, values = map_data(clusters, data)
    c_centers = centers(ids, values)
    closest = closes_cluster(c_centers)
    distances = compute_distances(values)
    a_values = compute_a_values(ids, distances)
    b_values = compute_b_values(ids, distances, closest)
    s_values = compute_s_values(a_values, b_values)
    print s_values.mean()