import numpy as np
import leidenalg
import igraph as ig
import sys
import json

def main():
    # Example adjacency matrix (replace with your data)
    adj_matrix = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]])

    # Check for command line argument for input file
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        adj_matrix = np.loadtxt(input_file)

    # Create a graph from the adjacency matrix
    graph = ig.Graph.Adjacency((adj_matrix > 0).tolist())
    graph.es['weight'] = adj_matrix[adj_matrix.nonzero()]

    # Perform Leiden clustering
    partition = leidenalg.find_partition(graph, leidenalg.ModularityVertexPartition)

    # Get the clustering result
    clusters = partition.membership

    # Save the result to a file
    with open('clusters.json', 'w') as f:
        json.dump(clusters, f)

    print("Leiden clustering completed. Results saved to clusters.json")

if __name__ == '__main__':
    main()
