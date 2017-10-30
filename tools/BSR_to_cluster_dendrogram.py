#!/usr/bin/env python3
import sys
try:
    import numpy as np
except:
    print("numpy must be installed")
    sys.exit()
import itertools

# Tree creation
try:
    import pandas as pd
except:
    print("pandas must be installed")
    sys.exit()
try:
    from scipy.cluster.hierarchy import weighted
except:
    print("scipy must be installed")
    sys.exit()
try:
    from skbio.tree import TreeNode
except:
    print("skbio must be installed")
    sys.exit()

def main(infile):
    fields, data = read_file(sys.argv[1])
    distances = calculate_distances(data)
    matrix = make_symmetric_matrix(distances)
    write_matrix(matrix, fields)
    write_tree()

def read_file(file):
    array = []
    with open(file, "rU") as infile:
        fields = infile.readline().split()
        for line in infile:
            array.append([float(x) for x in line.split()[1:]])
    return fields, np.array(array)

def calculate_distances(data):
    length = data.shape[1]
    out = {}
    for i in range(length):
        out[i] = [0]*(i+1)
    def _distance(cols):
        return np.mean(np.abs(np.diff(data[:, cols])))
    for combo in itertools.combinations(range(length), 2):
        out[combo[0]].append(_distance(combo))
    return out

def make_symmetric_matrix(distances):
    a = []
    i = 0
    while i in distances:
        a.append(distances[i])
        i += 1
    array = np.around(np.array(a), 2)
    return array + array.T - np.diag(array.diagonal())

def write_matrix(matrix, fields):
    with open("distance_matrix", "w") as outfile:
        outfile.write("\t" + "\t".join(fields))
        outfile.write("\n")
        for idx, row in enumerate(matrix):
            outfile.write(fields[idx] + "\t")
            outfile.write("\t".join(map(str,row)))
            outfile.write("\n")

def write_tree():
    dmx = pd.read_csv("distance_matrix", index_col=0, sep="\t")
    ids = dmx.index.tolist()
    triu = np.square(dmx.as_matrix())
    hclust = weighted(triu)
    t = TreeNode.from_linkage_matrix(hclust, ids)
    nw = t.__str__().replace("'", "")
    outfile = open("bsr_matrix.tree", "w")
    outfile.write(nw)
    outfile.close()

if __name__ == "__main__":
    try:
        main(sys.argv[1])
    except:
        print("usage: script bsr_matrix")
        sys.exit()

