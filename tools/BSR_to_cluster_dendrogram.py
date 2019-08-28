#!/usr/bin/env python3
import sys
import optparse
try:
    import numpy as np
except:
    print("numpy must be installed, try: conda install -c anaconda numpy")
    sys.exit()
import itertools

# Tree creation
try:
    import pandas as pd
except:
    print("pandas must be installed, try: conda install -c anaconda pandas")
    sys.exit()
try:
    from scipy.cluster.hierarchy import average
    from scipy.cluster.hierarchy import weighted
except:
    print("scipy must be installed, try: conda install -c anaconda scipy")
    sys.exit()
try:
    from skbio.tree import TreeNode
except:
    print("skbio must be installed, try: conda install -c anaconda scikit-bio")
    sys.exit()

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s file cannot be opened' % option)
        sys.exit()

def read_file(file):
    array = []
    with open(file) as infile:
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
    outfile.close()

def write_tree(cluster_method):
    import scipy.spatial.distance as ssd
    dmx = pd.read_csv("distance_matrix", index_col=0, sep="\t")
    ids = dmx.index.tolist()
    triu = np.square(dmx.values)
    distArray = ssd.squareform(triu)
    if cluster_method == "average":
        hclust = average(distArray)
    elif cluster_method == "weighted":
        hclust = weighted(distArray)
    else:
        print("invalid cluster method chosen")
        sys.exit()
    t = TreeNode.from_linkage_matrix(hclust, ids)
    nw = t.__str__().replace("'", "")
    outfile = open("bsr_matrix.tree", "w")
    outfile.write(nw)
    outfile.close()

def main(bsr_matrix,cluster_method):
    fields, data = read_file(bsr_matrix)
    distances = calculate_distances(data)
    matrix = make_symmetric_matrix(distances)
    write_matrix(matrix, fields)
    write_tree(cluster_method)

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-b", "--bsr_matrix", dest="bsr_matrix",
                      help="PATH to bsr matrix",
                      type="string", action="callback", callback=test_file)
    parser.add_option("-m", "--method", dest="cluster_method",
                      help="cluster method. Choose from weighted [default] or average",
                      type="string", action="store", default="weighted")
    options,args = parser.parse_args()
    mandatories = ["bsr_matrix"]
    for m in mandatories:
        if not getattr(options, m, None):
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)
    main(options.bsr_matrix,options.cluster_method)
