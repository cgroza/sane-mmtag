import cigar
import sys
import pickle

node_sizes_path = sys.argv[1]
pickle_in = sys.argv[2]
pickle_out = sys.argv[3]

def parse_path(s):
    node_pos = 1
    nodes = []
    while node_pos > -1:
        next_node_pos = s.find("s", node_pos + 1)
        node_name = None
        if next_node_pos == -1:
            node_name = s[node_pos:]
        else:
            node_name = s[node_pos:next_node_pos - 1]
        nodes.append((s[node_pos - 1], node_name))
        node_pos = next_node_pos
    return(nodes)

def increment_count(d, node, pos):
    # all elements should already be here
    d[node][pos] = d[node][pos] + 1

pb_aln = sys.stdin

# load node sizes
node_sizes = {}
with open(node_sizes_path) as node_sizes_file:
    for line  in node_sizes_file:
        node_name, node_size = line.split()
        node_sizes[node_name] = int(node_size)

node_bmod_count = pickle.load(open(pickle_in, "rb"))
node_cov_count = {}
for node in node_bmod_count:
    node_cov_count[node] = {}
    for offset in node_bmod_count[node]:
        node_cov_count[node][offset] = 0
        node_cov_count[node][-offset] = 0

nrecords = 0
for line in pb_aln:
    if nrecords % 1000 == 0:
        print("Parsed", nrecords, "records", file=sys.stderr)

    nrecords = nrecords + 1

    fields = line.split()

    assert fields[4] == "+"

    qlen = int(fields[1])
    qstart = int(fields[2])
    qend = int(fields[3])

    path = fields[5]
    plen = int(fields[6])
    pstart = int(fields[7])
    pend = int(fields[8])

    # cigar string
    assert fields[18].startswith("cg:Z:")
    cg = cigar.Cigar(fields[18][5:])

    parsed_path = parse_path(path)

    # list of intervals describing matches on path, left closed, right closed intervals
    pmatches = []
    # list of intervals describing matches on read, left closed, right closed intervals
    qmatches = []

    # record matches along path length
    # cursor to keep track of position in strings
    qi = qstart
    pi = pstart
    nmatch = 0
    nmismatch = 0
    for op in cg.items():
        match op[1]:
            case 'X' | '=':
                # count matches vs mismatches
                if op[1] == "X":
                    nmismatch = nmismatch + op[0]
                else:
                    nmatch = nmatch + op[0]
                pmatches.append((pi, pi + op[0]))
                qmatches.append((qi, qi + op[0]))
                # move the cursors
                qi = qi + op[0]
                pi = pi + op[0]
            # skip on query
            case 'I':
                qi = qi + op[0]
            # skip on path
            case 'D':
                pi = pi + op[0]
            case _:
                raise Exception("Unknown CIGAR operation " + op[1])

    node_breaks = []
    j = 0
    for node in parsed_path:
        node_breaks.append((j, j + node_sizes[node[1]], node))
        j = j + node_sizes[node[1]]


    # for every node covered by the path
    for nb in node_breaks:
        node = nb[2]
        # find the positions covered by a mC
        for offset in node_bmod_count[node[1]]:
            # check if strands match
            if offset >= 0 and node[0] == "<":
                continue

            if offset < 0 and node[0] == ">":
                continue

            # revert the mC offset back to reverse strand if node inversed
            if node[0] == "<":
                offset = node_sizes[node[1]] - abs(offset) - 1  # see offset in lift_5mC.py

            # get the position of this mC relative to the sequence of nodes
            pb = nb[0] + abs(offset)

            # now find if this offset is inside a match in the aligned path
            # pmatches are relative to the node strand
            for pm in pmatches:
                # match
                if pb >= pm[0] and pb <= pm[1]:
                    # increment node coverage
                    # need to go back to forward strand offset if pmatches was on reverse strand
                    if node[0] == "<":
                        offset = -(node_sizes[node[1]] - abs(offset) - 1)
                    increment_count(node_cov_count, node[1], offset)
                    break       # there is only one match on path for this basepair

with open(pickle_out, "wb") as f:
    pickle.dump(node_cov_count, f)
