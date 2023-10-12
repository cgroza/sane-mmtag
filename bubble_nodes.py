import sys

uniq_nodes = set()

bubbles_path = sys.argv[1]

chrom = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])

with open(sys.argv[1], "r") as bubbles:
    for bubble in bubbles:
        fields = bubble.split()
        bstart = int(fields[1])
        bend = int(fields[2])
        if not (fields[0] == chrom and bstart >= start and bend <= end):
            continue
        nodes = fields[11].split(",")
        for node in nodes:
            uniq_nodes.add(node)

for node in uniq_nodes:
    print(node)
