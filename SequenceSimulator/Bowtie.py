# Takes in our sam file and Bowtie's sam file and does a comparison, gives result in percent
def compare_sam_files(nmd_sam_file, bwa_sam_file):
    alignments = {}
    our_sam = open(nmd_sam_file, "r")
    correctly_aligned = calculate_correctly_aligned(alignments, our_sam, bwa_sam_file)
    results = correctly_aligned * 100 / len(alignments)
    result_file = open("resultsBowTie.txt", "w")
    result_file.write("BWA-MEM accuracy of alignment is: {}%.\n".format(results))

# Helping method for comparing SAMs, takes into account only alignment lines, skips those that start with a '@'
# which usually contain a name or some other identifier
def calculate_correctly_aligned(alignments, our_sam, bwa_sam_file):
    correctly_aligned = 0
    for line in our_sam.readlines():
        if line[0] == '@':
            continue
        split = line.split('\t')
        map_id = split[0] + split[2]
        alignments[map_id] = split[1]
    our_sam.close()
    bwa_sam = open(bwa_sam_file, "r")
    for line in bwa_sam.readlines():
        if line[0] == '@':
            continue
        split = line.split('\t')
        map_id = split[0] + split[9]
        if alignments.get(map_id) == split[3]:
            correctly_aligned += 1

    return correctly_aligned


compare_sam_files('Test.sam', 'bowtiesam')
