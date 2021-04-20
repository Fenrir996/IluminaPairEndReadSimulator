def compare_sam_files(nmd_sam_file, bwa_sam_file):
    alignments = {}
    our_sam = open(nmd_sam_file, "r")
    correctly_aligned = calculate_correctly_aligned(alignments, our_sam, bwa_sam_file)
    results = correctly_aligned * 100 / len(alignments)
    result_file = open("resultsBowTie.txt", "w")
    result_file.write("BWA-MEM accuracy of alignment is: {}%.\n".format(results))


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
