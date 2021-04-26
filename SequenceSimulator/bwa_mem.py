import subprocess
from main import sequence_simulator

# Takes in our sam file and Bowtie's sam file and does a comparison, gives result in percent
def compare_sam_files(nmd_sam_file, bwa_sam_file):
    alignments = {}
    our_sam = open(nmd_sam_file, "r")
    correctly_aligned = calculate_correctly_aligned(alignments, our_sam, bwa_sam_file)
    return correctly_aligned * 100 / len(alignments)

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


# add the name of refGenome file - without extension!

def simulate_bwa(file_name, delete_error_rate, insert_error_rate, snv_error_rate):
    result_file = open("results.txt", "w")
    sequence_simulator("{}.fa".format(file_name), 70, 4, 150, 500, delete_error_rate, insert_error_rate, snv_error_rate)

    subprocess.run(["bwa", "index", "{}.fa".format(file_name)])
    subprocess.run(["bwa", "mem", "{}.fa".format(file_name), "{}_read1.fastq".format(file_name),
                    "{}_read2.fastq".format(file_name)], stdout=open("{}_bwa.sam".format(file_name), "w"))

    result = compare_sam_files("{}.sam".format(file_name), "{}_bwa.sam".format(file_name))
    result_file.write(
        "Error rate for DELETE: {}, error rate for INSERT: {}, error rate for SNV: {}, "
        "BWA-MEM accuracy of alignment is: {}%.\n".format(delete_error_rate, insert_error_rate, snv_error_rate, result))
    return result

# simulate_bwa("Test")
