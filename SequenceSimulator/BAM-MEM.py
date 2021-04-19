import subprocess
from main import sequence_simulator


def compare_sam_files(nmd_sam_file, bwa_sam_file):
    alignments = {}
    our_sam = open(nmd_sam_file, "r")
    correctly_aligned = calculate_correctly_aligned(alignments, our_sam, bwa_sam_file)
    return correctly_aligned * 100 / len(alignments)


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


def simulate_bwa():
    file_name = "69820_ref_ASM270686v1_chr12"  # add the name of refGenome file - without extension!
    delete_error_rate = int(0)
    insert_error_rate = int(0)
    snv_error_rate = int(0)
    result_file = open("results.txt", "w")
    sequence_simulator("{}.fa".format(file_name), 70, 4, 150, 500, delete_error_rate, insert_error_rate, snv_error_rate)

    subprocess.run(["bwa", "index", "{}.fa".format(file_name)])
    subprocess.run(["bwa", "mem", "{}.fa".format(file_name), "read1.fastq",
                    "read2.fastq"], stdout=open("{}_bwa.sam".format(file_name), "w"))

    result = compare_sam_files("{}.sam".format(file_name), "{}_bwa.sam".format(file_name))
    result_file.write(
        "Error rate for DELETE: {}, error rate for INSERT: {}, error rate for SNV: {}, "
        "BWA-MEM accuracy of alignment is: {}%.\n".format(delete_error_rate, insert_error_rate, snv_error_rate, result))


simulate_bwa()
