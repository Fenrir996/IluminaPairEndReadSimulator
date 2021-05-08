import random
import time

import numpy

NORMAL_DIST_SIGMA = 2
LEFT = 1
RIGHT = 2
ERROR = -1

# Qualities go from 33 to 126 in ascii so we have them mapped here. In case you want to change the scale
# for example go for 1 to 91, you're free to do that here
qualities_in_ascii = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
# Bases were at first only ATCG, but there was a need for expanding so we had to delete a frankenstein method
# which used to convert non base nucleotides to base nucleotides
bases = ['A', 'T', 'C', 'G', 'N']

# @read_id#read_direction(1 - right, 2 - left)
# read
# +
# qualities
fastq_entry = "@{}/{}\n{}\n+\n{}\n"

# @SQ SN:sequence_name LN:sequence_length
sam_QNAME = "@SQ SN:{} LN:{}\n"

# read_id read_position read read_quality
sam_data = "{}\t{}\t{}\t{}\n"


# FASTA format
# >gi|186681228|ref|YP_001864424.1| phycoerythrobilin:ferredoxin oxidoreductase -> name
# MNSERSDVTLYQPFLDYAIAYMRSRLDLEPYPIPTGFESNSAVVGKGKNQEEVVTTSYAFQTAKLRQIRA -> nucleotides
# AHVQGGNSLQVLNFVIFPHLNYDLPFFGADLVTLPGGHLIALDMQPLFRDDSAYQAKYTEPILPIFHAHQ
# QHLSWGGDFPEEAQPFFSPAFLWTRPQETAVVETQVFAAFKDYLKAYLDFVEQAEAVTDSQNLVAIKQAQ
# LRYLRYRAEKDPARGMFKRFYGAEWTEEYIHGFLFDLERKLTVVK

# Method which creates reverse complements for a given genome as string
# This mimics what polymerase does when creating a complement on a parallel strand
def create_reverse_complement_genome(genome):
    complementary_bases = {
        "A": "T",
        "C": "G",
        "T": "A",
        "G": "C",
        "N": "N"
    }
    result = ''
    for nucleotide in genome:
        result = complementary_bases[nucleotide] + result
    return result


# Method for reading fasta file, file needs to be in root of the project
def read_genome_from_fasta_file(file_name):
    file = open(file_name, 'r')
    lines = file.readlines()
    genome = ''
    genome_name = ''

    for line in lines:
        #  in case of multiple sequences
        if line[0] != '>':
            genome += line.strip()
        else:
            genome_name = line[1:].rstrip('\n').split()[0]

    file.close()
    return [genome.upper(), genome_name]


# sigma should probably be 1.0, mean 60-80
# Given the average quality for the whole genome we had to give each nucleotide in the genome its' own read quality
# this method does exactly that, creating an array of qualities where the index of the quality array fits the
# index of the nucleotide in genome array
def create_qualities_by_normal_distribution(length, mean, sigma):
    quality_values = numpy.random.normal(mean, sigma, length)
    qualities = ''
    for value in quality_values:
        int_val = int(value)
        if chr(int_val) in qualities_in_ascii:
            quality = int_val
        else:
            # if by some chance it drops under lowest value give it the lowest value
            if int(33) > int_val > int(0):
                quality = int(33)
            else:
                # if it goes above the max value, give it max value
                if int_val > int(126):
                    quality = int(126)
                else:
                    raise ValueError("Quality can't be negative, something bad happened with Normal distribution")

        qualities += chr(quality)

    return qualities


# keep error rates under 0.01 eg. 0.001, 0.005, 0.009
# This method is responsible for mutating the genome as the name states. However down the line we figured that
# even if we mutate the whole genome we could still take pieces untouched by the mutation. So, yes
# this method can be used to mutate the whole genome, but we use it to mutate reads actually to make sure we work
# with errors.
# reference_genome presents the nucleotides, while gen_read_data is a structure which holds the number of errors
# we have to make of each type
def mutate_genome(reference_genome, gen_read_data):
    genome_length = len(reference_genome)
    number_of_insertions = gen_read_data["add_errors"]
    number_of_deletions = gen_read_data["del_errors"]
    number_of_variations = gen_read_data["snv_errors"]

    if number_of_variations <= 0 or number_of_insertions <= 0 or number_of_deletions <= 0:
        return reference_genome

    gen_read_data["add_errors"] -= 1
    gen_read_data["del_errors"] -= 1
    gen_read_data["snv_errors"] -= 1

    # introduce single nucleotide variation - SNV
    # at a random position change the nucleotide into another random one, not the same one
    while number_of_variations > 0:
        variation_index = random.randint(0, genome_length - 1)
        base_index = random.randint(0, 100) % 4
        while base_index == bases.index(reference_genome[variation_index]):
            base_index = random.randint(0, 100) % 4
        reference_genome_pre = reference_genome[0:variation_index]
        reference_genome_post = reference_genome[variation_index + 1:genome_length]
        reference_genome = reference_genome_pre + bases[base_index] + reference_genome_post
        number_of_variations -= 1
    # introduce insertion of singular nucleotides
    # at a random position insert a random nucleotide
    while number_of_insertions > 0:
        insertion_index = random.randint(0, genome_length - 1)
        base_index = random.randint(0, 100) % 4
        reference_genome_pre = reference_genome[0:insertion_index]
        reference_genome_post = reference_genome[insertion_index:genome_length]
        reference_genome = reference_genome_pre + bases[base_index] + reference_genome_post
        genome_length += 1
        number_of_insertions -= 1

    # introduce deletion of singular nucleotides
    # at a random position delete the specified nucleotide
    while number_of_deletions > 0:
        deletion_index = random.randint(0, genome_length - 1)
        reference_genome_pre = reference_genome[0:deletion_index]
        reference_genome_post = reference_genome[deletion_index + 1:genome_length]
        reference_genome = reference_genome_pre + reference_genome_post
        genome_length -= 1
        number_of_deletions -= 1

    return reference_genome


# Utility checks
def check_positive_value(value):
    if value < 0:
        return 0
    return 1


def check_coverage(coverage):
    if coverage < 0:
        return 0
    return 1


def check_errors_values(error):
    if error < 0 or error > 1:
        return 0
    return 1


def check_quality(average_quality):
    if average_quality < int(33) or average_quality > int(126):
        return 0
    return 1


# the simulator
# file is the name of the file in the root of the project
# average_quality is the average read quality of the whole genome, so a number
# coverage, read_size, insert_size, also numbers, where read_size, has to be larger than insert_size
# error rates are between 0 and 1 which later get transferred into percent

def sequence_simulator(file, average_quality, coverage, read_size, insert_size, delete_error_rate,
                       insert_error_rate, snv_error_rate):
    start_time = time.time()

    # block of checks
    if not file:
        print('Reference genome mast be defined!')
        return ERROR

    if check_quality(average_quality) == 0:
        print('Quality mast be between 33 and 126')
        return ERROR

    if check_coverage(coverage) == 0:
        print('Coverage mast be positive number')
        return ERROR

    if check_positive_value(read_size) == 0:
        print('Read size mast be positive number')
        return ERROR

    if check_positive_value(insert_size) == 0:
        print('Insert size mast be positive number')
        return ERROR

    if read_size > insert_size:
        print('Read size can not be longer than insert size (size of the fragment of the genome)')
        return ERROR

    if check_errors_values(delete_error_rate) == 0:
        print('Delete error mast be value between 0 and 1')
        return ERROR

    if check_errors_values(insert_error_rate) == 0:
        print('Insert error mast be value between 0 and 1')
        return ERROR

    if check_errors_values(snv_error_rate) == 0:
        print('SNV error mast be value between 0 and 1')
        return ERROR

    if delete_error_rate + insert_error_rate + snv_error_rate > 1:
        print('Sum of error can not be greater than 1')
        return ERROR

    # read the genome from fasta
    referent_genome = read_genome_from_fasta_file(file)

    gen_read_data = {
        "ref_genome": referent_genome[0],
        "ref_genome_name": referent_genome[1],
        "quality": average_quality,
        "read_size": read_size,
        "insert_size": insert_size,
        "snv_errors": int(snv_error_rate * len(referent_genome[0])),
        "add_errors": int(insert_error_rate * len(referent_genome[0])),
        "del_errors": int(delete_error_rate * len(referent_genome[0])),
        "num_of_reads": int(coverage * len(referent_genome[0]) / read_size),
        "file_name": file.split(".")[0],
        "ref_genome_size": len(referent_genome[0])
    }

    # Method which does all the reading and writing
    generate_reads(gen_read_data)

    end_time = time.time()
    print("Sequencing finished. Time elapsed: {}".format(end_time - start_time))


# Utility method which generates one read
# It takes the error information via gen_read_data, position of read start and end
# And the direction of reading which is backwards for read2
def generate_read(gen_read_data, read_start, read_end, direction):
    read = gen_read_data["ref_genome"][read_start:read_end]
    read = mutate_genome(read, gen_read_data)
    # For read2 direction is right
    if direction == RIGHT:
        read = create_reverse_complement_genome(read)
    return read


# Main method of the simulator
def generate_reads(gen_read_data):
    # Creates fastq files and a SAM file
    fastq1 = open("{}_read1.fastq".format(gen_read_data["file_name"]), "w")
    fastq2 = open("{}_read2.fastq".format(gen_read_data["file_name"]), "w")
    sam = open("{}.sam".format(gen_read_data["file_name"]), "w")

    sam.write(sam_QNAME.format(gen_read_data["ref_genome_name"], gen_read_data["ref_genome_size"]))
    for i in range(gen_read_data["num_of_reads"]):
        read_id = gen_read_data["ref_genome_name"] + "_" + str(i)

        # Pick a random place to start reading read1
        read1_start = random.randint(0, gen_read_data["ref_genome_size"] - gen_read_data["insert_size"])
        read1_finish = read1_start + gen_read_data["read_size"]
        # Generate the read
        read1 = generate_read(gen_read_data, read1_start, read1_finish, LEFT)
        # Generate the qualities, since there's no reason to create qualities for the whole genome
        read1_qualities = create_qualities_by_normal_distribution(len(read1), gen_read_data["quality"],
                                                                  NORMAL_DIST_SIGMA)
        fastq1.write(fastq_entry.format(read_id, LEFT, read1, read1_qualities))

        # Same for read2
        read2_end = read1_start + gen_read_data["insert_size"]
        read2_start = read2_end - gen_read_data["read_size"]
        read2 = generate_read(gen_read_data, read2_start, read2_end, RIGHT)
        read2_qualities = create_qualities_by_normal_distribution(len(read2), gen_read_data["quality"],
                                                                  NORMAL_DIST_SIGMA)
        fastq2.write(fastq_entry.format(read_id, RIGHT, read2, read2_qualities))

        orig_read2 = create_reverse_complement_genome(read2)
        #  Generate SAM file
        sam.write(sam_data.format(read_id, read1_start + 1, read1, read1_qualities))
        sam.write(sam_data.format(read_id, read2_start + 1, orig_read2, read2_qualities))

    fastq1.close()
    fastq2.close()
    sam.close()
