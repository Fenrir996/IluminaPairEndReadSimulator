import random
import numpy

qualities_in_ascii = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
bases = ['A', 'T', 'C', 'G']

#FASTA format
#>gi|186681228|ref|YP_001864424.1| phycoerythrobilin:ferredoxin oxidoreductase -> name
#MNSERSDVTLYQPFLDYAIAYMRSRLDLEPYPIPTGFESNSAVVGKGKNQEEVVTTSYAFQTAKLRQIRA -> nucleotides
#AHVQGGNSLQVLNFVIFPHLNYDLPFFGADLVTLPGGHLIALDMQPLFRDDSAYQAKYTEPILPIFHAHQ
#QHLSWGGDFPEEAQPFFSPAFLWTRPQETAVVETQVFAAFKDYLKAYLDFVEQAEAVTDSQNLVAIKQAQ
#LRYLRYRAEKDPARGMFKRFYGAEWTEEYIHGFLFDLERKLTVVK

def change_into_base(nucleotide):
    number = random.randint(0, 100)

    if nucleotide == 'R':
        rng = number % 2
        if rng == 0:
            return 'G'
        else:
            return 'A'

    if nucleotide == 'Y':
        rng = number % 2
        if rng == 0:
            return 'T'
        else:
            return 'C'

    if nucleotide == 'K':
        rng = number % 2
        if rng == 0:
            return 'G'
        else:
            return 'T'

    if nucleotide == 'M':
        rng = number % 2
        if rng == 0:
            return 'C'
        else:
            return 'A'

    if nucleotide == 'S':
        rng = number % 2
        if rng == 0:
            return 'G'
        else:
            return 'C'

    if nucleotide == 'W':
        rng = number % 2
        if rng == 0:
            return 'T'
        else:
            return 'A'

    if nucleotide == 'B':
        rng = number % 3
        if rng == 0:
            return 'G'
        elif rng == 1:
            return 'T'
        else:
            return 'C'

    if nucleotide == 'D':
        rng = number % 3
        if rng == 0:
            return 'G'
        elif rng == 1:
            return 'T'
        else:
            return 'A'

    if nucleotide == 'H':
        rng = number % 3
        if rng == 0:
            return 'A'
        elif rng == 1:
            return 'T'
        else:
            return 'C'

    if nucleotide == 'V':
        rng = number % 3
        if rng == 0:
            return 'G'
        elif rng == 1:
            return 'A'
        else:
            return 'C'

    if nucleotide == 'N':
        rng = number % 4
        if rng == 0:
            return 'G'
        elif rng == 1:
            return 'T'
        elif rng == 2:
            return 'C'
        else:
            return 'A'

def create_reverse_complement_genome(genome):
    complementary_bases = {
        "A": "T",
        "C": "G",
        "T": "A",
        "G": "C"
    }
    result = ''
    for nucleotide in genome:
        result = complementary_bases[nucleotide] + result
    return result

def read_genome_from_fasta_file(file_name):
    file = open(file_name, 'r')
    lines = file.readlines()
    genome = ''

    for line in lines:
        #  in case of multiple sequences
        if line[0] != '>':
            genome += line.strip()

    file.close()
    # need to make a base nucleotide if some of them aren't
    genome_length = len(genome)
    new_genome = ''
    for index in range(genome_length):
        if genome[index] not in bases:
            pre_genome = genome[0:index]
            post_genome = genome[index+1:genome_length]
            new_base = change_into_base(genome[index])
            genome = pre_genome + new_base + post_genome

    return genome.upper()

# sigma should probably be 1.0, mean 60-80
def create_qualities_by_normal_distribution(length, mean, sigma):
    quality_values = numpy.random.normal(mean, sigma, length)
    qualities = ''
    for value in quality_values:
        int_val = int(value)
        if chr(int_val) in qualities_in_ascii:
            quality = int_val
        else:
            if ord('!') > int_val > 0:
                quality = ord('!')
            else:
                if int_val > ord('~'):
                    quality = ord('~')
                else:
                    raise ValueError("Quality can't be negative, something bad happened with Normal distribution")

        qualities += chr(quality)

    return qualities


# keep error rates under 0.01 eg. 0.001, 0.005, 0.009


def mutate_genome(reference_genome, insert_error_rate, delete_error_rate, snv_error_rate):
    genome_length = len(reference_genome)
    number_of_insertions = round(insert_error_rate * genome_length)
    number_of_deletions = round(delete_error_rate * genome_length)
    number_of_variations = round(snv_error_rate * genome_length)

    # introduce single nucleotide variation - SNV

    while number_of_variations > 0:
        variation_index = random.randint(0, genome_length)
        print(variation_index)
        base_index = random.randint(0, 100) % 4
        while base_index == bases.index(reference_genome[variation_index]):
            base_index = random.randint(0, 100) % 4
        reference_genome_pre = reference_genome[0:variation_index]
        reference_genome_post = reference_genome[variation_index+1:genome_length]
        reference_genome = reference_genome_pre + bases[base_index] + reference_genome_post
        number_of_variations -= 1
    # introduce insertion of singular nucleotides

    while number_of_insertions > 0:
        insertion_index = random.randint(0, genome_length)
        print(insertion_index)
        base_index = random.randint(0, 100) % 4
        reference_genome_pre = reference_genome[0:insertion_index]
        reference_genome_post = reference_genome[insertion_index:genome_length]
        reference_genome = reference_genome_pre + bases[base_index] + reference_genome_post
        genome_length += 1
        number_of_insertions -= 1

    # introduce deletion of singular nucleotides
    while number_of_deletions > 0:
        deletion_index = random.randint(0, genome_length)
        print(deletion_index)
        reference_genome_pre = reference_genome[0:deletion_index]
        reference_genome_post = reference_genome[deletion_index+1:genome_length]
        reference_genome = reference_genome_pre + reference_genome_post
        genome_length -= 1
        number_of_deletions -= 1

    return reference_genome

