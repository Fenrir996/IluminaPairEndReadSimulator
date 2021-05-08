from main import create_reverse_complement_genome, create_qualities_by_normal_distribution, \
    read_genome_from_fasta_file, sequence_simulator

ERROR = -1


# Testing our utility methods
def test_methods():
    # Test create_reverse_complement_genome
    reverse = create_reverse_complement_genome("ATCG")
    assert reverse == "CGAT"

    # create_qualities_by_normal_distribution
    qualities = create_qualities_by_normal_distribution(20, 80, 1)
    for quality in qualities:
        assert 77 <= ord(quality) <= 83

    # Test function for reading genomes
    read = read_genome_from_fasta_file('SampleGenome.fa')
    assert read[0] == 'ATGCATAACGCGTTA'
    assert read[1] == 'genome1'

    # Test simulator reference genome
    code = sequence_simulator("", 70, 4, 150, 500, 0.001, 0, 0)
    assert code == ERROR

    # Test average quality
    code = sequence_simulator("SampleGenome.fa", 25, -1, 25, 5, 0.5, 0.5, 0.5)
    assert code == ERROR

    # Test coverage
    code = sequence_simulator("SampleGenome.fa", 70, -1, 25, 5, 0.5, 0.5, 0.5)
    assert code == ERROR

    # Test read size is positive
    code = sequence_simulator("SampleGenome.fa", 70, 5, -1, 5, 0.5, 0.5, 0.5)
    assert code == ERROR

    # Test insert size is positive
    code = sequence_simulator("SampleGenome.fa", 70, 5, 5, -2, 0.5, 0.5, 0.5)
    assert code == ERROR

    # Test read size > insert size
    code = sequence_simulator("SampleGenome.fa", 70, 5, 10, 5, 0.5, 0.5, 0.5)
    assert code == ERROR

    # Test delete error value
    code = sequence_simulator("SampleGenome.fa", 70, 5, 2, 5, 5, 0.5, 0.5)
    assert code == ERROR

    # Test insert error value
    code = sequence_simulator("SampleGenome.fa", 70, 5, 2, 5, 0.1, -1, 0.5)
    assert code == ERROR

    # Test snv error value
    code = sequence_simulator("SampleGenome.fa", 70, 5, 2, 5, 0.1, 0.8, 6)
    assert code == ERROR

    # Test errors sum value
    code = sequence_simulator("SampleGenome.fa", 70, 5, 2, 5, 0.1, 0.8, 0.6)
    assert code == ERROR


def main():
    test_methods()


if __name__ == "__main__":
    main()
