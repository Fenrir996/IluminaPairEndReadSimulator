from main import create_reverse_complement_genome, change_into_base, create_qualities_by_normal_distribution


def test_methods():
    # Test create_reverse_complement_genome
    reverse = create_reverse_complement_genome("ATCG")
    assert reverse == "CGAT"

    # Test change_into_base
    nucleotide = change_into_base('N')
    assert nucleotide == 'A' or nucleotide == 'T' or nucleotide == 'C' or nucleotide == 'G'

    # Test create_qualities_by_normal_distribution
    qualities = create_qualities_by_normal_distribution(20, 80, 1)
    for quality in qualities:
        assert 77 <= ord(quality) <= 83


def main():
    test_methods()


if __name__ == "__main__":
    main()
