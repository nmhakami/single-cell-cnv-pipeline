import pysam


OUTPUT = "tests/data/test.bam"

header = {
    "HD": {"VN": "1.6"},
    "SQ": [{"SN": "chr1", "LN": 5000}],
}

# Starting positions of our fake reads.
positions = [
    # Bin 1: 5 reads
    100, 200, 300, 400, 500,

    # Bin 2: 5 reads total,
    # but positions 1450 and 1550 should be excluded.
    1100, 1200, 1300, 1450, 1550,

    # Bin 3: 4 reads
    2100, 2200, 2300, 2400,
]


with pysam.AlignmentFile(OUTPUT, "wb", header=header) as bam:
    for i, position in enumerate(positions):
        read = pysam.AlignedSegment()

        read.query_name = f"read_{i}"
        read.query_sequence = "A" * 100

        read.flag = 0
        read.reference_id = 0
        read.reference_start = position
        read.mapping_quality = 60
        read.cigar = ((0, 100),)

        read.query_qualities = pysam.qualitystring_to_array("I" * 100)

        bam.write(read)


print(f"Created {OUTPUT}")