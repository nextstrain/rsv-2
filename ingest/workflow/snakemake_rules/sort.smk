"""
This part of the workflow handles sorting downloaded sequences and metadata
into a and b by aligning them to reference sequences.

It produces output files as

    metadata = "data/{type}/metadata.tsv"

    sequences = "data/{type}/sequences.fasta"

"""

TIME = ['1', '2','3']

rule align:
    input:
        sequences = rules.transform.output.sequences,
        reference = "config/{type}_{time}_reference.fasta"
    output:
        alignment = "data/{type}/{time}_sequences.aligned.fasta"
    threads: 2
    shell:
        """
        nextalign run -j {threads} --silent \
        --reference {input.reference} \
        --output-fasta {output.alignment} \
        {input.sequences}
        """

rule metadataandsequences:
    input:
        alignment = rules.align.output.alignment,
        metadata = rules.transform.output.metadata,
        sequences = rules.transform.output.sequences
    output:
        metadata = "data/{type}/{time}_metadata.tsv",
        sequences = "data/{type}/{time}_sequences.fasta"

    shell:
        """
        python bin/sequencesandmetadata.py \
        --sortedalignment {input.alignment} \
        --allmetadata {input.metadata} \
        --allsequences {input.sequences} \
        --metadata {output.metadata} \
        --sequences {output.sequences}
        """

rule sort:
    input:
        allsequences = "data/sequences.fasta",
        metadata = "data/metadata.tsv",
        alignment_a = expand("data/a/{time}_sequences.aligned.fasta", time=TIME),
        alignment_b = expand("data/b/{time}_sequences.aligned.fasta", time=TIME),
        reference_a = expand("config/a_{time}_reference.fasta", time=TIME),
        reference_b = expand("config/b_{time}_reference.fasta", time=TIME),
        metadata_b = expand("data/b/{time}_metadata.tsv", time=TIME),
        metadata_a = expand("data/a/{time}_metadata.tsv", time=TIME)
    output:
        sequences_a = "data/a/sequences.fasta",
        metadata_a = "data/a/metadata_sorted.tsv",
        sequences_b = "data/b/sequences.fasta",
        metadata_b = "data/b/metadata_sorted.tsv"
    shell:
        """
        python bin/sort.py
        """

rule coverage:
    input:
        alignment_a = expand("data/a/{time}_sequences.aligned.fasta", time=TIME),
        alignment_b = expand("data/b/{time}_sequences.aligned.fasta", time=TIME),
        metadata_b = expand("data/b/{time}_metadata.tsv", time=TIME),
        metadata_a = expand("data/a/{time}_metadata.tsv", time=TIME),
        sorted_metadata_a = "data/a/metadata_sorted.tsv",
        sorted_metadata_b = "data/b/metadata_sorted.tsv"
    output:
        metadata_a = "data/a/metadata.tsv",
        metadata_b = "data/b/metadata.tsv"
    shell:
        """
        python bin/gene-coverage.py
        """