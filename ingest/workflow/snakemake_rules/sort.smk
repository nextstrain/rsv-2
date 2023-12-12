"""
This part of the workflow handles sorting downloaded sequences and metadata
into a and b by aligning them to reference sequences.

It produces output files as

    metadata = "data/{type}/metadata.tsv"

    sequences = "data/{type}/sequences.fasta"

"""

TIME = ['1', '2','3']


rule sort:
    input:
        sequences = rules.transform.output.sequences
    output:
        "data/a/sequences.fasta",
        "data/b/sequences.fasta"
    shell:
        '''
        nextclade3 sort {input.sequences} --output-dir tmp
        seqkit rmdup tmp/nextstrain/rsv/b/sequences.fasta > data/b/sequences.fasta
        seqkit rmdup tmp/nextstrain/rsv/a/sequences.fasta > data/a/sequences.fasta
        rm -r tmp
        '''

rule metadata:
    input:
        metadata = rules.transform.output.metadata,
        sequences = "data/{type}/sequences.fasta"
    output:
        metadata = "data/{type}/metadata_raw.tsv"
    run:
        import pandas as pd
        from Bio import SeqIO

        strains = [s.id for s in SeqIO.parse(input.sequences, 'fasta')]
        d = pd.read_csv(input.metadata, sep='\t', index_col='accession').loc[strains].drop_duplicates()
        d.to_csv(output.metadata, sep='\t')

rule nextclade_dataset:
    output:
        ref_a = "rsv-a_nextclade/reference.fasta",
        ref_b = "rsv-b_nextclade/reference.fasta"
    params:
        dataset_a = "nextstrain/rsv/a/EPI_ISL_412866",
        dataset_b = "nextstrain/rsv/b/EPI_ISL_1653999"
    shell:
        """
        nextclade3 dataset get -n {params.dataset_a} --output-dir rsv-a_nextclade
        nextclade3 dataset get -n {params.dataset_b} --output-dir rsv-b_nextclade
        """

rule nextclade:
    input:
        sequences = "data/{type}/sequences.fasta",
        ref = "rsv-a_nextclade/reference.fasta"
    output:
        nextclade = "data/{type}/nextclade.tsv"
    params:
        dataset = "rsv-a_nextclade",
        output_columns = "clade qc.overallScore qc.overallStatus alignmentScore  alignmentStart  alignmentEnd  coverage dynamic"
    shell:
        """
        nextclade3 run -D {params.dataset}  \
                          --output-columns-selection {params.output_columns} \
                          --output-tsv {output.nextclade} \
                          {input.sequences}
        """

rule extend_metadata:
    input:
        nextclade = "data/{type}/nextclade.tsv",
        metadata = "data/{type}/metadata_raw.tsv"
    output:
        metadata = "data/{type}/metadata.tsv"
    shell:
        """
        python3 bin/extend-metadata.py --metadata {input.metadata} \
                                       --id-field accession \
                                       --virus-type {wildcards.type} \
                                       --nextclade {input.nextclade} \
                                       --output {output.metadata}
        """


# rule coverage:
#     input:
#         alignment_a = expand("data/a/{time}_sequences.aligned.fasta", time=TIME),
#         alignment_b = expand("data/b/{time}_sequences.aligned.fasta", time=TIME),
#         metadata_b = expand("data/b/{time}_metadata.tsv", time=TIME),
#         metadata_a = expand("data/a/{time}_metadata.tsv", time=TIME),
#         dedup_metadata_a = rules.deduplication.output.dedup_metadata_a,
#         dedup_metadata_b = rules.deduplication.output.dedup_metadata_b
#     output:
#         metadata_a = "data/a/metadata.tsv",
#         metadata_b = "data/b/metadata.tsv"
#     shell:
#         """
#         python bin/gene-coverage.py
#         """