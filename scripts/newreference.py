from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, Seq
import shutil
import argparse
import sys

def new_reference(referencefile, outgenbank, outfasta, gene, start, end):
    ref = SeqIO.read(referencefile, "genbank")
    startofgene = None
    endofgene = None
    for feature in ref.features:
        if feature.type == 'source':
            ref_source_feature = feature
        if feature.type =='gene' or feature.type == 'CDS':
            a = list(feature.qualifiers.items())[0][-1][0]
            if a == gene:
                startofgene = int(list(feature.location)[0])
                endofgene =  int(list(feature.location)[-1])+1

    # If user provides a --gene 'some name' that is not found, error out as this may indicate that
    # the gene name is misspelled or the user may be using the wrong GenBank file.
    if(gene is not None and startofgene is None and endofgene is None):
        print(f"ERROR: No '{gene}' was found under 'gene' or 'CDS' features in the GenBank file.", file=sys.stderr)
        sys.exit(1)

    record = ref[startofgene:endofgene]
    # Allows for subgenic regions
    if (gene is not None and start is not None and end is not None):
        record = record[int(start):int(end)]

    source_feature =  SeqFeature(FeatureLocation(start=0, end=len(record)), type='source',
                                 qualifiers=ref_source_feature.qualifiers)
    record.features.append(source_feature)

    # For subgenic regions, adds a CDS feature {gene}_{start}_{end}
    if(gene is not None and start is not None and end is not None):
        gene = gene + "_" + start + "_" + end
        record.features.append(SeqFeature(FeatureLocation(start=0, end=len(record)), type='CDS',
                                      qualifiers={'gene': gene}))

    SeqIO.write(record, outgenbank, 'genbank')
    SeqIO.write(record, outfasta, "fasta")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="make new reference depending on whether the entire genome or only part is to be used for the tree",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--reference", required=True, help="GenBank file with reference sequences")
    parser.add_argument("--output-fasta", required=True, help="GenBank new reference file")
    parser.add_argument("--output-genbank", required=True, help="GenBank new reference file")
    parser.add_argument("--gene", help="gene name or genome for entire genome")
    parser.add_argument("--start", help="Start 0-based position relative to the gene (requires --gene argument)")
    parser.add_argument("--end", help="End 0-based position relative to the gene (requires --gene argument)")
    args = parser.parse_args()

    if args.gene=='genome':
        # Check if start and end are specified here
        if (args.start is not None and args.end is not None):
            print(f"ERROR: --start '{args.start}' --end '{args.end}' is not supported for full genome.", file=sys.stderr)
            sys.exit(1)


        shutil.copy(args.reference, args.output_genbank)
        SeqIO.write(SeqIO.read(args.reference, 'genbank'), args.output_fasta, 'fasta')
    else:
        new_reference(args.reference, args.output_genbank, args.output_fasta, args.gene, args.start, args.end)

