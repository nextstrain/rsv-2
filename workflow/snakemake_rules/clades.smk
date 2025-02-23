rule clades_genome:
    message:
        "adding clades based on the entire genome"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = "config/clades_genome_{a_or_b}.tsv"
    output:
        node_data = build_dir + "/{a_or_b}/{build_name}/{resolution}/clades_genome.json"
    log:
        "logs/{a_or_b}/clades_genome_{build_name}_{resolution}.txt"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --membership-name genome_clade \
            --label-name genome_clade \
            --output-node-data  {output.node_data} 2>&1 | tee {log}
        """


rule clades_Goya:
    message: "Adding internal clade labels"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = "config/clades_G_{a_or_b}.tsv"
    output:
        node_data = build_dir + "/{a_or_b}/{build_name}/{resolution}/clades_G.json"
    log:
        "logs/{a_or_b}/clades_{build_name}_{resolution}.txt"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --membership-name G_clade \
            --label-name G_clade \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """

rule clades_consortium:
    message: "Adding internal clade labels"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = "config/clades_consortium_{a_or_b}.tsv"
    output:
        node_data = build_dir + "/{a_or_b}/{build_name}/{resolution}/clades_consortium.json"
    log:
        "logs/{a_or_b}/clades_{build_name}_{resolution}.txt"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.node_data} 2>&1 | tee {log}
        """

rule download_clades:
    output:
        clades = "config/clades_consortium_{a_or_b}.tsv"
    params:
        url = lambda w: f"https://raw.githubusercontent.com/rsv-lineages/lineage-designation-{w.a_or_b.upper()}/main/.auto-generated/lineages.tsv"
    shell:
        """
        curl {params.url} --output {output.clades}
        """
