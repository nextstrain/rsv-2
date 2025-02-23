configfile: "config/configfile.yaml"

wildcard_constraints:
    a_or_b = r"a|b"

build_dir = 'results'
auspice_dir = 'auspice'

rule all:
    input:
        expand("auspice/rsv_{subtype}_{build}_{resolution}.json",
               subtype = config.get("subtypes",['a']),
               build = config.get("builds_to_run", ['genome']),
               resolution = config.get("resolutions_to_run", ["all-time"])),
        expand("auspice/rsv_{subtype}_{build}_{resolution}_tip-frequencies.json",
               subtype = config.get("subtypes",['a']),
               build = config.get("builds_to_run", ['genome']),
               resolution = config.get("resolutions_to_run", ["all-time"])),

include: "workflow/snakemake_rules/chores.smk"

include: "workflow/snakemake_rules/core.smk"

include: "workflow/snakemake_rules/export.smk"

include: "workflow/snakemake_rules/download.smk"

include: "workflow/snakemake_rules/glycosylation.smk"

include: "workflow/snakemake_rules/clades.smk"


if "deploy_url" in config:
    include: "workflow/snakemake_rules/nextstrain_automation.smk"

rule clean:
    params:
        targets = ["auspice", "results"]
    shell:
        """
        rm -rf {params.targets}
        """

rule clobber:
    params:
        targets = ["data", "auspice", "results"]
    shell:
        """
        rm -rf {params.targets}
        rm config/clades*tsv
        """
