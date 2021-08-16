import os
import pandas as pd

LINEAGES = []
LINEAGES_df = pd.read_csv((os.getcwd()+"/"+(config["lineage_splits"])))
for i,row in LINEAGES_df.iterrows():
    LINEAGES.append(row["lineage"])


rule merge_and_create_new_redcap_lineages:
    input:
        rules.add_lin_to_annotations.output.traits
    params:
        script = os.path.join(workflow.current_basedir, "../utilities/curate_lineages.py"),
        country = config["country"],
        country_code = config["country_code"]
    output:
        config["output_path"] + "/5/updated_traits.csv"
    log:
        config["output_path"] + "/logs/5_merge_and_create_new_redcap_lineages.log"
    resources: 
        mem_per_cpu=10000
    shell:
        """
        python {params.script} {input} {params.country} {params.country_code} {output} &> {log}
        """

#issue with r function 'normalizePath'
#skipping for now
#rule step_5_generate_sankey_plot:
#    input:
#        old_traits = config["output_path"] + "/4/all_traits.csv",
#        new_traits = config["output_path"] + "/5/updated_traits.csv",
#    output:
#        links = config["output_path"] + "/5/sankey_links.txt",
#        plot = config["output_path"] + "/5/sankey.html"
#    params:
#        python_script = os.path.join(workflow.current_basedir, "../utilities/get_sankey_links.py"),
#        R_script = os.path.join(workflow.current_basedir, "../utilities/plot_sankey.R")
#    log:
#        config["output_path"] + "/logs/5_generate_sankey_plot.log"
#    resources: mem_per_cpu=10000
#    shell:
#        """
#        python {params.python_script} {input.old_traits} {input.new_traits} {output.links} &> {log}
#        Rscript {params.R_script} {output.links} {output.plot} &>> {log}
#        """


#currently no global lineages
#new_global_lineages is currently empty
#skipping for now
#rule five_update_global_lineage_metadata:
#    input:
#        metadata = config["output_path"] + "/3/cog_gisaid.lineages.expanded.csv",
#        global_lineages = config["global_lineages"],
#        new_global_lineages = config["output_path"] + "/2/normal_pangolin/lineage_report.csv"
#    output:
#        metadata_temp = temp(config["output_path"] + "/5/cog_gisaid.global.lineages.with_all_traits.temp.csv"),
#        metadata = config["output_path"] + "/5/cog_gisaid.global.lineages.with_all_traits.csv"
#    log:
#        config["output_path"] + "/logs/5_five_update_global_lineage_metadata.log"
#    shell:
#        """
#        fastafunk add_columns \
#          --in-metadata {input.metadata} \
#          --in-data {input.global_lineages} \
#          --index-column strain \
#          --join-on taxon \
#          --new-columns lineage lineage_support lineages_version \
#          --where-column lineage_support=UFbootstrap \
#          --out-metadata {output.metadata_temp} &> {log}
#
#        fastafunk add_columns \
#          --in-metadata {output.metadata_temp} \
#          --in-data {input.new_global_lineages} \
#          --index-column strain \
#          --join-on taxon \
#          --new-columns lineage lineage_support lineages_version \
#          --where-column lineage_support=UFbootstrap \
#          --out-metadata {output.metadata} &>> {log}
#        """


rule update_lineage_metadata:
    input:
        metadata = config["output_path"] + "/3/redcap_gisaid.lineages.expanded.csv",
#        metadata = rules.five_update_global_lineage_metadata.output.metadata,
        traits = rules.output_annotations.output.traits,
        updated_lineages = rules.merge_and_create_new_redcap_lineages.output
    params:
        country_code = config["country_code"]
    output:
        traits_metadata = temp(config["output_path"] + "/5/redcap_gisaid.lineages.with_traits.csv"),
        all_metadata = config["output_path"] + "/5/redcap_gisaid.lineages.with_all_traits.csv"
    log:
        config["output_path"] + "/logs/5_update_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.traits} \
          --index-column strain \
          --join-on taxon \
          --new-columns del_lineage del_introduction \
          --out-metadata {output.traits_metadata} &> {log} ;

        fastafunk add_columns \
          --in-metadata {output.traits_metadata} \
          --in-data {input.updated_lineages} \
          --index-column strain \
          --join-on taxon \
          --new-columns {params.country_code}_lineage microreact_lineage \
          --out-metadata {output.all_metadata} &> {log}
        """

# rule run_5_subroutine_on_lineages:
#     input:
#         metadata = rules.update_lineage_metadata.output.all_metadata,
#     params:
#         path_to_script = workflow.current_basedir,
#         output_path = config["output_path"],
#         publish_path = config["publish_path"],
#         export_path = config["export_path"],
#     output:
#         metadata = config["output_path"] + "/5/cog_gisaid.lineages.with_all_traits.with_phylotype_traits.csv",
#         full_tree = config["output_path"] + "/5/cog_gisaid_full.tree.nexus",
#     log:
#         config["output_path"] + "/logs/5_run_5_subroutine_on_lineages.log"
#     threads: 40
#     shell:
#         """
#         snakemake --nolock \
#           --snakefile {params.path_to_script}/5_subroutine/5_process_lineages.smk \
#           --cores {threads} \
#           --configfile {params.path_to_script}/5_subroutine/config.yaml \
#           --config \
#           output_path={params.output_path} \
#           publish_path={params.publish_path} \
#           export_path={params.export_path} \
#           metadata={input.metadata} &> {log}
#         """

################################################################################
rule step_5_annotate_tree:
    input:
        tree = rules.merge_sibling_del_introduction.output.tree,
        metadata = rules.update_lineage_metadata.output.all_metadata
    params:
        country_code = config["country_code"]
    output:
        tree=config["output_path"] + "/5/redcap_gisaid_grafted.annotated.tree"
    log:
        config["output_path"] + "/logs/5_annotate.log"
    resources: 
        mem_per_cpu=20000
    shell:
        """
        clusterfunk annotate_tips \
          --in-metadata {input.metadata} \
          --trait-columns {params.country_code}_lineage \
          --index-column strain \
          --input {input.tree} \
          --output {output.tree} &> {log}
        """


#need to understand what this is doing
#why -ge 3?
#          I=`echo ${{LINE}} | cut -d"K" -f2`
rule get_redcap_lineage_samples:
    input:
        metadata = rules.merge_and_create_new_redcap_lineages.output
    output:
        outdir = directory(config["output_path"] + "/5/samples/")
    log:
        config["output_path"] + "/logs/5_get_redcap_lineage_samples.log"
    shell:
        """
        mkdir -p {output.outdir} &> {log}

        tail -n+2 {input.metadata} | cut -d, -f2 | sort | uniq | while read LINE
        do
          if [ $(grep -c ",${{LINE}}," {input.metadata}) -ge 3 ]
          then
            grep ",${{LINE}}," {input.metadata} | cut -d"," -f1 > {output.outdir}/${{LINE}}.samples.txt
          fi
        done 2>> {log}
        """


rule dequote_tree:
    input:
        full_newick_tree = rules.sort_collapse.output.sorted_collapsed_tree
    output:
        tree = config["output_path"] + "/5/redcap_gisaid_full.tree.noquotes.newick"
    log:
        config["output_path"] + "/logs/5_dequote_tree.log"
    resources: 
        mem_per_cpu=10000
    shell:
        """
        sed "s/'//g" {input.full_newick_tree} > {output.tree} 2> {log}
        """


#country info not in tips but that's because this pipeline was originally dealing with all UK countries
#PH has only itself, so country info is unnecessary
#wrong, some redcap records have china as country
rule cut_out_trees:
    input:
        full_tree = rules.dequote_tree.output.tree,
        indir = config["output_path"] + "/5/samples"
    params:
        country_code = config["country_code"]
    output:
        outdir = directory(config["output_path"] + "/5/trees/")
    log:
        config["output_path"] + "/logs/5_cut_out_trees.log"
    resources: 
        mem_per_cpu=10000
    shell:
        """
        mkdir -p {output.outdir} 2> {log}

        for FILE in {input.indir}/{params.country_code}*.samples.txt
        do
            LIN=`echo ${{FILE}} | rev | cut -d"/" -f1 | rev | cut -d"." -f1`
            gotree prune -r \
                -i {input.full_tree} \
                --tipfile ${{FILE}} \
                -o {output.outdir}/${{LIN}}.tree
        done 2>> {log}
        """


#collapse parameter not used?
rule phylotype_cut_trees:
    input:
        treedir = config["output_path"] + "/5/trees/"
    output:
        phylotypedir = directory(config["output_path"] + "/5/phylotyped_trees/")
    params:
        collapse=5E-6,
        threshold=2E-5
    log:
        config["output_path"] + "/logs/5_phylotype_cut_trees.log"
    resources: 
        mem_per_cpu=20000
    shell:
        """
        mkdir -p {output.phylotypedir} 2> {log}

        for FILE in {input.treedir}/*.tree
        do
            LIN=`echo ${{FILE}} | rev | cut -d"/" -f1 | rev | cut -d"." -f1`
            clusterfunk phylotype \
                --threshold {params.threshold} \
                --prefix ${{LIN}}_1 \
                --input ${{FILE}} \
                --in-format newick \
                --output {output.phylotypedir}/${{LIN}}.tree
        done 2>> {log}
        """


rule get_redcap_phylotypes_csv:
    input:
        phylotypedir = config["output_path"] + "/5/phylotyped_trees/"
    output:
        csvdir = directory(config["output_path"] + "/5/phylotype_csvs/")
    log:
        config["output_path"] + "/logs/5_get_redcap_phylotypes_csv.log"
    resources: 
        mem_per_cpu=20000
    shell:
        """
        mkdir -p {output.csvdir} 2> {log}

        for FILE in {input.phylotypedir}/*.tree
        do
            LIN=`echo ${{FILE}} | rev | cut -d"/" -f1 | rev | cut -d"." -f1`
            clusterfunk extract_tip_annotations \
              --traits phylotype \
              --input ${{FILE}} \
              --output {output.csvdir}/${{LIN}}.csv
        done 2>> {log}
        """


# def aggregate_input_csv(wildcards):
#     checkpoint_output_directory = checkpoints.get_uk_lineage_samples.get(**wildcards).output[0]
#     required_files = expand( "%s/5/phylotyped_trees/uk_lineage_UK{i}.csv" %(config["output_path"]),
#                             i=glob_wildcards(os.path.join(checkpoint_output_directory, "UK{i}.samples.txt")).i)
#     return (required_files)
#
# def aggregate_input_trees(wildcards):
#     checkpoint_output_directory = checkpoints.cut_out_trees.get(**wildcards).output[0]
#     print(checkpoints.cut_out_trees.get(**wildcards).output[0])
#     lineage = wildcards.lineage
#     required_files = expand( "%s/5/%s/phylotyped_trees/uk_lineage_UK{i}.tree" %(config["output_path"],lineage),
#                             i=glob_wildcards(os.path.join(checkpoint_output_directory, "uk_lineage_UK{i}.tree")).i)
#     return (sorted(required_files))
#
# def aggregate_input_labels(wildcards):
#     checkpoint_output_directory = checkpoints.cut_out_trees.get(**wildcards).output[0]
#     print(checkpoints.cut_out_trees.get(**wildcards).output[0])
#     labels = expand( "UK{i}",i=glob_wildcards(os.path.join(checkpoint_output_directory, "uk_lineage_UK{i}.tree")).i)
#     return (sorted(labels))


rule combine_phylotypes_csv:
    input:
        csvdir = config["output_path"] + "/5/phylotype_csvs/"
    output:
        phylotype_csv = config["output_path"] + "/5/redcap_phylotypes.csv"
    log:
        config["output_path"] + "/logs/5_traits_combine_phylotype_csv.log"
    resources: 
        mem_per_cpu=20000
    run:
        import pandas as pd
        import glob
        import os

        mypath = str(input.csvdir)

        #print(mypath)

        files = glob.glob(os.path.join(mypath, "*csv"))
        
        #print(files)

        dfs = [pd.read_csv(x) for x in files]

        #print(dfs)

        result = pd.concat(dfs)
        result.to_csv(output[0], index=False)


#adding country lineage column to metadata from rule 2
#for import to redcap
rule add_country_lineage_to_redcap_metadata:
    input:
        metadata = rules.redcap_add_pangolin_lineages_to_metadata.output.metadata,
        lineage = rules.merge_and_create_new_redcap_lineages.output
    params:
        country_code = config["country_code"]
    output:
        metadata = config["output_path"] + "/5/redcap_metadata.with_country_lineage.csv"
    log:
        config["output_path"] + "/logs/5_add_phylotypes_to_redcap_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.lineage} \
          --index-column strain \
          --join-on taxon \
          --new-columns {params.country_code}_cluster \
          --where-column {params.country_code}_cluster={params.country_code}_lineage\
          --out-metadata {output.metadata} &>> {log}
        """


rule get_country_lineage_for_analysis_instrument_and_import_to_redcap:
    input:
        metadata = rules.add_country_lineage_to_redcap_metadata.output.metadata,
        redcap_db = config["redcap_access"]
    output:
        metadata = config["output_path"] + "/5/analysis_instrument_phylotype.csv"
    run:
        import redcap
        import pandas as pd

        with open(input.redcap_db, 'r') as f:
            read_file = f.read().strip('\n')
            url = read_file.split(',')[0]
            key = read_file.split(',')[1]
        f.close()

        proj = redcap.Project(url, key)
        df = pd.read_csv(input.metadata)

        #select columns necessary for import
        df = df.loc[:,['central_id', 'redcap_repeat_instance', 'ph_cluster']]

        #add column necessary for import
        df.insert(1, 'redcap_repeat_instrument', 'Analysis')
        df.loc[:,'redcap_repeat_instrument'] = df.loc[:,'redcap_repeat_instrument'].str.casefold()

        proj.import_records(df)
        df.to_csv(output.metadata, index=False)


rule merge_with_metadata:
    input:
        metadata = rules.update_lineage_metadata.output.all_metadata,
        traits = rules.combine_phylotypes_csv.output.phylotype_csv
    output:
        metadata = config["output_path"] + "/5/redcap_gisaid.lineages.with_all_traits.with_phylotype_traits.csv"
    log:
        config["output_path"] + "/logs/5_merge_with_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.traits} \
          --index-column strain \
          --join-on taxon \
          --new-columns phylotype \
          --out-metadata {output.metadata} &> {log}
        """


rule annotate_phylotypes:
    input:
        tree=rules.step_5_annotate_tree.output.tree,
        metadata = rules.merge_with_metadata.output.metadata
    output:
        annotated_tree = config["output_path"] + "/5/redcap_gisaid_full.tree.nexus"
    log:
        config["output_path"] + "/logs/5_annotate_phylotypes.log"
    resources: 
        mem_per_cpu=20000
    shell:
        """
        clusterfunk annotate_tips \
          --in-metadata {input.metadata} \
          --trait-columns phylotype \
          --index-column strain \
          --input {input.tree} \
          --output {output.annotated_tree} &> {log}
        """

################################################################################



#        echo '{{"text":"' > {params.json_path}/5b_data.json
#        echo "*Step 5: Generate {params.date} UK lineage trees is complete*\\n" >> {params.json_path}/5b_data.json
#        echo "> _Look at this Sankey plot_: {input.sankey_plot}\\n" >> {params.json_path}/5b_data.json
#        echo '"}}' >> {params.json_path}/5b_data.json
#        echo "webhook {params.grapevine_webhook}"
#        curl -X POST -H "Content-type: application/json" -d @{params.json_path}/5b_data.json {params.grapevine_webhook}
rule summarize_define_redcap_lineages_and_cut_out_trees:
    input:
        annotated_tree = rules.annotate_phylotypes.output.annotated_tree,
#        sankey_plot = rules.step_5_generate_sankey_plot.output.plot,
        import_country_lineage = rules.get_country_lineage_for_analysis_instrument_and_import_to_redcap.output.metadata
    params:
#        grapevine_webhook = config["grapevine_webhook"],
#        json_path = config["json_path"],  
#        date = config["date"]
    log:
        config["output_path"] + "/logs/5_summarize_define_redcap_lineages_and_cut_out_trees.log"
    shell:
        """
        echo "5_subroutine complete" &>> {log}
        """
