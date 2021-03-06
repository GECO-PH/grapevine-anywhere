#previous_stage defined but never used, commenting out for now
#only determines pangolin lineage for sequences which pass filters
#AND
#don't already have a pangolin lineage
#but might we want to allow for lineage re-assignment?
#lineage re-assignment would probably be done infrequently enough that
#it could just be done by manually making the lineage column blank
#when reading in the metadata from redcap
rule redcap_normal_pangolin:
    input:
#        previous_stage = config["output_path"] + "/logs/1_summarise_preprocess_uk.log",
        fasta = rules.redcap_add_dups_to_lineageless.output.fasta
    params:
        outdir = config["output_path"] + "/2/normal_pangolin",
        tmpdir = config["output_path"] + "/2/normal_pangolin/tmp"
    output:
        lineages = config["output_path"] + "/2/normal_pangolin/lineage_report.csv"
    log:
        config["output_path"] + "/logs/2_redcap_normal_pangolin.log"
    shell:
        """
        pangolin {input.fasta} \
        --outdir {params.outdir} \
        --tempdir {params.tmpdir}  >> {log} 2>&1
        """


rule redcap_filter_unassignable_lineage:
    input:
        lineages = rules.redcap_normal_pangolin.output.lineages
    output:
        lineages = config["output_path"] + "/2/normal_pangolin/lineage_report_filtered.csv", 
        previous_lineages = config["previous_outputs"] + "/most_recent/most_recent_redcap_lineages.csv",
        save_lineages = config["previous_outputs"] + "/" + config["date"] + "/past_redcap_pango_lin_report_filtered.csv"
    log:
        config["output_path"] + "/logs/2_redcap_filter_unassignable_lineage.log"
    run:
        import pandas as pd

        df = pd.read_csv(input.lineages)

        df_filtered = df.loc[df['lineage'].notnull()]
        df_unassigned = df.loc[df['lineage'].isnull()]

        df_filtered.to_csv(output.lineages, index=False)

        with open(str(log), "w") as log_out:
            log_out.write("The following sequences were not assigned a pangolin lineage: \n")
            [log_out.write(i + "\n") for i in df_unassigned['taxon']]
        log_out.close()

        shell("""
            cp {output.lineages} {output.previous_lineages}
            cp {output.lineages} {output.save_lineages}
            """)


#index column should be strain
#originally took rules.uk_add_previous_lineages_to_metadata.output.metadata from rule_1
#will append new lineages to pango column
rule redcap_add_pangolin_lineages_to_metadata:
    input:
        metadata = rules.redcap_add_del_finder_result_to_metadata.output.metadata,
        lineages = rules.redcap_filter_unassignable_lineage.output.lineages
    output:
        metadata = config["output_path"] + "/2/redcap.with_new_lineages.csv"
    log:
        config["output_path"] + "/logs/2_redcap_add_normal_pangolin_lineages_to_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.lineages} \
          --index-column strain \
          --join-on taxon \
          --new-columns pango lineage_support pango_version \
          --where-column pango=lineage lineage_support=probability pango_version=pangoLEARN_version \
          --out-metadata {output.metadata} &>> {log}
        """


rule get_filled_analysis_instrument:
    input:
        metadata = rules.redcap_add_pangolin_lineages_to_metadata.output.metadata
    output:
        metadata = config["output_path"] + "/2/filled_analysis_instrument.csv"
    run:
        import pandas as pd

        df = pd.read_csv(input.metadata)

        df.loc[:,'sequence_length'] = df.loc[:,'length']
        df.loc[:, 'gaps'] = df.loc[:, 'missing']
        df.loc[:, 'missing'] = df.loc[:, 'coverage']

        df = df.loc[:,['central_id', 'redcap_repeat_instance', \
                        'consensus', 'ave_depth', 'sequence_length', \
                        'missing', 'gaps', 'pango', 'lineage_support', \
                        'pango_version', 'p323l', 'd614g', \
                        'n439k', 'del_1605_3', 'epi_week']]

        df.to_csv(output.metadata, index=False)


rule import_analysis_instrument_to_redcap:
    input:
        metadata = rules.get_filled_analysis_instrument.output.metadata,
        redcap_db = config["redcap_access"]
    output:
        metadata = config["output_path"] + "/2/analysis_intrument_form_exported.csv"
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
        df.insert(1, 'redcap_repeat_instrument', 'Analysis')
        df.loc[:,'redcap_repeat_instrument'] = df.loc[:,'redcap_repeat_instrument'].str.casefold()

        #convert missing column to percentage
        df.loc[:,'missing'] = df.loc[:,'missing'].apply(lambda x:round(x*100,2))

        proj.import_records(df)
        df.to_csv(output.metadata, index=False)


#rename pango back to lineage
#uk_lineage/ph_cluster will be empty?
#need to use pandas to rename index column
#as the fetch command doesn't seem to support that
rule redcap_output_lineage_table:
    input:
        fasta = rules.redcap_filter_omitted_sequences.output.fasta,
        metadata = rules.redcap_add_pangolin_lineages_to_metadata.output.metadata
    params:
        country_code = config["country_code"]
    output:
        fasta = config["output_path"] + "/2/redcap.matched.fasta",
        metadata = config["output_path"] + "/2/redcap.matched.lineages.csv"
    log:
        config["output_path"] + "/logs/2_redcap_output_full_lineage_table.log"
    shell:
        """
        fastafunk fetch \
        --in-fasta {input.fasta} \
        --in-metadata {input.metadata} \
        --index-column strain \
        --filter-column strain country adm1 adm2 \
                        sample_date epi_week \
                        lineage {params.country_code}_lineage \
        --where-column country=adm0 lineage=pango {params.country_code}_lineage={params.country_code}_cluster\
        --out-fasta {output.fasta} \
        --out-metadata {output.metadata} \
        --log-file {log} \
        --low-memory \
        --restrict
         """


#strain is the gisaid equivalent to fasta_header
#will be necessary for combine step
#rule rename_fasta_header_to_strain:
#    input:
#        metadata = rules.uk_output_lineage_table.output.metadata
#    output:
#        metadata = config["output_path"] + "/2/RC_metadata.lineages.strain.csv"
#    run:
#        import pandas as pd
#
#        df = pd.read_csv(input.metadata)
#
#        df.rename(columns={'fasta_header':'strain'}, inplace=True)
#
#        df.to_csv(output.metadata, index=False)


#commented out as it doesn't seem immediately useful
#had to include 'analysis_form' input to upload to redcap
#there's probably a better way of doing that
rule summarise_pangolin_lineage_typing:
    input:
        fasta = rules.redcap_output_lineage_table.output.fasta,
        metadata = rules.redcap_output_lineage_table.output.metadata,
        analysis_form = rules.import_analysis_instrument_to_redcap.output.metadata
    params:
        grapevine_webhook = config["grapevine_webhook"],
        json_path = config["json_path"],
        date = config["date"]
    log:
        config["output_path"] + "/logs/2_summarise_pangolin_lineage_typing.log"
#    shell:
#        """
#        echo '{{"text":"' > {params.json_path}/2_data.json
#        echo "*Step 2: {params.date} COG-UK pangolin typing complete*\\n" >> {params.json_path}/2_data.json
#        echo '"}}' >> {params.json_path}/2_data.json
#        echo "webhook {params.grapevine_webhook}"
#        curl -X POST -H "Content-type: application/json" -d @{params.json_path}/2_data.json {params.grapevine_webhook}
#        """
