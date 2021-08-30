import os

#output.summary used until issue with liburl and chardet warning is resolved
#maybe summaries and tables should be under log rather than output
rule get_redcap_metadata:
    input:
        redcap_db = config["redcap_access"]
    params:
        script = os.path.join(workflow.current_basedir, "../utilities/get_redcap_metadata.py"),
        dag_dir = config["output_path"] + "/1/redcap_records_missing_data_by_dag/"
    output:
        metadata = config["output_path"] + "/1/redcap_metadata.csv",
        summary = config["output_path"] + "/1/get_redcap_metadata_summary.txt",
        no_consensus_table = config["output_path"] + "/1/redcap_records_without_consensus.csv",
        no_dates_table = config["output_path"] + "/1/redcap_records_without_dates.csv",
        dag_summary = config["output_path"] + "/1/redcap_dag_summary.txt"
    log:
        config["output_path"] + "/logs/1_get_redcap_metadata.log"
    shell:
        """
        mkdir -p {params.dag_dir}

        python {params.script} {input.redcap_db} {output.metadata} \
               {output.summary} {output.no_consensus_table} {output.no_dates_table} \
               {output.dag_summary} {params.dag_dir} &> {log}
        """


#filter out records if they exist in gisaid data
rule deduplicate_gisaid:
    input:
        redcap_metadata = rules.get_redcap_metadata.output.metadata,
        gisaid_metadata = rules.gisaid_output_all_matched_metadata.output.metadata
    output:
        metadata = config["output_path"] + "/1/redcap_metadata.gisaid_filtered.csv"
    log:
        config["output_path"] + "/logs/1_deduplicate_gisaid.log"
    run:
        import pandas as pd

        redcap_df = pd.read_csv(input.redcap_metadata)
        gisaid_df = pd.read_csv(input.gisaid_metadata)

        gisaid_names = [i.split('/')[2] for i in gisaid_df.loc[:,'strain']]
        gisaid_names_bool = redcap_df['gisaid_name'].isin(gisaid_names)
        gisaid_names_inverse_bool = ~redcap_df['gisaid_name'].isin(gisaid_names)

        redcap_filtered_df = redcap_df[gisaid_names_inverse_bool]
        gisaid_duplicates_df = redcap_df[gisaid_names_bool]

        redcap_filtered_df.to_csv(output.metadata, index=False)

        with open(str(log), 'w') as log_out:
            log_out.write("The following sequences were found in the input gisaid dataset and filtered out: \n")
            [log_out.write(i + "\n") for i in gisaid_duplicates_df['gisaid_name']]
        log_out.close()


rule get_redcap_fasta:
    input:
        redcap_db = config["redcap_access"],
        metadata = rules.deduplicate_gisaid.output.metadata
    output:
        fasta = config["output_path"] + "/1/redcap_fasta.fasta"
    run:
        import redcap
        import pandas as pd

        with open(input.redcap_db, 'r') as f:
            read_file = f.read().strip('\n')
            url = read_file.split(',')[0]
            key = read_file.split(',')[1]
        f.close()

        proj = redcap.Project(url, key)
        proj_df = pd.read_csv(input.metadata, index_col='central_id')

        with open(str(output.fasta), 'wb') as f:
            for i in proj_df.index:
                file_contents, headers = proj.export_file((str(i)), 'consensus')
                f.write(file_contents)        
        f.close()


#add sample date column based on collection date or received date
#originally took config["latest_uk_metadata"] as input
#there's already a function for this in datafunk/fastafunk, should probably use that
rule add_sample_date:
    input:
        metadata = rules.deduplicate_gisaid.output.metadata
    output:
        metadata = config["output_path"] + "/1/redcap_metadata.strain.sample_date.csv"
    log:
        config["output_path"] + "/logs/1_redcap_add_sample_date.log"
    resources: 
        mem_per_cpu=20000
    run:
        import pandas as pd

        df = pd.read_csv(input.metadata)

        sample_date = []

        for i,row in df.iterrows():

            if not pd.isnull(row['date_collected']) and row['date_collected'] != "None":
                sample_date.append(row['date_collected'])
            elif not pd.isnull(row['date_received']) and row['date_received'] != "None":
                sample_date.append(row['date_received'])
            else:
                sample_date.append("")

        df['sample_date'] = sample_date
        df.to_csv(output.metadata, index=False, sep = ",")


rule add_strain:
    input:
        metadata = rules.add_sample_date.output.metadata,
        fasta = rules.get_redcap_fasta.output.fasta
    output:
        metadata = config["output_path"] + "/1/redcap_metadata.strain.csv"
    run:
        import pandas as pd
        from Bio import SeqIO
        
        df = pd.read_csv(input.metadata)

        strain = []

        fasta_in = SeqIO.parse(str(input.fasta), "fasta")
        for record in fasta_in:
            strain.append(record.description)

        df['strain'] = strain
        df.to_csv(output.metadata, index=False, sep=",")


rule format_redcap_fasta_header_and_strain:
    input:
        fasta = rules.get_redcap_fasta.output.fasta,
        metadata = rules.add_strain.output.metadata
    params:
        country = config["country"].capitalize()
    output:
        fasta = config["output_path"] + "/1/redcap_formatted.fasta",
        metadata = config["output_path"] + "/1/redcap_metadata.strain.formatted.csv"
    log:
        config["output_path"] + "/logs/1_redcap_strip_header_digits.log"
    run:
        import pandas as pd
        from Bio import SeqIO

        df = pd.read_csv(input.metadata)
        fasta_in = SeqIO.index(str(input.fasta), "fasta")

        new_strain_col = []

        with open(str(output.fasta), 'w') as fasta_out:
            for i,row in df.iterrows():
                fasta_header = row["strain"].split(" ")[0] #split by whitespace since fasta headers may have them and then 'fasta_in[fasta_header]' will throw an KeyError
                gisaid_name = row["gisaid_name"]
                year = row["sample_date"].split('-')[0]
                new_header = "hCoV-19/" + params.country + "/" + gisaid_name + "/" + year
                new_strain_col.append(new_header)

                record = fasta_in[fasta_header]
                fasta_out.write(">" + new_header + "\n")
                fasta_out.write(str(record.seq) + "\n")

        df['strain'] = new_strain_col
        df.to_csv(output.metadata, index=False)


#not necessary
#rule uk_add_pillar_2:
#    input:
#        metadata = rules.uk_add_sample_date.output.metadata,
#    output:
#        metadata = config["output_path"] + "/1/uk_latest.add_pillar_2.csv",
#    log:
#        config["output_path"] + "/logs/1_uk_add_pillar_2.log"
#    resources: mem_per_cpu=20000
#    run:
#        df = pd.read_csv(input.metadata, sep = ",")
#
#        pillar_2 = []
#
#        for i,row in df.iterrows():
#
#            if row['collection_pillar'] == 2 or row['central_sample_id'][0:4] in ["ALDP", "CAMC", "MILK", "QEUH"]:
#                pillar_2.append(True)
#            else:
#                pillar_2.append(False)
#
#        df['pillar_2'] = pillar_2
#        df.to_csv(output.metadata, index=False, sep = ",")


#already have rule for formatting fasta header and strain column
#rule uk_make_sequence_name:
#    input:
#        metadata = rules.uk_add_pillar_2.output.metadata,
#    output:
#        metadata = config["output_path"] + "/1/uk_latest.add_sample_date.add_sequence_name.csv",
#    log:
#        config["output_path"] + "/logs/1_uk_make_sequence_name.log"
#    resources: mem_per_cpu=20000
#    run:
#        df = pd.read_csv(input.metadata)
#
#        adm1_to_country = {"UK-SCT": "Scotland",
#                           "UK-WLS": "Wales",
#                           "UK-ENG": "England",
#                           "UK-NIR": "Northern_Ireland"}
#
#        sequence_name = []
#
#        for i,row in df.iterrows():
#            country = adm1_to_country[row['adm1']]
#            id = row['central_sample_id']
#            year = str(row['sample_date']).split("-")[0]
#            name = country + "/" + id + "/" + year
#
#            sequence_name.append(name)
#
#        df['sequence_name'] = sequence_name
#        df.to_csv(output.metadata, index=False)


#not necessary, as gisaid ID is already in metadata
#rule uk_add_gisaid_accession:
#    input:
#        metadata = rules.uk_make_sequence_name.output.metadata,
#        accessions_table = config["latest_uk_accessions"],
#    output:
#        metadata = config["output_path"] + "/1/uk_latest.add_sample_date.add_sequence_name.accessions.csv"
#    log:
#        config["output_path"] + "/logs/1_uk_add_gisaid_accession.log"
#    resources: mem_per_cpu=20000
#    run:
#        logfile = open(str(log), "w")
#        df_acc = pd.read_csv(input.accessions_table, sep='\t')
#        accessions_dict = {}
#        for i,row in df_acc.iterrows():
#            central_sample_id = row["central_sample_id"]
#            run_name = row["run_name"]
#            gisaid_accession = row["gisaid.accession"]
#
#            if central_sample_id in accessions_dict:
#                if run_name in accessions_dict[central_sample_id]:
#                    logfile.write(f'duplicate central_sample_id * run_name in accessions list: {central_sample_id} {run_name}\n')
#                    continue
#                accessions_dict[central_sample_id][run_name] = gisaid_accession
#            else:
#                accessions_dict[central_sample_id] = {run_name: gisaid_accession}
#
#        df = pd.read_csv(input.metadata)
#        covv_accession_id = []
#        for i,row in df.iterrows():
#            acc = ""
#            if row["central_sample_id"] in accessions_dict:
#                if row["run_name"] in accessions_dict[row["central_sample_id"]]:
#                    acc = accessions_dict[row["central_sample_id"]][row["run_name"]]
#
#            covv_accession_id.append(acc)
#
#        df['covv_accession_id'] = covv_accession_id
#        df.to_csv(output.metadata, index=False)
#        logfile.close()


#not sure what duplicates are being removed here
#there shouldn't be duplicate fastas at this stage
#annotate adds length, missing and gaps columns
rule annotate_to_remove_duplicates:
    input:
        fasta = rules.format_redcap_fasta_header_and_strain.output.fasta,
        metadata = rules.format_redcap_fasta_header_and_strain.output.metadata
    output:
        metadata = config["output_path"] + "/1/redcap_latest.add_sample_date.add_sequence_name.accessions.annotated.csv"
    log:
        config["output_path"] + "/logs/1_redcap_annotate_to_remove_duplicates.log"
    resources: 
        mem_per_cpu=20000
    shell:
        """
        fastafunk annotate \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --index-column strain &> {log}
        """
          # --add-cov-id \


#need to make sure method of calculating coverage is okay
rule add_coverage_column:
    input:
        fasta = rules.format_redcap_fasta_header_and_strain.output.fasta,
        metadata = rules.annotate_to_remove_duplicates.output.metadata
    output:
        metadata = config["output_path"] + "/1/redcap_latest.add_sample_date.add_sequence_name.accessions.annotated.covg.csv"
    run:
        import pandas as pd
        from Bio import SeqIO

        df = pd.read_csv(input.metadata)
        coverage = []

        record_dict = SeqIO.index(input.fasta, "fasta")

        for record in df['strain']:
            record_seq = str(record_dict[record].seq)
            seq_sans_N = record_seq.replace("N", "")
            coverage.append(float(len(seq_sans_N)) / len(record_seq))
        
        df['coverage'] = coverage

        df.to_csv(output.metadata, index=False)


#fasta sequences deduplicated based on which has the fewest gaps
#we currently don't have a good way to deal with linking repeat instruments in redcap
#so at the moment there are no duplicate central IDs at this point
#because they are removed during the 'get_redcap_metadata' rule
#rule uk_remove_duplicates_COGID_by_gaps:
#    input:
#        fasta = rules.uk_strip_header_digits.output.fasta,
#        metadata = rules.uk_annotate_to_remove_duplicates.output.metadata
#    output:
#        fasta = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id.fasta",
#        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id.csv"
#    log:
#        config["output_path"] + "/logs/1_uk_filter_duplicates_bycovid.log"
#    resources: mem_per_cpu=20000
#    shell:
#        """
#        fastafunk subsample \
#          --in-fasta {input.fasta} \
#          --in-metadata {input.metadata} \
#          --group-column central_sample_id \
#          --index-column fasta_header \
#          --out-fasta {output.fasta} \
#          --out-metadata {output.metadata} \
#          --sample-size 1 \
#          --select-by-min-column gaps &> {log}
#        """


#not necessary, as metadata is always pulled directly from redcap
#rule uk_update_sample_dates:
#    input:
#        metadata = rules.uk_remove_duplicates_COGID_by_gaps.output.metadata,
#        updated_dates = config["uk_updated_dates"],
#    output:
#        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id.sample_date.updated_sample_date.csv",
#    log:
#        config["output_path"] + "/logs/1_uk_update_sample_dates.log"
#    resources: mem_per_cpu=20000
#    shell:
#        """
#        fastafunk add_columns \
#          --in-metadata {input.metadata} \
#          --in-data {input.updated_dates} \
#          --index-column central_sample_id \
#          --join-on central_sample_id \
#          --new-columns sample_date \
#          --out-metadata {output.metadata} &>> {log}
#        """


#previously took rules.uk_update_sample_dates.output.metadata as input
rule add_epi_week:
    input:
        metadata = rules.add_coverage_column.output.metadata
    output:
        metadata = config["output_path"] + "/1/redcap_latest.add_header.annotated.deduplicated_cov_id.sample_date.updated_sample_date.epi_week.csv"
    log:
        config["output_path"] + "/logs/1_redcap_add_epi_week.log"
    resources: 
        mem_per_cpu=20000
    shell:
        """
        datafunk add_epi_week \
        --input-metadata {input.metadata} \
        --output-metadata {output.metadata} \
        --date-column sample_date \
        --epi-week-column-name epi_week \
        --epi-day-column-name epi_day &> {log}
        """


#biosample not used in our data so following 3 rules are unnecessary
#rule uk_annotate_to_remove_duplicates_by_biosample:
#    input:
#        metadata = rules.uk_add_epi_week.output.metadata,
#    output:
#        metadata = config["output_path"] + "/1/uk_latest.epi_week.annotated2.csv",
#    log:
#        config["output_path"] + "/logs/1_uk_annotate_to_remove_duplicates_by_biosample.log",
#    resources: mem_per_cpu=20000
#    run:
#        df = pd.read_csv(input.metadata)
#
#        edin_biosample = []
#        edin_biosample_root = []
#
#        for i,row in df.iterrows():
#            if pd.isnull(row['biosample_source_id']):
#                edin_biosample.append(row['central_sample_id'])
#            else:
#                edin_biosample.append(row['biosample_source_id'])
#
#            if pd.isnull(row['root_biosample_source_id']):
#                edin_biosample_root.append(row['central_sample_id'])
#            else:
#                edin_biosample_root.append(row['root_biosample_source_id'])
#
#        df['edin_biosample'] = edin_biosample
#        df['edin_biosample_root'] = edin_biosample_root
#        df.to_csv(output.metadata, index=False)
#
#
#rule uk_remove_duplicates_biosamplesourceid_by_date:
#    input:
#        fasta = rules.uk_remove_duplicates_COGID_by_gaps.output.fasta,
#        metadata = rules.uk_annotate_to_remove_duplicates_by_biosample.output.metadata
#    output:
#        fasta = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id_biosample_source_id.fasta",
#        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_cov_id_biosample_source_id.csv"
#    log:
#        config["output_path"] + "/logs/1_uk_filter_duplicates_by_biosample.log"
#    resources: mem_per_cpu=20000
#    shell:
#        """
#        fastafunk subsample \
#          --in-fasta {input.fasta} \
#          --in-metadata {input.metadata} \
#          --group-column edin_biosample \
#          --index-column fasta_header \
#          --out-fasta {output.fasta} \
#          --out-metadata {output.metadata} \
#          --sample-size 1 \
#          --select-by-min-column edin_epi_day &> {log}
#        """
#
#
#rule uk_remove_duplicates_root_biosample_by_gaps:
#    input:
#        fasta = rules.uk_remove_duplicates_biosamplesourceid_by_date.output.fasta,
#        metadata = rules.uk_remove_duplicates_biosamplesourceid_by_date.output.metadata
#    output:
#        fasta = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_rootbiosample.fasta",
#        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated_rootbiosample.csv"
#    log:
#        config["output_path"] + "/logs/1_uk_filter_duplicates_root_biosample_by_gaps.log"
#    resources: mem_per_cpu=20000
#    shell:
#        """
#        fastafunk subsample \
#          --in-fasta {input.fasta} \
#          --in-metadata {input.metadata} \
#          --group-column edin_biosample_root \
#          --index-column fasta_header \
#          --out-fasta {output.fasta} \
#          --out-metadata {output.metadata} \
#          --sample-size 1 \
#          --select-by-min-column gaps &> {log}
#        """


#already have rule for formatting fasta header and strain columns
#rule uk_unify_headers:
#    input:
#        fasta = rules.uk_remove_duplicates_root_biosample_by_gaps.output.fasta,
#        metadata = rules.uk_remove_duplicates_root_biosample_by_gaps.output.metadata
#    output:
#        fasta = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated.unify_headers.fasta",
#        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated.unify_headers.csv"
#    log:
#        config["output_path"] + "/logs/1_uk_unify_headers.log"
#    resources: mem_per_cpu=20000
#    run:
#        df = pd.read_csv(input.metadata)
#        header_dict = {}
#
#        for i,row in df.iterrows():
#            header_dict[row['fasta_header']] = row['sequence_name']
#
#        fasta_in = SeqIO.parse(str(input.fasta), "fasta")
#        with open(str(output.fasta), 'w') as f:
#            for record in fasta_in:
#                new_ID = header_dict[record.description]
#                f.write(">" + new_ID + "\n")
#                f.write(str(record.seq) + "\n")
#
#        df.to_csv(output.metadata, index=False)


#not sure if I should adapt to PH or just leave it out
#skip for now
#rule uk_sed_United_Kingdom_to_UK:
#    input:
#        metadata = rules.uk_unify_headers.output.metadata
#    output:
#        metadata = config["output_path"] + "/1/uk_latest.add_header.annotated.deduplicated.unify_headers.uk_sedded.csv"
#    log:
#        config["output_path"] + "/logs/1_uk_sed_United_Kingdom_to_UK.log"
#    resources: mem_per_cpu=20000
#    shell:
#        """
#        sed 's/United Kingdom/UK/g' {input.metadata} > {output.metadata}
#        """


#filters fasta seqs by length
#threshold set in config
#current length threshold: 29000
rule redcap_filter_by_length:
    input:
        fasta = rules.format_redcap_fasta_header_and_strain.output.fasta
    params:
        min_length = config["min_length"]
    output:
        fasta = config["output_path"] + "/1/redcap.RD.UH.filter_length.fasta"
    log:
        config["output_path"] + "/logs/redcap_filter_by_length.log"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min-length {params.min_length} &> {log}
        """


#originally took rules.uk_unify_headers.output.fasta as input
rule redcap_minimap2_to_reference:
    input:
        fasta = rules.redcap_filter_by_length.output.fasta,
        reference = config["reference_fasta"]
    output:
        sam = config["output_path"] + "/1/redcap_latest.unify_headers.epi_week.deduplicated.mapped.sam"
    log:
        config["output_path"] + "/logs/1_redcap_minimap2_to_reference.log"
    threads: 16
    resources: 
        mem_per_cpu=2000
    shell:
        """
        minimap2 -t {threads} -a -x asm5 {input.reference} {input.fasta} > {output.sam} 2> {log}
        """


rule redcap_get_variants:
    input:
        sam = rules.redcap_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"],
        genbank_anno = config["reference_genbank_annotation"]
    output:
        variants = config["output_path"] + "/1/redcap.variants.csv"
    log:
        config["output_path"] + "/logs/1_redcap_get_variants.log"
    threads: 12
    shell:
        """
        gofasta sam variants -t {threads} \
            --samfile {input.sam} \
            --reference {input.reference} \
            --genbank {input.genbank_anno} \
            --outfile {output.variants} &>> {log}
        """


#        gofasta sam toMultiAlign \
#        -t {threads} \
#        -s {input.sam} \
#        -o {output.fasta} \
#        --trim \
#        --trimstart {params.trim_start} \
#        --trimend {params.trim_end} \
#        --pad &> {log}
#
#        mv insertions.txt {params.insertions}
#        mv deletions.txt {params.deletions}
#        cp {params.insertions} {output.insertions}
#        cp {params.deletions} {output.deletions}
rule redcap_remove_insertions_and_trim_and_pad:
    input:
        sam = rules.redcap_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"]
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
#        insertions = config["output_path"] + "/1/uk_insertions.txt",
#        deletions = config["output_path"] + "/1/uk_deletions.txt"
    output:
        fasta = config["output_path"] + "/1/redcap_latest.unify_headers.epi_week.deduplicated.alignment.trimmed.fasta",
#        insertions = config["export_path"] + "/metadata/uk_insertions.txt",
#        deletions = config["export_path"] + "/metadata/uk_deletions.txt"
    log:
        config["output_path"] + "/logs/1_redcap_remove_insertions_and_trim_and_pad.log"
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam} \
          -r {input.reference} \
          -o {output.fasta} \
          -t [{params.trim_start}:{params.trim_end}] \
          --pad \
          --log-inserts \
          --log-deletions &> {log}
        """


#mask file is currently same one used in gisaid preprocessing step
rule redcap_mask_1:
    input:
        fasta = rules.redcap_remove_insertions_and_trim_and_pad.output.fasta,
        mask = config["redcap_mask_file"]
    output:
        fasta = config["output_path"] + "/1/redcap_latest.unify_headers.epi_week.deduplicated.alignment.trimmed.masked.fasta",
    log:
        config["output_path"] + "/logs/1_redcap_mask_1.log"
    shell:
        """
        datafunk mask \
          --input-fasta {input.fasta} \
          --output-fasta {output.fasta} \
          --mask-file \"{input.mask}\" 2> {log}
        """


rule redcap_filter_low_coverage_sequences:
    input:
        fasta = rules.redcap_mask_1.output.fasta
    params:
        min_covg = config["min_covg"]
    output:
        fasta = config["output_path"] + "/1/redcap_latest.unify_headers.epi_week.deduplicated.trimmed.low_covg_filtered.fasta"
    log:
        config["output_path"] + "/logs/1_redcap_filter_low_coverage_sequences.log"
    resources: 
        mem_per_cpu=20000
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output.fasta} \
          --min-covg {params.min_covg} &> {log}
        """


#will need to get omissions list
#do we have or need one?
rule redcap_filter_omitted_sequences:
    input:
        fasta = rules.redcap_filter_low_coverage_sequences.output.fasta,
        omissions = config["redcap_omissions"]
    output:
        fasta = config["output_path"] + "/1/redcap_latest.unify_headers.epi_week.deduplicated.trimmed.low_covg_filtered.omissions_filtered.fasta"
    log:
        config["output_path"] + "/logs/1_redcap_filter_omitted_sequences.log"
    resources: 
        mem_per_cpu=20000
    shell:
        """
        datafunk remove_fasta \
          -i {input.fasta} \
          -f {input.omissions} \
          -o {output.fasta}  &> {log}
        """


rule redcap_full_untrimmed_alignment:
    input:
        sam = rules.redcap_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"],
    output:
        fasta = config["output_path"] + "/1/redcap_latest.unify_headers.epi_week.deduplicated.alignment.full.fasta"
    log:
        config["output_path"] + "/logs/1_redcap_full_untrimmed_alignment.log"
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam} \
          -r {input.reference} \
          -o {output.fasta} \
          &> {log}
        """


#using the same mask file
rule redcap_mask_2:
    input:
        fasta = rules.redcap_full_untrimmed_alignment.output.fasta,
        mask = config["redcap_mask_file"]
    output:
        fasta = config["output_path"] + "/1/redcap_latest.unify_headers.epi_week.deduplicated.alignment.full.masked.fasta"
    log:
        config["output_path"] + "/logs/1_redcap_mask_2.log"
    shell:
        """
        datafunk mask \
          --input-fasta {input.fasta} \
          --output-fasta {output.fasta} \
          --mask-file \"{input.mask}\" 2> {log}
        """


# get the same alignment as we use for tree building but with no mask
rule redcap_get_unmasked_alignment:
    input:
        fasta_template = rules.redcap_filter_omitted_sequences.output.fasta,
        fasta = rules.redcap_full_untrimmed_alignment.output.fasta
    output:
        fasta = config["output_path"] + "/1/redcap_latest.unify_headers.epi_week.deduplicated.alignment.full.unmasked.fasta"
    log:
        config["output_path"] + "/logs/1_redcap_get_unmasked_alignment.log"
    run:
        from Bio import SeqIO

        fasta_template = SeqIO.index(str(input.fasta_template), "fasta")
        fasta_in = SeqIO.index(str(input.fasta), "fasta")

        with open(str(output.fasta), 'w') as fasta_out:
            for record in fasta_in:
                if record in fasta_template:
                    fasta_out.write('>' + fasta_in[record].id + '\n')
                    fasta_out.write(str(fasta_in[record].seq) + '\n')


rule redcap_AA_finder:
    input:
        fasta = rules.redcap_full_untrimmed_alignment.output.fasta,
        AAs = config["AAs"]
    output:
        found = config["output_path"] + "/1/redcap.AA_finder.csv"
    log:
        config["output_path"] + "/logs/1_redcap_AA_finder.log"
    shell:
        """
        datafunk AA_finder -i {input.fasta} --codons-file {input.AAs} --genotypes-table {output.found} &> {log}
        """


rule add_AA_finder_result_to_metadata:
    input:
        AAs = config["AAs"],
        metadata = rules.add_epi_week.output.metadata,
        new_data = rules.redcap_AA_finder.output.found
    output:
        metadata = config["output_path"] + "/1/redcap_latest.unify_headers.epi_week.deduplicated.with_AA_finder.csv"
    log:
        config["output_path"] + "/logs/1_add_AA_finder_result_to_metadata.log"
    shell:
        """
        columns=$(head -n1 {input.new_data} | cut -d',' -f2- | tr ',' ' ')
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.new_data} \
          --index-column strain \
          --join-on sequence_name \
          --new-columns $columns \
          --out-metadata {output.metadata} &>> {log}
        """


rule redcap_del_finder:
    input:
        fasta = rules.redcap_full_untrimmed_alignment.output.fasta,
        dels = config["dels"]
    output:
        metadata = config["output_path"] + "/1/redcap.del_finder.csv"
    log:
        config["output_path"] + "/logs/1_redcap_del_finder.log"
    shell:
        """
        datafunk del_finder \
            -i {input.fasta} \
            --deletions-file {input.dels} \
            --genotypes-table {output.metadata} &> {log}
        """


#dels defined but not used?
rule redcap_add_del_finder_result_to_metadata:
    input:
        dels = config["dels"],
        metadata = rules.add_AA_finder_result_to_metadata.output.metadata,
        new_data = rules.redcap_del_finder.output.metadata
    output:
        metadata = config["output_path"] + "/1/redcap_latest.unify_headers.epi_week.deduplicated.with_AA_finder.with_del_finder.csv"
    log:
        config["output_path"] + "/logs/1_add_del_finder_result_to_metadata.log"
    shell:
        """
        columns=$(head -n1 {input.new_data} | cut -d',' -f2- | tr ',' ' ')
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.new_data} \
          --index-column strain \
          --join-on sequence_name \
          --new-columns $columns \
          --out-metadata {output.metadata} &>> {log}
        """


"""
Instead of new sequences (as determined by a date stamp), it might be more robust
to extract sequences for lineage typing that don't currently have an associated
lineage designation in the metadata file.
"""
#our pipeline isn't currently using previous metadata or global lineages
#skip for now
#rule uk_add_previous_lineages_to_metadata:
#    input:
#        metadata = rules.uk_add_del_finder_result_to_metadata.output.metadata,
#        previous_metadata = config["previous_uk_metadata"],
#        global_lineages = config["global_lineages"]
#    output:
#        metadata_temp = temp(config["output_path"] + "/1/uk.with_previous_lineages.temp.csv"),
#        metadata = config["output_path"] + "/1/uk.with_previous_lineages.csv",
#    log:
#        config["output_path"] + "/logs/1_uk_add_previous_lineages_to_metadata.log"
#    shell:
#        """
#        fastafunk add_columns \
#          --in-metadata {input.metadata} \
#          --in-data {input.previous_metadata} \
#          --index-column fasta_header \
#          --join-on sequence_name \
#          --new-columns uk_lineage edin_date_stamp \
#          --out-metadata {output.metadata_temp} &>> {log}
#
#        fastafunk add_columns \
#          --in-metadata {output.metadata_temp} \
#          --in-data {input.global_lineages} \
#          --index-column fasta_header \
#          --join-on taxon \
#          --new-columns lineage lineage_support lineages_version \
#          --where-column lineage_support=probability lineages_version=pangoLEARN_version \
#          --out-metadata {output.metadata} &>> {log}
#        """


#previously took rules.uk_add_previous_lineages_to_metadata.output.metadata as input
#takes records without pango lineage for assignment in later rule
rule redcap_extract_lineageless:
    input:
        fasta = rules.redcap_filter_omitted_sequences.output,
        metadata = rules.redcap_add_del_finder_result_to_metadata.output.metadata
    output:
        fasta = config["output_path"] + "/1/redcap.new.pangolin_lineages.fasta"
    log:
        config["output_path"] + "/logs/1_extract_lineageless.log"
    resources: 
        mem_per_cpu=20000
    run:
        from Bio import SeqIO
        import pandas as pd

        fasta_in = SeqIO.index(str(input.fasta), "fasta")
        df = pd.read_csv(input.metadata)

        sequence_record = []

        with open(str(output.fasta), 'w') as fasta_out:
            for i,row in df.iterrows():
                if pd.isnull(row['pango']) or row['pango']=='?':
                    sequence_name = row['strain']
                    if sequence_name in fasta_in:
                        if sequence_name not in sequence_record:
                            record = fasta_in[sequence_name]
                            fasta_out.write('>' + record.id + '\n')
                            fasta_out.write(str(record.seq) + '\n')
                            sequence_record.append(sequence_name)


#dup logs commented out for now
#previously took rules.uk_add_previous_lineages_to_metadata.output.metadata as input
#not doing anything at the moment since deduplicating rules aren't in use
rule redcap_add_dups_to_lineageless:
    input:
        master_fasta = rules.redcap_filter_omitted_sequences.output.fasta,
        lineageless_fasta = rules.redcap_extract_lineageless.output.fasta,
        metadata = rules.redcap_add_del_finder_result_to_metadata.output.metadata,
        # dup_cogid_log = rules.uk_remove_duplicates_COGID_by_gaps.log,
        # dup_biosample_log = rules.uk_remove_duplicates_biosamplesourceid_by_date.log,
        # dup_rootbio_log = rules.uk_remove_duplicates_root_biosample_by_gaps.log
    output:
        fasta = config["output_path"] + "/1/redcap.new.dedupes.pangolin_lineages.fasta"
    log:
        config["output_path"] + "/logs/1_redcap_add_dups_to_lineageless.log"
    resources: 
        mem_per_cpu=20000
    run:
        from Bio import SeqIO
        import pandas as pd

        master_fasta_in = SeqIO.index(str(input.master_fasta), "fasta")

        df = pd.read_csv(input.metadata)

        lineageless_sequence_record = set()
        dup_record = set()

        lineageless_fasta_in = open(str(input.lineageless_fasta), "r")
        pangolin_fasta_out = open(str(output.fasta), "w")

        for record in SeqIO.parse(lineageless_fasta_in, "fasta"):
            pangolin_fasta_out.write(">" + record.id + "\n")
            pangolin_fasta_out.write(str(record.seq) + "\n")
            lineageless_sequence_record.add(record.id)

        #with open(str(input.dup_cogid_log), "r") as cogid_log:
        #    for line in cogid_log:
        #        l = line.rstrip()
        #        if l.startswith("COGUK/"):
        #            dup_record.add(l)

        # # can uncomment these if necessary, but it shouldn't be
        # with open(str(input.dup_biosample_log), "r") as biosample_log:
        #     for line in dup_biosample_log:
        #         l = line.rstrip()
        #         if l.startswith("COGUK/"):
        #             dup_record.add(l)
        #
        # with open(str(input.dup_rootbio_log), "r") as rootbio_log:
        #     for line in rootbio_log:
        #         l = line.rstrip()
        #         if l.startswith("COGUK/"):
        #             dup_record.add(l)

        for i,row in df.iterrows():
            if row["strain"] in dup_record:
                if row["sequence_name"] in lineageless_sequence_record:
                    continue
                if row["sequence_name"] in master_fasta_in:
                    record = master_fasta_in[row["sequence_name"]]
                    pangolin_fasta_out.write(">" + record.id + "\n")
                    pangolin_fasta_out.write(str(record.seq) + "\n")

        lineageless_fasta_in.close()
        pangolin_fasta_out.close()


#following echo lines removed for now:
#        echo "> Number of sequences after deduplication by central_sample_id: $(cat {input.deduplicated_fasta_by_covid} | grep ">" | wc -l)\\n" &>> {log}
#        echo "> Number of sequences after deduplication by bio_sample_id: $(cat {input.deduplicated_fasta_by_biosampleid} | grep ">" | wc -l)\\n" &>> {log}
#        echo "> Number of sequences after deduplication by root_source_biosample_id: $(cat {input.deduplicated_fasta_by_rootbiosample} | grep ">" | wc -l)\\n" &>> {log}
#        echo "> Number of sequences after unifying headers: $(cat {input.unify_headers_fasta} | grep ">" | wc -l)\\n" &>> {log}
#        echo "> Number of sequences after mapping: $(cat {input.removed_low_covg_fasta} | grep ">" | wc -l)\\n" &>> {log}
#        echo "> Number of sequences after removing those in omissions file: $(cat {input.removed_omitted_fasta} | grep ">" | wc -l)\\n" &>> {log}
#
#        echo "> Number of new/deduped sequences passed to Pangolin for typing: $(cat {input.pangolin_fasta} | grep ">" | wc -l)\\n" &>> {log}
#
rule summarise_preprocess_redcap:
    input:
        webhook = config["webhook"],
        metadata_consensus_and_date_filter = rules.get_redcap_metadata.output.summary,
        dag_summary = rules.get_redcap_metadata.output.dag_summary,
        deduplicated_metadata_by_gisaid = rules.deduplicate_gisaid.output.metadata,
        raw_fasta = rules.format_redcap_fasta_header_and_strain.output.fasta,
        #deduplicated_fasta_by_covid = rules.uk_remove_duplicates_COGID_by_gaps.output.fasta,
        #deduplicated_fasta_by_biosampleid = rules.uk_remove_duplicates_biosamplesourceid_by_date.output.fasta,
        #deduplicated_fasta_by_rootbiosample = rules.uk_remove_duplicates_root_biosample_by_gaps.output.fasta,
        #unify_headers_fasta = rules.uk_unify_headers.output.fasta,
        low_length_fasta_filter = rules.redcap_filter_by_length.output.fasta,
        low_covg_fasta_filter = rules.redcap_filter_low_coverage_sequences.output.fasta,
        removed_omitted_fasta = rules.redcap_filter_omitted_sequences.output.fasta,
        full_unmasked_alignment = rules.redcap_full_untrimmed_alignment.output.fasta,
        full_metadata = rules.redcap_add_del_finder_result_to_metadata.output.metadata,
        #full_metadata = rules.uk_add_previous_lineages_to_metadata.output.metadata,
        #pangolin_fasta = rules.uk_add_dups_to_lineageless.output.fasta,
        variants = rules.redcap_get_variants.output.variants
    params:
        #grapevine_webhook = config["grapevine_webhook"],
        json_path = config["json_path"],
        date=config["date"]
    log:
        config["output_path"] + "/logs/1_summarise_preprocess_redcap.log"
    shell:
        """
        #creating log file
        cat {input.metadata_consensus_and_date_filter} &> {log}
        echo "Number of records after deduplicating based on gisaid name: $(cat {input.deduplicated_metadata_by_gisaid} | tail -n +2 | wc -l)\\n" &>> {log}
        echo "Number of sequences in redcap fasta: $(cat {input.raw_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "Number of sequences after filtering by length: $(cat {input.low_length_fasta_filter} | grep ">" | wc -l)\\n" &>> {log}
        echo "Number of sequences after filtering by coverage: $(cat {input.low_covg_fasta_filter} | grep ">" | wc -l)" &>> {log}

        #creating json to fit log into
        echo '{{ "attachments": [ {{ "color": "#d61c0f", "blocks": [ {{ "type" : "section", "text" : {{ "type" : "mrkdwn", "text" : "*Redcap preprocessing complete*\n" }} }}, {{ "type": "divider" }}, {{ "type": "section", "text": {{ "type": "mrkdwn", "text": "' > {params.json_path}/1_data.json
        cat {log} >> {params.json_path}/1_data.json
        echo '" }} }} ] }} ] }}' >> {params.json_path}/1_data.json
        curl -X POST -H "Content-type: application/json" -d @{params.json_path}/1_data.json $(cat {input.webhook} | xargs)

        #creating json to fit dag summary into
        echo '{{ "attachments": [ {{ "color": "#d61c0f", "blocks": [ {{ "type" : "section", "text" : {{ "type" : "mrkdwn", "text" : "*Redcap Records SNL Summary*\\n" }} }}, {{ "type": "divider" }}, {{ "type": "section", "text": {{ "type": "mrkdwn", "text": "```' > {params.json_path}/1_dag_summary.json
        cat {input.dag_summary} >> {params.json_path}/1_dag_summary.json
        echo '```"}} }} ] }} ] }}' >> {params.json_path}/1_dag_summary.json
        curl -X POST -H "Content-type: application/json" -d @{params.json_path}/1_dag_summary.json $(cat {input.webhook} | xargs)
        """
