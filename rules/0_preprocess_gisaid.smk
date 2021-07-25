#rule gisaid_process_json:
#    input:
#        json = config["latest_gisaid_json"],
#        omitted = config["previous_omitted_file"],
#    output:
#        fasta = config["output_path"] + "/0/gisaid.fasta",
#        metadata = config["output_path"] + "/0/gisaid.csv"
#    log:
#        config["output_path"] + "/logs/0_gisaid_process_json.log"
#    resources: mem_per_cpu=20000
#    shell:
#        """
#        datafunk process_gisaid_data \
#          --input-json {input.json} \
#          --input-metadata False \
#          --exclude-file {input.omitted} \
#          --output-fasta {output.fasta} \
#          --output-metadata {output.metadata} \
#          --exclude-undated &> {log}
#        """



#gisaid_unify_headers standardises sequence names in both fasta and metadata files
#this isn't necessary for current nextstrain data
#rule gisaid_unify_headers:
#    input:
#        fasta = config["gisaid_fasta"],
#        metadata = config["gisaid_meta"],
#        #fasta = rules.gisaid_process_json.output.fasta,
#        #metadata = rules.gisaid_process_json.output.metadata,
#    output:
#        fasta = config["output_path"] + "/0/gisaid.UH.fasta",
#        metadata = config["output_path"] + "/0/gisaid.UH.tsv",
#    log:
#        config["output_path"] + "/logs/0_gisaid_unify_headers.log"
#    run:
#        import pandas as pd
#        from Bio import SeqIO
#
#        fasta_in = SeqIO.index(str(input.fasta), "fasta")
#        df = pd.read_csv(input.metadata, sep='\t')
#
#        sequence_name = []
#
#        #start writing to output fasta
#        #
#        with open(str(output.fasta), 'w') as fasta_out:
#            for i,row in df.iterrows():
#                edin_header = row["edin_header"]
#                new_header = edin_header.split("|")[0]
#                sequence_name.append(new_header)
#
#                try:
#                    record = fasta_in[edin_header]
#                    fasta_out.write(">" + new_header + "\n")
#                    fasta_out.write(str(record.seq) + "\n")
#                except:
#                    continue
#
#        df['sequence_name'] = sequence_name
#        df.to_csv(output.metadata, index=False, sep = ",")


#in-data and in-metadata are assumed csv apparently
#index-column in metadata is shared by data
#join-on refers to column in data shared by metadata
#new-columns are columns in in-data
#rule not really necessary at the moment
#rule gisaid_add_previous_lineages:
#    input:
#        metadata = rules.gisaid_unify_headers.output.metadata,
#        previous_lineages = config["previous_gisaid_lineages"],
#    output:
#        metadata = config["output_path"] + "/0/gisaid.lineages.csv"
#    log:
#        config["output_path"] + "/logs/0_gisaid_add_previous_lineages.log"
#    shell:
#        """
#        fastafunk add_columns \
#          --in-metadata {input.metadata} \
#          --in-data {input.previous_lineages} \
#          --index-column sequence_name \
#          --join-on sequence_name \
#          --new-columns lineage lineage_support lineages_version \
#          --out-metadata {output.metadata} &> {log}
#        """


#reads metadata and fasta for lineage outgroups
#required for rerooting in treebuilding step
#will be removed in later rule if duplicate exists
rule add_gisaid_outgroups:
    input:
        metadata = config["gisaid_meta"],
        fasta = config["gisaid_fasta"],
        outgroups_metadata = config["outgroups_metadata"],
        outgroups_fasta = config["outgroups_fasta"]
    output:
        metadata = config["output_path"] + "/0/gisaid.with_outgroups.csv",
        fasta = config["output_path"] + "/0/gisaid.with_outgroups.fasta"
    run:
        import pandas as pd
        from Bio import SeqIO

        gisaid_meta_df = pd.read_csv(input.metadata, sep="\t")
        outgroup_meta_df = pd.read_csv(input.outgroups_metadata, sep="\t")

        combine = [gisaid_meta_df, outgroup_meta_df]
        combine_df = pd.concat(combine)
        combine_df.to_csv(output.metadata, index=False, sep = ",")

        fastas = [input.fasta, input.outgroups_fasta]

        with open((output.fasta), 'w') as out_handle:
            for fasta in fastas:
                for record in SeqIO.parse(fasta, "fasta"):
                    SeqIO.write(record, out_handle, "fasta")
        out_handle.close()


#gets difference between value in 'date' column and current date
#returns number of days which is added as a new column
#there is a datafunk/fastafunk function which does the same thing...
#change later
rule add_days_and_weeks_since_epi:
    input:
        metadata = rules.add_gisaid_outgroups.output.metadata
    output:
        metadata = config["output_path"] + "/0/gisaid.UH.epi_day.csv"
    run:
        import datetime
        import pandas as pd

        date_format = "%Y-%m-%d"

        df = pd.read_csv(input.metadata)

        df['date'] = pd.to_datetime(df['date'], format= date_format)
        df['epi'] = pd.to_datetime('2019-12-22', format= date_format)
        df['epi_day'] = df['date'] - df['epi']
        df['epi_day'] = [df['epi_day'][i].days for i in range(len(df))]
        df['epi_week'] = df['epi_day']//7

        temp_day = df['epi_day']
        temp_week = df['epi_week']
        df.drop(columns=['epi_day', 'epi_week'], inplace=True)
        df.insert(loc=5, column='epi_week', value=temp_week)
        df.insert(loc=5, column='epi_day', value=temp_day)

        df.drop(columns=['epi'], inplace=True)

        df.to_csv(output.metadata, index=False, sep = ",")


#groups together by sequence name and then takes one from each group,
#removing duplicates
#is it possible for fasta sequences with identical names to have non-identical sequences?
#note that 'sequence_name' column is now 'strain'
#of duplicates, most recent sequence is kept, given by 'select-by-min'
#when deduplicating redcap by date, the oldest is kept...
rule gisaid_remove_duplicates:
    input:
        fasta = rules.add_gisaid_outgroups.output.fasta,
        metadata = rules.add_days_and_weeks_since_epi.output.metadata,
#        fasta = rules.gisaid_unify_headers.output.fasta,
#        metadata = rules.gisaid_add_previous_lineages.output.metadata
    output:
        fasta = config["output_path"] + "/0/gisaid.UH.RD.fasta",
        metadata = config["output_path"] + "/0/gisaid.UH.RD.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_filter_duplicates.log"
    resources: mem_per_cpu=20000
    shell:
        """
        fastafunk subsample \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --group-column strain \
          --index-column strain \
          --out-fasta {output.fasta} \
          --sample-size 1 \
          --out-metadata {output.metadata} \
          --select-by-min-column 'epi_week' &> {log}
        """


#counts number of records by grouping on 'country' column
#was previously 'edin_admin_0' column
#not really a priority
#rule gisaid_counts_by_country:
#    input:
#        metadata = rules.gisaid_remove_duplicates.output.metadata
#    output:
#        counts = config["output_path"] + "/0/gisaid_counts_by_country.csv",
#    log:
#        config["output_path"] + "/logs/0_gisaid_counts_by_country.log"
#    resources: mem_per_cpu=20000
#    shell:
#        """
#        fastafunk count \
#          --in-metadata {input.metadata} \
#          --group-column country \
#          --log-file {output.counts} &> {log}
#        """


#filters fasta seqs by length
#threshold set in config
#current length threshold: 29000
rule gisaid_filter_1:
    input:
        fasta = rules.gisaid_remove_duplicates.output.fasta
    params:
        min_covg = config["min_covg"],
        min_length = config["min_length"]
    output:
        fasta = config["output_path"] + "/0/gisaid.RD.UH.filt1.fasta"
    log:
        config["output_path"] + "/logs/0_gisaid_filter_1.log"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output} \
          --min-length {params.min_length} &> {log}
        """


rule gisaid_minimap2_to_reference:
    input:
        fasta = rules.gisaid_filter_1.output.fasta,
        reference = config["reference_fasta"]
    output:
        sam = config["output_path"] + "/0/gisaid.mapped.sam"
    log:
        config["output_path"] + "/logs/0_gisaid_minimap2_to_reference.log"
    threads: 8
    shell:
        """
        minimap2 -t8 -a -x asm5 {input.reference} {input.fasta} -o {output} &> {log}
        """


#head line gets column names
#tail line gets lines (except for column names) containing 'Philippines'
rule gisaid_get_variants:
    input:
        sam = rules.gisaid_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"],
        genbank_anno = config["reference_genbank_annotation"],
    output:
        variants = config["output_path"] + "/0/gisaid.variants.csv",
        global_variants = config["output_path"] + "/0/gisaid.global.variants.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_get_variants.log"
    threads: 8
    shell:
        """
        gofasta sam variants -t {threads} \
            --samfile {input.sam} \
            --reference {input.reference} \
            --genbank {input.genbank_anno} \
            --outfile {output.variants} &>> {log}

        head -n1 {output.variants} > {output.global_variants}
        tail -n+2 {output.variants} | grep -v -E "^Philippines" >> {output.global_variants}
        """


#insertion and deletion params are found but don't appear to be used?
rule gisaid_remove_insertions_and_pad:
    input:
        sam = rules.gisaid_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"]
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
        insertions = config["output_path"] + "/0/gisaid_insertions.txt",
        deletions = config["output_path"] + "/0/gisaid_deletions.txt",
    output:
        fasta = config["output_path"] + "/0/gisaid.RD.UH.filt1.mapped.fasta"
    threads: 8
    log:
        config["output_path"] + "/logs/0_gisaid_remove_insertions_and_pad.log"
    shell:
        """
        gofasta sam toMultiAlign -t {threads} -s {input.sam} -o {output.fasta} --trim --trimstart {params.trim_start} --trimend {params.trim_end} --pad &> {log}
        """


#filters fasta seqs by coverage
#threshold set in config
#current coverage threshold: 93
rule gisaid_filter_2:
    input:
        fasta = rules.gisaid_remove_insertions_and_pad.output.fasta
    output:
        fasta = config["output_path"] + "/0/gisaid.RD.UH.filt.mapped.filt2.fasta"
    params:
        min_covg = config["min_covg"]
    log:
        config["output_path"] + "/logs/0_gisaid_filter_2.log"
    shell:
        """
        datafunk filter_fasta_by_covg_and_length \
          -i {input.fasta} \
          -o {output.fasta} \
          --min-covg {params.min_covg} &> {log}
        """


#current config only masks one position: 11083
rule gisaid_mask:
    input:
        fasta = rules.gisaid_filter_2.output.fasta,
        mask = config["gisaid_mask_file"]
    output:
        fasta = config["output_path"] + "/0/gisaid.RD.UH.filt.mapped.filt2.masked.fasta",
    shell:
        """
        datafunk mask \
          --input-fasta {input.fasta} \
          --output-fasta {output.fasta} \
          --mask-file \"{input.mask}\"
        """


# TODO: use Julia closest not Python here
# rule gisaid_distance_to_WH04:
#     input:
#     output:
#     log:
#     resources: mem_per_cpu=20000
#     shell:
#         """
#         """


#distance_to_root assumes column names in metadata for ommitting data
#may need to make sure these columns are created later
#skip for now
#rule gisaid_distance_QC:
#    input:
#        fasta = rules.gisaid_mask.output.fasta,
#        metadata = rules.gisaid_remove_duplicates.output.metadata
#    log:
#        config["output_path"] + "/logs/0_gisaid_distance_QC.log"
#    output:
#        table = config["output_path"] + "/0/QC_distances.tsv",
#    resources: mem_per_cpu=20000
#    shell:
#        """
#        datafunk distance_to_root \
#          --input-fasta {input.fasta} \
#          --input-metadata {input.metadata} &> {log}
#
#        mv distances.tsv {output.table}
#        """


#skip for now
#rule gisaid_filter_on_distance_to_WH04:
#    input:
#        fasta = rules.gisaid_mask.output.fasta,
#        table = rules.gisaid_distance_QC.output.table,
#    output:
#        fasta = config["output_path"] + "/0/gisaid.RD.UH.filt.mapped.filt2.masked.filt3.fasta",
#    log:
#        config["output_path"] + "/logs/0_gisaid_filter_on_distance_to_WH04.log"
#    resources: mem_per_cpu=20000
#    run:
#        from Bio import SeqIO
#        import pandas as pd
#
#        fasta_in = SeqIO.index(str(input.fasta), "fasta")
#        df = pd.read_csv(input.table, sep='\t')
#
#        with open(str(output.fasta), 'w') as fasta_out, open(str(log), 'w') as log_out:
#            for i,row in df.iterrows():
#                sequence_name = row['sequence_name']
#                distance = row['distance_stdevs']
#                if distance < 4.0:
#                    if sequence_name in fasta_in:
#                        record = fasta_in[sequence_name]
#                        fasta_out.write('>' + record.id + '\n')
#                        fasta_out.write(str(record.seq) + '\n')
#                else:
#                    log_out.write(sequence_name + ' was filtered for having too high a distance to WH04 (' + str(distance) + ' epi-week std devs)\n')


#Would normally take fasta output from gisaid_filter_on_distance_to_WH04
rule gisaid_AA_finder:
    input:
        fasta = rules.gisaid_mask.output.fasta,
        AAs = config["AAs"]
    output:
        found = config["output_path"] + "/0/gisaid.AA_finder.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_AA_finder.log"
    shell:
        """
        datafunk AA_finder -i {input.fasta} --codons-file {input.AAs} --genotypes-table {output.found} &> {log}
        """


#input snps defined but not used
rule gisaid_add_AA_finder_result_to_metadata:
    input:
        snps = config["snps"],
        metadata = rules.gisaid_remove_duplicates.output.metadata,
        new_data = rules.gisaid_AA_finder.output.found
    output:
        metadata = config["output_path"] + "/0/gisaid.RD.UH.AAfinder.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_add_AA_finder_result_to_metadata.log"
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


#Would normally take fasta output from gisaid_filter_on_distance_to_WH04
rule gisaid_extract_lineageless:
    input:
        fasta = rules.gisaid_mask.output.fasta,
        metadata = rules.gisaid_add_AA_finder_result_to_metadata.output.metadata,
    output:
        fasta = config["output_path"] + "/0/gisaid.new.pangolin_lineages.fasta",
    log:
        config["output_path"] + "/logs/0_extract_lineageless.log"
    resources: mem_per_cpu=20000
    run:
        from Bio import SeqIO
        import pandas as pd

        fasta_in = SeqIO.index(str(input.fasta), "fasta")
        df = pd.read_csv(input.metadata)

        sequence_record = []

        with open(str(output.fasta), 'w') as fasta_out:
            if 'pango_lineage' in df.columns:
                for i,row in df.iterrows():
                    if pd.isnull(row['pango_lineage']):
                        sequence_name = row['strain']
                        if sequence_name in fasta_in:
                            if sequence_name not in sequence_record:
                                record = fasta_in[sequence_name]
                                fasta_out.write('>' + record.id + '\n')
                                fasta_out.write(str(record.seq) + '\n')
                                sequence_record.append(sequence_name)
            else:
                for i,row in df.iterrows():
                    sequence_name = row['strain']
                    if sequence_name in fasta_in:
                        if sequence_name not in sequence_record:
                            record = fasta_in[sequence_name]
                            fasta_out.write('>' + record.id + '\n')
                            fasta_out.write(str(record.seq) + '\n')
                            sequence_record.append(sequence_name)


rule gisaid_normal_pangolin:
    input:
        fasta = rules.gisaid_extract_lineageless.output.fasta
    params:
        outdir = config["output_path"] + "/0/normal_pangolin",
        tmpdir = config["output_path"] + "/0/normal_pangolin/tmp"
    output:
        lineages = config["output_path"] + "/0/normal_pangolin/lineage_report.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_normal_pangolin.log"
    resources: mem_per_cpu=20000
    shell:
        """
        pangolin {input.fasta} --tempdir {params.tmpdir} --outdir {params.outdir} --verbose > {log} 2>&1
        """


rule gisaid_add_pangolin_lineages_to_metadata:
    input:
        metadata = rules.gisaid_add_AA_finder_result_to_metadata.output.metadata,
        normal_lineages = rules.gisaid_normal_pangolin.output.lineages
    output:
        metadata = config["output_path"] + "/0/gisaid.RD.UH.SNPfinder.lineages.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_add_pangolin_lineages_to_metadata.log"
    shell:
        """
        fastafunk add_columns \
          --in-metadata {input.metadata} \
          --in-data {input.normal_lineages} \
          --index-column strain \
          --join-on taxon \
          --new-columns pango_lineage lineage_support lineages_version \
          --where-column pango_lineage=lineage lineage_support=probability lineages_version=pangoLEARN_version \
          --out-metadata {output.metadata} &> {log}
        """


#Would normally take output from gisaid_filter_on_distance_to_WH04
rule gisaid_del_finder:
    input:
        fasta = rules.gisaid_mask.output.fasta,
        dels = config["dels"]
    output:
        metadata = config["output_path"] + "/0/gisaid.del_finder.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_del_finder.log"
    shell:
        """
        datafunk del_finder \
            -i {input.fasta} \
            --deletions-file {input.dels} \
            --genotypes-table {output.metadata} &> {log}
        """


rule gisaid_add_del_finder_result_to_metadata:
    input:
        dels = config["dels"],
        metadata = rules.gisaid_add_pangolin_lineages_to_metadata.output.metadata,
        new_data = rules.gisaid_del_finder.output.metadata
    output:
        metadata = config["output_path"] + "/0/gisaid.RD.UH.SNPfinder.lineages.del_finder.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_add_del_finder_result_to_metadata.log"
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

        sed -i.bak '/Italy\/ABR-IZSGC-TE5166\/2020/s/,B.1.5,/,B.1,/g' {output.metadata} 2> {log}
        """


#Would normally take output from gisaid_filter_on_distance_to_WH04
rule gisaid_output_all_matched_metadata:
    input:
        fasta = rules.gisaid_mask.output.fasta,
        metadata = rules.gisaid_add_del_finder_result_to_metadata.output.metadata
    output:
        fasta = config["output_path"] + "/0/gisaid.all.fasta",
        metadata = config["output_path"] + "/0/gisaid.all.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_output_all_matched_metadata.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column strain \
          --filter-column \
                         strain \
                         virus \
                         gisaid_epi_isl \
                         genbank_accession \
                         date \
                         epi_day \
                         epi_week \
                         region \
                         country \
                         division \
                         location \
                         region_exposure \
                         country_exposure \
                         division_exposure \
                         segment \
                         length \
                         host \
                         age \
                         sex \
                         Nextstrain_clade \
                         GISAID_clade \
                         originating_lab \
                         submitting_lab \
                         authors \
                         url \
                         title \
                         paper_url \
                         date_submitted \
                         purpose_of_sequencing \
                         variant \
                         subsample_omit \
                         p323l \
                         n439k \
                         d614g \
                         lineage \
                         lineage_support \
                         lineages_version \
                         del_1605_3 \
          --where-column lineage=pango_lineage \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --low-memory \
          --log-file {log} \
          --restrict
        """


#Would normally take fasta output from gisaid_filter_on_distance_to_WH04
#filter columns are retained
rule gisaid_output_global_matched_metadata:
    input:
        fasta = rules.gisaid_mask.output.fasta,
        metadata = rules.gisaid_add_del_finder_result_to_metadata.output.metadata
    output:
        fasta = temp(config["output_path"] + "/0/gisaid.global.temp.fasta"),
        metadata = config["output_path"] + "/0/gisaid.global.csv"
    log:
        config["output_path"] + "/logs/0_gisaid_output_global_matched_metadata.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column strain \
          --filter-column \
                         strain \
                         virus \
                         gisaid_epi_isl \
                         genbank_accession \
                         date \
                         epi_day \
                         epi_week \
                         region \
                         country \
                         division \
                         location \
                         region_exposure \
                         country_exposure \
                         division_exposure \
                         segment \
                         length \
                         host \
                         age \
                         sex \
                         Nextstrain_clade \
                         GISAID_clade \
                         originating_lab \
                         submitting_lab \
                         authors \
                         url \
                         title \
                         paper_url \
                         date_submitted \
                         purpose_of_sequencing \
                         variant \
                         subsample_omit \
                         p323l \
                         n439k \
                         d614g \
                         lineage \
                         lineage_support \
                         lineages_version \
                         del_1605_3 \
          --where-column lineage=pango_lineage \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --low-memory \
          --restrict
        """


#strips redundancy from sequences
#equivalent to hash_seqs in rule_3?
#requires julia and julialign
#can't get julia to install
#skipping for now
#rule gisaid_collapse:
#    input:
#        fasta = rules.gisaid_exclude_uk_seqs.output.fasta,
#    output:
#        fasta = config["output_path"] + "/0/gisaid.global.collapsed.fasta",
#        tip_to_redudants = config["output_path"] + "/0/tip_to_redundants.csv",
#        redundant_to_tips = config["output_path"] + "/0/redundant_to_tips.csv",
#    log:
#        config["output_path"] + "/logs/0_gisaid_collapse.log"
#    resources: mem_per_cpu=12000
#    threads: 8
#    shell:
#        """
#        /cephfs/covid/bham/climb-covid19-jacksonb/programs/julia-1.5.0/bin/julia --threads {threads} /cephfs/covid/bham/climb-covid19-jacksonb/programs/julialign/src/collapse.jl \
#        -i {input.fasta} \
#        -o {output.fasta} &> {log}
#
#            mv tip_to_redundants.csv {output.tip_to_redudants} &>> {log}
#            mv redundant_to_tips.csv {output.redundant_to_tips} &>> {log}
#        """


#takes redundant to tips output from collapse rule
#skipping for now
#rule gisaid_get_unique_redundants:
#    input:
#        fasta = rules.gisaid_exclude_uk_seqs.output.fasta,
#        redundant_to_tips = rules.gisaid_collapse.output.redundant_to_tips,
#    output:
#        fasta = config["output_path"] + "/0/unique_redundants.fasta",
#    log:
#        config["output_path"] + "/logs/0_gisaid_get_unique_redundants.log"
#    run:
#        from Bio import SeqIO
#
#        fasta_in = SeqIO.index(str(input.fasta), "fasta")
#
#        firstline = True
#        with open(str(input.redundant_to_tips), "r") as f_in, open(str(output.fasta), "w") as f_out:
#            for line in f_in:
#                if firstline:
#                    firstline = False
#                    continue
#
#                l = line.rstrip().split(",")
#                redundant = l[0]
#                possible_tips = l[1].split("|")
#                if len(possible_tips) > 1:
#                    continue
#
#                record = fasta_in[redundant]
#
#                f_out.write(">" + record.id + "\n")
#                f_out.write(str(record.seq) + "\n")


#fasta input originally from collapsed output
rule gisaid_get_collapsed_metadata:
    input:
        fasta = rules.gisaid_mask.output.fasta,
        metadata = rules.gisaid_add_del_finder_result_to_metadata.output.metadata,
    output:
        fasta = temp(config["output_path"] + "/0/gisaid.global.collapsed.temp.fasta"),
        metadata = config["output_path"] + "/0/gisaid.global.collapsed.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_get_collapsed_metadata.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column strain \
          --filter-column \
                         strain \
                         virus \
                         gisaid_epi_isl \
                         genbank_accession \
                         date \
                         epi_day \
                         epi_week \
                         region \
                         country \
                         division \
                         location \
                         region_exposure \
                         country_exposure \
                         division_exposure \
                         segment \
                         length \
                         host \
                         age \
                         sex \
                         Nextstrain_clade \
                         GISAID_clade \
                         originating_lab \
                         submitting_lab \
                         authors \
                         url \
                         title \
                         paper_url \
                         date_submitted \
                         purpose_of_sequencing \
                         variant \
                         subsample_omit \
                         p323l \
                         n439k \
                         d614g \
                         lineage \
                         lineage_support \
                         lineages_version \
                         del_1605_3 \
          --where-column lineage=pango_lineage \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --low-memory \
          --restrict
        """


#fasta input originally from collapsed output
#unique fasta input originally used in the cat line between input.collapsed and output.fasta
#unique_fasta = rules.gisaid_get_unique_redundants.output.fasta,
rule gisaid_get_collapsed_expanded_metadata:
    input:
        collapsed_fasta = rules.gisaid_mask.output,
        metadata = rules.gisaid_add_del_finder_result_to_metadata.output.metadata,
    output:
        fasta = config["output_path"] + "/0/gisaid.global.collapsed.unique_expanded.fasta",
        metadata = config["output_path"] + "/0/gisaid.global.collapsed.unique_expanded.csv",
    log:
        config["output_path"] + "/logs/0_gisaid_get_collapsed_expanded_metadata.log"
    shell:
        """
        cat {input.collapsed_fasta} > {output.fasta} 2> {log}

        fastafunk fetch \
          --in-fasta {output.fasta} \
          --in-metadata {input.metadata} \
          --index-column strain \
          --filter-column \
                         strain \
                         virus \
                         gisaid_epi_isl \
                         genbank_accession \
                         date \
                         epi_day \
                         epi_week \
                         region \
                         country \
                         division \
                         location \
                         region_exposure \
                         country_exposure \
                         division_exposure \
                         segment \
                         length \
                         host \
                         age \
                         sex \
                         Nextstrain_clade \
                         GISAID_clade \
                         originating_lab \
                         submitting_lab \
                         authors \
                         url \
                         title \
                         paper_url \
                         date_submitted \
                         purpose_of_sequencing \
                         variant \
                         subsample_omit \
                         p323l \
                         n439k \
                         d614g \
                         lineage \
                         lineage_support \
                         lineages_version \
                         del_1605_3 \
          --where-column lineage=pango_lineage \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --low-memory \
          --restrict 2>> {log}
        """


rule check_root_pangolin_lineages:
    input:
        metadata = rules.gisaid_add_del_finder_result_to_metadata.output.metadata,
        lineage_splits = config["lineage_splits"],
    output:
        metadata = config["output_path"] + "/0/root_pangolin_lineages.txt",
    log:
        config["output_path"] + "/0/check_root_pangolin_lineages.log"
    resources: mem_per_cpu=20000
    run:
        import pandas as pd

        lineage_splits = pd.read_csv(input.lineage_splits)

        roots = []
        for i,row in lineage_splits.iterrows():
            roots.append(row["outgroup"])

        df = pd.read_csv(input.metadata)

        with open(str(output.metadata), "w") as f:
            f.write("root\tlineage\n")
            for i,row in df.iterrows():
                if row["strain"] in roots:
                    f.write(row["strain"] + "\t" + row["lineage"] + "\n")


#cp {input.counts} {params.published_counts} <- was causing issues
#commenting out doesn't work, need to remove?
#for all removed 'echo' lines, see previous version of rule 0 file
#echo "Number of non-UK sequences: $(cat {input.global_fasta} | grep '>' | wc -l)\\n" >> {log}
rule summarize_preprocess_gisaid:
    input:
        latest_fasta = rules.gisaid_remove_duplicates.input.fasta,
        deduplicated_fasta = rules.gisaid_remove_duplicates.output.fasta,
        #unify_headers_fasta = rules.gisaid_unify_headers.output.fasta,
        removed_short_fasta = rules.gisaid_filter_1.output,
        removed_low_covg_fasta = rules.gisaid_filter_2.output.fasta,
        #removed_distance_to_root_fasta = rules.gisaid_filter_on_distance_to_WH04.output.fasta,
        all_fasta = rules.gisaid_output_all_matched_metadata.output.fasta,
        all_metadata = rules.gisaid_output_all_matched_metadata.output.metadata,
        #global_fasta = rules.gisaid_exclude_uk_seqs.output.fasta,
        global_metadata = rules.gisaid_output_global_matched_metadata.output.metadata,
        #collapsed_fasta = rules.gisaid_collapse.output.fasta,
        collapsed_metadata = rules.gisaid_get_collapsed_metadata.output.metadata,
        collapsed_expanded_metadata = rules.gisaid_get_collapsed_expanded_metadata.output.metadata,
        #counts = rules.gisaid_counts_by_country.output.counts,
        variants = rules.gisaid_get_variants.output.variants,
        root_lineages = rules.check_root_pangolin_lineages.output.metadata,
    params:
        publish_path = config["publish_path"] + "/GISAID/",
        #published_counts = config["publish_path"] + "/GISAID/gisaid_counts_by_country.csv",
        published_all_fasta = config["publish_path"] + "/GISAID/gisaid.all.fasta",
        published_all_metadata = config["publish_path"] + "/GISAID/gisaid.all.csv",
        published_variants = config["publish_path"] + "/GISAID/gisaid.all.variants",
        grapevine_webhook = config["grapevine_webhook"],
        date = config["date"],
    log:
        config["output_path"] + "/logs/0_summarize_preprocess_gisaid.log"
    shell:
        """
        mkdir -p {params.publish_path}

        cp {input.all_fasta} {params.published_all_fasta}
        cp {input.all_metadata} {params.published_all_metadata}
        cp {input.variants} {params.published_variants}

        echo "Number of sequences in latest GISAID download: $(cat {input.latest_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequence after deduplicating: $(cat {input.deduplicated_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequences after removing sequences <29000bps: $(cat {input.removed_short_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "Number of sequences after mapping and removing those with <95% coverage: $(cat {input.removed_low_covg_fasta} | grep '>' | wc -l)\\n" >> {log}
        echo "> \\n" >> {log}
        echo "pangolin lineages for roots: \\n"
        cat {input.root_lineages} >> {log}
        echo "> \\n" >> {log}
        echo "> Full alignment published to {params.published_all_fasta}\\n" >> {log}
        echo "> Matched metadata published to {params.published_all_metadata}\\n" >> {log}
        echo "> Variants published to {params.published_variants}\\n" >> {log}
        """

#rule alert_sam:
#    input: rules.summarize_preprocess_gisaid.log,
#    log: config["output_path"] + "/logs/0_alert_sam.log"
#    params:
#        date = config["date"],
#    shell:
#        """
#        ln -sfn /cephfs/covid/bham/raccoon-dog/{params.date}_gisaid /cephfs/covid/bham/raccoon-dog/gisaid-latest 2> {log}
#
#        ~/.conda/envs/ben-ipc/bin/python /cephfs/covid/software/sam/public/mqtt-message.py -t 'COGUK/infrastructure/housekeeping/gisaid/status' --attr status finished &>> {log}
#        """
