#todays_date = datetime.datetime.strptime(params.date, '%Y-%m-%d').date()
#used if seq_rec.id in outgroups from just 'seq_rec',
#as comparisons using a SeqRecord object aren't implemented
#
#doesn't filter outgroup sequences
#should we be excluding sequences that are "too old"?
#current window paramater is 180 days
#might set to arbitrarily large number since we have a relatively low number of sequences
rule filter_by_date:
    input:
        fasta = rules.redcap_output_lineage_table.output.fasta,
        metadata = rules.redcap_output_lineage_table.output.metadata,
        lineage_splits = config["lineage_splits"]
    params:
        date = config["date"],
        time_window = config["time_window"]
    output:
        fasta = config["output_path"] + "/3/redcap.filter_by_date.fasta"
    log:
        config["output_path"] + "/logs/3_filter_by_date.log"
    run:
        import datetime
        from Bio import SeqIO
        import csv

        outgroups = []
        with open(str(input.lineage_splits), "r") as outgroup_handle:
            line = outgroup_handle.readline()
            while line:
                try:
                    outgroup = line.strip().split(",")[-1]
                    outgroups.append(outgroup)
                except:
                    continue
                line = outgroup_handle.readline()
                

        indexed_fasta = SeqIO.index(str(input.fasta), "fasta")
        
        window = datetime.timedelta(int(params.time_window))
        todays_date = datetime.datetime.strptime(str(datetime.date.today()), '%Y-%m-%d').date()
        
        with open(str(input.metadata), 'r', newline = '') as csv_in, \
             open(str(output.fasta), "w") as fasta_out:

            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")

            for row in reader:
                if row["strain"] not in indexed_fasta:
                    print("%s not in fasta" %row["strain"])
                    continue
                seq_rec = indexed_fasta[row["strain"]]
                if seq_rec.id in outgroups:
                    fasta_out.write(">" + seq_rec.id + "\n")
                    fasta_out.write(str(seq_rec.seq) + "\n")
                    continue

                try:
                    date = datetime.datetime.strptime(row["sample_date"], '%Y-%m-%d').date()
                except:
                    continue

                if (todays_date - window) > date:
                    continue

                fasta_out.write(">" + seq_rec.id + "\n")
                fasta_out.write(str(seq_rec.seq) + "\n")


#removes identical sequences from fastas
#why would this not be done sooner?
#hashmap would contain names of fastas for identical sequences
#the sequence name kept appears arbitrary; the first encountered
#outgroup fastas are retained
#if there is more than one redundant name, they will be '|' delimited
rule redcap_hash_seqs:
    input:
        fasta = rules.filter_by_date.output.fasta,
        lineage_splits = config["lineage_splits"]
    output:
        fasta = config["output_path"] + "/3/redcap.hashed.fasta",
        metadata = config["output_path"] + "/3/redcap.hashmap.csv",
    log:
        config["output_path"] + "/logs/3_redcap_hash_seqs.log",
    resources: 
        mem_per_cpu=8000
    run:
        from Bio import SeqIO

        outgroups = []
        with open(input.lineage_splits, "r") as outgroup_handle:
            line = outgroup_handle.readline()
            while line:
                try:
                    outgroup = line.strip().split(",")[-1]
                    outgroups.append(outgroup)
                except:
                    continue
                line = outgroup_handle.readline()

        input_fasta = SeqIO.index(str(input.fasta), "fasta")

        hash_dict = {}

        with open(str(input.fasta), "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                id = record.id
                seq = str(record.seq)

                if seq in hash_dict:
                    hash_dict[seq] = hash_dict[seq] + [id]
                else:
                    hash_dict[seq] = [id]

        with open(str(output.fasta), "w") as fasta, open(str(output.metadata), "w") as metadata:
            metadata.write("tip,redundant\n")

            for key, value in hash_dict.items():
                if len(value) == 1:

                    r = input_fasta[value[0]]

                    fasta.write(">" + r.id + "\n")
                    fasta.write(str(r.seq) + "\n")

                elif len(value) > 1:
                    r = None

                    for id in value:
                        if id in outgroups:
                            r = input_fasta[id]
                            fasta.write(">" + r.id + "\n")
                            fasta.write(str(r.seq) + "\n")
                            value.remove(id)

                    if not r:
                        r = input_fasta[value[0]]
                        fasta.write(">" + r.id + "\n")
                        fasta.write(str(r.seq) + "\n")
                        value.remove(value[0])

                    metadata.write(r.id + ",")
                    metadata.write("|".join(value) + "\n")


#where-column to match with gisaid metadata
rule redcap_output_hashed_lineage_table:
    input:
        fasta = rules.redcap_hash_seqs.output.fasta,
        metadata = rules.redcap_output_lineage_table.output.metadata
    params:
        country_code = config["country_code"]
    output:
        fasta = temp(config["output_path"] + "/3/redcap.hashed.temp.fasta"),
        metadata = config["output_path"] + "/3/redcap.hashed.lineages.csv"
    log:
        config["output_path"] + "/logs/3_redcap_output_hashed_lineage_table.log"
    shell:
          """
          fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column strain \
          --filter-column strain country adm1 adm2 \
                          sample_date epi_week \
                          lineage {params.country_code}_lineage \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --low-memory \
          --restrict
          """


#changed the where-column option a bit
#input was originally config["GISAID_background_fasta"] and
#config["GISAID_background_metadata"].
#not sure what is meant by 'background' data.
rule gisaid_output_lineage_table:
    input:
        fasta = rules.gisaid_output_all_matched_metadata.output.fasta,
        metadata = rules.gisaid_output_all_matched_metadata.output.metadata
    params:
        country_code = config["country_code"]
    output:
        fasta = config["output_path"] + "/3/gisaid.matched.fasta",
        metadata = config["output_path"] + "/3/gisaid.matched.lineages.csv"
    log:
        config["output_path"] + "/logs/3_gisaid_output_lineage_table.log"
    shell:
        """
        fastafunk fetch \
          --in-fasta {input.fasta} \
          --in-metadata {input.metadata} \
          --index-column strain \
          --filter-column strain country adm1 adm2 \
                          sample_date epi_week \
                          lineage {params.country_code}_lineage \
          --where-column adm1=division adm2=location sample_date=date \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --log-file {log} \
          --low-memory \
          --restrict
        """


#previous stage define but not used
#uses hashed fastas
#passed to split_by_lineages
rule combine_gisaid_and_redcap:
    input:
        previous_stage = config["output_path"] + "/logs/2_summarise_pangolin_lineage_typing.log",
        gisaid_fasta = rules.gisaid_output_lineage_table.output.fasta,
        gisaid_metadata = rules.gisaid_output_lineage_table.output.metadata,
        redcap_fasta = rules.redcap_hash_seqs.output.fasta,
        redcap_metadata = rules.redcap_output_hashed_lineage_table.output.metadata
    output:
        fasta = config["output_path"] + "/3/redcap_gisaid.fasta",
        metadata = config["output_path"] + "/3/redcap_gisaid.lineages.csv"
    log:
        config["output_path"] + "/logs/3_combine_gisaid_and_redcap.log"
    shell:
        """
        fastafunk merge \
          --in-fasta {input.gisaid_fasta} {input.redcap_fasta} \
          --in-metadata {input.gisaid_metadata} {input.redcap_metadata} \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --index-column strain \
          --low-memory \
          --log-file {log}
        """


#previous stage define but not used
#uses unhashed fastas
rule combine_gisaid_and_redcap_expanded:
    input:
        previous_stage = config["output_path"] + "/logs/2_summarise_pangolin_lineage_typing.log",
        gisaid_fasta = rules.gisaid_output_lineage_table.output.fasta,
        gisaid_metadata = rules.gisaid_output_lineage_table.output.metadata,
        redcap_fasta = rules.redcap_output_lineage_table.output.fasta,
        redcap_metadata = rules.redcap_output_lineage_table.output.metadata
    output:
        fasta = config["output_path"] + "/3/redcap_gisaid.expanded.fasta",
        metadata = config["output_path"] + "/3/redcap_gisaid.lineages.expanded.csv"
    log:
        config["output_path"] + "/logs/3_combine_gisaid_and_redcap_expanded.log"
    shell:
        """
        fastafunk merge \
          --in-fasta {input.gisaid_fasta} {input.redcap_fasta} \
          --in-metadata {input.gisaid_metadata} {input.redcap_metadata} \
          --out-fasta {output.fasta} \
          --out-metadata {output.metadata} \
          --index-column strain \
          --low-memory \
          --log-file {log}
        """


#removed the following echo lines:
#        echo '{{"text":"' > {params.json_path}/3_data.json
#        echo "*Step 3: Combine {params.date} COG-UK and GISAID data complete*\\n" >> {params.json_path}/3_data.json
#        cat {log} >> {params.json_path}/3_data.json
#        echo '"}}' >> {params.json_path}/3_data.json
#        echo "webhook {params.grapevine_webhook}"
#        curl -X POST -H "Content-type: application/json" -d @{params.json_path}/3_data.json {params.grapevine_webhook}
#as they weren't necessary
rule summarise_combine_gisaid_and_redcap:
    input:
        redcap_hashed_fasta = rules.redcap_hash_seqs.output.fasta,
        fasta = rules.combine_gisaid_and_redcap.output.fasta,
        metadata = rules.combine_gisaid_and_redcap.output.metadata,
        full_fasta = rules.combine_gisaid_and_redcap_expanded.output.fasta,
        full_metadata = rules.combine_gisaid_and_redcap_expanded.output.metadata
    params:
        grapevine_webhook = config["grapevine_webhook"],
        json_path = config["json_path"],
        date = config["date"]
    log:
        config["output_path"] + "/logs/3_summarise_combine_gisaid_and_redcap.log"
    shell:
        """
        echo "> Number of sequences in total REDCap and GISAID matched files for later steps: $(cat {input.full_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences in reduced REDCap fasta: $(cat {input.redcap_hashed_fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> Number of sequences in collapsed REDCap and GISAID matched files for tree building: $(cat {input.fasta} | grep ">" | wc -l)\\n" &>> {log}
        echo "> \\n" &>> {log}
        """
