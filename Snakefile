##### Import #####
import datetime
import pytz

##### Configuration #####
configfile: workflow.current_basedir + "/config.yaml"

if config.get("output_path"):
    config["output_path"] = config["output_path"].rstrip('/')

if config.get("publish_path"):
    config["publish_path"] = config["publish_path"].rstrip('/')
config["publish_path"] = os.path.abspath(config["publish_path"])

if config.get("export_path"):
    config["export_path"] = config["export_path"].rstrip('/')
config["export_path"] = os.path.abspath(config["export_path"])

if config.get("json_path"):
    config["json_path"] = config["json_path"].rstrip('/')
config["json_path"] = os.path.abspath(config["json_path"])

if not config.get("date"):
    utc_now = pytz.utc.localize(datetime.datetime.utcnow())
    config["date"] = utc_now.astimezone(pytz.timezone("Europe/London")).strftime('%d-%m-%Y')

if not config.get("time"):
    utc_now = pytz.utc.localize(datetime.datetime.utcnow())
    config["time"] = utc_now.astimezone(pytz.timezone("Europe/London")).strftime('%H:%M %Z')

##### Target rules #####

rule all:
    input:
#        config["output_path"] + "/logs/6_summarize_publish.log",
        config["output_path"] + "/logs/5_summarise_define_redcap_lineages_and_cut_out_trees.log",
        config["output_path"] + "/logs/4_summarise_make_trees.log",
        config["output_path"] + "/logs/3_summarise_combine_gisaid_and_redcap.log",
        config["output_path"] + "/logs/2_summarise_pangolin_lineage_typing.log",
        config["output_path"] + "/logs/1_summarise_preprocess_redcap.log",
        config["output_path"] + "/logs/0_summarise_preprocess_gisaid.log",
#        config["output_path"] + "/snakejunk/all"

#clean_up doesn't appear to be used at this stage
#rule clean_up:
#    input:
#        #config["output_path"] + "/logs/6_summarize_publish.log",
#        config["output_path"] + "/logs/0_summarize_preprocess_gisaid.log",
#    output:
#        config["output_path"] + "/snakejunk/all"
#    shell:
#        """
#        mkdir -p {output}
#        mv slurm-*.out *_data.json {output}/
#        for file in pre trace default.profraw
#        do
#          if [ -f "$file" ]
#          then
#            rm $file
#          fi
#        done
#        """

##### Modules #####
include: "rules/0_preprocess_gisaid.smk"
include: "rules/1_preprocess_redcap.smk"
include: "rules/2_pangolin_lineage_typing.smk"
include: "rules/3_combine_gisaid_and_redcap.smk"
include: "rules/4_make_trees.smk"
include: "rules/5_define_redcap_lineages_and_cut_out_trees.smk"
#include: "rules/6_treetime.smk"
#include: "rules/7_publish.smk"
