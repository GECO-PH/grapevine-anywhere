# grapevine-anywhere
The original grapevine pipeline (https://github.com/COG-UK/grapevine), was created for the CoG-UK project for the processing and analysis of UK SARS-COV-2 sequence data.
Here, it has been adapted for the use in the GECO project (https://www.geco-seqlab.org) for Philippines SARS-COV-2 sequence data.
grapevine-anywhere can also be used for other countries besides the Philippines. This can be set in 'config.yaml'.

## Install
```
git clone https://github.com/josephhughes/grapevine-anywhere.git
cd grapevine-anywhere
mamba env create -f environment.yml
conda activate grapevine
```
Note that, while conda can be used to create the environment, mamba is currently much faster.

## Configure

### `Required config keys`
The following are keys in the config.yaml file which require values for the pipeline to run:

`country`
> *Country of interest. Where the fastas stored in REDCap database were sequenced. E.g. philippines*

`country_code`
> *Abreviated form of country. E.g. ph*

`gisaid_meta`
> *Path to input gisaid metadata.*

`gisaid_fasta`
> *Path to input gisaid fasta sequences.*

`gisaid_mask_file`
> *Path to gisaid mask file. Used to mask sites in gisaid fasta sequences.*

`redcap_access`
> *Path to file containing REDCap database API URL and a redcap API token for that database.*
> *The format of the file should be such that it contains only the URL and token on a single line, separated by a comma.*
> *E.g. https://geco.ritm-edc.net/redcap/api/,TOKEN*

`redcap_mask_file`
> *Path to REDCap maskfile. Used to mask sites in REDCap fasta sequences.*

`redcap_omissions`
> *Path to REDCap omissions file. Used to remove specific fasta sequences.*

`output_path`
> *Path for outputs.*

`previous_outputs`
> *Path to outputs from previous runs.*

`publish_path`

`export_path`

`previous_gisaid_lineages`
> *Path to the gisaid pangolin lineage report from the previous run. If pipeline has not been run yet, this will be a dummy file.*

`global_lineages`

`min_covg`
> *Minimum coverage threshold represented by percentage. Fasta sequences will be filtered if they fall below this threshold.*
> *Default is set to 93%.*

`min_length`
> *Minimum length threshold in base pairs. Fasta sequences will be filtered if they fall below this threshold.*
> *Default is set to 29000bp.*

`reference_fasta`
> *Path to reference genome used for mapping.*

`reference_genbank_annotation`
> *Path to genbank annotation used for determining variants.*

`trim_start`
`trim_end`
> *Parameters used to set start and end co-ordinates for trimming.*
> *Default is set to 265bp and 29674bp, respectively.*

`snps`
> *Path to file used for snp_finder. Identifies SNPs of interest from input fastas.*

`dels`
> *Path to file used for del_finder. Identifies deletions of interest from input fastas.*

`AAs`
> *Path to file used for AA_finder. Identifies amino acid sites of interest from input fasta.*

`lineage_splits`
> *Path to file used to split sequences into groups to create subtrees.*

`lineage_aliases`
> *Path to file which defines aliases to use for lineages.*

`outgroups_metadata`
> *Path to outgroups metadata. Outgroups are added to gisaid data.*

`outgroups_fasta`
> *Path to outgroups fasta sequences. Outgroups are added to gisaid data.*

`guide_tree`
> *Path to guide tree file used to graft subtrees together. Guide tree consists of outgroups.*

`date`
> *Current date. Can be left blank and it will be determined in Snakefile. Used for summaries sent to slack via webhook.*

`time`
> *Current time. Can be left blank and it wil be determined in Snakefile. Used for summaries sent to slack via webhook.*

`time_window`
> *Parameter for filtering sequences by date. Window should be in number of days. If date of sequence is older than number of days in time window, it will be filtered.*
> *Set to an arbitrarily large value '99999' by default so that sequences will not be filtered by date.*

`webhook`
> *Path to file containing webhook URL. Used for sending summaries to slack.*
> *The format of the file should be such that it only contains the URL, which should look like: https://hooks.slack.com/services/T00000000/B00000000/XXXXXXXXXXXXXXXXXXXXXXXX.*
> *A good tutorial for this can be found here: https://api.slack.com/messaging/webhooks*


## Pipeline

### `Overview`

0) Process GISAID data.

1) Process REDCap data.

2) Clean and filter.

3) Align to reference.

4) Lineage type using PANGOLIN.

5) Merge GISAID and REDCap alignments.

6) Split into lineages based on Pangolin typing to relieve some of the burden on tree building. (optional)

7) Build a Fasttree maximum likelihood tree (for each lineage if so defined):

```
FastTreeMP -nosupport -nt <sequences.fasta> > <output.tree>
```

8) If split, graft together the subtrees into a complete tree.

9) Define and extract clusters from country of interest (e.g. Philippines) using `cluster-funk`.

10) Build reports.


### `Rules`

`0 - preprocess_gisaid`
> Input GISAID data is cleaned, aligned and lineage typed (for sequences which don't already have a lineage designation)..

`1 - preprocess_redcap`
> Input REDCap data is cleaned and aligned.

`2 - pangolin_lineage_typing`
> Cleaned and aligned REDCap data is lineage typed.

`3 - combine_gisaid_and_redcap`
> Processed GISAID and REDCap data are combined.

`4 - make_trees`
> Sequences are split into outgroups, if they are provided, and a phylogenetic tree is generated.
> Ancestral state reconstruction is also carried out.

`5 - define_redcap_lineages_and_cut_out_trees`
> Tree is annotated with clusters from ancestral state reconstruction.
> Sub-trees with >=3 sequences in a cluster are cut out and phylotyped.
> Full tree is annotated with phylotypes.


## Outputs


## Acknowledgements

Grapevine uses the following tools:

```
bioconda
conda
biopython
minimap2
python
snakemake
pandoc
fasttree
gotree
TreeTime (https://github.com/neherlab/treetime)
Pangolin (https://github.com/hCoV-2019/pangolin)
```
