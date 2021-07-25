import datetime
import redcap
import pandas as pd
import math
import sys


date_format = "%Y-%m-%d"
today = pd.to_datetime(datetime.datetime.today().strftime(date_format))


def parse_redcap_access(access_path):
    with open(access_path, 'r') as f:
        read_file = f.read().strip('\n')
        url = read_file.split(',')[0]
        key = read_file.split(',')[1]
    f.close()
    return url,key


def get_redcap_metadata(url, key, outpath):

    proj = redcap.Project(url, key)
    proj_df = proj.export_records(format='df', forms=['case','analysis'], raw_or_label='label')

    #first filter out rows without consensus
    for i in set(proj_df.index):
        if type(proj_df.consensus[i]) == float: #a lone row not containing a fasta will be filtered
            proj_df.drop(i, inplace=True)
            continue
        if not any(~proj_df.consensus[i].isna()): #IDs containing no fastas will be filtered
            proj_df.drop(i, inplace=True)

    #set multiindex to create unique keys
    proj_df.set_index(['redcap_repeat_instrument', 'redcap_repeat_instance'], append=True, inplace=True)

    #filter repeat instances with no dates
    #loop through central IDs
    for i in proj_df.index.levels[0]:
        temp_df = proj_df.loc[(i,'Case',slice(None)),:] #case specified since number of repeat instances may differ between instruments
        temp_index = set(temp_df.index.get_level_values('redcap_repeat_instance'))
        #loop through repeat instances
        for j in temp_index:
            dates_df = proj_df.loc[(i,'Case',j),:]
            if not any(~dates_df[['date_collected', 'date_received']].isna()): #if both date columns are NaN
                proj_df.drop(proj_df.loc[(i,slice(None),j),:].index, inplace=True) #drop case and analysis repeat instance

    #filter repeat instances based on date
    #loop through central IDs
    for i in set([proj_df.index[i][0] for i in range(len(proj_df.index))]): #using this ugly code because proj_df.index.levels[0] gives an outdated central ID set for some reason:
        temp_df = proj_df.loc[(i,'Case',slice(None)),:] #case specified
        temp_index = set(temp_df.index.get_level_values('redcap_repeat_instance'))
        #check if there are repeat instances to filter
        if len(temp_index) > 1:
            delta_max = datetime.timedelta(days=0) #initial timedelta
            inst_to_keep = 1
            for j in temp_index:
                rep_inst_df = proj_df.loc[(i,'Case',j),:]
                if math.isnan(rep_inst_df['date_collected']): #if date_collected is NaN, then date_received is used
                    date = pd.to_datetime(rep_inst_df['date_received'], format=date_format)
                    new_delta = today-date
                    if new_delta>delta_max: #if date of repeat is older, it will be kept
                        delta_max = new_delta
                        inst_to_keep = j
                else:
                    date = pd.to_datetime(rep_inst_df['date_collected'], format=date_format)
                    new_delta = today-date
                    if new_delta>delta_max: #if date of repeat is older, it will be kept
                        delta_max = new_delta
                        inst_to_keep = j
            for j in temp_index: #loop through repeat instances and drop rows that aren't the instance to keep
                if j == inst_to_keep:
                    continue
                else:
                    proj_df.drop(proj_df.loc[(i,slice(None),j),:].index, inplace=True) #drop case and analysis repeat instance

    #split dataframe into case and analysis instruments and remove repeat_instrument columns
    case_df = proj_df.loc[(slice(None),'Case',slice(None)),:'case_complete'].droplevel('redcap_repeat_instrument')
    analysis_df = proj_df.loc[(slice(None),'Analysis',slice(None)),'consensus':].droplevel('redcap_repeat_instrument')
    #merge into new dataframe such that case and analysis info is on one row
    merged_df = case_df.join(analysis_df)

    merged_df.to_csv(outpath, sep=',')


if __name__ == "__main__":
    input = sys.argv[1]
    output_file = sys.argv[2]
    url,key = parse_redcap_access(input)
    get_redcap_metadata(url, key, output_file)
