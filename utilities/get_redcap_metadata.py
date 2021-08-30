import datetime
import redcap
import pandas as pd
import math
import sys
import os
from tabulate import tabulate


date_format = "%Y-%m-%d"
today = pd.to_datetime(datetime.datetime.today().strftime(date_format))


def parse_redcap_access(access_path):
    with open(access_path, 'r') as f:
        read_file = f.read().strip('\n')
        url = read_file.split(',')[0]
        key = read_file.split(',')[1]
    f.close()
    return url,key


def get_redcap_metadata(url, key, outpath, summary, consensus, dates, dag_table, dag_folder):

    proj = redcap.Project(url, key)
    proj_df = proj.export_records(format='df', forms=['case','analysis'], export_data_access_groups=True, raw_or_label='label')

    #init variables for logging
    init_cent_id_count = len(set(proj_df.index))
    no_consensus = {}
    no_dates = {}
    missing_data_by_dag = {}

    #create dataframe for splitting counts by data access group
    dag_df = pd.DataFrame(columns=['Total Records', 'Passing Records', 'Records Missing Data'], index=(proj_df['redcap_data_access_group'].unique()), dtype=object)
    for dag in dag_df.index:
        missing_data_by_dag[dag] = []   #use dags as keys in dictionary
    dag_df = dag_df.reindex(index=dag_df.index.dropna()) #in case of nans in dag column
    dag_totals = []
    for dag in dag_df.index:
        dag_totals.append(len(proj_df.loc[(proj_df['redcap_repeat_instrument'] == 'Case') & (proj_df['redcap_data_access_group'] == dag),]))
    dag_df['Total Records'] = dag_totals

    #first filter out rows without consensus
    for i in set(proj_df.index):
        if type(proj_df.consensus[i]) == float: #a lone row not containing a fasta will be filtered
            no_consensus[i] = proj_df.loc[i,'local_id']
            missing_data_by_dag[(proj_df.loc[str(i),'redcap_data_access_group'])].append([str(i),str((proj_df.loc[str(i),'local_id']))])   #add central id and local id to dag key in dict
            proj_df.drop(i, inplace=True)
            continue
        if not any(~proj_df.consensus[i].isna()): #IDs containing no fastas will be filtered
            no_consensus[i] = proj_df.loc[i,'local_id']
            missing_data_by_dag[(proj_df.loc[str(i),'redcap_data_access_group'][0])].append([str(i),(proj_df.loc[str(i),'local_id'][0])])
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
                no_dates[i] = proj_df.loc[(i,'Case',j),'local_id']
                missing_data_by_dag[(proj_df.loc[str(i),'redcap_data_access_group'][0])].append([str(i),(proj_df.loc[str(i),'local_id'][0])])
                proj_df.drop(proj_df.loc[(i,slice(None),j),:].index, inplace=True) #drop case and analysis repeat instance

    #finish filling in dag_df
    temp_df = proj_df.reset_index()
    dag_filled = []
    for dag in dag_df.index:
        dag_filled.append(len(temp_df.loc[(temp_df['redcap_repeat_instrument'] == 'Case') & (temp_df['redcap_data_access_group'] == dag),]))
    dag_df['Passing Records'] = dag_filled
    dag_df['Records Missing Data'] = dag_df['Total Records'] - dag_df['Passing Records']
    dag_df.loc['Totals'] = dag_df.sum()

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

    #logging
    with open(summary, 'w') as sum_handle:
        sum_handle.write("Number of unique central IDs from case and analysis instruments in redcap database: " + str(init_cent_id_count) + "\\n")
        sum_handle.write("Number of records after filtering based on consensus sequence: " + str(init_cent_id_count - len(no_consensus)) + "\\n")
        sum_handle.write("Number of records after filtering based on date columns: " + str(init_cent_id_count - len(no_consensus) - len(no_dates)) + "\\n")
    sum_handle.close()

    with open(consensus, 'w') as con_handle:
        con_handle.write("central_id,local_id\n")
        for key,val in no_consensus.items():
            cent_id = key
            local_id = val
            if type(val) == pd.core.series.Series: #there is currently an odd case where central id 11-2 has multiple analysis repeats without a consensus seq, probably for testing
                for i in val:
                    if type(i) == str:
                        local_id = i
                        break
                    else:
                        continue
            con_handle.write(str(cent_id) + "," + str(local_id) + "\n")
    con_handle.close()

    with open(dates, 'w') as dates_handle:
        dates_handle.write("central_id,local_id\n")
        for key,val in no_dates.items():
            cent_id = key
            local_id = val
            dates_handle.write(str(cent_id) + "," + str(local_id) + "\n")
    dates_handle.close()

    #save dag_df as ascii table to be sent to slack via webhook
    with open(dag_table, 'w') as dag_handle:
        str_table = tabulate(dag_df, headers='keys', tablefmt='psql').replace('\n', '\\n\n') #replace needed for formatting for sending table to slack via webhook
        dag_handle.write(str_table)

    #write records with missing data to files based on dag
    folder = dag_folder
    if not os.path.exists(os.path.dirname(folder)): #make directory if it doesn't exist
        os.makedirs(os.path.dirname(folder))
    for dag in missing_data_by_dag.keys():
        with open((folder + str(dag) + '.csv'), 'w') as dag_handle:
            dag_handle.write('central_id,local_id\n')
            dag_handle.writelines([str(i) + "," + str(j) + "\n" for i,j in missing_data_by_dag[dag]])


if __name__ == "__main__":
    input = sys.argv[1]
    output_file = sys.argv[2]
    output_summary = sys.argv[3]
    output_no_consensus_table = sys.argv[4]
    output_no_dates_table = sys.argv[5]
    output_dag_summary = sys.argv[6]
    output_dag_directory = sys.argv[7]
    url,key = parse_redcap_access(input)
    get_redcap_metadata(url, key, output_file, output_summary, output_no_consensus_table, output_no_dates_table, output_dag_summary, output_dag_directory)
