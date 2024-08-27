
# coding: utf-8

# In[1]:

# Based on Cluster_Code_ReleaseQuality/originalPA-Java/Biomine_beforeJavaPA_FullStagesTesting_part1.py
# This file include:
# Mapping string node name to integer and output csv and pickle file

import pandas as pd
import pickle
from datetime import datetime
import sys
import os

# In[38]:


# General environmental settings
# Toggle local/cluster codes
Local_flag = False
# Toggle version ID
version_ID = '0.0.1'
# placeholder function
def no_job_done():
  print("No job done here.")
groupType = ["COPD", "Control"]

# In[ ]:


if __name__ == '__main__':
    print('Data pre-processing started:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '--------------------------------------------------------------------------------------')
    for curr_groupType in groupType:
        # load original dataset
        my_data_name = f"{sys.argv[1]}_{curr_groupType}_MAGIC_processed_correlationMTX_toList.txt"
        print('txt file to be read in:', os.path.join(sys.argv[2], my_data_name))
        my_data = pd.read_csv(os.path.join(sys.argv[2], my_data_name), sep="\t", header=0)

        print('Shape of the input txt:')
        print(my_data.shape)
        print('\nHead 20 lines of the input txt:')
        print(my_data.head(20))
        print('\nInfo of the input txt:')
        print(my_data.info())
        print('\nTail 20 lines of the input txt:')
        print(my_data.tail(20))
        print('\nSummary of the prob column in input txt:')
        print(my_data.Freq.describe())

        # create a header cache df for later to check node name mapping
        sanity_check_df = my_data.head(20)

        print('Mapping my_data string nodes with number:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        ## map my_data string nodes with number
        node_name_list = pd.concat([my_data['Var1'], my_data['Var2']]).unique()
        node_name_mapping = { node_name_list[i]: i for i in range(len(node_name_list)) }

        print('len of node_name_mapping:', len(node_name_mapping))

        print('my_data string nodes mapped with number, saving pickle to:', f'node_name_mapping_{sys.argv[1]}_{curr_groupType}_expression_network.pickle', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        with open(f'node_name_mapping_{sys.argv[1]}_{curr_groupType}_expression_network.pickle', 'wb') as f4:
            pickle.dump(node_name_mapping, f4)

        print('replace my_data string nodes first column with number:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        my_data['Var1'] = my_data['Var1'].map(node_name_mapping)
        print('replace my_data string nodes second column with number:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        my_data['Var2'] = my_data['Var2'].map(node_name_mapping)
        
        print('\nHead 20 lines of the current df:')
        print(my_data.head(20))

        print("Begin sanity checking...")
        print(f"For {sanity_check_df.iloc[0].Var1}, the integer node name is {node_name_mapping[sanity_check_df.iloc[0].Var1]}, current df has value {my_data.head(20).iloc[0].Var1}, matches: {node_name_mapping[sanity_check_df.iloc[0].Var1] == my_data.head(20).iloc[0].Var1}")
        print(f"For {sanity_check_df.iloc[0].Var2}, the integer node name is {node_name_mapping[sanity_check_df.iloc[0].Var2]}, current df has value {my_data.head(20).iloc[0].Var2}, matches: {node_name_mapping[sanity_check_df.iloc[0].Var2] == my_data.head(20).iloc[0].Var2}")
        print(f"For {sanity_check_df.iloc[1].Var1}, the integer node name is {node_name_mapping[sanity_check_df.iloc[1].Var1]}, current df has value {my_data.head(20).iloc[1].Var1}, matches: {node_name_mapping[sanity_check_df.iloc[1].Var1] == my_data.head(20).iloc[1].Var1}")
        print(f"For {sanity_check_df.iloc[1].Var2}, the integer node name is {node_name_mapping[sanity_check_df.iloc[1].Var2]}, current df has value {my_data.head(20).iloc[1].Var2}, matches: {node_name_mapping[sanity_check_df.iloc[1].Var2] == my_data.head(20).iloc[1].Var2}")
        print(f"For {sanity_check_df.iloc[2].Var1}, the integer node name is {node_name_mapping[sanity_check_df.iloc[2].Var1]}, current df has value {my_data.head(20).iloc[2].Var1}, matches: {node_name_mapping[sanity_check_df.iloc[2].Var1] == my_data.head(20).iloc[2].Var1}")
        print(f"For {sanity_check_df.iloc[3].Var1}, the integer node name is {node_name_mapping[sanity_check_df.iloc[3].Var1]}, current df has value {my_data.head(20).iloc[3].Var1}, matches: {node_name_mapping[sanity_check_df.iloc[3].Var1] == my_data.head(20).iloc[3].Var1}")
        print(f"For {sanity_check_df.iloc[4].Var2}, the integer node name is {node_name_mapping[sanity_check_df.iloc[4].Var2]}, current df has value {my_data.head(20).iloc[4].Var2}, matches: {node_name_mapping[sanity_check_df.iloc[4].Var2] == my_data.head(20).iloc[4].Var2}")

        # print('my_data string nodes replaced with number, saving pickle:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        # with open(sys.argv[1] + '_expression_network.pickle', 'wb') as f2:
        #     pickle.dump(my_data, f2) # don't need this since 1) replace str with int is very quick and 2) by default pickle cannot serialize object larger than 4 GiB
        
        print('Saving final csv:', f'{sys.argv[1]}_{curr_groupType}_expression_network.txt', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        my_data[['Var1','Var2','Freq']].to_csv(f'{sys.argv[1]}_{curr_groupType}_expression_network.txt', index=False, sep='\t', header=False, float_format='%.18f')

    print('Data pre-processing finished:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '--------------------------------------------------------------------------------------')