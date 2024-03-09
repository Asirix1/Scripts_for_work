import os 
import pandas as pd


i=os.environ['i']

data=pd.read_csv('/beegfs/data/hpcws/ws1/popov-GENA/Insulation_with_degrons_Human/merged_nodups/Tab_delim_'+i, sep='\t',header=None)

data = data[data[1].notna()]
data = data[data[2].notna()]
data = data[data[3].notna()]
data = data[data[4].notna()]
data.to_csv('/beegfs/data/hpcws/ws1/popov-GENA/Insulation_with_degrons_Human/merged_nodups/Tab_delim_subsamp_'+i+'',header=None, index=False, sep='\t')