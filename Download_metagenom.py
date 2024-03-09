import pandas as pd
import numpy as np
import os
import wget
import json


list_species_trash=pd.read_csv('/beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/List_of_species.txt',header=None,sep=' ')
list_species=list_species_trash[46]
for i in list_species:
   if not os.path.exists('/beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/Species/'+i):
    os.makedirs('/beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/Species/'+i)
   try:
    wget_json='https://dnazoo.s3.wasabisys.com/'+i+'README.json'
    wget_json_out='/beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/Species/'+i
    os.system(f'wget {wget_json} -O {wget_json_out}')
   except Exception:
     continue
   try: 
    with open('/beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/Species/'+i+'README.json') as json_file:
     data = json.load(json_file)
     file=data['chromlengthAssembly']['name']
     i=i[:-1]
     wget_hic='https://dnazoo.s3.wasabisys.com/'+i+'/'+file+'.hic'
     wget_hic_out='/beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/hic/'+i+'.hic'
     os.system(f'wget {wget_hic} -O {wget_hic_out}')
     wget_fasta='https://dnazoo.s3.wasabisys.com/'+i+'/'+file+'.fasta.gz'
     wget_fasta_out='/beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/genomes/'+i+'.fasta.gz'
     os.system(f'wget {wget_fasta} -O {wget_fasta_out}')
   except Exception:
      continue