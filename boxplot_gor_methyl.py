import pandas as pd 
import numpy as np 
import random
import tqdm 

window=10000

methyl_gpc=pd.read_csv('GalGal_T2T_CpG.txt',sep='\t',header=None,names=['chr','pos','orient','meth','unmeth','din','trip'])
methyl_cpg=pd.read_csv('GalGal_T2T_GpC.txt',sep='\t',header=None,names=['chr','pos','orient','meth','unmeth','din','trip'])
big_tads=pd.read_csv('TADs_over_2MB.bed',sep='\t',header=None,names=['chr','start','end'])
small_tads=pd.read_csv('TADs_less_2MB.bed',sep='\t',header=None,names=['chr','start','end'])

prop_of_meth_gpc=[]
prop_of_meth_cpg=[]
for i, reg in big_tads.iterrows():
    print(reg)
    print(reg['end'])
    center=(reg['start']+reg['end'])/2
    reads_in_reg_gpc=methyl_gpc[(methyl_gpc['chr']==reg['chr']) & (methyl_gpc['pos'] >= reg['start']) & (methyl_gpc['pos'] <= reg['end'])]
    reads_in_reg_cpg=methyl_gpc[(methyl_cpg['chr']==reg['chr']) & (methyl_cpg['pos'] >= reg['start']) & (methyl_cpg['pos'] <= reg['end'])]
    perc_of_meth_gpc=reads_in_reg_gpc['meth'].sum()/(reads_in_reg_gpc['meth'].sum()+reads_in_reg_gpc['unmeth'].sum())
    perc_of_meth_cpg=reads_in_reg_cpg['meth'].sum()/(reads_in_reg_cpg['meth'].sum()+reads_in_reg_cpg['unmeth'].sum())
    '''
    if pd.isnull(perc_of_meth_gpc):
        print(reads_in_reg_gpc)
        perc_of_meth_gpc=0
    if pd.isnull(perc_of_meth_cpg):
        print(reads_in_reg_cpg)
        perc_of_meth_cpg=0
    '''
    #assert(perc_of_meth_gpc <= 1)
    #assert(perc_of_meth_cpg <= 1)
    prop_of_meth_gpc.append(perc_of_meth_gpc)
    prop_of_meth_cpg.append(perc_of_meth_cpg)
print(prop_of_meth_gpc)
print(prop_of_meth_cpg)
prop_of_meth_gpc_small=[]
prop_of_meth_cpg_small=[]
for i,reg in small_tads.iterrows():
    center=(reg['start']+reg['end'])/2
    reads_in_reg_gpc=methyl_gpc[(methyl_gpc['chr']==reg['chr']) & (methyl_gpc['pos'] >= reg['start']) & (methyl_gpc['pos'] <= reg['end'])]
    reads_in_reg_cpg=methyl_gpc[(methyl_cpg['chr']==reg['chr']) & (methyl_cpg['pos'] >= reg['start']) & (methyl_cpg['pos'] <= reg['end'])]
    perc_of_meth_gpc=reads_in_reg_gpc['meth'].sum()/(reads_in_reg_gpc['meth'].sum()+reads_in_reg_gpc['unmeth'].sum())
    perc_of_meth_cpg=reads_in_reg_cpg['meth'].sum()/(reads_in_reg_cpg['meth'].sum()+reads_in_reg_cpg['unmeth'].sum())
    '''
    if pd.isnull(perc_of_meth_gpc):
        print(reads_in_reg_gpc)
        perc_of_meth_gpc=0
    if pd.isnull(perc_of_meth_cpg):
        print(reads_in_reg_cpg)
        perc_of_meth_cpg=0
    '''
    #assert(perc_of_meth_gpc <= 1)
    #assert(perc_of_meth_cpg <= 1)
    prop_of_meth_gpc_small.append(perc_of_meth_gpc)
    prop_of_meth_cpg_small.append(perc_of_meth_cpg)   
print(prop_of_meth_gpc_small)
print(prop_of_meth_cpg_small)

data=pd.DataFrame({'Experiment_over_2MB_TADS_GpC':prop_of_meth_gpc,'Experiment_over_2MB_TADS_CpG':prop_of_meth_cpg})
data_1=pd.DataFrame({'Experiment_under_2MB_TADS_GpC':prop_of_meth_gpc_small,'Experiment_under_2MB_TADS_CpG':prop_of_meth_cpg_small})
data.to_csv('Data_for_methylation_in_regions_big_TADs.txt',sep='\t',header=True,index=False)
data_1.to_csv('Data_for_methylation_in_regions_small_TADS.txt',sep='\t',header=True,index=False)

def calc_random_control(bed):
 chroms = pd.unique(bed['chr'])
 random_control = pd.DataFrame(columns=["chr", "start", "end"])
 chr_size = pd.read_csv('/beegfs/data/hpcws/ws1/popov-popov_gnomad/Chicken_oocytes/GalGal_T2T.chrom.sizes', sep='\t',header=None)
 for chrom in chroms:  
    starts = bed.query('chr==@chrom')["start"].values
    ends = bed.query('chr==@chrom')["end"].values
    region_lengths = ends - starts
    new_starts=np.random.randint(low=np.absolute(min(region_lengths)),high=(max(region_lengths)+chr_size[chr_size[0]==chrom][1]),size=starts.size)# well done!
    np.random.shuffle(region_lengths)

    new_ends=new_starts+region_lengths

     
    #assert inter_region_lengths.min() >= 0 # good point!
    '''
    random.shuffle(region_lengths)
    random.shuffle(inter_region_lengths)
    
    new_starts = np.insert(np.cumsum(region_lengths[:-1] + inter_region_lengths)+starts.min(), 0, starts.min())

    #assert (new_starts[1:] - new_starts[:-1]).min() >= 0
    
    new_ends = new_starts +  region_lengths
    '''
    chroms = [chrom] * len(new_starts)
    
 

    random_for_chrom = pd.DataFrame({'chr': chroms, 'start': new_starts, 'end': new_ends})
    random_control = pd.concat([random_control, random_for_chrom], axis=0, ignore_index=True)
    random_control = random_control.sort_values(by=['chr', 'start'])
 return random_control
control_big_TADS=calc_random_control(big_tads)
control_small_TADS=calc_random_control(small_tads)

prop_of_meth_gpc=[]
prop_of_meth_cpg=[]

print(len(control_big_TADS))
print(len(control_small_TADS))
for i, reg in tqdm.tqdm(control_big_TADS.iterrows()):
    center=(reg['start']+reg['end'])/2
    reads_in_reg_gpc=methyl_gpc[(methyl_gpc['chr']==reg['chr']) & (methyl_gpc['pos'] >= reg['start']) & (methyl_gpc['pos'] <= reg['end'])]
    reads_in_reg_cpg=methyl_gpc[(methyl_cpg['chr']==reg['chr']) & (methyl_cpg['pos'] >= reg['start']) & (methyl_cpg['pos'] <= reg['end'])]
    perc_of_meth_gpc=reads_in_reg_gpc['meth'].sum()/(reads_in_reg_gpc['meth'].sum()+reads_in_reg_gpc['unmeth'].sum())
    perc_of_meth_cpg=reads_in_reg_cpg['meth'].sum()/(reads_in_reg_cpg['meth'].sum()+reads_in_reg_cpg['unmeth'].sum())
    '''
    if pd.isnull(perc_of_meth_gpc):
        print(reads_in_reg_gpc)
        perc_of_meth_gpc=0
    if pd.isnull(perc_of_meth_cpg):
        print(reads_in_reg_cpg)
        perc_of_meth_cpg=0
    '''
    #assert(perc_of_meth_gpc <= 1)
    #assert(perc_of_meth_cpg <= 1)
    prop_of_meth_gpc.append(perc_of_meth_gpc)
    prop_of_meth_cpg.append(perc_of_meth_cpg)
print(prop_of_meth_gpc)
print(prop_of_meth_cpg)

data=pd.DataFrame({'Control_over_2MB_TADS_GpC':prop_of_meth_gpc,'Control_over_2MB_TADS_CpG':prop_of_meth_cpg})
data.to_csv('Control_data_for_methylation_in_regions_big_TADs.txt',sep='\t',header=True,index=False)

prop_of_meth_gpc_small=[]
prop_of_meth_cpg_small=[]

for i,reg in tqdm.tqdm(control_small_TADS.iterrows()):
    center=(reg['start']+reg['end'])/2
    reads_in_reg_gpc=methyl_gpc[(methyl_gpc['chr']==reg['chr']) & (methyl_gpc['pos'] >= reg['start']) & (methyl_gpc['pos'] <= reg['end'])]
    reads_in_reg_cpg=methyl_gpc[(methyl_cpg['chr']==reg['chr']) & (methyl_cpg['pos'] >= reg['start']) & (methyl_cpg['pos'] <= reg['end'])]

    perc_of_meth_gpc=reads_in_reg_gpc['meth'].sum()/(reads_in_reg_gpc['meth'].sum()+reads_in_reg_gpc['unmeth'].sum())
    perc_of_meth_cpg=reads_in_reg_cpg['meth'].sum()/(reads_in_reg_cpg['meth'].sum()+reads_in_reg_cpg['unmeth'].sum())
    '''
    if pd.isnull(perc_of_meth_gpc):
        print(reads_in_reg_gpc)
        perc_of_meth_gpc=0
    if pd.isnull(perc_of_meth_cpg):
        print(reads_in_reg_cpg)
        perc_of_meth_cpg=0
    '''
    #assert(perc_of_meth_gpc <= 1)
    #assert(perc_of_meth_cpg <= 1)
    prop_of_meth_gpc_small.append(perc_of_meth_gpc)
    prop_of_meth_cpg_small.append(perc_of_meth_cpg)   
print(prop_of_meth_gpc_small)
print(prop_of_meth_cpg_small)

data_1=pd.DataFrame({'Control_under_2MB_TADS_GpC':prop_of_meth_gpc_small,'Control_under_2MB_TADS_CpG':prop_of_meth_cpg_small})
data_1.to_csv('Control_data_for_methylation_in_regions_small_TADS.txt',sep='\t',header=True,index=False)



