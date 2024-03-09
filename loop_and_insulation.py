import pandas as pd
import cooler
import numpy as np
import seaborn as sns
from coolpuppy import coolpup
import coolpuppy.lib.numutils as num
from tqdm import tqdm
from coolpuppy import plotpup
import matplotlib.pyplot as plt


def making_random_pairs(path_random_cool,random_pairs_path):
    c = cooler.Cooler(path_random_cool)
    df = c.pixels(join=True)[:]
    pairs_df = df[['chrom1', 'start1', 'chrom2', 'start2']]
    count_vector = df['count']

    pairs = pairs_df.loc[pairs_df.index.repeat(count_vector)].reset_index(drop=True)

    

    pairs.sort_values(['chrom1', 'chrom2'],inplace=True)
    print(pairs)

    pairs.to_csv(random_pairs_path, sep='\t', header=None, index=False)  



def making_pairs(path_random_cool):
    c = cooler.Cooler(path_random_cool)
    df = c.pixels()[:]
    pairs_df = df[['bin1_id', 'bin2_id']]
    count_vector = df['count']

    pairs = pairs_df.loc[pairs_df.index.repeat(count_vector)].reset_index(drop=True)

    pairs.sort_values(['bin1_id', 'bin2_id'],inplace=True)
    
    return pairs
        

def making_bad_cool(normal_cool, random_cool, df_normal, df_random, new_cool_name, percent_of_noise):
    c_normal=cooler.Cooler(normal_cool)
    c_random=cooler.Cooler(random_cool)
    
    bins_normal=c_normal.bins()[:]
    bins_random=c_random.bins()[:]

    num_random_rows = int(len(df_normal) * percent_of_noise / 100)
    random_df = df_random.sample(n=num_random_rows)
    normal_df = df_normal.sample(n=(len(df_normal) - num_random_rows))
    
    df = pd.concat([normal_df, random_df]).reset_index()
    df.sort_values(['bin1_id', 'bin2_id'],inplace=True)
    
    df['count']=df.groupby(['bin1_id','bin2_id'],as_index=False).transform('count')
    df = df.drop('index', axis=1).reset_index(drop=True)  
    df = df.drop_duplicates()
    print(df)
    cooler.create_cooler(new_cool_name,bins_normal,df)

def calculate_cis_trans(cool_with_noise):
    c_noisy=cooler.Cooler(cool_with_noise)
    df=c_noisy.pixels(join=True)[:]
    filtered_df = df[df['chrom1'] == df['chrom2']]
    count_matches = filtered_df.shape[0]
    percentage_matches = (count_matches / (df.shape[0]-count_matches)) * 100
    print(f'Cis/trans = {round(percentage_matches, 2)}')
    
def calc_insulation_and_loop(cool_path,ctcf_bed_path, ctcf_pairs_bedpe):
    clr = cooler.Cooler(cool_path)
    ctcf=pd.read_csv(ctcf_bed_path,sep='\t')
    insulation = coolpup.pileup(clr, ctcf, features_format='bed',clr_weight_name=False,coverage_norm=False, ooe=False, local=True,flank=200000, min_diag=2,nproc=50, seed=1, nshifts=0)
    plotpup.plot(insulation,sym=False,score=True)
    plt.savefig(f'{cool_path}_insulation.png')
    map_for_ins=insulation['data'].values
    score_insulation=num.get_insulation_strength(map_for_ins[0],ignore_central=3)
    ctcf_pairs=pd.read_csv(ctcf_pairs_bedpe, sep='\t')
    loop = coolpup.pileup(clr, ctcf_pairs, features_format='bedpe', clr_weight_name=False, coverage_norm=False, ooe=False, flank=200000, min_diag=2,nproc=50, seed=1, nshifts=0)
    plotpup.plot(loop, cmap='coolwarm', sym=False,score=True)
    plt.savefig(f'{cool_path}_loop.png')
    map_for_loop=loop['data'].values
    score_loop=num.get_enrichment(map_for_loop[0], 3)

    return score_insulation,score_loop

ins_score_list=[]
loop_score_list=[]

df_normal=making_pairs('HCT116_deep_seq_10e8.cool')
df_random=making_pairs('Random_HCT116_10e8.cool')
df_normal=df_normal.sample(n=len(df_random)).reset_index(drop=True)
print(df_normal)
print(df_random)
for i in tqdm(np.arange(0,120,20)):
    making_bad_cool('HCT116_deep_seq_10e8.cool','Random_HCT116_10e8.cool', df_normal, df_random, f'Cool_with_noise_{i}.cool', i)
    #calculate_cis_trans(f'Cool_with_noise_{i}.cool')
    #cooler balance --force Cool_with_noise_{i}.cool
    ins,loop=calc_insulation_and_loop('Cool_with_noise_'+str(i)+'.cool','/beegfs/data/hpcws/ws1/popov-popov_gnomad/Random/Tracks/TADs_boundaries_HCT116_hg38.bed','/beegfs/data/hpcws/ws1/popov-popov_gnomad/Random/Tracks/Loop_boundaries_HCT116_hg38.bedpe')
    #!coolpup.py /beegfs/data/hpcws/ws1/popov-popov_gnomad/Random/Cool/Cool_with_noise.cool /beegfs/data/hpcws/ws1/popov-popov_gnomad/Random/Tracks/Loop_boundaries_HCT116_hg38.bedpe --flank 200000 --nproc 0 -o /beegfs/data/hpcws/ws1/popov-popov_gnomad/Random/Loop/Loop_with_noise_percent_{i}.clpy
    #!plotpup.py --dpi 600  --input_pups /beegfs/data/hpcws/ws1/popov-popov_gnomad/Random/Loop/Loop_with_noise_percent_{i}.clpy --output /beegfs/data/hpcws/ws1/popov-popov_gnomad/Random/Figures/Loop_with_noise_percent_{i}.png
    #!coolpup.py /beegfs/data/hpcws/ws1/popov-popov_gnomad/Random/Cool/Cool_with_noise.cool /beegfs/data/hpcws/ws1/popov-popov_gnomad/Random/Tracks/TADs_boundaries_HCT116_hg38.bed --local --flank 200000 --nproc 0 -o /beegfs/data/hpcws/ws1/popov-popov_gnomad/Random/Insulation/Insulation_with_noise_percent_{i}.clpy
    #!plotpup.py --dpi 600 --input_pups /beegfs/data/hpcws/ws1/popov-popov_gnomad/Random/Insulation/Insulation_with_noise_percent_{i}.clpy --output /beegfs/data/hpcws/ws1/popov-popov_gnomad/Random/Figures/Insulation_with_noise_percent_{i}.png
    ins_score_list.append(ins)
    loop_score_list.append(loop)
    #making_random_pairs(f'Cool_with_noise_{i}.cool',f'Pairs_with_noise_{i}.pairs')
    #!java -Xmx20G -jar /home/aapopov/tool/juicer1.6_compact-main/scripts/common/juicer_tools.jar pre  Pairs_with_noise_{i}.pairs Hic_with_noise_{i}.hic --threads 20 /beegfs/data/hpcws/ws1/popov-popov_gnomad/Random/reference.chromsizes.for_cooler


df=pd.DataFrame({'noise percent':np.arange(0,120,20),'loop score':loop_score_list,'insulation score':ins_score_list})
plt.figure(figsize=(10, 5))
sns.lineplot(x='noise percent', y='loop score', data=df)
plt.title('Figure loop score od noise percent')
plt.show()
plt.savefig('loop_noise.png')

plt.figure(figsize=(10, 5))
sns.lineplot(x='noise percent', y='insulation score', data=df)
plt.title('Figure insulation score of noise percent')
plt.show()
plt.savefig('insulation_noise.png')