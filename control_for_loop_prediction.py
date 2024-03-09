import cooler
import cooltools
import pandas as pd
import numpy as np
import gzip
import shutil
import tqdm
import os
from pysam import FastaFile
import h5py

#species = 'Muntiacus_muntjak+Colobus_angolensis+Dugong_dugon+Xerus_rutilus+Leopardus_tigrinus+Herpailurus_yagouaroundi+Martes_foina+Martes_martes+Gymnobelideus_leadbeateri+Myocastor_coypus+Phascolarctos_cinereus+Prionailurus_bengalensis+Manis_javanica+Martes_zibellina+Otocolobus_manul+Pusa_sibirica+Urocyon_littoralis_catalinae+Loxodonta_africana+Potorous_gilbertii+Meriones_libycus+Peromyscus_crinitus+Osphranter_robustus+Rhincodon_typus+Felis_nigripes+Peromyscus_eremicus+Bassariscus_sumichrasti+Odocoileus_hemionus+Phascogale_tapoatafa+Halichoerus_grypus+Arctocephalus_forsteri+Myrmecobius_fasciatus+Theropithecus_gelada+Sylvilagus_bachmani+Balaenoptera_borealis+Uromys_caudimaculatus+Manis_pentadactyla+Martes_flavigula+Papio_papio+Dendrolagus_matschiei+Leptonychotes_weddellii+Miopithecus_talapoin+Otolemur_crassicaudatus+Eubalaena_australis'
#species = species.split('+')
species=['Homo_sapiens']
for i in species:
    ctct_pairs = pd.read_csv(f'/beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/ctcf_pairs/{i}_CTCF_pairs.bedpe', sep='\t', header=None, names=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'])
    cool_path = f'/beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/cool/{i}.cool'
    c = cooler.Cooler(cool_path)
    df = c.pixels(join=True)[:]
    expected = cooltools.api.expected.expected_cis(c)
    file_path = f'/beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/genomes/{i}.fasta.gz'
    '''
    with gzip.open(file_path, 'rb') as f_in:
        with open(file_path.replace('.gz', ''), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    '''
    genome = file_path[:-3]
    sequences_object = FastaFile(genome)
    file = h5py.File(f'/beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/for_GENA/{i}_all_pairs_random_control.hdf5', 'w')
    chrom_size_path = f'/beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/genomes/{i}.fasta.fai'
    chrom_size = pd.read_csv(chrom_size_path, usecols=[0, 1],names=['chr','value'],sep='\t',header=None)
    random_control=pd.DataFrame(columns=['chrom1', 'start1', 'end1', 'chrom2','start2', 'end2'])
    for chrom in ctct_pairs['chrom1'].unique():
        subset_df = ctct_pairs[ctct_pairs['chrom1'] == chrom]
        value = chrom_size.loc[chrom_size['chr'] == chrom, 'value'].values - max(subset_df['end2']-subset_df['start1'])
        rand_start1 = np.random.randint(0,value,size=len(subset_df))
        rand_reg1=subset_df['end1'].values - subset_df['start1'].values
        rand_inter_reg=subset_df['start2'].values - subset_df['end1'].values
        rand_reg2=subset_df['end2'].values - subset_df['start2'].values

        np.random.shuffle(rand_reg1)
        np.random.shuffle(rand_inter_reg)
        np.random.shuffle(rand_reg2)

        rand_end_1 = rand_start1 + rand_reg1
        rand_start2 = rand_end_1 + rand_inter_reg
        rand_end2 = rand_start2 + rand_reg2
        chrom_list = [chrom]*len(subset_df)
        random_for_chrom=pd.DataFrame({'chrom1':chrom_list,'start1':rand_start1,'end1':rand_end_1,'chrom2':chrom_list,'start2':rand_start2,'end2':rand_end2})
        random_control = pd.concat([random_control, random_for_chrom], axis=0, ignore_index=True)
    random_control.to_csv(f'/beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/ctcf_pairs/{i}_CTCF_all_pairs_random_control.bedpe', sep='\t',index=False,header=None)
    for num, row in tqdm.tqdm(random_control.iterrows()):
        observed = df[(df['chrom1'] == row['chrom1']) & (df['start1'] <= row['start1']) & (df['end1'] >= row['end1']) &
                      (df['chrom2'] == row['chrom2']) & (df['start2'] <= row['start2']) & (df['end2'] >= row['end2'])]
        if observed.empty:
            continue
        else:
            bin_1_id = observed['start1'].values[0] / 10000
            bin_2_id = observed['start2'].values[0] / 10000
            dist = abs(bin_2_id - bin_1_id)
            chrom = observed['chrom1'].values[0]
            expected_reg = expected[(expected['region1'] == chrom) & (expected['dist'] == int(dist))]
            obs_vs_exp=observed['count'].sum()/expected_reg['count.avg'].values
            sequence_1=sequences_object.fetch(observed['chrom1'].values[0], observed['start1'].values[0], observed['end1'].values[0])
            sequence_2=sequences_object.fetch(observed['chrom2'].values[0], observed['start2'].values[0], observed['end2'].values[0])
            file.create_dataset(f'pair_{num}_sequence1', data=sequence_1)
            file.create_dataset(f'pair_{num}_sequence2', data=sequence_2)
            file.create_dataset(f'pair_{num}_value', data=obs_vs_exp[0])
    file.close()