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



def find_multiple(start, end):
    start_multiple = start + (10000 - start % 10000) if start % 10000 != 0 else start
    end_multiple = end - (end % 10000) if end % 10000 != 0 else end

    return start_multiple, end_multiple



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
    genome = file_path[:-3]
    sequences_object = FastaFile(genome)
    file = h5py.File(f'/beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/for_GENA/{i}_seq_btw_pairs.hdf5', 'w')
    for num, row in tqdm.tqdm(ctct_pairs.iterrows()):
        observed = df[(df['chrom1'] == row['chrom1']) & (df['start1'] <= row['start1']) & (df['end1'] >= row['end1']) &
                      (df['chrom2'] == row['chrom2']) & (df['start2'] <= row['start2']) & (df['end2'] >= row['end2'])]

        if observed.empty:
            continue
        else:
            observed_btw = df[(df['chrom1'] == row['chrom1']) & (df['end1'] >= row['end1']) & 
                          (df['chrom2'] == row['chrom2']) & (df['start2'] <= row['start2'])].copy(deep=True).reset_index()
            if observed_btw.empty:
                obs_vs_expected=0
            else:
                bin_1_id = observed_btw['start1'].values / 10000
                bin_2_id = observed_btw['start2'].values / 10000
                observed_btw['dist'] = abs(bin_2_id - bin_1_id)
                chrom = observed_btw['chrom1'].values[0]
                expected_reg = expected[(expected['region1'] == chrom) & (expected['dist'] < int(observed_btw['dist'].max()))]
                average_count = observed_btw.groupby('dist')['count'].mean()
                result = pd.merge(expected_reg, average_count, left_on='dist', right_index=True)
                result = result[~result['dist'].isin([0, 1])]
                if result.empty:
                    obs_vs_expected=0
                else:
                    obs_vs_expected=np.nanmean(result['count'] / result['count.avg'])
                
            start_multiple, end_multiple = find_multiple(observed['end1'].values[0], observed['start2'].values)
            sequence_1=sequences_object.fetch(observed['chrom1'].values[0], start_multiple, end_multiple)
            
            file.create_dataset(f'seq_btw_pair_{num}', data=sequence_1)
            file.create_dataset(f'pair_{num}_value', data=obs_vs_expected)
    file.close()

