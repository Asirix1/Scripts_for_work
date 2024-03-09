#!/bin/sh
#SBATCH --job-name=making_train_set     
#SBATCH --error=random_contacts-%j.err        
#SBATCH --output=random_contacts-%j.log        
#SBATCH --time=71:00:00 
#SBATCH --partition=operation
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G



source /home/popov/miniconda3/bin/activate
conda activate cooltools
#for i in Muntiacus_muntjak Colobus_angolensis Dugong_dugon Xerus_rutilus Leopardus_tigrinus Herpailurus_yagouaroundi Martes_foina Martes_martes Gymnobelideus_leadbeateri Myocastor_coypus Phascolarctos_cinereus Prionailurus_bengalensis Manis_javanica Martes_zibellina Otocolobus_manul Pusa_sibirica Urocyon_littoralis_catalinae Loxodonta_africana Potorous_gilbertii Meriones_libycus Peromyscus_crinitus Osphranter_robustus Rhincodon_typus Felis_nigripes Peromyscus_eremicus Bassariscus_sumichrasti Odocoileus_hemionus Phascogale_tapoatafa Halichoerus_grypus Arctocephalus_forsteri Myrmecobius_fasciatus Theropithecus_gelada Sylvilagus_bachmani Balaenoptera_borealis Uromys_caudimaculatus Manis_pentadactyla Martes_flavigula Papio_papio Dendrolagus_matschiei Leptonychotes_weddellii Miopithecus_talapoin Otolemur_crassicaudatus Eubalaena_australis
#do 
#hic2cool convert -r 10000 -p 20 /beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/hic/${i}.hic /beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/cool/${i}.cool
#i='Homo_sapiens'
#cooler balance -p 50 /beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/cool/${i}.cool
#done 
#cd /beegfs/data/hpcws/ws1/popov-popov_gnomad/GENA_loop_predict/ctcf_pairs
#shuf -n 50000 ./Homo_sapiens_CTCF_pairs.bedpe -o ./Homo_sapiens_CTCF_pairs_1.bedpe
#mv Homo_sapiens_CTCF_pairs_1.bedpe Homo_sapiens_CTCF_pairs.bedpe
#cd ..
python gena_loop_predict.py

#python tads_obs_exp.py