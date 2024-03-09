import os
import scipy as sc
import nucleotide_transformer
import h5py

import haiku as hk
import jax
import jax.numpy as jnp
from nucleotide_transformer.pretrained import get_pretrained_model

i=os.environ['i'] 
number_sequences=int(os.environ['number_sequences'])
layer=int(os.environ['layer'])

model_name = '2B5_multi_species'

parameters, forward_fn, tokenizer, config = get_pretrained_model(
    model_name=model_name,
    mixed_precision=False,
    embeddings_layers_to_save=(layer,),
    attention_maps_to_save=((1, 4), (7, 18)),
    max_positions=100   
)
forward_fn = hk.transform(forward_fn)

sequences=[]
line_in_file=[]
t=0
with open('/beegfs/data/hpcws/ws1/popov-GENA/GENA/genomic_reg_fa/'+i+'_new_coordinates.bed.fa') as f:
 for line in f:
    t+=1 
    if len(line)<=6000 and 'N' not in line:
        sequences+=[line.rstrip()]
        line_in_file+=t
    if len(sequences)>=number_sequences:
        print(len(sequences))
        break
f.close()

tokens_ids = [b[1] for b in tokenizer.batch_tokenize(sequences)]
tokens_str = [b[0] for b in tokenizer.batch_tokenize(sequences)]
tokens = jnp.asarray(tokens_ids, dtype=jnp.int32)

random_key = jax.random.PRNGKey(0)

outs = forward_fn.apply(parameters, random_key, tokens)

embeddings = outs["embeddings_"+ str(layer)][:, 1:, :]  # removing CLS token
print(embeddings.shape)

six_1000_bed_regions=pd.read_csv('/beegfs/data/hpcws/ws1/popov-GENA/GENA/genomic_reg_fa/'+i+'_new_coordinates.bed', sep='\t', header=None)
bed_regions=pd.read_csv('/beegfs/data/hpcws/ws1/popov-GENA/GENA/genomic_reg_fa/'+i+'.bed', sep='\t', header=None)

six_1000_bed_regions[3]=bed_regions[1]-six_1000_bed_regions[1]
six_1000_bed_regions[4]=six_1000_bed_regions[2]-bed_regions[2]
six_1000_bed_regions[3]=(six_1000_bed_regions[3]/6).astype(int) #number of left tokens for deleting
six_1000_bed_regions[4]=(six_1000_bed_regions[4]/6).astype(int) #number of right tokens for deleting
count=0
lengths=[]
line_in_file=np.arange(0,10,1)
for i in line_in_file:     #deleting informative tokens
    embeddings[count,:six_1000_bed_regions[3][i],:]=0
    embeddings[count,six_1000_bed_regions[4][i]:,:]=0
    lengths.append(bed_regions[2][i]-bed_regions[1][i])     #collecting lengths of sequences 
    count+=1    
lengths=np.array(lengths)
mean_embeddings= np.true_divide(np.sum(embeddings, axis=1),lengths[:, np.newaxis])

with h5py.File('mean_embedings_'+str(layer)+'_'+i+'_number_sequences_'+str(number_sequences)+'.hdf5', 'w') as f:
    dset = f.create_dataset('mean_embedings_'+str(layer), data=mean_embeddings)
f.close()