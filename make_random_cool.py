import cooler
import cooltools
import pandas as pd
import numpy as np
from scipy.sparse import dok_matrix
from tqdm import tqdm
from sklearn. preprocessing import normalize
from math import ceil


def making_random_contacts(cool_path,path_random_cool,depth_of_seq=10**7,reduction_bin=100):
    c = cooler.Cooler(cool_path)
    df = c.pixels()[:]
    print(df)
    cov=cooltools.coverage(c, ignore_diags=0)
    bins=c.bins()[:]
    df_only_main=df.copy(deep=True)
    df_only_main.loc[df_only_main['bin1_id'] != df_only_main['bin2_id'],['count']]=0
    cov_main_diag=np.bincount(df_only_main["bin1_id"], weights=df_only_main['count'], minlength=len(bins['chrom']))    
    cov=np.array((cov-cov_main_diag),dtype=int)
    coverage_per_bin=cov[1,:]
    reduced_coverage_per_bin = np.add.reduceat(coverage_per_bin, np.arange(0, len(coverage_per_bin), reduction_bin))
    np.random.seed(1)
    reduced_matrix = np.outer(reduced_coverage_per_bin, reduced_coverage_per_bin) * 2
    reduced_matrix = np.triu(reduced_matrix, k=0)
    reduced_matrix_flat=reduced_matrix.reshape(1, -1)
    probabil=normalize(reduced_matrix_flat, norm="l1")
    final_matrix = np.random.multinomial(depth_of_seq, pvals=probabil.ravel()).reshape(reduced_matrix.shape) 
    print(final_matrix)
    adding_for_not_main = ceil(np.trace(final_matrix)/2)
    print(adding_for_not_main)
    print(reduced_matrix)
    reduced_matrix = np.triu(reduced_matrix, k=1)
    print(reduced_matrix)
    probabil=normalize(reduced_matrix_flat, norm="l1")
    adding_for_not_main_matrix=np.random.multinomial(adding_for_not_main, pvals=probabil.ravel()).reshape(reduced_matrix.shape) 
    final_matrix=final_matrix+adding_for_not_main_matrix
    print(final_matrix)
    print(final_matrix.sum())
    
    
    sparse_matrix = dok_matrix(final_matrix)
    
    bin1_id_list=[]
    bin2_id_list=[]
    count_list=[]
    for i in tqdm(sparse_matrix.keys()):
        rows=coverage_per_bin[i[0]*reduction_bin:i[0]*reduction_bin+reduction_bin]
        cols=coverage_per_bin[i[1]*reduction_bin:i[1]*reduction_bin+reduction_bin]
        matrix_matrix=np.outer(rows,cols)*2
        number_to_random_add=sparse_matrix[i]
        if i[0]==i[1]:
            print(i)
            print(sparse_matrix[i])
            number_to_random_add=ceil(sparse_matrix[i]/2)
            print(number_to_random_add)
            matrix_matrix=np.triu(matrix_matrix,k=0)
        if matrix_matrix.sum() == 0:
            continue
        matrix_matrix_flat=matrix_matrix.reshape(1, -1)
        probabil_bin=normalize(matrix_matrix_flat, norm="l1")
        random_matrix=np.random.multinomial(number_to_random_add, pvals=probabil_bin.ravel()).reshape(matrix_matrix.shape)
        sparse_random_matrix = dok_matrix(random_matrix)
        
        for g in sparse_random_matrix.keys():
            new_df=pd.DataFrame(columns=['bin1_id', 'bin2_id', 'count'])
            coordinate1=g[0]+i[0]*reduction_bin
            coordinate2=g[1]+i[1]*reduction_bin
            bin1_id_list.append(coordinate1)
            bin2_id_list.append(coordinate2)
            count_list.append(sparse_random_matrix[g])
        lines=len(bin1_id_list)
    final_random_df=pd.DataFrame({'bin1_id': bin1_id_list,'bin2_id':bin2_id_list,'count':count_list})
    final_random_df = final_random_df.astype(int)
    cooler.create_cooler(path_random_cool,bins,final_random_df)



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
    df = c.pixels(join=True)[:]
    pairs_df = df[['chrom1', 'start1', 'chrom2', 'start2']]
    count_vector = df['count']

    pairs = pairs_df.loc[pairs_df.index.repeat(count_vector)].reset_index(drop=True)

    pairs.sort_values(['chrom1', 'chrom2'],inplace=True)
    
    return pairs

    
#making_random_contacts('HCT116_deep_seq.cool','Random_HCT116_10e8.cool',depth_of_seq=10**8,reduction_bin=100)
#making_random_pairs('Random_HCT116_10e8.cool','Random_with_no_dots.pairs')

df_pairs=making_pairs('HCT116_deep_seq.cool')
df_pairs_subsamp=df_pairs.sample(n=100000740).reset_index(drop=True)
print(df_pairs_subsamp)
df_pairs_subsamp.to_csv('HCT116_deep_seq_10e8.pairs', header=None, sep='\t', index=False)
