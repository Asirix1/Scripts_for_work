def making_random_contacts(cool_path,depth_of_seq=10**7,path_random_cool):
	c = cooler.Cooler(cool_path)
	cov=cooltools.coverage(c, ignore_diags=0)
	df = c.pixels()[:]
	bins=c.bins()[:]
	df_only_main=df.copy(deep=True)
	df_only_main.loc[df_only_main['bin1_id'] != df_only_main['bin2_id'],['count']]=0
	cov_main_diag=np.bincount(df_only_main["bin1_id"], weights=df_only_main['count'], minlength=len(bins['chrom']))    
	cov=np.array((cov-cov_main_diag),dtype=int)



	coverage_per_bin=cov[1,:]

	all_coverage_bins=coverage_per_bin.sum()

	bin_ids=np.arange(0,len(coverage_per_bin))

	dict_bin_coverage=dict(zip(bin_ids,coverage_per_bin))
	coverage_pairs_df=df[['bin1_id','bin2_id']]
	coverage_pairs_df['bin1_id'] = coverage_pairs_df['bin1_id'].map(dict_bin_coverage)
	coverage_pairs_df['bin2_id'] = coverage_pairs_df['bin2_id'].map(dict_bin_coverage)


	random_contacts=2*coverage_pairs_df['bin1_id']*coverage_pairs_df['bin2_id']

	probability_vector=normalize(random_contacts.values.reshape(1, -1), norm='l1')

	probability_vector=probability_vector.flatten()

	vector_count=np.array(df['count'])

	random_vector=np.random.multinomial(depth_of_seq,pvals=probability_vector)
	df['count']=random_vector


	cooler.create_cooler(path_random_cool,bins,df)
    return df
def making_random_pairs(path_random_cool,random_pairs_path):
    c = cooler.Cooler('Random_contacts.cool')
    df = c.pixels(join=True)[:]
    pairs_df = df[['chrom1', 'start1', 'chrom2', 'start2']]
    count_vector = df['count']

    pairs = pairs_df.loc[pairs_df.index.repeat(count_vector)].reset_index(drop=True)
    print(pairs)
    pairs.to_csv('random.pairs', sep='\t', header=None, index=False)
    return pairs
making_random_contacts(