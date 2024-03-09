import pandas as pd
import numpy as np
file='/storage2/asirix/branch/Danio/test.bedgraph'
bed_graph = pd.read_csv(file, sep='\t', header=None)    
coordinates = np.concatenate([np.arange(start, end) for start, end in bed_graph[[1, 2]].values])
chr_list = np.concatenate([np.full(int(start), chromosome) for chromosome, start, end in bed_graph[[0, 1, 2]].values])
score_list = np.concatenate([np.full(int(start), score) for score, start, end in bed_graph[[3, 1, 2]].values])  
bed_graph_converted = pd.DataFrame({"chromosome": chr_list, "start": coordinates, "end": coordinates + 1, "score": score_list})
bed_graph_converted.to_csv(file + "_converted.bedgraph", sep='\t', header=False, index=False)