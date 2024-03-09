def bed_graph_converter(file, chunksize=5*10**5, tempdir="./tmp"):
    import os
    import pandas as pd
    import numpy as np
    # Create temporary directory if it does not exist
    os.makedirs(tempdir, exist_ok=True)
    
    with open(file) as f:
        num_rows = sum(1 for line in f)

    num_chunks = int(np.ceil(num_rows/chunksize))
    chunks = pd.read_csv(file, sep='\t', header=None, chunksize=chunksize)
    
    # Write each chunk to a temporary file
    for i, chunk in enumerate(chunks):
        coordinates = np.concatenate([np.arange(start, end) for start, end in chunk[[1, 2]].values])
        chr_score_values = chunk[[0, 3]].values
        chr_list = np.repeat(chr_score_values[:, 0], chunk[2]-chunk[1])
        score_list = np.repeat(chr_score_values[:, 1], chunk[2]-chunk[1])
        chunk_converted = pd.DataFrame({"chromosome": chr_list, "start": coordinates, "score": score_list})
        
        chunk_file = os.path.join(tempdir, f"{i}.csv")
        chunk_converted.to_csv(chunk_file, index=False)
    
    # Combine temporary files into final output
    results = []
    for i in range(num_chunks):
        chunk_file = os.path.join(tempdir, f"{i}.csv")
        chunk_converted = pd.read_csv(chunk_file)
        results.append(chunk_converted)
        os.remove(chunk_file)
    
    bed_graph_converted = pd.concat(results).sort_values(["chromosome", "start"]).reset_index(drop=True)
    
    return bed_graph_converted