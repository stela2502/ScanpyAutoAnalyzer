

def addSampleDataRust( adata, sname, sampleInfo ):
    """Add the sample information obtained by either running demux10x or quantifyRhapsody to an anndata object"""
    adata.obs['sname'] = sname
    adata.obs_names = [ sname +"."+ n for n in adata.obs_names]
    sampleInfo = pd.read_csv( sampleInfo, sep="\t",index_col= 0)
    sampleInfo['CellIDs'] = [ sname +"."+n+"-1" for n in sampleInfo.index.values]
    sampleInfo.set_index('CellIDs', inplace=True)
    adata.obs['CellIDs'] = adata.obs_names
    obs = adata.obs.merge( sampleInfo, left_index=True, right_index=True)
    print(obs.shape)
    obs.drop_duplicates()
    obs.drop( ['CellIDs'], axis=1, inplace=True)
    print(obs.shape)
    adata = adata[ obs.index.values]
    adata.obs = obs
    ADT = pd.read_csv( 'Cell2Sample.ADT'+sname+'.R1_001.fastq.gz.tsv', sep="\t",index_col= 0)
    ADT['CellIDs'] = [ sname+"."+n+"-1" for n in ADT.index.values]
    ADT.set_index( 'CellIDs' , inplace=True)
    
    obs = adata.obs.merge( ADT, left_index=True, right_index=True )
    #print(obs)
    #print (adata)
    adata = adata[ obs.index.values]
    #print(adata)
    adata.obs = obs
    return (adata)

def testRP(x, RP):
    r= True
    if RP.match(x):
        r = False
    return (r)

def dropRP ( adata, drop = True ):
    """ remove all genes matching to /^R[pP][SsLl]/ """
    RP = re.compile('^R[pP][SsLl]')
    RPgenes = [ x  for x in adata.var.index if not testRP(x, RP)]
    print ( Counter(RPgenes) )
    if  len(RPgenes) > 0:
        adata.obs['RP[LS]sum'] = adata[:,RPgenes].X.sum(axis=1)
    else:
        adata.obs['RP[LS]sum'] = 0
    if drop and len(RPgenes) > 0:
        OK = [ testRP(x, RP) for x in adata.var.index ]
        print(Counter(OK))
        adata._inplace_subset_var( np.array(OK) )
    adata

def testMT(x, MT):
    r= True
    if MT.match(x):
        r = False
    return (r)

def dropMT (adata, drop = True ) :
    """ remove all genes matching to /^[Mm][Tt]-/ """
    MT = re.compile('^[Mm][Tt]-')
    MTgenes = [ x  for x in adata.var.index if not testMT(x, MT)]
    print ( Counter(MTgenes) )
    if len(MTgenes) > 0:
        adata.obs['MTsum'] = adata[:,MTgenes].X.sum(axis=1)
    else:
        adata.obs['MTsum'] = 0
    if drop and len(MTgenes) > 0:
        OK = [ testMT(x, MT) for x in adata.var.index ]
        adata._inplace_subset_var( np.array(OK) )


def testGene(x, MT, RP):
    r= True
    if MT.match(x) or RP.match(x):
        r = False
    return (r)

def write_top_genes_per_cluster(adata, stats_name, n_top_genes=20, output_file="top_genes_per_cluster.csv"):
    """
    Extracts the top `n_top_genes` highly expressed genes per cluster from a Scanpy `AnnData` object
    and writes them into a CSV file. Each row in the file contains the top genes for a cluster,
    comma-separated.

    Parameters:
        adata (AnnData): The annotated data matrix from Scanpy containing gene expression data.
        n_top_genes (int): The number of top genes to extract per cluster. Default is 20.
        output_file (str): The output CSV file where top genes per cluster will be saved.
    """


    # Check if the `rank_genes_groups` function has been run
    if stats_name not in adata.uns:
        raise ValueError(f"No {stats_name} results found in the AnnData object. Please run `sc.tl.rank_genes_groups()` first.")

    # Get the ranked genes for each group
    ranked_genes = pd.DataFrame(adata.uns[ stats_name ]['names']).head(n_top_genes)

    with open(output_file, 'w') as f:
        # Write top genes for each cluster
        for cluster in ranked_genes.columns:
            top_genes = ranked_genes[cluster].tolist()
            f.write(f"{cluster}: {','.join(top_genes)}\n")

    print(f"Top {n_top_genes} genes per cluster have been written to {output_file}")

# Example usage:
# Assuming you have an AnnData object `adata` and have run `sc.tl.rank_genes_groups` already:
# write_top_genes_per_cluster(adata, n_top_genes=20, output_file="top_genes.csv")



def write_stats_tables( adata, key_added,cn ="louvain" ):
    """
    Creates the stats using the 'wilcoxon' method, plots the top genes using scanpy.pl.rank_genes_groups
    and finally writes tab separated tables to a <key_added> foilder.
    """ 

    scanpy.tl.rank_genes_groups(
        adata, 
        groupby   = cn,
        key_added = key_added,
        method    = 'wilcoxon',
    )

    scanpy.pl.rank_genes_groups(adata, key = key_added, save= key_added+"_overview")

    diff_results = adata.uns[key_added]
    columns = ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']
    try: 
        os.mkdir(f"{key_added}") 
    except OSError as error: 
        print(error)  

    for i in adata.obs[cn].unique():
        table = {}
        for j in columns:
            #print(f"Analyszing column {diff_results[j]} {i}")
            table[j] = pd.DataFrame(diff_results[j])[str(i)]
        table = pd.DataFrame(table).to_csv(f"{key_added}/{cn}_cluster_{i}.csv")



def remove_gzipped_files(path):
    # Ensure the directory exists
    if os.path.exists(path):
        for file_name in os.listdir(path):
            # Check if the file is a gzipped file
            if file_name.endswith(".gz"):
                gzipped_file_path = os.path.join(path, file_name)
                os.remove(gzipped_file_path)
                print(f"Removed gzipped file: {gzipped_file_path}")
    else:
        print(f"The directory {path} does not exist.")

def gzip_files_in_directory(path):
    # Ensure the directory exists
    if os.path.exists(path):
        # Run the gzip command on all files in the directory
        subprocess.run(f"gzip {path}/*", shell=True, check=True)
        print(f"All files in {path} have been gzipped.")
    else:
        print(f"The directory {path} does not exist.")

def write_data_to_directory(adata, output_dir):
    """ 
    exports the anndata.X sparse matrix as 10X formated matrix.mtx barcodes.tsv and features.tsv set 
    """
    # Ensure the directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Remove gzipped files before writing new ones
    remove_gzipped_files(output_dir)

    # Write the matrix data to a .mtx file
    from scipy.io import mmwrite
    mmwrite(os.path.join(output_dir, 'matrix.mtx'), adata.X.T )

    # Write the barcodes (cell identifiers) to a .tsv file
    with open(os.path.join(output_dir, 'barcodes.tsv'), 'w') as f:
        for barcode in adata.obs_names.values:
            f.write(f"{barcode}\n")

    # Write the features (gene identifiers) to a .tsv file
    # Create the dataframe for the features.tsv.gz file
    features_df = pd.DataFrame({
        'gene_id': adata.var.index,  # Gene identifiers (e.g., Ensembl IDs)
        'gene_name': adata.var['gene_name'] if 'gene_name' in adata.var.columns else adata.var.index,  # Gene names (or use gene_id if unavailable)
        'gene_type': adata.var['feature_types'] if 'feature_types' in adata.var.columns else  ['Gene Expression'] * len(adata.var)  # Gene type, usually 'Gene Expression'
    })
    features_df.to_csv(os.path.join(output_dir, 'features.tsv'), sep='\t', index=False, header=False )

    # Optionally, gzip the files in the directory
    gzip_files_in_directory(output_dir)


def impute_expression_data( adata, output_dir ):
    """
    Uses scvi.model.SCVI to impute the expression data and adds that (model.posterior_predictive_sample) as the adata.X sparse matrix.
    In addition the imputed data is saved as 10X matrix.mtx, features.tsv and barcodes.tsv file triplet.
    """

    # Step 1: Get the number of CPU cores
    num_cores = os.cpu_count()
    
    # Step 2: Set num_workers to the number of cores + 1
    num_workers = num_cores + 1
    
    # Convert to CSR format
    adata.X = scipy.sparse.csr_matrix(adata.X)
    # Set up the model
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata)
    # Train the model
    model.train()
    # Get the imputed expression values for cells
    adata.X = model.posterior_predictive_sample(adata, n_samples=1).to_scipy_sparse()
    write_data_to_directory( adata, output_dir )


def grouped_obs_mean(adata, group_key, layer=None, gene_symbols=None):
    """
    calculate the mean expression for all genes in all groups
    """
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = adata.var[idx]
    else:
        new_idx = adata.var_names

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))
    return out


def mergeClosest( adata, group = 'louvain11_merged.9', cut = 0.9 ):
    """ 
    create mean expression per group, correlate all mean groups and scan all groups for corr > cut.
    The first group with cor > cut will merge with the best correlating groups available.
    Only one group will be merged at a time.
    """
    print (str(len(Counter(adata.obs[group])))+ " -", end=" ")
    df = grouped_obs_mean( adata, group )
    df = df.corr()
    for i in range(0,df.shape[0]):
        df.at[df.columns[i],df.columns[i]] = 0
    #print("processing groups")
    #print (df)
    for i in range( 0, len(df.columns)):
        i = df.columns[i]
        col = df[i]
        #print ( "max value="+ str( max( col )) )
        if max( col ) >= cut:
            gr = np.array(adata.obs[group])
            m = [df.columns[a] for a in range(0, len(col)) if col[a] == max(col) ]
            m = m[0]
            #print ("max value="+str(max[col])+ " merge " + str(m) + " to "+ str(i)+".")
            for x in range(0,len(adata.obs[group])):
                    if str(gr[x]) == i:
                        #print ( "changed" + str(x)+ " change " + gr[x]  + " to "+ str(m))
                        gr[x] = m
            #type(gr)
            #print("finished")
            adata.obs[group] = gr
            return True
    print ( "no cluster did match to any other cluster with cor > "+ str(cut) )
    return False


def reIDgroup( adata, group ):
    """
    Re-set the cluster ids to be in the tange of 0-max cluster id
    """
    gr = np.array(adata.obs[group])
    n = 0
    for i in set(adata.obs[group]):
        for x in range(0,len(adata.obs[group])):
            if adata.obs[group][x] == i:
                gr[x] = n
        n = n +1
    adata.obs[group+'_renamed'] = gr
