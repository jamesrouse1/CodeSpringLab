##jpplot version 0.8
## February 4, 2021

## This is a repository of plotting and preprocessing functions to be used in conjunction with Scanpy
import pandas as pd
import numpy as np
import scanpy as sc

def get_markers(adata, n_markers=20):
	##updated for gprofiler version 1.0.0
    import pandas as pd
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    marker_dataframe=pd.DataFrame(
    {group + '_' + key[:7]: result[key][group]
    for group in groups for key in ['names', 'logfoldchanges','pvals']})
    return marker_dataframe[:n_markers]
		
def kneeplot(adata, vline=None, min_genes=1, min_cells=1):
    import matplotlib
    import matplotlib.pyplot as pl
    import scanpy as sc
    import numpy as np
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    if 'n_counts' not in adata.obs.columns:
        print('Saving n_counts to .obs dataframe')
        adata.obs['n_counts'] = np.sum(adata.X, axis=1).A1
    pl.plot(sorted(np.array(adata.obs['n_genes']), reverse=True), linewidth=3, color='Blue', label='Unique Genes')
    pl.plot(sorted(np.array(adata.obs['n_counts']), reverse=True), linewidth=3, color='Red', label='UMIs')
    pl.yscale('log')
    pl.xscale('log')
    pl.legend()
    pl.xlabel('Barcodes')
    pl.ylabel('Genes Detected')
    if vline:
        pl.vlines(ymin=1, ymax=10000, x=vline)
        
def soup_kneeplot_from_raw(adata, vline=None, min_genes=1, min_cells=1, return_soup=False, Sample_Name=None, obs_key='n_genes', plot_cells=10000):
    from sklearn.mixture import GaussianMixture as GMM
    import matplotlib.pyplot as pl
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import jpplot
    
    adata.var_names_make_unique()
    
    if not Sample_Name:
        if 'Cellranger_Library_ID' in adata.uns.keys():
            Sample_Name = adata.uns['Cellranger_Library_ID']
            print('Detected sample name as:',Sample_Name[0])
    else:
        Sample_Name = 'Soup_Counts'
        
    if 'n_counts' not in adata.obs.columns:
        print('Saving n_counts to .obs dataframe')
        adata.obs['n_counts'] = np.sum(adata.X, axis=1).A1
    if 'n_genes' not in adata.obs.columns:
        print('Saving n_genes to .obs dataframe')
        adata.obs['n_genes'] = np.sum(adata.X > 0, axis=1)
        
    assert obs_key in adata.obs.select_dtypes(include=np.number), 'Warning, obs_key must be a numerical data key in adata.obs'
            
    #sc.pp.filter_genes(adata, min_cells=min_cells)
    #sc.pp.filter_cells(adata, min_genes=min_genes)

        
    X = adata.obs[obs_key]
    X = np.array(X).reshape(-1, 1)
    gmm = GMM(n_components=3).fit(X)
    labels = gmm.predict(X)
    adata.obs['GMM'] = labels
    
    df = adata.obs.loc[:,['n_genes','GMM']]
    df = df.sort_values(by='n_genes', ascending=False)
    df['order'] = np.log10(np.arange(len(df)))
    compartments = ['Empty','Soup','Cells']
    groups = df.groupby('GMM').agg({'n_genes':'median'}).sort_values(by='n_genes').index.tolist()
    groups = [str(i) for i in groups]
    mydict = dict(zip(groups,compartments))
    df['Compartment'] = df['GMM'].copy().astype('str')
    df = df.replace({'Compartment':mydict})
    compartment_counts = df['Compartment'].value_counts()
    

    adata.obs['Compartment'] = adata.obs['GMM'].copy().astype('str')
    adata.obs = adata.obs.replace({'Compartment':mydict})
    
    #subsample to speed up plotting:
    #df = df.sample(plot_cells)
    
	# Plot
    df_1 = df[df['Compartment'] == 'Empty'].sample(plot_cells)
    df_2 = df[df['Compartment'] == 'Soup']
    df_3 = df[df['Compartment'] == 'Cells']
    
    pl.scatter(x=df_1['order'], y=df_1['n_genes'],c="red", label='Empty')
    pl.scatter(x=df_2['order'], y=df_2['n_genes'],c="grey", label='Soup')
    pl.scatter(x=df_3['order'], y=df_3['n_genes'],c="black", label='Cells')
    
    
    pl.yscale('log')
    pl.legend()
    pl.grid(None)
    pl.xlabel('Barcodes')
    pl.ylabel('Genes Detected')
    if vline:
        pl.vlines(ymin=1, ymax=10000, x=vline)
    print(compartment_counts)
    pl.show()
    
    # Return a dataframe with soup counts
    if return_soup:
        soup = adata[adata.obs['Compartment'] == 'Soup']
        soup_df = pd.DataFrame(jpplot.matrix_to_df(soup).sum(0), columns = Sample_Name)
        print('Returning dataframe of soup counts...')
        return soup_df
	
def topdiff(adata, groupby='Cluster', n_genes=10, pval_adj_thresh=0.01):
    import scanpy as sc
    import pandas as pd
    topnames = []
    
    # Re-do rank_genes_groups if the correct data isn't current stored in .uns
    current_diffex_key = adata.uns['rank_genes_groups']['params']['groupby']
    if current_diffex_key != groupby:
        sc.tl.rank_genes_groups(adata, groupby=groupby)
        
    for cluster in adata.obs[groupby].unique():
        df = sc.get.rank_genes_groups_df(adata, group=cluster)
        df = df[df['pvals_adj'] < pval_adj_thresh]
        df = df.nlargest(columns='logfoldchanges', n=n_genes)
        topnames += df['names'].tolist()
        
    return(topnames)
    
def umap_topdiff(adata, n_genes=4, color_map='Reds', pvals_adj_thresh=1e-2, vmin=None, vmax=None, layer=None, groups=None):
    import scanpy as sc
    import pandas as pd
	
    if layer:
        assert layer in adata.layers.keys(), 'layer must be in adata.layers'
		
    topnames = []
    groupby=adata.uns['rank_genes_groups']['params']['groupby']
    
    if groups:
        if not isinstance(groups, list):
            groups = [groups]
            groups = [name for name in groups if groups in adata.obs[groupby].unique()]
    else:
        groups = adata.obs[groupby].unique()
        
    for g in groups:
        g_top = sc.get.rank_genes_groups_df(adata, group=g)
        g_top = g_top[g_top['pvals_adj'] < pvals_adj_thresh]
        g_top = g_top.sort_values(by='logfoldchanges', ascending=False)['names'][:n_genes].tolist()
        topnames += g_top

    sc.pl.umap(adata, color=topnames, color_map=color_map, vmin=vmin, vmax=vmax, layer=layer)
    
def umap_topscore(adata, n_genes=4, color_map='Reds', vmin=None, vmax=None):
    import scanpy as sc
    import pandas as pd
	
    topnames = []
    groupby=adata.uns['rank_genes_groups']['params']['groupby']
    for g in adata.obs[groupby].unique():
        g_top = sc.get.rank_genes_groups_df(adata, group=g)['names'][:n_genes].tolist()
        topnames += g_top

    sc.pl.umap(adata, color=topnames, color_map=color_map, vmin=vmin, vmax=vmax)
    
def matrix_toppval(adata, groupby=None, n_genes=4, standard_scale=None, color_map="Reds"):
    import scanpy as sc
    df = get_markers(adata, len(adata.obs[groupby].unique())*100)
    cluster_name_headers = [name for name in df.columns if name.endswith('names')]
    cluster_logfold_headers = [name for name in df.columns if name.endswith('logfold')]
    cluster_logfold_headers = [name for name in df.columns if name.endswith('pvals')]
    cluster_column_headers = [name for name in df.columns if name.endswith('names')]
    cluster_names = [word.replace('_names','') for word in cluster_name_headers]


    topnames = []
    for i in cluster_names:
        df1 = df.loc[:,[i+'_names',i+'_logfold',i+'_pvals']]
        topnames = topnames + df1.sort_values(by=i+'_pvals',ascending=True)[:n_genes][i+'_names'].tolist()

    sc.settings.set_figure_params(color_map=color_map, dpi=80)
    sc.pl.matrixplot(adata, groupby=groupby, var_names=topnames, standard_scale=standard_scale)
    sc.settings.set_figure_params(color_map="Reds", dpi=80)
    
def matrix_topscore(adata, groupby, n_genes=4, kind='matrixplot', standard_scale=None, color_map="Reds"):
    import scanpy as sc
    df = get_markers(adata, len(adata.obs[groupby].unique())*100)
    cluster_name_headers = [name for name in df.columns if name.endswith('names')]
    cluster_logfold_headers = [name for name in df.columns if name.endswith('logfold')]
    cluster_logfold_headers = [name for name in df.columns if name.endswith('pvals')]
    cluster_column_headers = [name for name in df.columns if name.endswith('names')]
    cluster_names = [word.replace('_names','') for word in cluster_name_headers]


    topnames = []
    for i in cluster_names:
        df1 = df.loc[:,[i+'_names',i+'_logfold',i+'_pvals']]
        topnames = topnames + df1.sort_values(by=i+'_pvals',ascending=True)[:n_genes][i+'_names'].tolist()

    sc.settings.set_figure_params(color_map=color_map, dpi=80)
    if kind == 'dotplot':
        sc.pl.dotplot(adata, groupby=groupby, var_names=topnames, standard_scale=standard_scale, cmap=color_map)
    elif kind == 'matrixplot':
        sc.pl.matrixplot(adata, groupby=groupby, var_names=topnames, standard_scale=standard_scale, cmap=color_map)
    else:
        print('Warning, you did not choose a valid plot type.  Must be matrixplot or dotplot.')
    sc.settings.set_figure_params(color_map="Reds", dpi=80)
    
def matrix_topdiff(adata, groupby=None, n_genes=4, kind='matrixplot', standard_scale='var', color_map="Reds", pval_thresh=1e-6, min_fold_change=None):
    import scanpy as sc
    import numpy as np
    
    assert kind in ['matrixplot','dotplot'], 'kind must be one of \[\'matrixplot\',\'dotplot\'\]'
    
    print('Showing genes with adj_pval <',pval_thresh)
    if min_fold_change:
        print('Showing genes with fold change >',min_fold_change)
    
    valid_groupby_keys = adata.obs.select_dtypes('category').columns
    existing_groupby = adata.uns['rank_genes_groups']['params']['groupby']
    
    if not groupby and existing_groupby:
    	groupby = existing_groupby
    	print('Grouping by previously used key for differential expression:', existing_groupby)
    	
    
    if not groupby == existing_groupby:
    	if groupby in valid_groupby_keys:
    		sc.tl.rank_genes_groups(adata, groupby=groupby)
    	else:
    		print('Warning: ' +str(groupby) +' not detected as a categorical type in obs keys' )
    		
    df = get_markers(adata, len(adata.obs[groupby].unique())*100)
    cluster_name_headers = [name for name in df.columns if name.endswith('names')]
    cluster_logfold_headers = [name for name in df.columns if name.endswith('logfold')]
    cluster_logfold_headers = [name for name in df.columns if name.endswith('pvals')]
    cluster_column_headers = [name for name in df.columns if name.endswith('names')]
    cluster_names = [word.replace('_names','') for word in cluster_name_headers]


    topnames = []
    
    for i in cluster_names:
        df1 = df.loc[:,[i+'_names',i+'_logfold',i+'_pvals']]
        if pval_thresh:
            assert type(pval_thresh) == float
            df1 = df1[df1[i+'_pvals'] < pval_thresh]
            
        if min_fold_change:
            assert isinstance(min_fold_change, (int, float, complex))
            df1 = df1[df1[i+'_logfold'] > min_fold_change]
    		
        topnames += df1.sort_values(by=i+'_logfold',ascending=False)[:n_genes][i+'_names'].tolist()
        
    sc.settings.set_figure_params(color_map=color_map, dpi=80)
    if kind == 'dotplot':
        sc.pl.dotplot(adata, groupby=groupby, var_names=topnames, standard_scale=standard_scale, cmap=color_map)
    elif kind == 'matrixplot':
        sc.pl.matrixplot(adata, groupby=groupby, var_names=topnames, standard_scale=standard_scale, cmap=color_map)
    else:
        print('Warning, you did not choose a valid plot type.  Must be matrixplot or dotplot.')
    sc.settings.set_figure_params(color_map="Reds", dpi=80)
		    
def toppval(adata, n_genes=10):
    import scanpy as sc
    df = get_markers(adata, len(list(adata.uns['rank_genes_groups']['names'].dtype.names))*100)
    cluster_name_headers = [name for name in df.columns if name.endswith('names')]
    cluster_logfold_headers = [name for name in df.columns if name.endswith('logfold')]
    cluster_logfold_headers = [name for name in df.columns if name.endswith('pvals')]
    cluster_column_headers = [name for name in df.columns if name.endswith('names')]
    cluster_names = [word.replace('_names','') for word in cluster_name_headers]


    topnames = []
    for i in cluster_names:
        df1 = df.loc[:,[i+'_names',i+'_logfold',i+'_pvals']]
        topnames = topnames + df1.sort_values(by=i+'_pvals',ascending=True)[:n_genes][i+'_names'].tolist()
    return(topnames)
    
#def corrplot(adata, n_genes, pixels=2000):
#    import scanpy as sc
#    import jpplot
#    import pandas as pd
#    from yellowbrick.features import Rank2D

#    topn = jpplot.topdiff(adata, n_genes=20)
#    df = pd.DataFrame(adata[:,topn].X.todense())
#    df.columns = topn

#    visualizer = Rank2D(algorithm="pearson", size=(pixels,pixels))
#   visualizer.fit_transform(df)
#   visualizer.poof()


def topscore(adata, n_genes=10):
    import scanpy as sc
    import pandas as pd
    df = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    cluster_names = df.columns
    
    topnames = []
    for i in cluster_names:
    	topnames = topnames + df[:n_genes][i].tolist()
    return(topnames)


### gProfiler for Gene Ontology scoring
def gprofiler_mouse(genelist, organism = 'mmusculus', cmap="rainbow"):
	from gprofiler import GProfiler
	import pandas as pd
	import numpy as np
	import matplotlib.pyplot as pl
	from matplotlib import rcParams
	from matplotlib import colors
	import scanpy as sc
	import seaborn as sns
	
	#mm10_genes=pd.read_csv('/Users/jonpreall/Data/Gene_lists/mm10_genes.csv')
	mm10_genes=pd.read_csv('https://www.dropbox.com/s/j0kqojg4ukihet6/mm10_genes.csv?dl=1')
	bg_list = mm10_genes['Gene_name'].tolist()
	
	gp = GProfiler(return_dataframe=True)
	
	df = gp.profile(genelist, organism=organism, background=bg_list, no_evidences=False, combined=False, ordered=True)
	print(df['source'].value_counts())
	return df
	
def gprofiler_human(genelist, organism = 'hsapiens', cmap="rainbow"):
	from gprofiler import GProfiler
	import pandas as pd
	import numpy as np
	import matplotlib.pyplot as pl
	from matplotlib import rcParams
	from matplotlib import colors
	import scanpy as sc
	import seaborn as sns
	
	#mm10_genes=pd.read_csv('/Users/jonpreall/Data/Gene_lists/mm10_genes.csv')
	hg38_genes=pd.read_csv('https://www.dropbox.com/s/7id4mjsc8f44uu5/GRCh38_genes.csv?dl=1')
	bg_list = hg38_genes['Gene_name'].tolist()
	
	gp = GProfiler(return_dataframe=True)
	
	df = gp.profile(genelist, organism=organism, background=bg_list, no_evidences=False, combined=False, ordered=True)
	print(df['source'].value_counts())
	return df
	
def gprofiler_plot(genelist, organism = 'mmusculus', cmap="rainbow", domain='BP'):
	from gprofiler import GProfiler
	import pandas as pd
	import numpy as np
	import matplotlib.pyplot as pl
	from matplotlib import rcParams
	from matplotlib import colors
	import scanpy as sc
	import seaborn as sns
	
	#mm10_genes=pd.read_csv('/Users/jonpreall/Data/Gene_lists/mm10_genes.csv')
	mm10_genes=pd.read_csv('https://www.dropbox.com/s/j0kqojg4ukihet6/mm10_genes.csv?dl=1')
	bg_list = mm10_genes['Gene_name'].tolist()
	
	gp = GProfiler(return_dataframe=True)
	
	df = gp.profile(genelist, organism=organism, background=bg_list, no_evidences=False, combined=False, ordered=True)
	#df = df.sort_values('p_value').iloc[:,[2,3,5,6,11]]


	#df['-logp'] = -np.log(df['p_value'])

	data_to_plot = df.copy()


	#df[df['term.size'] < 25].nlargest(columns='-logp', n=10).loc[:,['name','-logp']].set_index('term.name').plot(kind='barh')


	data_to_plot['go.id'] = data_to_plot.index
	data_to_plot = data_to_plot.nsmallest(columns='p_value', n=20)

	norm = colors.LogNorm(data_to_plot['p_value'].min()/10, data_to_plot['p_value'].max()*10)
	sm = pl.cm.ScalarMappable(cmap=cmap, norm=norm)
	sm.set_array([])
	
	pl.rcParams['figure.figsize'] = [4.0, 8.0]
	rcParams.update({'font.size': 14, 'font.weight': 'bold'})

	sns.set(style="whitegrid")

	path = pl.scatter(x='recall', y="name", c='p_value', cmap=cmap, norm=colors.LogNorm(data_to_plot['p_value'].min(), data_to_plot['p_value'].max()), data=data_to_plot, linewidth=1, edgecolor="grey", s=[(i+10)**1.5 for i in data_to_plot['intersection_size']])
	ax = pl.gca()
	ax.invert_yaxis()

	ax.set_ylabel('')
	ax.set_xlabel('Frac. of GO Genes recalled', fontsize=14, fontweight='bold')
	ax.xaxis.grid(False)
	ax.yaxis.grid(True)

	# Shrink current axis by 20%
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

	#Colorbar
	fig = pl.gcf()
	cbaxes = fig.add_axes([0.8, 0.15, 0.03, 0.4]) 
	cbar = ax.figure.colorbar(sm, ticks=[float(1e-14), float(1e-12), float(1e-10), float(1e-8), float(1e-6), float(1e-4)], shrink=0.5, anchor=(0,0.1), cax=cbaxes)
	cbar.ax.set_yticklabels(['$10^{-14}$', '$10^{-12}$', '$10^{-10}$', '$10^{-8}$', '$10^{-6}$', '$10^{-4}$'])
	cbar.set_label("Adjusted p_value", fontsize=14, fontweight='bold')

	#Size legend
	l1 = pl.scatter([],[], s=(5+10)**1.5, edgecolors='none', color='black')
	l2 = pl.scatter([],[], s=(25+10)**1.5, edgecolors='none', color='black')
	l3 = pl.scatter([],[], s=(50+10)**1.5, edgecolors='none', color='black')
	l4 = pl.scatter([],[], s=(75+10)**1.5, edgecolors='none', color='black')

	labels = ["5", "25", "50", "75"]

	leg = pl.legend([l1, l2, l3, l4], labels, ncol=1, frameon=False, fontsize=12,
	handlelength=1, loc = 'center left', borderpad = 1, labelspacing = 1.4,
	handletextpad=2, title='Gene overlap', scatterpoints = 1,  bbox_to_anchor=(-2, 1.5), facecolor='black')

	#pl.savefig('./figures/Haber_paneth_cell_GO_BP_enrichment.pdf', dpi=300, format='pdf')

	pl.show()
	pl.rcParams['figure.figsize'] = [4.0, 8.0]
	sc.settings.set_figure_params(dpi=80, color_map="Reds")



### gProfiler for TF enrichment:
#def gprofiler(genelist, organism = 'mmusculus', cmap="rainbow", domain='BP'):
def tf_gprofiler(adata, obs_key, cluster, n_genes=100, organism='mmusculus'):
	from gprofiler import GProfiler
	import jpplot
	import scanpy as sc
	import pandas as pd
	import numpy as np


	mm10_genes = pd.read_csv('https://www.dropbox.com/s/j0kqojg4ukihet6/mm10_genes.csv?dl=1')['Gene_name'].tolist()
	n_clust = len(adata.obs[obs_key].unique())

	sc.tl.rank_genes_groups(adata, groupby=obs_key)
	markerstemp = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
	markers_cols = markerstemp.columns.tolist()
	myclust_num = markers_cols.index(cluster)

	topdiff = jpplot.topdiff(adata, n_genes=n_genes)

	mylist = topdiff[myclust_num*n_genes:myclust_num*n_genes+n_genes]

	gp = GProfiler(return_dataframe=True)
	
	#gprofiler.gprofiler(mylist, organism='hsapiens', custom_bg=GRCh38_genes, correction_method='fdr', src_filter=['GO:tf'], ordered_query=True)
	#df = gprofiler.gprofiler(mylist, organism='mmusculus', custom_bg=mm10_genes, correction_method='fdr', ordered_query=True)
	df = gp.profile(mylist, organism=organism, background=bg_list, no_evidences=False, combined=False, ordered=True)
	tf_df = df[df['source'] == 'TF']
	tf_df = tf_df.sort_values('p_value', ascending=True)
	return tf_df

def import_degradation(adata, output_csv):
	import pandas as pd
	import numpy as np
	import scanpy as sc
	
	if adata.obs.index.str.replace('[ATGC]*','').unique()[0] == '':
		#print("Detected STARsolo style cell names, adding '-1' to match Cellranger output...")
		#adata.obs.index = [name + '-1' for name in adata.obs.index]
		numbered_adata = True
    
	sc.pp.filter_cells(adata, min_genes=1)
	sc.pp.filter_cells(adata, min_counts=1)
	
	degrade_df = pd.read_csv(output_csv, index_col=0)
	degrade_df.columns = ['total_alignment','degraded','percent_degraded']
	
	if degrade_df.index.str.replace('[ATGC]*','').unique()[0] == '':
		#print("Detected STARsolo style cell names, adding '-1' to match Cellranger output...")
		#adata.obs.index = [name + '-1' for name in adata.obs.index]
		numbered_degradation = True
	
	CR_filt_bc = [name for name in degrade_df.index if name in adata.obs.index]
	degrade_df = degrade_df.loc[CR_filt_bc,:]
	adata.obs = adata.obs.merge(degrade_df, left_index=True, right_index=True, how='left')
	sc.pl.scatter(adata, x='n_counts', y='n_genes', color='percent_degraded', color_map='jet')
	
def getsoup(adata, soup_counts=1e6):
	import pandas as pd
	import numpy as np
	import scanpy as sc
	
	adata.var_names_make_unique()
	adata.obs['n_counts'] = adata.X.sum(1)
	cumsum = lambda x,adata: adata[adata.obs['n_counts'] < x].X.sum()
	result_array = np.array([])

	for i in np.arange(0,50):
		result = cumsum(i,adata)
		result_array = np.append(result_array, result)
	z = pd.DataFrame(np.arange(0,50), columns=['n_counts'])
	z['cumulative_counts'] = result_array
	soup_thresh = z[z['cumulative_counts'] < soup_counts]['n_counts'].idxmax()
    
	soup = adata[adata.obs['n_counts'] < soup_thresh]
	print("Estimating soup using droplets with n_counts <",soup_thresh)
	print("Total reads in soup:", np.round(cumsum(soup_thresh, adata), 0).astype('int64'))
	return(soup)
	

	
def set_genome_build(org, use_local=None, local_gene_list_dir=None):
	import pandas as pd
	import numpy as np
	import os
	
	assert org in ['mm10', 'hg38','GRCh38','GRCm38'], 'organism must be one of: mm10, hg38, GRCh38, GRCm38'


	print('Example Usage: ribo_genes, mito_genes, S, G2M, IEG, cell_cycle = jpplot.set_genome_build(\'mm10\')')
	print('Currently available builds: mm10 and hg38')
	
	if not local_gene_list_dir:
		local_gene_list_dir='/Users/jpreall/Dropbox/Preall_Lab/Preall/Gene_Lists/'
		
	if not use_local:
		if os.path.exists(local_gene_list_dir):
			use_local=True
	
	if org in ['mm10','GRCm38']:
		if use_local:
			mm10_genes=pd.read_csv(local_gene_list_dir + '/mm10_genes.csv')['Gene_name'].tolist()
			IEG = pd.read_csv(local_gene_list_dir + 'IEG_mouse.txt')['Gene_name'].tolist()
			IEG = [name for name in IEG if name in mm10_genes]
			mouse_cc = pd.read_csv(local_gene_list_dir + 'regev_cc_mouse.csv')
			cell_cycle = mouse_cc['Gene_Name'].tolist()
			S = mouse_cc['Gene_Name'][mouse_cc['Stage'] == 'S'] .tolist()
			G2M = mouse_cc['Gene_Name'][mouse_cc['Stage']== 'G2M'] .tolist()

			mito_genes = [name for name in mm10_genes if name.startswith('mt-')]
			Rpl_genes = [name for name in mm10_genes if name.startswith('Rpl')]
			Rps_genes = [name for name in mm10_genes if name.startswith('Rps')]
			ribo_genes = Rpl_genes + Rps_genes
			cell_cycle = S + G2M
			return ribo_genes, mito_genes, S, G2M, IEG, cell_cycle
		else:
			mm10_genes=pd.read_csv('https://www.dropbox.com/s/j0kqojg4ukihet6/mm10_genes.csv?dl=1')['Gene_name'].tolist()
			IEG = pd.read_csv('https://www.dropbox.com/s/dhy1h8r8fp23ps0/IEG_mouse.txt?dl=1')['Gene_name'].tolist()
			IEG = [name for name in IEG if name in mm10_genes]
			mouse_cc = pd.read_csv('https://www.dropbox.com/s/mlfaltik2b3p16q/regev_cc_mouse.csv?dl=1')
			cell_cycle = mouse_cc['Gene_Name'].tolist()
			S = mouse_cc['Gene_Name'][mouse_cc['Stage'] == 'S'] .tolist()
			G2M = mouse_cc['Gene_Name'][mouse_cc['Stage']== 'G2M'] .tolist()

			mito_genes = [name for name in mm10_genes if name.startswith('mt-')]
			Rpl_genes = [name for name in mm10_genes if name.startswith('Rpl')]
			Rps_genes = [name for name in mm10_genes if name.startswith('Rps')]
			ribo_genes = Rpl_genes + Rps_genes
			cell_cycle = S + G2M
			return ribo_genes, mito_genes, S, G2M, IEG, cell_cycle
        
	elif org in ['hg38','GRCh38']:
		if use_local:
			GRCh38_genes=pd.read_csv(local_gene_list_dir + '/GRCh38_genes.csv')['Gene_name'].tolist()
			IEG = pd.read_csv(local_gene_list_dir + 'IEG_human.csv')['Gene_name'].tolist()
			IEG = [name for name in IEG if name in GRCh38_genes]
			human_cc = pd.read_csv(local_gene_list_dir + 'regev_lab_cell_cycle_genes.csv')
			cell_cycle = human_cc['Gene_Name'].tolist()
			S = human_cc['Gene_Name'][human_cc['Stage'] == 'S'] .tolist()
			G2M = human_cc['Gene_Name'][human_cc['Stage']== 'G2M'] .tolist()

			mito_genes = [name for name in GRCh38_genes if name.startswith('MT-')]
			Rpl_genes = [name for name in GRCh38_genes if name.startswith('RPL')]
			Rps_genes = [name for name in GRCh38_genes if name.startswith('RPS')]
			ribo_genes = Rpl_genes + Rps_genes
			cell_cycle = S + G2M
			return ribo_genes, mito_genes, S, G2M, IEG, cell_cycle
		else:
			GRCh38_genes=pd.read_csv('https://www.dropbox.com/s/7id4mjsc8f44uu5/GRCh38_genes.csv?dl=1')['Gene_name'].tolist()
			human_cc = pd.read_csv('https://www.dropbox.com/s/ipkg33kqlm0hm7x/regev_lab_cell_cycle_genes.csv?dl=1')
			IEG = pd.read_csv('https://www.dropbox.com/s/zpkre5bq8gkfi1f/IEG_human.csv?dl=1')['Gene_name'].tolist()
			IEG = [name for name in IEG if name in GRCh38_genes]
			S = human_cc['Gene_Name'][human_cc['Stage'] == 'S'] .tolist()
			G2M = human_cc['Gene_Name'][human_cc['Stage']== 'G2M'] .tolist()
			mito_genes = [name for name in GRCh38_genes if name.startswith('MT-')]
			Rpl_genes = [name for name in GRCh38_genes if name.startswith('RPL')]
			Rps_genes = [name for name in GRCh38_genes if name.startswith('RPS')]
			ribo_genes = Rpl_genes + Rps_genes
			cell_cycle = S + G2M
			return ribo_genes, mito_genes, S, G2M, IEG, cell_cycle
	
def split_draw_graph(adata, gene, s=100, cmap="Reds", fig_scale=1, frameon=True, use_raw=True, batch_key='Sample', layer=None, vmin=None):
	import pandas as pd
	import numpy as np
	import matplotlib.pyplot as pl
	from matplotlib import rcParams
	from matplotlib import colors
	import scanpy as sc
	import seaborn as sns
	
		
	samples = adata.obs[batch_key].unique()
	assert len(samples) == 2, 'Split UMAP function only works when comparing 2 samples'
    
	tmp1 = adata[adata.obs[batch_key] == samples[0]]
	tmp2 = adata[adata.obs[batch_key] == samples[1]]    

	if gene in adata.raw.var_names:
		genenum = adata.raw.var_names.get_loc(gene)

		v1 = tmp1.raw.obs_vector(k=genenum).max()
		v2 = tmp2.raw.obs_vector(k=genenum).max()
		both = np.asarray((v1,v2))
		vmax = both.max()
        
		fig, ax = pl.subplots(1,2,figsize=(6*fig_scale,2.75*fig_scale))
		sc.pl.draw_graph(tmp1, color=gene, s=s, title=str(gene) + ': ' + str(samples[0]), ax=ax[0], vmax=vmax, vmin=vmin, cmap=cmap, show=False, frameon=frameon, use_raw=use_raw)
		sc.pl.draw_graph(tmp2, color=gene, s=s, title=str(gene) + ': ' + str(samples[1]), ax=ax[1], vmax=vmax, vmin=vmin, cmap=cmap, show=False, frameon=frameon, use_raw=use_raw)
		pl.tight_layout()
    
	else:
		if adata.obs[gene].dtype.kind in 'if': # check if the key is numerical
			fig, ax = pl.subplots(1,2,figsize=(6*fig_scale,2.75*fig_scale))
			v1 = tmp1.obs[gene].max()
			v2 = tmp2.obs[gene].max()
			both = np.asarray((v1,v2))
			vmax = both.max()
        
			sc.pl.draw_graph(tmp1, color=gene, s=s, title=str(gene) + ': ' + str(samples[0]), ax=ax[0], vmax=vmax, vmin=vmin, cmap=cmap, show=False, legend_loc='None', frameon=frameon, use_raw=use_raw)
			sc.pl.draw_graph(tmp2, color=gene, s=s, title=str(gene) + ': ' + str(samples[1]), ax=ax[1], vmax=vmax, vmin=vmin, cmap=cmap, show=False, frameon=frameon, use_raw=use_raw)
			pl.tight_layout()
            
		else: #if it is not numerical, plot as categorical
			fig, ax = pl.subplots(1,2,figsize=(6*fig_scale,2.75*fig_scale))
			sc.pl.draw_graph(tmp1, color=gene, s=s, title=str(gene) + ': ' + str(samples[0]), ax=ax[0], cmap=cmap, show=False, legend_loc='None', frameon=frameon, use_raw=use_raw)
			sc.pl.draw_graph(tmp2, color=gene, s=s, title=str(gene) + ': ' + str(samples[1]), ax=ax[1], cmap=cmap, show=False, frameon=frameon, use_raw=use_raw)
			pl.tight_layout()
        
	del(tmp1)
	del(tmp2)
		
def split_umap(adata, gene, s=100, cmap="Reds", fig_scale=1, frameon=True, use_raw=True, batch_key='Sample', layer=None, vmin=None):
	import pandas as pd
	import numpy as np
	import matplotlib.pyplot as pl
	from matplotlib import rcParams
	from matplotlib import colors
	import scanpy as sc
	import seaborn as sns
	
		
	samples = adata.obs[batch_key].unique()
	assert len(samples) == 2, 'Split UMAP function only works when comparing 2 samples'
    
	tmp1 = adata[adata.obs[batch_key] == samples[0]]
	tmp2 = adata[adata.obs[batch_key] == samples[1]]    

	if gene in adata.raw.var_names:
		genenum = adata.raw.var_names.get_loc(gene)

		v1 = tmp1.raw.obs_vector(k=genenum).max()
		v2 = tmp2.raw.obs_vector(k=genenum).max()
		both = np.asarray((v1,v2))
		vmax = both.max()
        
		fig, ax = pl.subplots(1,2,figsize=(6*fig_scale,2.75*fig_scale))
		sc.pl.umap(tmp1, color=gene, s=s, title=str(gene) + ': ' + str(samples[0]), ax=ax[0], vmax=vmax, vmin=vmin, cmap=cmap, show=False, frameon=frameon, use_raw=use_raw)
		sc.pl.umap(tmp2, color=gene, s=s, title=str(gene) + ': ' + str(samples[1]), ax=ax[1], vmax=vmax, vmin=vmin, cmap=cmap, show=False, frameon=frameon, use_raw=use_raw)
		pl.tight_layout()
    
	else:
		if adata.obs[gene].dtype.kind in 'if': # check if the key is numerical
			fig, ax = pl.subplots(1,2,figsize=(6*fig_scale,2.75*fig_scale))
			v1 = tmp1.obs[gene].max()
			v2 = tmp2.obs[gene].max()
			both = np.asarray((v1,v2))
			vmax = both.max()
        
			sc.pl.umap(tmp1, color=gene, s=s, title=str(gene) + ': ' + str(samples[0]), ax=ax[0], vmax=vmax, vmin=vmin, cmap=cmap, show=False, legend_loc='None', frameon=frameon, use_raw=use_raw)
			sc.pl.umap(tmp2, color=gene, s=s, title=str(gene) + ': ' + str(samples[1]), ax=ax[1], vmax=vmax, vmin=vmin, cmap=cmap, show=False, frameon=frameon, use_raw=use_raw)
			pl.tight_layout()
            
		else: #if it is not numerical, plot as categorical
			fig, ax = pl.subplots(1,2,figsize=(6*fig_scale,2.75*fig_scale))
			sc.pl.umap(tmp1, color=gene, s=s, title=str(gene) + ': ' + str(samples[0]), ax=ax[0], cmap=cmap, show=False, legend_loc='None', frameon=frameon, use_raw=use_raw)
			sc.pl.umap(tmp2, color=gene, s=s, title=str(gene) + ': ' + str(samples[1]), ax=ax[1], cmap=cmap, show=False, frameon=frameon, use_raw=use_raw)
			pl.tight_layout()
        
	del(tmp1)
	del(tmp2)

def diffexplot(adata, genes, use_raw=False, groupby='Sample', sig_thresh=1, cmap='agsunset', fig_height=300):
    import scanpy as sc
    import numpy as np
    import matplotlib.pyplot as pl
    import pandas as pd
    import plotly.express as px
    import jpplot
    
    samples = adata.obs[groupby].unique()
    assert len(samples) == 2, 'Diffexplot function only works when comparing 2 samples'
    
    sc.tl.rank_genes_groups(adata, groupby=groupby, use_raw=use_raw, n_genes=len(adata.var_names))
    df = jpplot.get_markers(adata, n_markers=len(adata.var_names))
    sampname1 = [name for name in df.columns if name.endswith('_names')][0].replace('_names','')
    sampname2 = [name for name in df.columns if name.endswith('_names')][1].replace('_names','')
    df = df.iloc[:,:3]

    name_col = [name for name in df.columns if name.endswith('_names')]
    logfold_col = [name for name in df.columns if name.endswith('_logfold')]
    pval_col = [name for name in df.columns if name.endswith('_pvals')]

    df.set_index(name_col, inplace=True)

    df = df.loc[genes]
    df['gene'] = df.index
    df = df.sort_values(by=logfold_col)
    df['significant'] = df[pval_col] < sig_thresh
    df = df[df['significant'] == True]
    df['-logP'] = -np.log10(df[pval_col])

    #sigmap = {True:'red',False:'lightgrey'}
    #colors = np.array(df['significant'].map(sigmap)).tolist()
    
    fig = px.bar(df, y='gene', x=logfold_col[0],
                 hover_data=[pval_col[0], 'significant'], 
                 color='-logP',
                 labels={logfold_col[0]:'log2('+sampname1+'/'+sampname2+')'}, 
                 #height=150+15*len(genes),
                 height=fig_height,
                 width=300,
                 color_continuous_scale=cmap,
                 orientation='h'
                )
    fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)',
                       'paper_bgcolor': 'rgba(0, 0, 0, 0)',
                      })
    fig.show()

def corrplot(adata, corr_thresh=0.6, tree_trim=3, min_frac_cells=0, font_scale=1, max_matrix_size=3e6):
    import pandas as pd
    import scanpy as sc
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as pl
    import scipy.cluster as cluster
    from scipy import sparse
    import sys
    
    print('To return cut tree and correlation dataframe, run as \'cutree, corrdf = jpplot.corrplot(adata...)\'')
    
    hvgdf = adata.var[adata.var['highly_variable'] == True]
    
    #Filter HVGs to only include those that appear in at least (CUTOFF%) cells
    hvgdf = hvgdf[hvgdf['n_cells'] > min_frac_cells*len(adata.obs)]
    hivar = hvgdf.index.tolist()
    
    bdata = adata.copy()
    bdata.X = sparse.csr_matrix(adata.X)
    df = pd.DataFrame(bdata[:,hivar].X.todense())
    df.index = bdata.obs.index
    df.columns = bdata[:,hivar].var_names
    del(bdata)
    
    matrix_size=df.shape[0]*df.shape[1]
    print('matrix size:', matrix_size)
    if matrix_size > 3000000:
        print("Warning, matrix is really big and might take a while...")
    if matrix_size > max_matrix_size:
        print("Matrix too big, aborting...")
        sys.exit("Goodbye.")
    else:
        print("Matrix is ok sized, attempting...")
        corrdf = df.corr()
        corrdf = corrdf.replace(1,0)
        keep_genes = corrdf[corrdf.max(1) > corr_thresh].index.tolist()
        corrdf = corrdf.loc[keep_genes,:]
        corrdf = corrdf.loc[:,keep_genes]
        Z = cluster.hierarchy.ward(corrdf)
        
        cutree = cluster.hierarchy.cut_tree(Z, n_clusters=[tree_trim])
        cutree = pd.DataFrame(cutree)
        cutree.index = corrdf.index
        
        #Create first layer of labels
        labels = np.array(cutree[0])
        lut = dict(zip(set(labels), sns.hls_palette(len(set(labels)), l=0.5, s=0.8)))
        row_colors = pd.DataFrame(labels)[0].map(lut)
        
        sns.set(font_scale=font_scale)
        sns.clustermap(corrdf, col_cluster=True, cmap='coolwarm', row_colors=[row_colors])
        sc.settings.set_figure_params(dpi=80, color_map="Reds")
        
        return cutree, corrdf
        pl.show()
        
def plot_3d(adata, use_rep='tsne', key='Cluster', n_pcs=50, size=10, color_map="Reds", marker_line_color="black", marker_line_width=1, components=(0,1,2), layer='raw'):
    import plotly.express as px
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as pl
    import scanpy as sc
    import sys

    projections = [name.replace('X_','') for name in list(adata.obsm.keys())]
    if use_rep not in projections:
        print('Warning, use_rep must specify one of ',projections)
        sys.exit()
        
    if adata.obsm['X_'+str(use_rep)].shape[1] < 3:
        print('N '+str(use_rep) + ' dimensions is not 3.  Re-running with 3 dimensions...')
        
        if use_rep == 'tsne':
            sc.tl.tsne(adata, n_pcs=n_pcs, n_components=3)
        
        elif use_rep == 'umap':
            sc.tl.umap(adata, n_components=3)

        elif use_rep == 'diffmap':
            sc.tl.diffmap(adata)
            adata.obsm['X_diffmap'] = adata.obsm['X_diffmap'][:,1:]
	        
    #else:
        
    coords = pd.DataFrame(adata.obsm['X_'+str(use_rep)])
    coords.columns = ['dim' + str(name) for name in coords.columns]
    coords['size'] = size

    if key in adata.obs.columns:
        coords[key] = np.array(adata.obs[key])
    	
    elif key in adata.var_names:
        if layer == 'raw':
            from scipy.sparse import issparse
            if issparse(adata.raw.X):
                coords[key] = np.array(adata.raw[:,key].X.toarray())
            else:
                coords[key] = np.array(adata.raw[:,key].X)
                
        elif layer in adata.layers.keys():
            from scipy.sparse import issparse
            if issparse(adata.layers[layer]):
                coords[key] = np.array(adata[:,key].layers[layer].toarray())
            else:
                coords[key] = np.array(adata[:,key].layers[layer])
                
        elif layer == 'X':
            from scipy.sparse import issparse
            if issparse(adata.X):
                coords[key] = np.array(adata[:,key].X.toarray())
            else:
                coords[key] = np.array(adata[:,key].X)
    	
    import plotly.express as px
    
    ## Color map
    def vdir(obj):
        return [x for x in dir(obj) if not x.startswith('_')]
        
    if color_map not in vdir(px.colors.sequential):
        print('color_map must be one of:', vdir(px.colors.sequential))
        
        
    colors = getattr(px.colors.sequential, color_map)
        
    #elif color_map not in vdir(px.colors.diverging):
    #    print('Invalid color map')
    #    print('color_map must be one of:', vdir(px.colors.diverging))

    xdim = 'dim' + str(components[0])
    ydim = 'dim' + str(components[1])
    zdim = 'dim' + str(components[2])
	
    fig = px.scatter_3d(coords, x=xdim, y=ydim, z=zdim, 
    	size='size',
    	size_max=size,
    	color=key,
    	opacity=1,
    	color_continuous_scale=colors,
    	)
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
    fig.update_traces(dict(marker_line_width=marker_line_width, marker_line_color=marker_line_color))

    fig.show()

def preprocess_annotate(adata, genome='', exclude_highly_expressed=True, min_counts=1, min_cells=1, min_genes=1, exclude_genes=[], score_only=False, use_local=False):
    import jpplot
    import pandas as pd
    import numpy as np
    import scanpy as sc
    import sys
    
    ## Detect genome
    if genome == '':
        jpplot.detect_genome(adata)
        genome = adata.uns['genome']
		#print('Detected genome: ',genome)
    else:
        print('Using specified genome: ',genome)
    
    if len(adata.obs) > 1e6:
    	print('Warning, it appears you are preprocessing a raw, unfiltered data matrix.')
    	
    	if min_counts > 200:
    		print('You are filtering with a UMI cutoff >= 200. This should adequately trim the matrix.  Proceeding...')
    		#continue
    	
    	else:
    		print('Try re-running with min_counts=200 or greater')
    		exit()
    		
    assert genome in ['mm10', 'hg38', 'GRCh38','GRCm38'], 'organism must be one of: mm10, hg38'
    assert exclude_highly_expressed in [True, False], 'exclude_highly_expressed must be one of: True, False'

    #print('Example Usage: ribo_genes, mito_genes, S, G2M, IEG, cell_cycle = jpplot.set_genome_build(\'mm10\')')
    #print('Currently available builds: mm10 and hg38')
    
    if adata.var_names.duplicated().any():
        print('Duplicate gene names detected. Running adata.var_names_make_unique() to deduplicate...')
        adata.var_names_make_unique()
	
    ribo_genes, mito_genes, S, G2M, IEG, cell_cycle = jpplot.set_genome_build(genome, use_local=use_local)
	
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    adata.obs['log10GenesPerUMI'] = np.log10(adata.obs['n_genes']) / np.log10(adata.obs['n_counts'])
    
    num_orig_mito_genes = len(mito_genes)
    mito_genes = list(set(mito_genes).intersection(adata.var_names))
    num_filt_mito_genes = len(mito_genes)
    
    num_orig_ribo_genes = len(ribo_genes)
    ribo_genes = list(set(ribo_genes).intersection(adata.var_names))
    num_filt_ribo_genes = len(ribo_genes)
    
    IEG = list(set(IEG).intersection(adata.var_names))
    cell_cycle = list(set(cell_cycle).intersection(adata.var_names))

	
    ##Count percent mitochondrial genes and calculate percentage per cell
    if num_filt_mito_genes > 2:
        adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
    else:
        print("Mitochondrial genes not detected in matrix.  Skipping percent mito calculation.")
        
	## add the total counts per cell as observations-annotation to adata
	## Legacy code - can remove later
	## this should already be taken care of by the Scanpy filtering steps above
	#adata.obs['n_counts'] = np.sum(adata.X, axis=1)
	
	##Count IEG genes and calculate percentage per cell
    #filt_IEG = [name for name in IEG if name in adata.var_names]
    #adata.obs['percent_IEG'] = np.sum(adata[:, filt_IEG].X, axis=1) / np.sum(adata.X, axis=1)
    sc.tl.score_genes(adata, gene_list=IEG, score_name='IEG_score')	
    
	##Calculate the percent Ribo genes per cell
    if num_filt_ribo_genes > 2:
        adata.obs['percent_ribo'] = np.sum(adata[:, ribo_genes].X, axis=1) / np.sum(adata.X, axis=1)
    else:
        print("Ribosomal genes not detected in matrix.  Skipping percent ribo calculation.")
    	
	#Calculate the percent cell cycle genes per cell
    #adata.obs['percent_cell_cycle'] = np.sum(adata[:, cell_cycle].X, axis=1) / np.sum(adata.X, axis=1)

	
    ## Save the raw counts to the layer adata.layers['counts']   
    print ('Saving raw count data to .layers[\'counts\']')
    counts = adata.X.astype(int)
    adata.layers['counts'] = counts
    
    if score_only == False:
        ## Optional: filter away a list of genes from the raw matrix.
        exclude_genes = list(set(exclude_genes).intersection(adata.var_names))
        print('Removing genes from matrix prior to normalization:', len(exclude_genes), 'including ',exclude_genes[:4])
        keep_genes = list(adata.var_names.difference(set(exclude_genes)))
        adata._inplace_subset_var(keep_genes)
        print('current data shape:',adata.shape)

	
        ## Normalize, and log transform:

        print ('Normalizing total, target sum = 1e6, exclude_highly_expressed = ' + str(exclude_highly_expressed))
        sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=exclude_highly_expressed, key_added='Norm_Factor')
        print ('Log transforming...')
        sc.pp.log1p(adata)
        print ('Saving normalized log-transformed data to .raw')
        adata.raw = adata.copy()
        print('Done.')
        return(adata)
		
    else:
        print('Skipping log-normalization and generation of raw data layer')
    
    #Annotate genes with chromosome start and stop coordinates
    if genome == 'mm10':
        if use_local == True:
            var_coords = pd.read_csv('/Users/jpreall/Dropbox/Preall_Lab/Preall/Gene_Lists/mm10_genes_with_coords.csv', index_col=0)
        else:
            var_coords = pd.read_csv('https://www.dropbox.com/s/alli12d3jp4dz09/mm10_genes_with_coords.csv?dl=1', index_col=0)
        shared_genes = set(adata.var.index).intersection(var_coords.index)
        adata.var = adata.var.merge(var_coords.loc[shared_genes, ['Chromosome','Start','End','Strand']], left_index=True, right_index=True, how='left')


def gprofiler_summary(adata, groupby, source='KEGG', species='human', n_per_cluster=2, n_genes=1000, sig_thresh=1e-3, min_fold_change=1.5):
    import jpplot
    import pandas as pd
    import numpy as np
    import scanpy as sc
    import os
    
    print('Usage: gprofiler_summary, gprofiler_summary_filt = jpplot.gprofiler_summary(adata, groupby=\'Cluster\'...')
    
    exclude_terms = ['KEGG root term', 'REACTOME root term','CORUM root','WIKIPATHWAYS']

    output = pd.DataFrame()
    sc.tl.rank_genes_groups(adata, groupby=groupby, n_genes=n_genes)
	
    if species == 'human':
	
        for i in adata.obs[groupby].unique():
            print('Running gprofiler on', str(i) + '...')
            df = sc.get.rank_genes_groups_df(adata, group=i)
            df = df[df['pvals_adj'] < sig_thresh]
            df = df[df['logfoldchanges'] > min_fold_change]
            mygenes = df['names'].tolist()
            #mygenes = sc.get.rank_genes_groups_df(adata, group=i)['names'].tolist()
            gp = jpplot.gprofiler_human(mygenes)
            gp[groupby] = str(i)
            output = pd.concat([output,gp])
    
    elif species == 'mouse':
	
        for i in adata.obs[groupby].unique():
            print('Running gprofiler on', str(i) + '...')
            df = sc.get.rank_genes_groups_df(adata, group=i)
            df = df[df['pvals_adj'] < sig_thresh]
            df = df[df['logfoldchanges'] > min_fold_change]
            mygenes = df['names'].tolist()
            #mygenes = sc.get.rank_genes_groups_df(adata, group=i)['names'].tolist()
            gp = jpplot.gprofiler_mouse(mygenes)
            gp[groupby] = str(i)
            output = pd.concat([output,gp])
	
    else:
        print('Error: unrecognized species.  Use human or mouse')

    name_values = output['name'].unique()
    keep_names = [i for i in name_values if i not in exclude_terms]
    keep_names = [i for i in keep_names if 'root term' not in i]
    output = output[output['name'].isin(keep_names)]

    print('Total term counts:', output['source'].value_counts())

    ## filter output
    output_filt = pd.DataFrame()

    for i in output[groupby].unique():
        out_split = output[output[groupby] == i]
        out_split = out_split[out_split['source'] == source][:n_per_cluster]
        output_filt = pd.concat([output_filt,out_split]).loc[:,[groupby,'source','name','p_value','intersections']]
     
    return output, output_filt
    output_filt
    
def filter_gprofiler_summary(gprofiler_summary, source='KEGG', n_per_cluster=2):
    import jpplot
    import pandas as pd
    
    
    print('Filter profiler summary table by source.')  
    print('First, run:\' gprofiler_summary = jpplot.gprofiler_summary(adata, groupby=\'...\'')
    print('Options: source=', gprofiler_summary['source'].unique().tolist())
    print()

    print('Total term counts:', gprofiler_summary['source'].value_counts())

    ## filter output
    output_filt = pd.DataFrame()
    groupby = gprofiler_summary.columns[-1]
	
    for i in gprofiler_summary[groupby].unique():
        out_split = gprofiler_summary[gprofiler_summary[groupby] == i]
        out_split = out_split[out_split['source'] == source][:n_per_cluster]
        output_filt = pd.concat([output_filt,out_split]).loc[:,[groupby,'source','name','p_value','intersections']]
     
    return output_filt
    
def corr_with_gene(adata, key, n_genes=None, layer='raw'):
	import jpplot
	import pandas as pd
	import numpy as np
	
	from scipy.sparse import issparse
	
	assert layer in list(adata.layers.keys()) + ['X','raw'],'layer must be a key in adata.layers, X, or raw. Default: raw'
	if layer == 'raw':
		if not adata.raw:
			adata.raw = adata
			
		if issparse(adata.raw.X):
			dense_matrix = pd.DataFrame(adata.raw.X.todense())
			dense_matrix.columns = adata.raw.var_names
		else:
			dense_matrix = pd.DataFrame(adata.raw.X)
			dense_matrix.columns = adata.raw.var_names

	elif layer == 'X':
		if issparse(adata.X):
			dense_matrix = pd.DataFrame(adata.X.todense())
			dense_matrix.columns = adata.var_names
		else:
			dense_matrix = pd.DataFrame(adata.X)
			dense_matrix.columns = adata.var_names
				
	else:
		if issparse(adata.layers[layer]):
			dense_matrix = pd.DataFrame(adata.layers[layer].todense())
			dense_matrix.columns = adata.var_names
			
		else:
			dense_matrix = pd.DataFrame(adata.layers[layer])
			dense_matrix.columns = adata.var_names
	
	if not n_genes:
			n_genes = dense_matrix.shape[1]
	
	valid_keys = adata.var_names.tolist() + adata.obs.select_dtypes('number').columns.tolist()		
	assert key in valid_keys, 'key must be a gene in adata.var_names or a numerical column in adata.obs'
	
	if key in adata.obs.select_dtypes('number').columns.tolist():
		dense_matrix.loc[:,key] = np.array(adata.obs.loc[:,key])
	
	correlation = dense_matrix.corrwith(dense_matrix.loc[:,key])
	correlation = correlation.nlargest(n_genes)
	return correlation
	

def pyScTransform(adata, output_file=None):
    import scanpy as sc
    from scipy.sparse import issparse

    """
    Function to call scTransform from Python
    """
    import rpy2.robjects as ro
    import anndata2ri

    ro.r('library(Seurat)')
    ro.r('library(scater)')
    anndata2ri.activate()

    sc.pp.filter_genes(adata, min_cells=5)
    
    if issparse(adata.X):
        if not adata.X.has_sorted_indices:
            adata.X.sort_indices()

    for key in adata.layers:
        if issparse(adata.layers[key]):
            if not adata.layers[key].has_sorted_indices:
                adata.layers[key].sort_indices()

    ro.globalenv['tmp'] = adata

    ro.r('seurat_obj = as.Seurat(adata, counts="X", data = NULL)')

    ro.r('res <- SCTransform(object=seurat_obj, return.only.var.genes = FALSE, do.correct.umi = FALSE)')

    norm_x = ro.r('res@assays$SCT@scale.data').T

    adata.layers['normalized'] = norm_x

    if output_file:
        adata.write(output_file)
        
def pyScTransform2(adata, output_file=None):
    import scanpy as sc
    from scipy.sparse import issparse

    """
    Function to call scTransform from Python
    """
    import rpy2.robjects as ro
    import anndata2ri
    

    ro.r('library(sctransform)')
    #ro.r('library(scater)')
    anndata2ri.activate()

    sc.pp.filter_genes(adata, min_cells=5)
    tmp=adata.copy()
    
    if issparse(tmp):
        if not tmp.X.has_sorted_indices:
            tmp.X.sort_indices()

    for key in tmp.layers:
        if issparse(tmp.layers[key]):
            if not tmp.layers[key].has_sorted_indices:
                tmp.layers[key].sort_indices()

    
    ro.globalenv['tmp'] = tmp
    ro.globalenv['float_counts'] = tmp.layers['counts'].astype('float')
    ro.r('data_matrix = as.matrix(float_counts)')
    ro.r('head(data_matrix)')
    
    #ro.r('seurat_obj = as.Seurat(adata, counts="X", data = NULL)')
    #

    #ro.r('data.sct <- sctransform::vst(float_counts, latent_var = c("log_umi"), return_gene_attr = TRUE, return_cell_attr = TRUE)')

   # sct_data = ro.r('as.matrix(yh@assays$SCT@data)').T

    #adata.layers['SCT'] = sct_data

    if output_file:
        adata.write(output_file)
        
def convert_mex_to_h5(datadir, genome=None, file_prefix='filtered_feature_bc_matrix'):
    import os
    import numpy as np
    import pandas as pd
    import h5py
    import h5sparse
    from scipy.io import mmread
    from scipy.sparse import csr_matrix
    
    if not genome:
        genome = 'custom_genome'
        
    output_folder = os.getcwd()
    #file_prefix = 'filtered_feature_bc_matrix'
    outfile=os.path.join(output_folder,file_prefix +'.h5')
    dtype = 'float32'
    matrixfile = [name for name in os.listdir(datadir) if name.startswith('matrix')][0]
    barcodesfile = [name for name in os.listdir(datadir) if name.startswith('barcodes')][0]
    featuresfile = [name for name in os.listdir(datadir) if name.startswith(('features','genes'))][0]

    X = mmread(matrixfile).astype(dtype)
    MATRIX = csr_matrix(X.T)
    BCS = np.array(pd.read_table(barcodesfile, header=None)[0]).astype('S')
    feature_df = pd.read_table(featuresfile, header=None)
    FEATURES = np.array(feature_df[1]).astype('S')
    FEATURE_IDS = np.array(feature_df[0]).astype('S')
    GENOME=genome
    all_tag_keys = np.array([b'genome'])
    LIBRARY_IDS = np.array([file_prefix.encode()])
    ORIG_GEM_GROUPS = np.array([1])
    SHAPE = MATRIX.T.shape

    with h5sparse.File(outfile, 'w') as h5f:
        h5f.create_dataset('matrix/', data=MATRIX, compression="gzip")
        h5f.close()
    
    with h5py.File(os.path.join(output_folder,file_prefix +'.h5'), 'r+') as f:
        f.create_dataset('matrix/barcodes', data=BCS)
        f.create_dataset('matrix/shape', (2,),dtype='int32', data=SHAPE)
        features = f.create_group('matrix/features')
        features.create_dataset('_all_tag_keys', (1,),'S6', data=all_tag_keys)  
        features.create_dataset('feature_type', data=np.array([b'Gene Expression'] * SHAPE[0]))
        features.create_dataset('genome', data=np.array([GENOME.encode()] * SHAPE[0]))
        features.create_dataset('id', data=FEATURE_IDS)
        features.create_dataset('name', data=FEATURES)
        
        f.attrs['chemistry_description'] = b'Single Cell 3\' v3'
        f.attrs['filetype'] = 'matrix'
        f.attrs['library_ids'] = LIBRARY_IDS
        f.attrs['original_gem_groups'] = ORIG_GEM_GROUPS
        f.attrs['version'] = 2
        f.close()
        
def mex_to_h5_old(matrixdir, chemistry_version='3', genome='', input_format='CellrangerV3', output_prefix=''):
	import os
	import pandas as pd
	import gzip
	
	
	tmp_genes = matrixdir + "features.tsv.gz"
	tmp_bar = matrixdir + "barcodes.tsv.gz"
	tmp_mtx = matrixdir + "matrix.mtx.gz"
	if output_prefix:
		output_prefix = output_prefix + '_'
    
	#matrixfile = matrixdir + '.mtx.gz'
	if  os.path.exists(matrixdir + "features.tsv.gz"):
		print('Detected features.tsv.gz, this seems like a Cellranger v3 output folder...')
		
	if not os.path.exists(matrixdir + "features.tsv.gz"):
	
		## convert v2 or STARsolo genes.tsv files to features.tsv.gz files
		if os.path.exists(matrixdir + "genes.tsv"):
			print('Detected genes.tsv, converting to v3 features.tsv.gz format...')
			genes = pd.read_table(matrixdir + 'genes.tsv', header=None)
			genes[2] = 'Gene Expression'
			genes = genes.iloc[:,:3]
			genes.to_csv(tmp_genes, sep="\t", compression='gzip', header=None, index=None)
		
		## Compress any uncompressed matrix file:
		if not os.path.exists(matrixdir + "matrix.mtx.gz"):
			if os.path.exists(matrixdir + "matrix.mtx"):
				print('Detected matrix.mtx, converting to v3 matrix.mtx.gz format...')
		
				with open(matrixdir + "matrix.mtx", 'rb') as orig_file:
					with gzip.open(matrixdir + "matrix.mtx.gz", 'wb') as zipped_file:
						zipped_file.writelines(orig_file)
		
		## Compress any uncompressed barcodes file:
		if not os.path.exists(matrixdir + "barcodes.tsv"):
			if os.path.exists(matrixdir + "barcodes.tsv"):
				print('Detected barcodes.tsv, converting to v3 barcodes.tsv.gz format...')
				
				import gzip
				f_in = open(matrixdir + "barcodes.tsv", 'rb')
				f_out = gzip.open(tmp_bar, 'wb')
				f_out.writelines(f_in)
				f_out.close()
				f_in.close()

	
	#outdir = os.path.dirname(os.path.dirname(matrixdir)) +'/'
	outdir = os.path.dirname(matrixdir) +'/'
	
	assert chemistry_version in ['2','3'], 'Chemistry version must be one of: 2, 3. Default: 3'
	assert input_format in ['CellrangerV2','CellrangerV3','STARSolo'], 'Chemistry version must be one of: CellrangerV2, CellrangerV3, STARSolo. Default: CellrangerV3'
	assert genome in ['mm10', 'hg38','GRCh38','GRCm38'], 'organism must be one of: mm10, hg38, GRCh38, GRCm38'
	
	if 'raw' in os.path.basename(matrixdir.split('/')[-2]):
		outfile = 'raw_feature_bc_matrices_h5.h5'
	elif 'filtered' in os.path.basename(matrixdir.split('/')[-2]):
		outfile = 'filtered_feature_bc_matrices_h5.h5'
	else:
		outfile = 'gene_bc_matrix.h5'

	print('Reading in matrix from ', matrixdir)
	print('Assuming version ',chemistry_version,'chemistry...')
	
	import rpy2.robjects as ro
    #import anndata2ri

	ro.globalenv['matrixdir'] = matrixdir
	ro.globalenv['genome'] = genome
	ro.globalenv['input_format'] = input_format
	ro.globalenv['chemistry_version'] = chemistry_version
	ro.globalenv['outdir'] = outdir
	ro.globalenv['outfile'] = outfile


	ro.r('library(DropletUtils)')
	ro.r('library(Matrix)')
    #anndata2ri.activate()
	
	ro.r('matrix <- read10xCounts(matrixdir, col.names=FALSE, type="sparse", version="3", genome=NULL)')
    
	ro.r('output_file = file.path(outdir,outfile)')
	ro.r('write10xCounts(path = output_file, x = counts(matrix), barcodes = colData(matrix)$Barcode, gene.id = rowData(matrix)$ID, gene.symbol = rowData(matrix)$Symbol, gene.type = "Gene Expression", overwrite = TRUE, type = "HDF5", genome = genome, version = chemistry_version)')
	
	os.rename(outdir + outfile, outdir + str(output_prefix) + outfile)
	
	#os.remove(tmp_genes)
	#os.remove(tmp_mtx)
	#os.remove(tmp_bar)

def marker_gene_overlap(adata, species='mm10', groupby='Cluster', method='overlap_coef', cmap='viridis'):
	print('To return overlap dataframe, run as overlap = jpplot.marker_gene_overlap(adata, ...)')
	import scanpy as sc
	import matplotlib.pyplot as pl
	import seaborn as sns
	import pandas as pd
	import numpy as np
	
	if species in ['mm10','GRCm38']:
		markers=pd.read_csv('https://www.dropbox.com/s/o30c5bfyg2p182s/Jon_classes_v2.csv?dl=1')
		marker_dict = {Class: Gene_name.tolist() for Class, Gene_name in markers.groupby("Class")["Gene_name"]}
	else:
		print('Warning: no marker list compiled yet for genome', genome)
	
	assert method in ['overlap_count', 'overlap_coef', 'jaccard'], 'method must be one of: overlap_count, overlap_coef, jaccard. Default: overlap_coef' 
	
	sc.tl.rank_genes_groups(adata, groupby=groupby)
	overlap = sc.tl.marker_gene_overlap(adata, reference_markers=marker_dict, method=method, normalize=None)
	overlap = np.round(overlap, 3)
	
	pl.figure(figsize = (6,6))
	sns.heatmap(overlap, cmap=cmap)
	
	return overlap
	
def scrublet(adata, n_pcs=30, min_gene_variability_pctl=85, min_cells=3, min_counts=2):
	import scanpy as sc
	import matplotlib.pyplot as pl
	import seaborn as sns
	import pandas as pd
	import numpy as np
	import scrublet as scr
	from scipy.sparse import issparse
	JP_key = 'stringent_doublets'
	
	if 'counts' in adata.layers.keys():
		
		print('using raw counts from adata.layers[\'counts\']')
		if issparse(adata.layers['counts']):
			counts_matrix = adata.layers['counts'].todense()
		else:
			counts_matrix = adata.layers['counts']

	else:
		print('Using adata.X')
		if issparse(adata.X):
			counts_matrix = adata.X.todense()
		else:
			counts_matrix = adata.X
	
	genes = adata.var_names

	scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.009*len(adata)/1000)
	
	doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=min_counts, min_cells=min_cells, min_gene_variability_pctl=min_gene_variability_pctl, n_prin_comps=n_pcs)
	
	if isinstance(predicted_doublets, np.ndarray):
		adata.obs['predicted_doublets'] = predicted_doublets.astype('str')
		
	else:
		adata.obs['predicted_doublets'] = 'False'
		
		
		
	adata.obs['doublet_scores'] = doublet_scores
	adata.uns['predicted_doublets_colors'] = ['#EEEEEE','#EE0000']
    
	#Stringent doublet threshold:
	predicted_doublet_fraction = 1 - (0.009 *(len(adata)/1000))
	JP_doublet_threshold = adata.obs['doublet_scores'].quantile(predicted_doublet_fraction)
	adata.obs[JP_key] = (adata.obs['doublet_scores'] > JP_doublet_threshold).astype('str')
	adata.uns[JP_key+'_colors'] = ['#EEEEEE','#EE0000']

	print('# doublets detected:', len(adata.obs[adata.obs['predicted_doublets'] == 'True']))
	print('Predicting doublet threshold based on cell #')
	print('# doublets based on total cell number and threshold:', len(adata.obs[adata.obs[JP_key] == True]))
	sns.distplot(adata.obs['doublet_scores'])
	pl.axvline(JP_doublet_threshold, linestyle=':', color='black')
	pl.show()
    
	if 'X_umap' in adata.obsm.keys():
		sc.pl.umap(adata, color=['predicted_doublets','doublet_scores',JP_key])
	
def quick_pipe(adata, genome='', n_pcs=50, min_genes=1, min_counts=1, min_cells=1, cluster_resolution=0.5, n_top_genes=4000, use_highly_variable=True, hvg_batch_key=None, exclude_genes=[], vars_to_regress='', scale=True, max_scale_value=10, use_local=False):
	import scanpy as sc
	import matplotlib.pyplot as pl
	import seaborn as sns
	import pandas as pd
	import numpy as np
	import scrublet as scr
	from scipy.sparse import issparse
	import jpplot
	
	
	assert type(scale) == bool, 'scale must be True or False'
	
	adata.var_names_make_unique()
	
	if genome == '':
		jpplot.detect_genome(adata)
		genome = adata.uns['genome']
		#print('Detected genome: ',genome)
	else:
		print('Using specified genome: ',genome)
	
	print("Preprocessing...")
	print("Using arguments:")
	print("n_pcs = ",n_pcs)
	print("min_genes = ",min_genes)
	print("min_counts = ",min_counts)
	print("cluster_resolution = ",cluster_resolution)
	print("n_top_genes = ",n_top_genes)
	print("use_highly_variable = ",use_highly_variable)
	print("scale = ",scale)
	if scale==True:
		print("Max scale value = ",max_scale_value)
	print()
	
	ribo_genes, mito_genes, S, G2M, IEG, cell_cycle = jpplot.set_genome_build(genome, use_local=use_local)
	jpplot.preprocess_annotate(adata, genome=genome, min_genes=min_genes, min_counts=min_counts, min_cells=min_cells, exclude_genes=exclude_genes, use_local=use_local)
	
	print("Scoring genes for cell cycle...")	
	sc.tl.score_genes_cell_cycle(adata, g2m_genes=G2M, s_genes=S)
	print()
	print("Computing Highly Variable Genes with n_top_genes = ",n_top_genes,"...")
	sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=hvg_batch_key)
	
	if vars_to_regress:
		scale=True
		if not isinstance(vars_to_regress, list): 
			vars_to_regress = [vars_to_regress]
		
		for j in vars_to_regress:
			sc.pp.regress_out(adata, keys=j, n_jobs=4)
			print('Scaling to',max_scale_value)
			sc.pp.scale(adata, max_value=max_scale_value)
	
	#if vars_to_regress:
	#	if vars_to_regress in adata.obs.columns:
	#		sc.pp.regress_out(adata, keys=vars_to_regress)
	#		sc.pp.scale(adata, max_value=10)
	#	else:
	#		print('Warning, you are attempting to regress out a variable not in the columns. Skipping.')
	else:
		if scale==True:
			print()
			print("Scaling .X layer to Z-scores with a max value of " + str(max_scale_value))
			sc.pp.scale(adata, max_value=max_scale_value)
		

	print()	
	print("Masking ribosomal and mitochondrial genes from highly variable...")	
	adata.var.highly_variable = adata.var.highly_variable.mask(adata.var.index.isin(ribo_genes), False)
	adata.var.highly_variable = adata.var.highly_variable.mask(adata.var.index.isin(mito_genes), False)

	print("PCA...")
	sc.pp.pca(adata, use_highly_variable=use_highly_variable)
	sc.pl.pca_variance_ratio(adata, log=True)
	
	print("Computing KNN graph with n_pcs = ",n_pcs,"...")
	sc.pp.neighbors(adata, n_pcs=n_pcs)

	#print("Computing UMAP...")
	sc.tl.umap(adata)
	
	print("Clustering by Leiden Modularity Optimization, cluster_resolution = ",cluster_resolution,"...")
	sc.tl.leiden(adata,resolution=cluster_resolution, key_added='Cluster')
	
	sc.pl.umap(adata, color='Cluster')
	
def annotate_clusters(adata, genome='', input_key='Cluster', output_key='Cell_Type', custom_list=None, top_n_markers=50, boost_score=None, cell_type_key='Class'):
	import scanpy as sc
	import matplotlib.pyplot as pl
	import seaborn as sns
	import pandas as pd
	import numpy as np
	import jpplot
	
	if boost_score:
		if isinstance(boost_score, list):
			for el in boost_score:
				assert isinstance(el[0], (str)), 'boost_score must be a tuple of the form (\'Pericyte\',1.5)'
				assert isinstance(el[1], (int,float)), 'boost_score must be a tuple of the form (\'Pericyte\',1.5)'
		else:
			assert isinstance(boost_score[0], (str)), 'boost_score must be a tuple of the form (\'Pericyte\',1.5)'
			assert isinstance(boost_score[1], (int,float)), 'boost_score must be a tuple of the form (\'Pericyte\',1.5)'
	
	if adata.uns['genome']:
		print('Using genome detected in adata.uns:',adata.uns['genome'])
		genome = adata.uns['genome']
	else:
  	  print('No genome detected in adata.uns')
  	  jpplot.detect_genome(adata, use_local=use_local)
  	
  	## Read in table of markers  
	if custom_list:
		markers=pd.read_csv(custom_list)
	
	else:
	##Specify Genome and read in marker table from Dropbox:
		assert genome in ['mm10', 'hg38','GRCh38','GRCm38'], 'organism must be one of: mm10, hg38, GRCh38, GRCm38'
		if genome in ['mm10','GRCm38']:
			markers=pd.read_csv('https://www.dropbox.com/s/o30c5bfyg2p182s/Jon_classes_v2.csv?dl=1')
		elif genome in ['hg38','GRCh38']:
			markers=pd.read_csv('https://www.dropbox.com/s/geq43fu77osxilv/Jon_classes_human.csv?dl=1')
		else:
			print('Warning: unrecognized genome specified')

	markers.set_index('Gene_name', inplace=True)
	markers = markers.loc[[name for name in markers.index if name in adata.var_names]]
	markers['Gene_name'] = markers.index
	marker_dict = {Class: Gene_name.tolist() for Class, Gene_name in markers.groupby(cell_type_key)["Gene_name"]}

	## Compute cluster marker genes, quantify overlap:
	methods = ['overlap_count', 'overlap_coef', 'jaccard']
	normalize = ['reference', 'data', None]

	if 'rank_genes_groups' in adata.uns.keys():
		if adata.uns['rank_genes_groups']['params']['groupby'] != input_key:
			sc.tl.rank_genes_groups(adata, groupby=input_key)
	else:
		print('No differential expression data detected. Running sc.tl.rank_genes_groups with groupby=' + str(input_key))
		sc.tl.rank_genes_groups(adata, groupby=input_key)
	
	#Plot a heatmap of the overlap to assess the extent of clustering
	
	pl.figure(figsize = (10,8))

	df = sc.tl.marker_gene_overlap(adata, reference_markers=marker_dict, key='rank_genes_groups', method=methods[0], normalize=normalize[1], top_n_markers=top_n_markers)
	
	if boost_score:
		def boost(gene, factor):
			if gene in df.index:
				print('Scaling ' + str(gene) + ' by factor: ' + str(factor))
				df.loc[gene] = df.loc[gene] * factor
			else:
				print('Warning, score ' + str(gene) + ' not detected in dataframe.')
		
		if isinstance(boost_score, list):
			for el in boost_score:
				boost(el[0], el[1])
				
		else:
			boost(boost_score[0], boost_score[1])
	
	sns.heatmap(df)
	
	output_key_dict = dict(zip(df.columns.tolist(), df.idxmax().tolist()))
	adata.obs[output_key] = adata.obs[input_key].map(output_key_dict).astype('category')
	sc.pl.umap(adata, color=[input_key,output_key])
		
def nmf(adata, n_components=50, random_state=0, layer='raw'):
	import scanpy as sc
	import matplotlib.pyplot as pl
	import seaborn as sns
	import pandas as pd
	import numpy as np
	import jpplot
	
	
	from sklearn.decomposition import non_negative_factorization
	
	if layer == 'raw':
		data_values  = adata.raw[:,adata.var_names].X
	elif layer == 'X':
		data_values  = adata.X
	elif layer in adata.layers.keys():
		data_values  = adata.layers[layer].X
			
	
	W, H, n_iter = non_negative_factorization(data_values, n_components=n_components, init='random', random_state=random_state)

	adata.obsm['X_nmf'] = W
	adata.varm['NMF'] = H.T
	
def ica(adata, n_components=50, random_state=42, max_matrix_size=1e9, use_highly_variable=True):
	import sklearn as skl
	import time
	import numpy as np
	import pandas as pd
	
	start_time = time.time()
	
	def run_ica(X):
		ica = skl.decomposition.FastICA(n_components=n_components, random_state=random_state)
		print('Performing ICA reduction using n_components = ',n_components)
		print('using random_state = ',random_state)
		
		ica.fit(X)
		ica_dims = ica.fit_transform(X)
		return ica, ica_dims
	
	if adata.shape[0]*adata.shape[1] > max_matrix_size:
		print('Matrix is larger than',max_matrix_size,'aborting...')
	else:
		## Define the array on which to perform ICA
		if use_highly_variable:
			print('Using highly variable genes in adata.var')
			mask = adata.var['highly_variable']
			X = adata[:,mask].X.toarray()
		else:
			print('Using full data matrix adata.X')
			X = adata.X.toarray()
		
		## Run ICA
		print('Running ICA on matrix of dimensions:', X.shape)
		ica, ica_dims = run_ica(X)
		
		## Write the ICA loadings back into adata.varm
		print('Saving ICA loadings to adata.varm[\'ICA_components\']')
		if use_highly_variable:
			not_hvg = adata.var_names[~mask]
			hvg = adata.var_names[mask]
			ica_dims_df = pd.DataFrame(ica.components_.T)
			ica_dims_df.index = hvg
			nvg_df = pd.DataFrame(np.zeros([len(not_hvg),n_components]))
			nvg_df.index = not_hvg
			ica_dims_df = pd.concat([nvg_df,ica_dims_df])
			adata.varm['ICA_components'] = ica_dims_df.values
			
		else:
			adata.varm['ICA_components'] = ica.components_.T
		
		## Write the ICA dimensional reduction  into adata.obsm
		print('Saving ICA representation to adata.obsm[\'X_ica\']')
		adata.obsm['X_ica'] = ica_dims
		
		## Time it
		end_time = time.time()
		elapsed_time = end_time - start_time
		print('Elapsed time:', np.round(elapsed_time,1),'seconds')

def SoupX_quick(CellRanger_path, output_suffix='_SoupX_Strained', genome='mm10', chemistry_version='3'):
	# must be a valid cellranger 3.0+ path
	import os
	import warnings
	warnings.filterwarnings('ignore')
	import rpy2.robjects as ro

	assert genome in ['mm10', 'hg38','GRCh38','GRCm38'], 'organism must be one of: mm10, hg38, GRCh38, GRCm38.  Default: mm10'
	outfile = os.path.join(CellRanger_path,'SoupX_strained.h5')

	ro.globalenv['CellRanger_path'] = CellRanger_path
	ro.globalenv['genome'] = genome
	ro.globalenv['chemistry_version'] = chemistry_version
	ro.globalenv['outfile'] = outfile
	ro.globalenv['genome'] = genome

	ro.r('library(SoupX)')
	ro.r('library(DropletUtils)')
	ro.r('library(Matrix)')
	
	ro.r('sc = load10X(CellRanger_path, keepDroplets=TRUE)')
	ro.r('sc = autoEstCont(sc)')
	ro.r('out = adjustCounts(sc)')

	ro.r('genes <-read.table(file.path(CellRanger_path,"filtered_feature_bc_matrix/features.tsv.gz"),sep = "\t")')

	ro.r('DropletUtils:::write10xCounts(outfile, out, \
	gene.id = genes$V1, \
	gene.symbol = genes$V2, \
	gene.type = "Gene Expression", \
	overwrite = TRUE, \
	type = "HDF5", \
	genome = genome, \
	version = "3")')


def filter_cells_by_marker_gene(adata, gene, apply_filter=False, count_threshold=3):
	import seaborn as sns
	import matplotlib.pyplot as pl
	import scanpy as sc
	import numpy as np
	from scipy.sparse import issparse

	if issparse(adata.X):

		sns.distplot(adata[:,gene].X.todense())
		pl.yscale('log')
		sc.pl.umap(adata, color=gene)

	else:
		sns.distplot(adata[:,gene].X)
		pl.yscale('log')
		sc.pl.umap(adata, color=gene)    
	
	if issparse(adata.X):
		adata.obs[gene + '_filter'] = adata[:, gene].X.todense() > count_threshold
	else:
		adata.obs[gene + '_filter'] = adata[:, gene].X > count_threshold
		
	adata.obs[gene + '_filter'] = adata.obs[gene + '_filter'].astype('str')
	adata.uns[gene + '_filter_colors'] = ['#DDDDDD','#EE3300']
	sc.pl.umap(adata, color=gene + '_filter')
	print(str(len(adata[adata.obs[gene + '_filter'] == 'True'])) + ' cells will be filtered out')
    
	if apply_filter:
		print('Applying filter...')
		adata = adata[adata.obs[gene + '_filter'] == 'False']
		return adata

def import_cellranger(cellranger_outs_folder):
	print("Usage: adata = import_cellranger('/path/to/cellranger/sample/outs/')")
	import numpy as np
	import scanpy as sc
	import pandas as pd
	import os
	import matplotlib.pyplot as pl
	import bbknn
	import jpplot 
	import seaborn as sns
	adata = sc.read_10x_h5(cellranger_outs_folder + 'filtered_feature_bc_matrix.h5')
	cellranger_outs_folder + '/analysis/tsne/2_components/projection.csv'
	cellranger_outs_folder + '/analysis/tsne/2_components/projection.csv'
	adata.obsm['X_tsne'] = np.array(pd.read_csv(cellranger_outs_folder + '/analysis/tsne/2_components/projection.csv', index_col=0))
	adata.obsm['X_umap'] = np.array(pd.read_csv(cellranger_outs_folder + '/analysis/umap/2_components/projection.csv', index_col=0))
	
	adata.obsm['X_tsne_cellranger'] = adata.obsm['X_tsne'].copy()
	adata.obsm['X_umap_cellranger'] = adata.obsm['X_umap'].copy()
	for i in os.listdir(cellranger_outs_folder + '/analysis/clustering/'):
		adata.obs[i] = np.array(pd.read_csv(cellranger_outs_folder + '/analysis/clustering/' + str(i) +'/clusters.csv', index_col=0).astype('str'))
	return adata
	
def redo_pca_cluster(adata, genome='hg38', batch_key=None, key_added='Subset_Cluster', n_pcs=50, min_genes=1, min_counts=1, min_cells=4, cluster_resolution=0.5, n_top_genes=4000, use_highly_variable=True, exclude_highly_expressed=True, exclude_genes=[], scale=False, max_scale_value=10, vars_to_regress=None, use_local=False):
	import scanpy as sc
	import matplotlib.pyplot as pl
	import seaborn as sns
	import pandas as pd
	import numpy as np
	#import scrublet as scr
	#from scipy.sparse import issparse
	import jpplot
	
	##Specify Genome and read in marker table from Dropbox:
	if 'genome' in adata.uns.keys():
		genome = adata.uns['genome']
	else:
		jpplot.detect_genome(adata)
		genome = adata.uns['genome']
	assert genome in ['mm10', 'hg38','GRCh38','GRCm38'], 'organism must be one of: mm10, hg38, GRCh38, GRCm38'
	ribo_genes, mito_genes, S, G2M, IEG, cell_cycle = jpplot.set_genome_build(genome, use_local=use_local)
	
	adata.var_names_make_unique()
	
	print("Reclustering...")
	print("Using arguments:")
	print("n_pcs = ",n_pcs)
	#print("min_genes = ",min_genes)
	#print("min_counts = ",min_counts)
	print("min_cells = ",min_cells)
	print("cluster_resolution = ",cluster_resolution)
	print("n_top_genes = ",n_top_genes)
	print("use_highly_variable = ",use_highly_variable)
	
	if 'counts' in adata.layers.keys():
		print('Restoring raw counts from adata.layers[\'counts\']')
		adata.X = adata.layers['counts']
		sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=exclude_highly_expressed, key_added='Norm_Factor')
		print ('Log transforming...')
		sc.pp.log1p(adata)
		print ('Saving normalized log-transformed data to .raw')
		adata.raw = adata.copy()
		print('Done.')
		
	elif 'raw' in dir(adata):
		print('Restoring raw data')
		adata.X =  adata.raw[:,adata.var_names].X.copy()
	
	print("Filtering out genes expressed in fewer than ",min_cells," cells...")
	sc.pp.filter_genes(adata, min_cells=min_cells)	
	
	print("Computing Highly Variable Genes with n_top_genes = ",n_top_genes,"...")
	sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_key)
	
	print("Masking ribosomal and mitochondrial genes from highly variable...")	
	adata.var.highly_variable = adata.var.highly_variable.mask(adata.var.index.isin(ribo_genes), False)
	adata.var.highly_variable = adata.var.highly_variable.mask(adata.var.index.isin(mito_genes), False)
	
	print("Scoring genes for cell cycle...")	
	sc.tl.score_genes_cell_cycle(adata, g2m_genes=G2M, s_genes=S)
	
	if vars_to_regress:
		scale=True
		if not isinstance(vars_to_regress, list): 
			vars_to_regress = [vars_to_regress]
		
		for j in vars_to_regress:
			sc.pp.regress_out(adata, keys=j, n_jobs=4)
			print('Scaling to',max_scale_value)
			sc.pp.scale(adata, max_value=max_scale_value)
	
	#if vars_to_regress:
	#	if vars_to_regress in adata.obs.columns:
	#		sc.pp.regress_out(adata, keys=vars_to_regress)
	#		sc.pp.scale(adata, max_value=10)
	#	else:
	#		print('Warning, you are attempting to regress out a variable not in the columns. Skipping.')
	else:
		if scale==True:
			print()
			print("Scaling .X layer to Z-scores with a max value of " + str(max_scale_value))
			sc.pp.scale(adata, max_value=max_scale_value)

	print("PCA...")
	sc.pp.pca(adata, use_highly_variable=use_highly_variable)
	
	print("Computing KNN graph with n_pcs = ",n_pcs,"...")
	sc.pp.neighbors(adata, n_pcs=n_pcs)

	#print("Computing UMAP...")
	sc.tl.umap(adata)
	
	print("Clustering by Leiden Modularity Optimization, cluster_resolution = ",cluster_resolution,"...", "new clusters will be named:", key_added)
	sc.tl.leiden(adata,resolution=cluster_resolution, key_added=key_added)
	
	sc.pl.umap(adata, color=key_added)

def import_cellranger_aggrcsv(adata, aggr_csv='aggregation.csv', key_added='Sample'):
	import pandas as pd
	import numpy as np

	sample_numbers = adata.obs.index.str.replace('[ATGC]*-','')
	num_samps_in_matrix = len(sample_numbers.unique())

	aggr = pd.read_csv(aggr_csv)
	sample_id_column = aggr.columns[0]
	aggr['Sample_Num'] = aggr.index + 1
	aggr['Sample_Num'] = aggr['Sample_Num'].astype(str)
	num_samps_in_aggrcsv = len(aggr['Sample_Num'].unique())
    
	if num_samps_in_matrix == num_samps_in_aggrcsv:
    
		sampdict = dict(zip(aggr['Sample_Num'],aggr[sample_id_column]))

		adata.obs[key_added] = sample_numbers
		adata.obs[key_added] = adata.obs[key_added].astype(str)
		adata.obs = adata.obs.replace({key_added:sampdict})
		adata.obs[key_added].value_counts()
        
		print('Detected ' + str(num_samps_in_matrix) + ' GEM groups. Adding key "' + str(key_added) + '" to obs dataframe.')
        
	else: 
		print('Warning: Number of samples detected in your data matrix differs from the number in the provided aggregation.csv file')

def detect_genome(adata, genome=''):
    import pandas as pd
    if genome == '':
        if 'genome' in adata.uns.keys():
            print('Using genome detected in .uns:',adata.uns['genome'])
        elif 'genome' in adata.var.columns:
            if len(adata.var['genome'] == 1):
                genome = adata.var['genome'].unique()[0]
                adata.uns['genome'] = genome
                print('Adding',genome,'to adata.uns[\'genome\']')
            elif len(adata.var['genome'] == 1):
                print('Detected more than one genome in \'genome\' column of adata.var')
                genomes = "_".join(adata.var['genome'].unique())
                print('Adding',genomes,'to adata.uns[\'genome\']')
                adata.uns['genome'] = genomes
        else:
            print('No genome detected in adata.uns')
            var_names = pd.Series(adata.var_names.tolist() + ['FAKEUPPER','fakelower'])
            count_upper_lower = var_names.str.isupper()
            num_upper = count_upper_lower.value_counts()[True]
            num_lower = count_upper_lower.value_counts()[False]
            
            if num_upper > num_lower:
                print('This appears to be a human dataset, using genome hg38')
                genome = 'hg38'
                adata.uns['genome'] = genome
            else:
                print('This appears to be a mouse dataset, using genome mm10')
                genome = 'mm10'
                adata.uns['genome'] = genome
    else:
        print('No genome specified, defaulting to hg38')
        genome = 'hg38'
        adata.uns['genome'] = genome
        
def harmony(adata, batch_key='Sample', theta=4, max_iter_harmony=15, random_state=42, use_gpu=False, cluster_and_umap=False, key_added='Cluster_harmony', n_pcs=None):
    #requires pip install harmony-pytorch
    from harmony import harmonize
    import scanpy as sc
    import pandas as pd
    
    if n_pcs:
        X = adata.obsm['X_pca'][:,:n_pcs]
    else: 
        X = adata.obsm['X_pca']
    print('Computing harmonized PCA using ' + str(X.shape[1]) + ' PCs...')
    Z = harmonize(X, adata.obs, batch_key = batch_key, max_iter_harmony=max_iter_harmony, random_state=random_state, use_gpu=use_gpu)
    adata.obsm['X_harmony'] = Z
    print('Storing harmonized PCA in adata.obsm[\'X_harmony\']')
    
    if cluster_and_umap:
        sc.pp.neighbors(adata, use_rep='X_harmony')
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=0.5, key_added=key_added)
        sc.pl.umap(adata, color=key_added)

#def harmony(adata, batch_key='Sample', theta=4):
#	import os
#	import pandas as pd
#	import warnings
#	import numpy as np
	
#	warnings.filterwarnings('ignore')

#	import rpy2.robjects as ro
#	from rpy2.robjects import numpy2ri
#	ro.conversion.py2ri = numpy2ri

#	ro.globalenv['pca'] = adata.obsm['X_pca']
#	ro.globalenv['hem'] = hem
#	ro.globalenv['theta'] = theta
	
#	ro.r('library(harmony)')
#	ro.r('library(magrittr)')
	
#	hem <- HarmonyMatrix(pca, batch, theta=4)
#	hem = data.frame(hem)
	
#	adata.obsm['X_pca_harmony'] = hem.values

def random_color_palette(n_colors=8):
	import random
	color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]) for i in range(n_colors)]
	return(color)
	
def write_mtx(adata, use_raw_counts=True, destination_folder='', h5=False, output_prefix='', detect_genome=True):
	import os
	import pandas as pd
	import warnings
	import numpy as np
	import scanpy as sc
	import gzip
	from scipy.io import mmread, mmwrite
	import shutil
	import jpplot
	
	########  Assign which data layer to use
	if use_raw_counts:
		print('Exporting raw digital counts layer')
		data=adata.copy()
		#matrix=data.layers['counts'].T

	else:
		print('Exporting .raw layer.  This is log-normalized.  If running Cellranger Reanalyze, make sure to exclude a normalization step.')
		data = adata.raw.copy()
		#matrix = data.X.T
		
	##Detect Genome
	if detect_genome:
		jpplot.detect_genome(adata)
	genome = adata.uns['genome']
	
	#If necessary, import ensembl ids:
	if genome in ['hg38','GRCh38']:
		vardf = pd.read_csv('/Users/jonpreall/Dropbox/Preall_Lab/Preall/Gene_Lists/GRCh38_genes_CR31.csv', index_col=0)
		
	elif genome in ['hg38','GRCm38']:
		vardf = pd.read_csv('/Users/jonpreall/Dropbox/Preall_Lab/Preall/Gene_Lists/mm10_genes_CR31.csv', index_col=0)
	
	vars_in_data = list(set(data.var.index).intersection(set(vardf.index)))
	vardf = vardf.loc[vars_in_data]
	
	## Create Paths
	if os.path.exists(destination_folder):
		destination_folder = (os.path.join(destination_folder, output_prefix + '/gene_bc_matrix/'))
		if os.path.exists(destination_folder):
			print('Removing old output folder')
			shutil.rmtree(destination_folder)

		os.makedirs(destination_folder, exist_ok=True)
		print('Saving to ' + destination_folder)
	else:
		destination_folder = (os.path.join(os.getcwd(), output_prefix + '/gene_bc_matrix/'))
		os.makedirs(destination_folder, exist_ok=True)
		print('Cannot locate destination folder, saving to ' + destination_folder)
		

	
	######## Write metadata.tsv.gz
	obs = adata.obs
	
	## add TSNE, UMAP, DIFFMAP, FA projections to metadata dataframe
	
	dimreds = [name for name in adata.obsm.keys() if name.startswith('X_')]

	possible_dimreds = ['umap','tsne','diffmap','fa']
	search_dimreds = possible_dimreds + [name.upper() for name in possible_dimreds] + [name.capitalize() for name in possible_dimreds]

	found_dimreds = []
	for red in search_dimreds:
	    found_dimreds += [name for name in adata.obsm.keys() if red in name]

	df_dims = pd.DataFrame()
	for red in found_dimreds:
		dim_name = red.replace('X_','')
		df = pd.DataFrame(adata.obsm[red])
		df.columns = [str(dim_name) + "_" + str(name) for name in df.columns]
		df.index = df.index.rename('Barcode')
		
		#Write projection to its own .csv
		writedf = df.copy()
		writedf.index = obs.index
		writedf.to_csv(os.path.join(destination_folder, dim_name + ".csv"), index_label='Barcode',sep = ",", index = True)

		df_dims = pd.concat([df_dims,df], axis=1)
		
	df_dims.index = obs.index

	obs = obs.merge(df_dims, left_index=True, right_index=True, how='left')
	obs.to_csv(os.path.join(destination_folder, "metadata.tsv.gz"), sep = "\t", index = True, compression='gzip')

	#pd.DataFrame(adata.var.index).to_csv(os.path.join(destination_folder, "genes.tsv" ),   sep = "\t", index = False)
	
	########  Write features.tsv.gz
	#vars_in_data = data.var.index.tolist()
	#else:
	#import re
	#dfstr = adata.var.select_dtypes(include=['O'])
	#matching_columns = [x for x in dfstr.columns.tolist() if re.match('^ens', x, re.IGNORECASE)]
	#ensembl_column = ''
	#if matching_columns == []:
	#	print('Warning, no EnsemblID column detected')
	#else:
	#	ensembl_column = matching_columns[0]
	
	#dfstr = data.var.select_dtypes(include=['O'])
	ensembl_column = vardf.columns[vardf.apply(lambda x: x.str.match('ENS')).all()][0]
	
	vardf['Gene_name'] = vardf.index
	vardf = vardf.loc[:,[ensembl_column,'Gene_name']]
	vardf['Assay'] = 'Gene Expression'
	vardf.to_csv(os.path.join(destination_folder, "features.tsv.gz" ),   sep = "\t", header=None, index=None, compression='gzip')
	
	######## Write barcodes.tsv.gz
	pd.DataFrame(obs.index).to_csv(os.path.join(destination_folder, "barcodes.tsv.gz"), header=None, sep = "\t", index = False, compression='gzip')
	
	########  Write mtx.tsv
	if use_raw_counts:
		#print('Exporting raw digital counts layer')
		#data=adata[:,vars_in_data]
		matrix=adata[:,vars_in_data].layers['counts'].T

	else:
		#print('Exporting .raw layer.  This is log-normalized.  If running Cellranger Reanalyze, make sure to exclude a normalization step.')
		#data = 
		matrix = adata.raw[:,vars_in_data].X.T
			
	mmwrite(os.path.join(destination_folder, "matrix.mtx"), matrix)
		
	## Compress any uncompressed matrix file:
	if not os.path.exists(destination_folder + "matrix.mtx.gz"):
		if os.path.exists(destination_folder + "matrix.mtx"):
			#print('Detected matrix.mtx, converting to v3 matrix.mtx.gz format...')
		
			with open(destination_folder + "matrix.mtx", 'rb') as orig_file:
				with gzip.open(destination_folder + "matrix.mtx.gz", 'wb') as zipped_file:
					zipped_file.writelines(orig_file)
			os.remove(destination_folder + "matrix.mtx") 
	
	if h5:
		jpplot.mex_to_h5(matrixdir=destination_folder, output_prefix=output_prefix, genome=genome)
		## from jpplot.mex_to_h5
		#matrixdir = destination_folder
		#outdir = os.path.dirname(matrixdir) +'/'
		
		#if output_prefix:
		#	output_prefix = output_prefix + '_'
		
		#if 'raw' in os.path.basename(matrixdir.split('/')[-2]):
		#	outfile = 'raw_feature_bc_matrices_h5.h5'
		#elif 'filtered' in os.path.basename(matrixdir.split('/')[-2]):
		#	outfile = 'filtered_feature_bc_matrices_h5.h5'
		#else:
		#	outfile = 'gene_bc_matrix.h5'
			
		#h5file = outdir + str(output_prefix) + outfile
		#print('Writing Cellranger compatible .h5 file to ', h5file)
		
		#print("ATTEMPTING TO WRITE TO:", h5file)
		
		
		#import h5py
		#strList=['genome']
		#asciiList = [n.encode("ascii", "ignore") for n in strList]

		#openhdf5 = h5py.File(h5file,'a')
		#data = openhdf5['matrix/data']
		#X1 = data.value.astype('int32')
		#del openhdf5['matrix/data']
		#openhdf5.create_dataset('matrix/data', data=X1, compression="gzip", compression_opts=4) #, chunks=(80000,))	
		#openhdf5['matrix']['features'].create_dataset('_all_tag_keys', (1,),'S6', asciiList)
		#openhdf5.close()
	
		
		
def export_coordinates(adata, destination_folder='', output_prefix=''):
	import os
	import pandas as pd
	import warnings
	import numpy as np
	import gzip
	import shutil
	import jpplot
	
	if os.path.exists(destination_folder):
		destination_folder = os.path.join(destination_folder, output_prefix)
		if os.path.exists(destination_folder):
			print('Removing old output folder')
			shutil.rmtree(destination_folder)

		os.makedirs(destination_folder, exist_ok=True)
		print('Saving to ' + destination_folder)
	else:
		destination_folder = os.path.join(os.getcwd(), output_prefix)
		os.makedirs(destination_folder, exist_ok=True)
		print('Cannot locate destination folder, saving to ' + destination_folder)
	
	obs = adata.obs
	
	## find any low-dimensionality projections that are used for visualization (eg, TSNE, UMAP, FA, etc)
	
	for dimred in adata.obsm.keys():
		n_dims = adata.obsm[dimred].shape[1]
    	
		if n_dims < 4:
			dimred_name = dimred.strip('X_').upper()
			df = pd.DataFrame(adata.obsm[dimred], index=adata.obs_names, columns = [dimred_name + str(n) for n in range(0,n_dims)])

			df.to_csv(os.path.join(destination_folder, dimred_name + ".csv"), index_label='Barcode',sep = ",", index = True)
		
def write_cellranger_h5(adata, output_folder, file_prefix, layer='counts', aggr_key=None, genome=None, feature_types='feature_types', force_library_id=None):
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import h5py
    import sys, os
    import scipy
    from scipy.sparse import csc_matrix
    from scipy.sparse import csr_matrix
    import h5sparse
    import re
    
    if layer in adata.layers.keys():
        MATRIX=csr_matrix(adata.layers[layer])
    elif layer == 'X':
    	MATRIX=csr_matrix(adata.X)
        #MATRIX=csc_matrix(adata.X)
    else:
        print("Error: specified matrix layer can't be found in your data")
        
    ## identify, if it exists, the column in .var containing ensembl gene ids:
    import re
    dfstr = adata.var.select_dtypes(include=['O'])
    ens_columns = []
    for col in dfstr.columns:
        mostly_ENS_names = len([x for x in dfstr[col] if re.match('^ens', x, re.IGNORECASE)]) / len(dfstr) > 0.9
        if mostly_ENS_names:
            ens_columns += [col]
    if len(ens_columns) > 1:
        print('Warning, more than one potential ENSEMBL name column detected.  Defaulting to',ens_columns[0])
        ensembl_column = ens_columns[0]
        gene_ids = adata.var[ensembl_column]
         
    elif len(ens_columns) == 1:
        print('Using column',ens_columns[0],'as ENSEMBL ids')
        ensembl_column = ens_columns[0]
        gene_ids = adata.var[ensembl_column]
         
    else:
        print('Warning, no ENSEMBL gene id column detected.  Duplicating contents of adata.var.index')
        gene_ids = np.array(adata.var.index)

    	
    ## if not already present, add 'feature_types' column to .var
    if 'feature_types' not in adata.var.columns:
        adata.var['feature_types'] = feature_types
    FEATURE_TYPE = np.array(adata.var['feature_types']).astype('S')
        
    BCS = np.array(adata.obs.index).astype('S')
    
    FEATURES = np.array(adata.var.index).astype('S')
    
    ## Detect chieck genome
    
    if genome:
        GENOME=genome
    else:
        if 'genome' in adata.uns.keys():
            GENOME=adata.uns['genome']
        else:
            if detect_genome:
                detect_genome(adata)
                GENOME=adata.uns['genome']
            else:
                print('No genome specified, writing attribute as unspecified_genome')
                GENOME='unspecified_genome'
        
    FEATURE_IDS = np.array(gene_ids).astype('S') 
       
    all_tag_keys = np.array([b'genome'])
    
    ## Generate list of  Library_IDs based on barcode suffixes:
    if force_library_id == None:
        sampnames = np.array(['Sample_' + suffix for suffix in adata.obs.index.str.replace('[ATGC]+-','').unique()]).astype('S')      
        LIBRARY_IDS = sampnames
    else:
        if type(force_library_id) == str:
            LIBRARY_IDS = np.array([force_library_id]).astype('S')   
        elif type(force_library_id) == list:
            LIBRARY_IDS = np.array(force_library_id).astype('S')   
        else:
            print('Warning, force_library_id must be a str or list.')
        
    # LIBRARY_IDS = np.array([sample.encode() for sample in sampnames])
    
    ORIG_GEM_GROUPS = np.array([1])

    SHAPE=adata.T.shape
    #SHAPE=adata.shape
    
    outfile=os.path.join(output_folder,file_prefix +'.h5')
    
    if aggr_key:
        ## Check if aggr_key is a categorical dtype:
        categorical_obs_cols = adata.obs.select_dtypes(include=['category']).keys()
        assert aggr_key in categorical_obs_cols, 'aggr_key must be a key corresponding to a categorical observation in .obs dataframe. Eg. \'Sample\''
        if aggr_key in categorical_obs_cols:
            obs = pd.DataFrame(adata.obs.loc[:,aggr_key])
        else:
            print('Try assigning',aggr_key,'to a categorical data type.')
            
    ## Test if barcodes are 10X Genomics compatible:
        BARCODES_IN_10X_FORMAT = len([x for x in obs.index.tolist() if re.match('[ATGCN]{16}-[0-9]+$', x, re.IGNORECASE)]) == len(adata.obs)
        
        ## Export an aggr.csv file to serve as a key for sample names:
        if BARCODES_IN_10X_FORMAT:
            obs = pd.DataFrame(adata.obs.loc[:,aggr_key])

        else:
            print("Warning, barcodes not in 10X format.  Attempting to convert.")

    	## Convert barcodes to a cellranger compatible format:
        sampnames = obs[aggr_key].unique().tolist()
        samp_num_dict = {sample:sampnames.index(sample) + 1 for sample in sampnames}
        
        ## Extract the 16 base cell barcode
        import re
        current_barcodes = obs.index.tolist()
        cell_barcode_seq = [re.search('[ATGCN]{16}', name, re.IGNORECASE).group(0).upper() for name in current_barcodes] 
        obs['cell_barcode'] = cell_barcode_seq
        obs['aggr_index'] =  obs[aggr_key].copy()
        obs = obs.replace({'aggr_index':samp_num_dict})
        obs['aggr_index'] = obs['cell_barcode'] + '-' + obs['aggr_index'].astype('str')
        obs.set_index('aggr_index', inplace=True)
        obs = pd.DataFrame(obs.loc[:,aggr_key])
        BCS = np.array(obs.index).astype('S')
        #### BELOW IS WRONG
        #LIBRARY_IDS = np.array([sample.encode() for sample in sampnames])
        #ORIG_GEM_GROUPS = np.array([1]).repeat(len(sampnames))
        
        def write_aggr_key(obs):
            aggr_filename = os.path.join(output_folder,file_prefix+'_aggr.csv')
            print('Writing aggregation key file:',aggr_filename,'...')
            barcode_suffix = obs.index.str.replace('[ATGC]+-','')
            sampnames = obs[aggr_key]
            sample_dict = sorted(dict(zip(barcode_suffix,sampnames)).items())
            aggr_csv = pd.DataFrame(sample_dict)
            aggr_csv.columns = ['Barcode_Suffix','library_id']
            aggr_csv.to_csv(aggr_filename, index=None)
            print('Done')
        

        
        write_aggr_key(obs)
    print('Writing to',outfile)
    
    with h5sparse.File(outfile, 'w') as h5f:
        h5f.create_dataset('matrix/', data=MATRIX, compression="gzip")
        h5f.close()
    
    with h5py.File(outfile, 'r+') as f:
        f.create_dataset('matrix/barcodes', data=BCS)
        f.create_dataset('matrix/shape', (2,),dtype='int32', data=SHAPE)
        features = f.create_group('matrix/features')
        features.create_dataset('_all_tag_keys', (1,),'S6', data=all_tag_keys)  
        features.create_dataset('feature_type', data=np.array([b'Gene Expression'] * SHAPE[0]), compression="gzip")
        features.create_dataset('genome', data=np.array([GENOME.encode()] * SHAPE[0]), compression="gzip")
        features.create_dataset('id', data=FEATURE_IDS, compression="gzip")
        features.create_dataset('name', data=FEATURES, compression="gzip")
        
        f.attrs['chemistry_description'] = b'Single Cell 3\' v3'
        f.attrs['filetype'] = 'matrix'
        f.attrs['library_ids'] = LIBRARY_IDS
        f.attrs['original_gem_groups'] = ORIG_GEM_GROUPS
        f.attrs['version'] = 2
        f.close()
        
def scrublet_aggregated(adata):
    import pandas as pd
    import jpplot
    import scanpy as sc
    JP_key = 'stringent_doublets'
    
    def predict_doublet_fraction(adata):
    	predicted_doublet_fraction = 1 - (0.009 *(len(adata)/1000))
    	return predicted_doublet_fraction
    
    adata.obs['Barcode_Suffix'] = adata.obs.index.str.replace('[ATGC]+-','')
    predicted_doublets=pd.DataFrame()
    for i in adata.obs['Barcode_Suffix'].unique():
        print()
        print('Slicing out sample ' + str(i) + ' and predicting doublets...')
        tmp = adata[adata.obs['Barcode_Suffix'] == i]
        predicted_doublet_fraction = predict_doublet_fraction(tmp)
        jpplot.scrublet(tmp)
        JP_doublet_threshold = tmp.obs['doublet_scores'].quantile(predicted_doublet_fraction)
        tmp.obs[JP_key] = tmp.obs['doublet_scores'] > JP_doublet_threshold
        predicted_doublets = pd.concat([predicted_doublets,tmp.obs.loc[:,['predicted_doublets','doublet_scores',JP_key]]], axis=0)
        predicted_doublets['predicted_doublets'] = predicted_doublets['predicted_doublets'].astype('str')
        #
        predicted_doublets[JP_key] = predicted_doublets[JP_key].astype('str')
        print('# doublets detected:', len(predicted_doublets[predicted_doublets['predicted_doublets'] == 'True']))
        print('Predicting doublet threshold based on cell #')
        print('# doublets based on total cell number and threshold:', len(tmp.obs[tmp.obs[JP_key] == True]))
    adata.obs = adata.obs.merge(predicted_doublets, left_index=True, right_index=True, how='left')
    adata.uns['predicted_doublets_colors'] = ['#EEEEEE','#EE0000']
    adata.uns[JP_key+'_colors'] = ['#EEEEEE','#EE0000']
    #return predicted_doublets
    
    if 'X_umap' in adata.obsm.keys():
        sc.pl.umap(adata, color=['predicted_doublets',JP_key])


def gini(array):
    import pandas as pd
    import numpy as np
    
    """Calculate the Gini coefficient of a numpy array."""
    # based on bottom eq: http://www.statsdirect.com/help/content/image/stat0206_wmf.gif
    # from: http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    array_type = type(array)
    
    assert array_type in [pd.core.frame.DataFrame, pd.core.series.Series, np.ndarray], 'Input must be a Pandas dataframe or a numpy array'
    
    if type(array) is pd.core.frame.DataFrame:
        array = np.array(array)
    
    if type(array) is pd.core.series.Series:
        array = np.array(array)
    
    array = array.flatten() #all values are treated equally, arrays must be 1d
    if np.amin(array) < 0:
        array -= np.amin(array) #values cannot be negative
    array += 0.0000001 #values cannot be 0
    array = np.sort(array) #values must be sorted
    index = np.arange(1,array.shape[0]+1) #index per array element
    n = array.shape[0]#number of array elements
    return ((np.sum((2 * index - n  - 1) * array)) / (n * np.sum(array))) #Gini coefficient
    
def gini_all_genes(adata, layer='counts'):
    import jpplot
    import numpy as np
    import pandas as pd
    import scanpy as sc
    
    assert layer in adata.layers.keys(), 'layer must be in adata.layers'
    tmpdf = pd.DataFrame(adata.layers[layer].todense()).astype('float')
    tmpdf.columns = adata.var_names
    tmpdf.index = adata.obs.index
    adata.var['gini_index'] = tmpdf.apply(jpplot.gini, axis=0)
    del(tmpdf)

def set_computer(force=None):
    import os
    import socket
    print('Usage: dropbox_dir = jpplot.set_computer()')

    computername = socket.gethostname()
    dropbox_loc_dict = {
        'Ebelskiver.local':'/Users/jpreall/Dropbox/',
        'Zeppole.local':'/Users/jonpreall/Dropbox/',
        }
    jons_computers = list(dropbox_loc_dict.keys())
    print('Detected hostname',computername)
    if computername not in jons_computers: 
        print('computer name not recognized as one of ' + str(jons_computers))
        if 'jpreall' in os.listdir('/Users/'):
            computername = 'Ebelskiver.local'
            print('Inferring that this is Ebelskiver')
        elif 'jonpreall' in os.listdir('/Users/'):
            computername = 'Zeppole.local'
            print('Inferring that this is Zeppole')
    
    dropbox_dir = dropbox_loc_dict[computername]
    print('Setting dropbox directory to ','\''+dropbox_dir+'\'')
    return dropbox_dir

def markers(adata, groupby=None, n_genes=None, use_raw=True, kind='dotplot', standard_scale=None, cmap=None, vmin=None, vmax=None, swap_axes=False, min_logfoldchange=None):
    import jpplot
    import numpy as np
    import pandas as pd
    import scanpy as sc
    
    assert kind in ['dotplot','matrixplot'], 'kind must be one of: dotplot, matrixplot'
    
    if not groupby:
        #Pick a reasonable sounding cluster with the most unique keys
        acceptable_groups = ['Cell_Type','Cluster','Subset_Cluster','cluster']
        acceptable_groups = [g for g in acceptable_groups if g in adata.obs.columns]
        result = {}
        for g in acceptable_groups:
            result[g] = len(adata.obs[g].unique())
        groupby = max(result, key=result.get)
        print('Grouping by key:',groupby)
        
    n_groups = len(adata.obs[groupby].unique())
    
    if not n_genes:
        n_genes = np.floor(40 / len(adata.obs[groupby].unique())).astype('int')
    
    sc.tl.rank_genes_groups(adata, groupby=groupby, use_raw=use_raw)
    sc.tl.dendrogram(adata, groupby=groupby)
    
    if kind == 'dotplot':
        if not cmap:
            cmap = 'Reds'	
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=n_genes, standard_scale=standard_scale, cmap=cmap, vmin=vmin, vmax=vmax, swap_axes=swap_axes)
        
    if kind == 'matrixplot':
        if not cmap:
            cmap = 'cividis'	
        sc.pl.rank_genes_groups_matrixplot(adata, n_genes=n_genes, standard_scale=standard_scale, cmap=cmap, vmin=vmin, vmax=vmax, swap_axes=swap_axes, min_logfoldchange=min_logfoldchange)
    
def export_cellranger_reana(adata, file_prefix, output_folder='./', aggr_key=None, obs_keys=[], layer='counts', genome=None, feature_types='feature_types'):
    import jpplot
    import os
    import pandas as pd
    homedir = os.getcwd()
    
    obs_keys =  obs_keys if isinstance(obs_keys, list) else [obs_keys]
    
    output_folder = os.path.join(output_folder,file_prefix)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    jpplot.write_cellranger_h5(adata, aggr_key=aggr_key, output_folder=output_folder, file_prefix=file_prefix, layer=layer, feature_types=feature_types, genome=None)
    
    
    
    for obs_key in obs_keys:
        adata.obs.loc[:,obs_key].to_csv(output_folder + '/' + obs_key + '.csv', index_label='Barcode', header=True)
    #if not os.path.exists(os.path.join(homedir+file_prefix)):
    #    os.makedirs(os.path.join(homedir+file_prefix))
    jpplot.export_coordinates(adata, output_prefix='projections', destination_folder=output_folder)

def matrix_to_df(adata, layer='X', gene_name_key=None):
    import jpplot
    import os
    import pandas as pd
    import scanpy as sc
    from scipy.sparse import issparse
    import numpy as np
    
    valid_layers = ['X','raw'] + list(adata.layers.keys())
    assert layer in valid_layers, 'layer must be either: X, raw, or a layer in adata.layers'
    
    if layer == 'X':
        print('Using scaled data in adata.X')
        if issparse(adata.X):
            counts = pd.DataFrame(adata.X.toarray())
        else: 
            counts = pd.DataFrame(adata.X)
        
    elif layer == 'raw':
        print('Using raw lognorm data in adata.raw')
        if issparse(adata.raw.X):
            counts = pd.DataFrame(adata.raw.X.toarray())
        else: 
            counts = pd.DataFrame(adata.X)
            
    elif layer == 'counts':
        if issparse(adata.layers[layer]):
            counts = pd.DataFrame(np.floor(adata.layers[layer].toarray()).astype('float64'))
        else: 
            counts = pd.DataFrame(np.floor(adata.layers[layer]).astype('float64'))
    
    elif layer in adata.layers:
        if issparse(adata.layers[layer]):
            counts = pd.DataFrame(adata.layers[layer].toarray())
        else: 
            counts = pd.DataFrame(adata.layers[layer])
    
    else:
        print('Warning, layer',layer,'not found in adata.layers')
    
    if str(counts.__class__) == "<class 'pandas.core.frame.DataFrame'>":
        if layer != 'raw':
            if not gene_name_key:
                counts.columns = adata.var_names
            else:
                counts.columns = np.array(adata.var[gene_name_key])
            counts.index = adata.obs.index
        else:
            if not gene_name_key:
            	counts.columns = adata.raw.var_names
            	shared_genes = adata.var_names.intersection(adata.raw.var_names)
            else:
                counts.columns = adata.raw.var[gene_name_key]
                shared_genes = list(set(adata.var[gene_name_key]).intersection(set(adata.raw.var[gene_name_key])))
            counts.index = adata.raw.obs_names
        	#shared_genes = adata.var_names.intersection(adata.raw.var_names)
        	#print('shared genes:',len(shared_genes))
            shared_cells = adata.obs_names.intersection(adata.raw.obs_names)
        	#print('shared cells:',len(shared_cells))
            counts = counts.loc[shared_cells,shared_genes]
        
    return counts
    
def decontX(adata, layer='counts', output_file=None):
    import scanpy as sc
    from scipy.sparse import issparse
    """
    Function to call celda::decontX from Python
    """
    print('Estimating soup contamination with celda::decontX')
    
    print()
    
    import rpy2.robjects as ro
    import anndata2ri
    import sys

    ro.r('library(celda)')
    ro.r('library(DropletUtils)')
    anndata2ri.activate()
    
    tmp = adata.copy()
    
    if issparse(tmp.X):
        if not tmp.X.has_sorted_indices:
            tmp.X.sort_indices()

    for key in tmp.layers:
        if issparse(tmp.layers[key]):
            if not tmp.layers[key].has_sorted_indices:
                tmp.layers[key].sort_indices()
    
    if layer != 'X':
        
        if 'layers' in tmp.__dir__():
            if layer not in tmp.layers.keys():
                print('layer \''+str(layer)+'\' not detected. If adata.X is the raw digital counts matrix, either add counts layer like this: adata.layers[\'counts\'] = adata.X, or specify layer=\'X\'')
                sys.exit('Aborting')
            else:
                print('Using counts from adata.layers[\''+str(layer)+'\']')
                tmp.X = tmp.layers[layer]
    else:
        print('Using counts from adata.X')
        tmp.layers['counts'] = tmp.X
    
    ##Convert to R object
    ro.globalenv['tmp'] = tmp

    ## Run decontX
    ro.r('tmp <- decontX(x=tmp)')

    print('Saving adjusted counts to adata.layers[\'decontX\']')
    adata.layers['decontX'] = ro.r('tmp@assays@data$decontXcounts').T
    return adata
    
    if output_file:
        print('Saving to',output_file,'...')
        adata.write(output_file)
    print('Done')
    
    
def read_cellranger_h5(cellranger_h5):
    import scanpy as sc
    import h5py
    import re
    import numpy as np
    
    adata = sc.read_10x_h5(cellranger_h5)
    
    adata.var_names_make_unique()
    print('Making var names unique...')
    
    with h5py.File(cellranger_h5, 'r+') as f:
        ATTRIBUTE_KEYS = list(f.attrs.keys())
        print('Keys detected in hdf5 file:', ATTRIBUTE_KEYS)
        print()
        
        #### library id
        if 'library_ids' in ATTRIBUTE_KEYS:
            LIBRARY_ID = [name.decode('utf-8') for name in f.attrs['library_ids']]
            #LIBRARY_ID = [name for name in [f.attrs['library_ids'].decode('utf-8')]]
            print('Detected Library_ID as:',LIBRARY_ID,'... Saving to adata.uns[\'Cellranger_Library_ID\']')
            adata.uns['Cellranger_Library_ID'] = LIBRARY_ID
            
            #save dictionary of barcode suffixes and library_ids for posterity
            BARCODE_SUFFIX_DICT = {'-'+str(k+1):v for k,v in enumerate(LIBRARY_ID)}
            print()
            print('Saving dictionary of {barcode suffixes:library_ids} to adata.uns[\'Barcode_suffixes\']')
            adata.uns['Barcode_suffixes'] = BARCODE_SUFFIX_DICT
            print()
            
        #### chemistry
        if 'chemistry_description' in ATTRIBUTE_KEYS:
            if isinstance(f.attrs['chemistry_description'], (np.ndarray) ):   
                CHEMISTRY = [name.decode('utf-8') for name in f.attrs['chemistry_description']]
                #CHEMISTRY = f.attrs['chemistry_description']

                if len(CHEMISTRY) == 1:
                     CHEMISTRY = CHEMISTRY[0]
            elif type(f.attrs['chemistry_description'])  == bytes:
                CHEMISTRY = f.attrs['chemistry_description'].decode('utf-8')
            else:
                CHEMISTRY = f.attrs['chemistry_description'].decode('utf-8')
                
            print('Detected Chemistry as:',CHEMISTRY, '... Saving to adata.uns[\'Chemistry\']')
            adata.uns['Chemistry'] = CHEMISTRY
            print()
            
        ## cellranger version
        if 'software_version' in ATTRIBUTE_KEYS:
            CELLRANGER_VERSION = f.attrs['software_version']
            print('Detected Cellranger software version as:',CELLRANGER_VERSION, '... Saving to adata.uns[\'Cellranger_Version\']')
            adata.uns['Cellranger_Version'] = CELLRANGER_VERSION
            print()
            
        # genome
        GENOME = None
        TOP_LEVEL_H5_KEY = list(f.keys())[0]
        if 'features' in f[TOP_LEVEL_H5_KEY].keys():
            GENOME = list(set(f[TOP_LEVEL_H5_KEY]['features']['genome'][:].astype('str').tolist()))
            if len(GENOME) > 1:
                print('Warning, multiple genomes detected in hdf5 dataset. Choosing: ',GENOME[0])
            GENOME = GENOME[0]
            print('Detected genome as: ',GENOME)
            adata.uns['genome'] = GENOME
        if GENOME == None:
            GENOME == 'Unknown_Genome'
        f.close()
        
    ## Save library ID to adata.obs
    ## Test if barcodes are in the correct 10X Genomics format:
    
    BARCODES_IN_10X_FORMAT = len([x for x in adata.obs.index.tolist() if re.match('[ATGCN]{16}-[0-9]+$', x, re.IGNORECASE)]) == len(adata.obs)

    if BARCODES_IN_10X_FORMAT:
        print('Barcodes are in 10X format. Adding Library_ID to adata.obs')
        adata.obs['Library_ID'] = adata.obs.index.str.replace('^[ATCG]+-','-')
    
        if 'Barcode_suffixes' in adata.uns.keys():
            adata.obs = adata.obs.replace({'Library_ID':adata.uns['Barcode_suffixes']})
    else:
        print('Barcodes are NOT in 10X format')
    print()
    print(adata)
    return adata
    
    
def predict_sex_aggregated(adata, sex_marker_gene='Xist', sample_key='Sample'):
    
    from sklearn.cluster import KMeans
    import numpy as np
    import scanpy as sc
    import pandas as pd
    import jpplot

    if not sex_marker_gene in adata.var_names:
        print('Marker gene ' + str(sex_marker_gene) + ' not detected in var_names')
        if sex_marker_gene.upper() in adata.var_names:
            sex_marker_gene = sex_marker_gene.upper()
            print('However, ' + str(sex_marker_gene) + ' was detected. Using ' + str(sex_marker_gene))
        elif sex_marker_gene.lower().capitalize() in adata.var_names:
            sex_marker_gene = sex_marker_gene.lower().capitalize()
            print('However, ' + str(sex_marker_gene) + ' was detected. Using ' + str(sex_marker_gene))
            
    
    df = jpplot.matrix_to_df(adata[:,sex_marker_gene])
    df[sample_key] = adata.obs[sample_key]
    df = df.groupby(sample_key).agg('mean')
    X = np.array(df[sex_marker_gene]).reshape(-1,1)
    kmeans = KMeans(n_clusters=2, random_state=0).fit(X)
    df['Sex'] = kmeans.labels_.astype('int')
    mydict = {df.groupby('Sex').agg('mean')[sex_marker_gene].argmin():'M',df.groupby('Sex').agg('mean')[sex_marker_gene].argmax():'F'}
    df = df.replace({'Sex':mydict})
    obs_dict = df['Sex'].to_dict()
    print(df)
    adata.obs['Sex'] = adata.obs[sample_key].copy()
    adata.obs = adata.obs.replace({'Sex':obs_dict})
    adata.uns['Sex_colors'] = ['#d1474e','#6f96e8']


def convert_ensembl_to_gene_name(adata, genome=None):
    import pandas as pd
    import jpplot
	
    if not genome:
        jpplot.detect_genome(adata)
        genome = adata.uns['genome']
    if genome == 'mm10':    
        ens = pd.read_csv('https://www.dropbox.com/s/tzarco19zucgfwo/mm10_genes_CR31.csv?dl=1')
    elif genome in ['hg38','GRCh38']:
        ens = pd.read_csv('https://www.dropbox.com/s/xay40cv0mcb7y86/GRCh38_genes_CR31.csv?dl=1')
        
    ens.set_index('gene_ids', inplace=True)
    ens = ens.loc[adata.var_names.intersection(ens.index),['Gene_name']]
    
    #Duplicate the gene_ids column for later
    adata.var['gene_ids'] = adata.var.index.copy()
    
    #Import Gene_names
    adata.var = adata.var.merge(ens, left_index=True, right_index=True, how='left')
    
    #If no gene short name is in my reference, just use the existing gene_id
    mymask = adata.var.Gene_name.isna()
    adata.var.Gene_name = adata.var.Gene_name.mask(mymask,adata.var.gene_ids)
    
    #Reassign the index to Gene_name
    adata.var.set_index('Gene_name', inplace=True)
    
def TF_IDF(adata, layer='counts'):
    import jpplot
    import numpy as np
    from scipy import sparse
    print('Calculating TF-IDF on layer:',layer)
    
    matrix_df = jpplot.matrix_to_df(adata, layer=layer)

    term_freq_df = matrix_df.div(matrix_df.sum(axis=1), axis=0)
    N = len(adata)
    
    df_i = np.sum(matrix_df > 0, axis=0)
    IDF_df = np.log(N/df_i)
    
    TF_IDF_df = term_freq_df * IDF_df
    print('Adding layer TF_IDF to adata.layers')
    adata.layers['TF_IDF'] = sparse.csr_matrix(TF_IDF_df.values)
    print('Adding column \'IDF\' to adata.var')
    adata.var['IDF'] = np.array(IDF_df)
    
def smooth_gex(adata, neighbors_key=None):
    import pandas as pd
    import numpy as np
    import scanpy as sc
    from scipy import sparse
    

    if neighbors_key:
        neighbors_key = neighbors_key +'_'
        con_key = neighbors_key + '_connectivities'
        dis_key = neighbors_key + '_distances'
        
    else:
        con_key = 'connectivities'
        dis_key = 'distances'
        
    con = pd.DataFrame(adata.obsp[con_key].toarray())
    con.index = adata.obs_names
    con.columns = adata.obs_names

    dis = pd.DataFrame(adata.obsp[dis_key].toarray())
    dis.index = adata.obs_names
    dis.columns = adata.obs_names
    
    if sparse.issparse(adata.raw.X):
        matrix = pd.DataFrame(adata.raw[:,adata.raw.var_names].X.toarray())	
    else:
        matrix = pd.DataFrame(adata.raw[:,adata.raw.var_names].X)	
        
    matrix.index = adata.obs_names
    matrix.columns = adata.raw.var_names


    matrix_smooth = matrix.copy()
    for x in dis.index:
        matrix_smooth.loc[x] = matrix.loc[dis.index[np.where(dis.loc[x] > 0)],:].sum(0)

    print('Adding layer \'smooth\' to adata.layers')
    adata.layers['smooth'] = sparse.csr_matrix(matrix_smooth.loc[:,adata.var_names].values)

def obskey_zscore_by_cluster(adata, key='n_counts', log_transform=True, groupby='Cluster', layer='X', vmin=-2, vmax=3, cmap='RdYlBu_r'):
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy.stats import zscore
    from scipy.sparse import issparse
    
    assert key in adata.obs.select_dtypes('number').columns.tolist() + adata.var_names.tolist(), 'key must be a numerical observation in adata.obs' 
	
    if key in adata.obs.select_dtypes('number'):
        df = adata.obs.loc[:,[key,groupby]]
    elif key in adata.var_names:
        print('Detected',key,'in adata.var_names')
        assert layer in ['X','raw'] + list(adata.layers.as_dict().keys()), 'layer must be: \'X\', \'raw\', or a key in adata.layers'
        df = pd.DataFrame(adata.obs.loc[:,groupby])
    	
        if layer == 'X':
            if issparse(adata.X):
                df[key] = adata[:,key].X.toarray()
            else:
                df[key] = adata[:,key].X
    			
        elif layer == 'raw':
            if issparse(adata.X):
                df[key] = adata.raw[:,key].X.toarray()
            else:
                df[key] = adata.raw[:,key].X
    			
        else:
            if issparse(adata.layers[layer]):
                df[key] = adata[:,key].layers[layer].toarray()
            else:
                df[key] = adata[:,key].layers[layer]
        print(df.columns)
    
    
    if log_transform:
        if layer != 'raw':
            df[key] = np.log10(df[key])
        elif layer == 'raw':
            df[key] = df[key]
        

    out = pd.DataFrame()
    if log_transform:
        key_new_name = 'log_' + str(key) + '_grouped_by_' + str(groupby) + ': Zscore'
    else:
        key_new_name = str(key) + '_grouped_by_' + str(groupby) + ': Zscore'
        
    for c in df[groupby].unique():
        tmp = df[df[groupby] == c]
        tmp[key_new_name] = zscore(tmp[key])
        out = pd.concat([out,tmp])
    del(tmp, df)
    
    if key_new_name in adata.obs.columns:
    	adata.obs = adata.obs.drop(columns=key_new_name)
    	
    adata.obs[key_new_name] = out[key_new_name]
    adata.obs[key_new_name] = adata.obs[key_new_name].fillna(0)
    sc.pl.umap(adata, color=key_new_name, cmap=cmap, vmin=vmin, vmax=vmax)
    

def read_cellranger_atac(h5_file, peak_annotation_tsv=None):
    from scipy.sparse import csr_matrix
    import h5py
    import anndata as ad
    import os

    with h5py.File(h5_file,'r') as f:
        M, N = f['matrix/shape'][()]
        data = f['matrix/data'][()]
        matrix = csr_matrix(
            (data, f['matrix/indices'], f['matrix/indptr']),
            shape=(N, M),
        )


        #BARCODES = [name.decode() for name in f['matrix/barcodes']]
        BARCODES = f['matrix/barcodes'][()].astype('str')
        VAR_NAMES=f['matrix/features/name'][()].astype(str)
        GENE_IDS=f['matrix/features/id'][()].astype(str),
        FEATURE_TYPES=f['matrix/features/feature_type'][()].astype(str)
        GENOME=f['matrix/features/genome'][()].astype(str)
        SOFTWARE_VERSION = f.attrs['software_version']

    obs_dict = dict(obs_names=BARCODES)
    var_dict = dict(var_names=VAR_NAMES,gene_ids=GENE_IDS,feature_types=FEATURE_TYPES,genome=GENOME)

    adata = ad.AnnData(matrix, obs_dict, var_dict)
    adata.uns['sotware_version'] = SOFTWARE_VERSION
    if len(adata.var['genome'].unique() == 1):
        adata.uns['genome'] = adata.var['genome'].unique()[0]
        
    if peak_annotation_tsv:
        assert os.path.exists(peak_annotation_tsv), 'Warning, invalid path to peak_annotation.tsv file'
        import pandas as pd
        anno = pd.read_table(peak_annotation_tsv, index_col=0)
        anno['renamed'] = [str(g[0]+':'+g[1]+'-'+g[2]) for g in anno.index.str.split('_')]
        anno = anno.set_index('renamed')
        shared_peaks = list(adata.var_names.intersection(anno.index))
        diff_peaks = anno.index.intersection(adata.var_names)
        print('Excluding',len(diff_peaks),'peaks from adata with no annotations')
        adata = adata[:,shared_peaks].copy()
        adata.var = adata.var.merge(anno.loc[shared_peaks], left_index=True, right_index=True, how='left')
		
    print(adata)
    return adata

def import_gene_coords(adata, genome=None, use_local=True):

    
    if 'genome' in list(adata.uns.keys()):
        genome=adata.uns['genome']
    
    print('Using genome:',genome)
    assert genome, 'No genome specified'
    if genome == 'mm10':
        if use_local == True:
            var_coords = pd.read_csv('/Users/jpreall/Dropbox/Preall_Lab/Preall/Gene_Lists/mm10_genes_with_coords.csv', index_col='Accession')
        else:
            var_coords = pd.read_csv('https://www.dropbox.com/s/alli12d3jp4dz09/mm10_genes_with_coords.csv?dl=1', index_col='Accession')
            shared_genes = set(adata.var.index).intersection(var_coords.index)
            #adata.var = adata.var.merge(var_coords.loc[shared_genes, ['Chromosome','Start','End','Strand']], left_index=True, right_index=True, how='left')
            
    elif genome in ['hg38','GRCh38']:
        if use_local == True:
            var_coords = pd.read_csv('/Users/jpreall/Dropbox/Preall_Lab/Preall/Gene_Lists/GTF/GRCh38/1.2.0/gene_coords.csv.gz', index_col='gene_ids')
        else:
            var_coords = pd.read_csv('https://www.dropbox.com/s/fdvvbkdo21fapji/gene_coords.csv.gz?dl=1', index_col='gene_ids')
    
    ## Merge
    adata_obj_columns = adata.var.select_dtypes('object')
    ensembl_column = adata_obj_columns.columns[adata_obj_columns.apply(lambda x: x.str.match('ENS')).all()][0]
    shared_genes = set(adata.var[ensembl_column]).intersection(var_coords.index)
    
    var = pd.DataFrame(adata.var.loc[:,ensembl_column]).set_index(ensembl_column)
    var['gene_ids'] = var.index.copy()
    var = var.merge(var_coords, left_index=True, right_index=True, how='left').set_index('Gene_name')
    new_cols = [name for name in var.columns if name not in adata.var.columns]
    
    if len(new_cols) > 0:
        print('Adding columns',new_cols,'to adata.var')
        for col in new_cols:
            adata.var[col] = np.array(var[col])
    elif len(set(var.columns).intersection(set(adata.var.columns))) > 0:
        print('Already detected columns in adata.var:',list(set(var.columns).intersection(set(adata.var.columns))))
    else:
        print('Something seems to have gone wrong')
    
def transfer_cell_types_svm(adata_ref, adata_target, cell_type_key, n_features=50, pval_adj_thresh=1e-6, kernel='rbf', output_key='Predicted_Cell_Type', groupby_cluster=None):
    
    # optimize n_features:
    import pandas as pd
    import numpy as np
    import scanpy
    from sklearn.model_selection import train_test_split
    from sklearn import svm
    import jpplot
    
    train = adata_ref.copy()
    train.var_names_make_unique()
    target = adata_target.copy()
    
    #assertions
    assert cell_type_key in train.obs.select_dtypes('category').columns, 'cell_type_key must be a categorical column in adata.obs.  Available options: ' + str(train.obs.select_dtypes('category').columns.tolist())
    assert kernel in ['linear','poly','rbf'], 'kernel must be one of [\'linear\',\'poly\',\'rbf\']'
    if groupby_cluster:
        assert groupby_cluster in target.obs.select_dtypes('category').columns, 'groupby_cluster must be a categorical column in adata.obs.  Available options: ' + str(train.obs.select_dtypes('category').columns.tolist())
    
    kernel_text_dict = {'linear':' a linear kernel','poly':' a polynomial kernel','rbf':' a gaussian kernel'} 
    print('Generating an SVM using' + kernel_text_dict[kernel] + '...')
    
    
    print('Generating log(CPM + 1) matrix...')
    
    if train.raw:
        print('Using raw data in adata_train.raw')
        genes_exlusive_to_raw = train.var_names.difference(train.raw.var_names)
        if len(genes_exlusive_to_raw) > 0:
            print('Excluding genes not appearing in processed data frame:',genes_exlusive_to_raw)
            train = train[:,genes_exlusive_to_raw].copy()
            train.X = train.raw.X
            
        else:
            print("adata.raw.var_names matches adata.var_names")
    elif 'counts' in train.layers:
        train.X = train.layers['counts']
        sc.pp.normalize_total(train, target_sum=1e4)
        sc.pp.log1p(train)

    ## instantiate results dataframe
    results_total = pd.DataFrame(columns=['genes_per_group','accuracy','precision','recall'])
    results_total.loc[len(results_total)] = np.array([0,0,0,0])

    sc.tl.rank_genes_groups(train, groupby=cell_type_key)
    empirical_markers = jpplot.topdiff(train, groupby=cell_type_key, n_genes=n_features, pval_adj_thresh=pval_adj_thresh)
    #deduplicate marker genes in test data
    empirical_markers = list(set([name for name in empirical_markers if name in train.var_names]))
    
    # exclude marker genes not occurring in target data
    filtered_markers = list(set(empirical_markers).intersection(target.var_names))
    
    print('Excluding',len(empirical_markers) - len(filtered_markers),'marker genes that do not appear in target dataset')
    print('Excluded genes:', list(set(empirical_markers).difference(set(filtered_markers))))
    
          
    empirical_cell_types = train.obs[cell_type_key]
    
    print("Number of marker genes used:",len(filtered_markers))
    df = jpplot.matrix_to_df(train[:,filtered_markers], layer='X')
    
    # Split dataset into training set and test set
    X_train, X_test, y_train, y_test = train_test_split(df, empirical_cell_types, test_size=0.3,random_state=42) # 70% training and 30% test

    #Create a svm Classifier
    #clf = svm.SVC(kernel='linear') # Linear Kernel
    #clf = svm.SVC(kernel='poly') # Polynomial Kernel
    #clf = svm.SVC(kernel='rbf') # Gaussian Kernel
    clf = svm.SVC(kernel=kernel) 

    #Train the model using the training sets
    print('Training model using n_features = ',n_features,'per cell type')
    clf.fit(X_train, y_train)

    #Predict the response for test dataset
    y_pred = clf.predict(X_test)

    #Import scikit-learn metrics module for accuracy calculation
    from sklearn import metrics

    # Model Accuracy: how often is the classifier correct?
    acc = metrics.accuracy_score(y_test, y_pred)

    # Model Precision: what percentage of positive tuples are labeled as such?
    prec = metrics.precision_score(y_test, y_pred, average='weighted')

    # Model Recall: what percentage of positive tuples are labelled as such?
    rec = metrics.recall_score(y_test, y_pred, average='weighted')
    
    print('Accuracy:', acc)
    print('Precision:',prec)
    print('Recall:',rec)

    print('Predicting class in adata_target...')
    ## Add some code here to allow for use of a UMI counts layer instead...
    target_data = jpplot.matrix_to_df(target[:,filtered_markers], layer='raw')
    SVM_prediction = clf.predict(target_data)
    
    print("Adding predicted cell type to adata.obs[\'" + str(output_key) + "\']")
    adata_target.obs[output_key] = SVM_prediction
    
    if groupby_cluster:
        df = pd.crosstab(adata_target.obs[groupby_cluster],adata_target.obs[output_key])
        cluster_dict = dict(zip(df.index.tolist(), df.idxmax(1).tolist()))
        adata_target.obs[output_key + '_grouped'] = adata_target.obs[groupby_cluster].copy()
        adata_target.obs = adata_target.obs.replace({output_key + '_grouped':cluster_dict})
        
    # Visualize the results with an available embedding:
    if 'X_umap' in adata_target.obsm.keys():
        sc.pl.umap(adata_target, color=output_key)
        if groupby_cluster:
            sc.pl.umap(adata_target, color=output_key + '_grouped')
        
    elif 'X_tsne' in adata_target.obsm.keys():
        sc.pl.tsne(adata_target, color=output_key)
        if groupby_cluster:
            sc.pl.tsne(adata_target, color=output_key + '_grouped')

def cmapjp():
	## Usage: jpcolors = cmapjp()
	## This produces a truncated version of the RdYlBu_r color map where the lowest value is light blue
	## Then when plotting, you can use color_map=jpcolors or cmap=jpcolors
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as pl
    from matplotlib import cm
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap

    # ``matplotlib.cm.get_cmap`` was removed in Matplotlib 3.10.  Use the
    # current colormap registry, while retaining compatibility with older
    # Matplotlib versions that may still be used outside the Scanpy container.
    if hasattr(matplotlib, "colormaps"):
        RdYlBu_r = matplotlib.colormaps.get_cmap('RdYlBu_r').resampled(256)
    else:
        RdYlBu_r = cm.get_cmap('RdYlBu_r', 256)
    jpcolors = ListedColormap(RdYlBu_r(np.linspace(.25,1,256)))
    print('Elegant colormap returned as: jpcolors.  Use this in plotting.')
    return jpcolors 
    
    
def volcano_plot(adata, groupby, group, pval_thresh=1e-20, labelgenes=None, s=5, alpha=None):
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as pl
    import seaborn as sns
    
    if adata.uns['rank_genes_groups']:
        if adata.uns['rank_genes_groups']['params']['groupby'] != groupby:
            sc.tl.rank_genes_groups(adata, groupby=groupby)
    assert group in adata.obs[groupby].cat.categories, 'group must be a valid category in the column being grouped'
    
    gex = sc.get.rank_genes_groups_df(adata, group=group).set_index('names')
    gex['-logP'] = -np.log10(gex['pvals_adj'])
    sns.scatterplot(data=gex, x='logfoldchanges', y='-logP', linewidth=0, s=s, alpha=alpha, color='lightgrey')
    
    siggex = gex[gex['pvals_adj'] < pval_thresh]
    sns.scatterplot(data=siggex, x='logfoldchanges', y='-logP', linewidth=0, s=s, alpha=alpha, color='Red')
    
    # plot text labels:
    if labelgenes:
        for gene in labelgenes:
            x = gex.loc[gene]['logfoldchanges']
            y = gex.loc[gene]['-logP']
            pl.text(x*3,y,gene)
    pl.grid(None)
    

def prep_mouse_cellphonedb_df(adata, labels, output_folder=None, file_prefix=None, layer='raw', use_local=True):
    import scanpy
    import jpplot
    import numpy as np
    import pandas as pd
    import os, sys
    
    from datetime import datetime
    now = datetime.now().strftime("%Y%m%d.%H%M")

    if not output_folder:
        output_folder = os.getcwd()
        
    if not file_prefix:
        print('No file prefix specified. Defaulting to ',now)
        file_prefix = now
        
    assert labels in adata.obs.select_dtypes('category').columns, 'labels must be a categorial column in adata.obs'
    assert layer in adata.layers.keys(), 'layer must be a key in adata.layers'
    
    orthofile_loc = '/Users/jpreall/Dropbox/Preall_Lab/Preall/Gene_Lists/biomaRt_ortholog_mouse_human.csv'
    orthofile_db = 'https://www.dropbox.com/s/atg8tks7djp41kd/biomaRt_ortholog_mouse_human.csv?dl=1'
    
    if use_local:
        if os.path.exists(orthofile_loc):
            ortho = pd.read_csv(orthofile_loc)
        else:
            print('Local copy of orthologs file not detected.  This must not be Jon\'s computer...')
            print('Reading orthologs file from Dropbox.')
            ortho = pd.read_csv(orthofile_db)
    else:
        ortho = pd.read_csv(orthofile_db)
    
    hum_genes_ens_loc = '/Users/jpreall/Dropbox/Preall_Lab/Preall/Gene_Lists/GRCh38_genes.csv'
    hum_genes_ens_db = 'https://www.dropbox.com/s/7id4mjsc8f44uu5/GRCh38_genes.csv?dl=1'
    
    if use_local:
        if os.path.exists(hum_genes_ens_loc):
            humgenes = pd.read_csv(hum_genes_ens_loc)
        else:
            print('Local copy of gene_name -> Ensembl file not detected.  This must not be Jon\'s computer...')
            print('Reading file from Dropbox.')
            humgenes = pd.read_csv(hum_genes_ens_db)
    else:
        humgenes = pd.read_csv(hum_genes_ens_db)
    
    ortho.index = ortho['HGNC.symbol']
    humgenes.index = humgenes['Gene_name']
    ortho = ortho.merge(humgenes.loc[:,'EnsemblID'], left_index=True, right_index=True, how='left')
    ortho = ortho[~ortho.index.duplicated()].sort_index()
    ortho.index = ortho['MGI.symbol']
    ortho = ortho[~ortho.index.duplicated()].sort_index()
    ortho = ortho.loc[[name for name in ortho.index if name in adata.var_names]]
    gex = jpplot.matrix_to_df(adata[:,ortho.index], layer=layer).T
    gex.index = ortho['EnsemblID']
    
    labels_output_file = os.path.join(output_folder, file_prefix + '_' + labels + '_cellphonedb.txt')
    print('Writing labels to ' + labels_output_file)
    adata.obs.loc[:,labels].to_csv(labels_output_file, header=True,index_label='Cell', sep='\t')
    
    gex_output_file = os.path.join(output_folder, file_prefix + '_gex_cellphonedb.txt')
    print('Writing GEX to ' + gex_output_file, 'using layer:',layer)
    gex.to_csv(gex_output_file, header=True,index_label='Gene', sep='\t')
    
def score_jplate_sigs(adata, species='Mouse'):
    import glob
    import os
    import scanpy as sc
    import numpy as np
    import pandas as pd
    siglist = glob.glob('/Users/jpreall/Dropbox/Preall_Lab/Preall/Gene_Lists/JPlate_sigs/' + '*pos_markers.csv')
    
    for sig in siglist:
        SIGNAME = os.path.basename(sig).replace('_pos_markers.csv','')
        #SIGGENES = [name for name in pd.read_csv('/Users/jpreall/Dropbox/Preall_Lab/Preall/Gene_Lists/JPlate_sigs/IFNg_pos_markers.csv')['Gene_Name'] if name in adata.var_names]
        SIGGENES = pd.read_csv(sig)['Gene_Name'].tolist()
        sc.tl.score_genes(adata, score_name=SIGNAME+'_score', gene_list=SIGGENES)
        
def isIntegers(array_like):
	from scipy.sparse import issparse
	import numpy
	elements_to_check = 100000
	if issparse(array_like):
		data = array_like.data[:elements_to_check]
	else:
		rows_to_check = np.ceil(elements_to_check/array_like.shape[1]).astype('int')
		data = array_like[:rows_to_check,:]
	result = np.all(np.mod(data, 1) == 0)
	if result:
		print('This array contains only integers')
	else:
		print('This array contains non-integers')
	return result
	
def umap(adata, color, color_map=None, layer=None, vmax=None, vmin=None, frameon=None, title=None, legend_loc=None, legend_fontsize=None):
	import scanpy as sc
	import jpplot
	
	if not color_map:
		color_map = jpplot.cmapjp()
	
	sc.pl.umap(adata, 
		color=color, 
		color_map=color_map, 
		layer=layer,
		vmax=vmax,
		vmin=vmin,
		frameon=frameon,
		title=title,
		legend_loc=legend_loc,
		legend_fontsize=legend_fontsize)

def magic(adata, layer='raw'):
    import scanpy as sc
    tmp = adata
    available_layers = ['X','raw'] + list(adata.layers)
    assert layer in available_layers, 'layer must be one of ' + str(available_layers)
    
    if layer == 'raw':
        tmp.X = tmp.raw[:,tmp.var_names].X
    elif layer == 'X':
        pass
    else:
        tmp.X = tmp.layers[layer]
        
    sc.external.pp.magic(tmp)
    adata.layers['MAGIC'] = tmp.X.copy()
    
def topdiff_df(adata, group, groupby=None, pvals_adj_thresh=1e-12):
    
    if not groupby:
        if adata.uns['rank_genes_groups']:
            groupby = adata.uns['rank_genes_groups']['params']['groupby']
            print('Using groupby = ',groupby)
        if groupby != adata.uns['rank_genes_groups']['params']['groupby']:
            sc.tl.rank_genes_groups(adata, groupby=groupby)
            
    df = sc.get.rank_genes_groups_df(adata, group=group)

    df = df[df['pvals_adj'] < pvals_adj_thresh]
    df = df.sort_values(by='logfoldchanges', ascending=False)
    return df
