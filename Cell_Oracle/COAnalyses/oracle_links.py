import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import celloracle as co

class oracle_links:
    def __init__(self, file_name):
        """Load the links object from the cell oracle. The links object should be ending with xx.celloracle.links.

        Args:
            file_name (str): the path to the object
        """
        self.links = co.load_hdf5(file_name)

    # filter the links by p value first.
    def filter_pval(self, p_thresh = 10**-5):
        """Before doing futhur filtering, first delete low quality links by using p value.
            This method will change the links_dict in the object.
        Args:
            p_thresh (int, optional): the p value threshold. Defaults to 10**-5.
        """
        for c in range(len(self.links.links_dict.keys())):
            print(f"Before Filtering: cluster {str(c)} has {len(self.links.links_dict[str(c)])} links.")
            post_links = pd.DataFrame(self.links.links_dict[str(c)])
            post_links = post_links[post_links['p'] < p_thresh]
            print(f"After Filtering: cluster {str(c)} has {len(post_links)} links.")
            self.links.links_dict[str(c)] = post_links

    def calc_weighted_logp(self):
        """Add a column of weighted logp to each dataframe for each cluster in the \
            links_dict in the object.
        """
        for c in range(len(self.links.links_dict.keys())):
           df =  pd.DataFrame(self.links.links_dict[str(c)])
           df['weighted_logp'] = df['-logp'] * df['coef_abs']
           df = df.sort_values(by=['weighted_logp'], ascending=False)
           self.links.links_dict[str(c)] = df

    def filter_human_TF(self, human_TF):
        """Filter out links(edges) that neither source or target node is \
           in the human_TF dataframe
        Args:
            human_TF (data frame): a dataframe with one column containing human TF names.
        """
        for c in range(len(self.links.links_dict.keys())):
            df = self.links.links_dict[str(c)]
            filtered_df = df[df['source'].isin(human_TF[0]) \
                               | df['target'].isin(human_TF[0])]
            self.links.links_dict[str(c)] = filtered_df
            print(f"Cluster {str(c)} has {len(df) - len(filtered_df)} links that are not found in database.")


    def get_top_links(self, percentile = 0.1):
        """For each cluster and each TF, get the top percentile links based on weighted logp.\
           Save the degree dataframe in self.TF_degree_dict.
           Save the top links in filtered_links_dict

        Args:
            percentile (float, optional): the percentage of top links. Defaults to 0.1.
        """
        TF_degree_dict = {}
        filtered_links_dict = {}
        for c in range(len(self.links.links_dict.keys())):
            df = self.links.links_dict[str(c)]
            print(f"Before filtering: Cluster {str(c)} has {len(df)} links.")
            # get the list of all TFs in the cluster
            list_TFs = list(set(df['source']).union(set(df['target'])))
            
            # create a dataframe to store the degree of each TF
            out_counts = df['source'].value_counts()
            in_counts = df['target'].value_counts()
            degree_df = pd.DataFrame(index = list_TFs)
            degree_df['out_degree'] = degree_df.index.map(out_counts).fillna(0)
            degree_df['in_degree'] = degree_df.index.map(in_counts).fillna(0)
            degree_df['degree'] = degree_df['out_degree'] + degree_df['in_degree']
            degree_df['filtered_out_degree'] = 0
            degree_df['filtered_in_degree'] = 0
            degree_df['filtered_degree'] = 0
            
            filtered_edges = pd.DataFrame()
            for TF in list_TFs:
                temp_df = df[(df['source'] == TF) | (df['target'] == TF)]
                temp_df = temp_df.sort_values(by=['weighted_logp'], ascending=False)
                temp_df = temp_df.iloc[0: int(len(temp_df) * percentile)]
                filtered_edges = filtered_edges.append(temp_df)
                degree_df.at[TF, 'filtered_out_degree'] = len(temp_df[temp_df['source'] == TF])
                degree_df.at[TF, 'filtered_in_degree'] = len(temp_df[temp_df['target'] == TF])
                degree_df.at[TF, 'filtered_degree'] = len(temp_df)

            print(f"After filtering: Cluster {str(c)} has {len(filtered_edges)} links.")
            TF_degree_dict[str(c)] = degree_df.sort_values(by=['filtered_degree'], ascending=False)
            filtered_links_dict[str(c)] = filtered_edges
        
        # save the 2 dicts to object
        self.TF_degree_dict = TF_degree_dict
        self.filtered_links_dict = filtered_links_dict
    
    def find_TF_overlap_SLIDE(self, latent_factors, human_TF):
        """Pull out TFs discovered in SLIDE results. Store the result in a dict in \
        oracle_links.degree_1_overlap.
        Args:
            latent_factors (dict): a dictionary of dataframes, \
                where each dataframe is a latent factor.
        """
        degree_1_overlap = dict()
        base_GRN = co.data.load_human_promoter_base_GRN()
        for c in range(len(self.filtered_links_dict.keys())):
            df = self.filtered_links_dict[str(c)]

            # get all the node names for that cluster
            cluster_TF = list(set(df['source']).union(set(df['target'])))
            print(f"Cluster {str(c)} has {len(cluster_TF)} nodes.")
            # get the overlap between cluster node names and human TFs
            cluster_TF = [x for x in cluster_TF if x in list(human_TF[0])]
            print(f"Cluster {str(c)} has {len(cluster_TF)} TFs.")

            overlap_df = pd.DataFrame()
            for lf in latent_factors.keys():
                lf_df = latent_factors[lf]
                lf_TF = list(lf_df['names'])
                overlap = list(set(cluster_TF).intersection(set(lf_TF)))
                temp_df = pd.DataFrame({'cluster':str(c), 'latent_factor': lf, 'overlap': overlap})
                overlap_df = overlap_df.append(temp_df)
            degree_1_overlap[str(c)] = overlap_df
        self.degree_1_overlap = degree_1_overlap

    def find_gene_TF_SLIDE(self, latent_factors, human_TF):
        """Find the genes in SLIDE latent factors that are connected to a TF in the context specific GRN.
            The resulting dictionary is saved in the object at self.gene_TF_SLIDE.
            The list of the linked TF is saved at self.linked_TF_SLIDE.
        Args:
            latent_factors (dict): a dictionary of dataframes, \
                where each dataframe is a latent factor.
        """
        gene_TF = dict() # this is a dict of dict, where the first key is the cluster, the second key is the latent factor
        linked_TF = dict()
        # for each louvain cluster
        for c in range(len(self.filtered_links_dict.keys())):
            df = self.filtered_links_dict[str(c)]
            
            # for each latent factor
            gene_TF_lf = dict()
            linked_TF_df = pd.DataFrame()
            for lf in latent_factors.keys():
                lf_df = latent_factors[lf]
                lf_genes = lf_df['names'].tolist()
                # find the genes in the latent factor that are not TFs
                lf_genes = [x for x in lf_genes if x not in human_TF[0].tolist()]
                # find the edges of the latent factor that are also in the GRN
                overlap_df = df[df['source'].isin(lf_genes) | df['target'].isin(lf_genes)]
                assert overlap_df['target'].isin(lf_genes).sum() == len(overlap_df), "not all target nodes are genes in SLIDE lfs."
                
                # add a column to count for the degree of the source node
                # this is to get the information on what are the master TFs for this lf
                degree = overlap_df['source'].value_counts()
                overlap_df['source_degree'] = overlap_df['source'].map(degree)
                overlap_df = overlap_df.sort_values(by=['source_degree'], ascending=False)
                gene_TF_lf[lf] = overlap_df

                # pull out the TF names that link to genes in the LF
                linked_TF_list = list(set(overlap_df['source']).union(set(overlap_df['target'])))
                linked_TF_list = [x for x in linked_TF_list if x in human_TF[0].tolist()]
                temp_df = pd.DataFrame({'cluster':str(c), 'latent_factor': lf, 'linked_TF': linked_TF_list})
                linked_TF_df = linked_TF_df.append(temp_df)
            
            gene_TF[str(c)] = gene_TF_lf
            linked_TF[str(c)] = linked_TF_df
        
        self.gene_TF_SLIDE = gene_TF
        self.linked_TF_SLIDE = linked_TF

    