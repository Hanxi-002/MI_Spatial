#%%
def _adata_to_matrix(adata, layer_name, transpose=True):
    """
    Extract an numpy array from adata and returns as numpy matrix.

    Args:
        adata (anndata): anndata

        layer_name (str): name of layer in anndata

        trabspose (bool) : if True, it returns transposed array.

    Returns:
        2d numpy array: numpy array
    """
    if isinstance(adata.layers[layer_name], np.ndarray):
        matrix = adata.layers[layer_name].copy()
    else:
        matrix = adata.layers[layer_name].todense().A.copy()

    if transpose:
        matrix = matrix.transpose()

    return matrix.copy(order="C")

def update_adata(adata):
    # Update Anndata
    # Anndata generated with Scanpy 1.4 or less should be updated with this function
    # This function will be depricated in the future.

    try:
        lo = adata.uns['draw_graph']['params']['layout']
        if isinstance(lo, np.ndarray):
            lo = lo[0]
        adata.uns['draw_graph']['params']['layout'] = lo
    except:
        pass

def _check_color_information_and_create_if_not_found(adata, cluster_column_name, embedding_name):
    if f"{cluster_column_name}_colors" in adata.uns.keys():
        pass
    else:
        message = f"Color information for the {cluster_column_name} is not found in the anndata. CellOracle is plotting the clustering data, {cluster_column_name}, to create color data."
        warnings.warn(message, UserWarning)
        sc.pl.embedding(adata, basis=embedding_name, color=cluster_column_name)

def _get_clustercolor_from_anndata(adata, cluster_name, return_as):
    """
    Extract clor information from adata and returns as palette (pandas data frame) or dictionary.

    Args:
        adata (anndata): anndata

        cluster_name (str): cluster name in anndata.obs

        return_as (str) : "palette" or "dict"

    Returns:
        2d numpy array: numpy array
    """
        # return_as: "palette" or "dict"
    def float2rgb8bit(x):
        x = (x*255).astype("int")
        x = tuple(x)

        return x

    def rgb2hex(rgb):
        return '#%02x%02x%02x' % rgb

    def float2hex(x):
        x = float2rgb8bit(x)
        x = rgb2hex(x)
        return x

    def hex2rgb(c):
        return (int(c[1:3],16),int(c[3:5],16),int(c[5:7],16), 255)

    pal = get_palette(adata, cluster_name)
    if return_as=="palette":
        return pal
    elif return_as=="dict":
        col_dict = {}
        for i in pal.index:
            col_dict[i] = np.array(hex2rgb(pal.loc[i, "palette"]))/255
        return col_dict
    else:
        raise ValueErroe("return_as")
    return 0
def get_palette(adata, cname):
    c = [i.upper() for i in adata.uns[f"{cname}_colors"]]
    #c = sns.cubehelix_palette(24)
    """
    col = adata.obs[cname].unique()
    col = list(col)
    col.sort()
    """
    try:
        col = adata.obs[cname].cat.categories
        pal = pd.DataFrame({"palette": c}, index=col)
    except:
        col = adata.obs[cname].cat.categories
        c = c[:len(col)]
        pal = pd.DataFrame({"palette": c}, index=col)
    return pal
#%%
from sklearn.svm import SVR
winsorize = False
min_expr_avg=0
max_expr_avg=20
min_expr_cells=2
N=3000
svr_gamma=None

#%%
from typing import Dict, Any, List, Union, Tuple
def score_cv_vs_mean(self, N: int=3000, min_expr_cells: int=2, max_expr_avg: float=20, min_expr_avg: int=0, svr_gamma: float=None,
                        winsorize: bool=False, winsor_perc: Tuple[float, float]=(1, 99.5), sort_inverse: bool=False, plot: bool=False) -> np.ndarray:

    X = _adata_to_matrix(self.adata, "raw_count")
    print("I'm here")
    if winsorize:
        if min_expr_cells <= ((100 - winsor_perc[1]) * X.shape[1] * 0.01):
            min_expr_cells = int(np.ceil((100 - winsor_perc[1]) * X.shape[0] * 0.01)) + 2
            logging.debug(f"min_expr_cells is too low for winsorization with upper_perc ={winsor_perc[1]}, upgrading to min_expr_cells ={min_expr_cells}")

    detected_bool = ((X > 0).sum(1) > min_expr_cells) & (X.mean(1) < max_expr_avg) & (X.mean(1) > min_expr_avg)
    Sf = X[detected_bool, :]
    if winsorize:
        down, up = np.percentile(Sf, winsor_perc, 1)
        Sfw = np.clip(Sf, down[:, None], up[:, None])
        mu = Sfw.mean(1)
        sigma = Sfw.std(1, ddof=1)
    else:
        mu = Sf.mean(1)
        sigma = Sf.std(1, ddof=1)

    cv = sigma / mu
    log_m = np.log2(mu)
    log_cv = np.log2(cv)

    if svr_gamma is None:
        svr_gamma = 150. / len(mu)
        logging.debug(f"svr_gamma set to {svr_gamma}")
    # Fit the Support Vector Regression
    clf = SVR(gamma=svr_gamma)
    clf.fit(log_m[:, None], log_cv)
    fitted_fun = clf.predict
    ff = fitted_fun(log_m[:, None])
    score = log_cv - ff
    if sort_inverse:
        score = - score
    nth_score = np.sort(score)[::-1][N]
    if plot:
        scatter_viz(log_m[score > nth_score], log_cv[score > nth_score], s=3, alpha=0.4, c="tab:red")
        scatter_viz(log_m[score <= nth_score], log_cv[score <= nth_score], s=3, alpha=0.4, c="tab:blue")
        mu_linspace = np.linspace(np.min(log_m), np.max(log_m))
        plt.plot(mu_linspace, fitted_fun(mu_linspace[:, None]), c="k")
        plt.xlabel("log2 mean S")
        plt.ylabel("log2 CV S")
    self.cv_mean_score = np.zeros(detected_bool.shape)
    self.cv_mean_score[~detected_bool] = np.min(score) - 1e-16
    self.cv_mean_score[detected_bool] = score
    self.cv_mean_selected = self.cv_mean_score >= nth_score
    self.cv_mean_selected_genes = self.adata.var.index[self.cv_mean_selected].values
# %%
embedding_name="X_umap"
cluster_column_name="louvain"
#oracle = co.Oracle()

oracle.adata = adata.copy()
update_adata(oracle.adata)
oracle.cluster_column_name = cluster_column_name
oracle.embedding_name = embedding_name
oracle.embedding = oracle.adata.obsm[embedding_name].copy()
oracle.adata.layers["raw_count"] = oracle.adata.X.copy()

sc.pp.log1p(oracle.adata)

oracle.adata.layers["normalized_count"] = oracle.adata.X.copy()


_check_color_information_and_create_if_not_found(adata=oracle.adata,
                                                    cluster_column_name=cluster_column_name,
                                                    embedding_name=embedding_name)


col_dict = _get_clustercolor_from_anndata(adata=oracle.adata,
                                            cluster_name=oracle.cluster_column_name,
                                            return_as="dict")

oracle.colorandum = np.array([col_dict[i] for i in oracle.adata.obs[oracle.cluster_column_name]])

n = min(adata.shape[1], 1000) - 1

oracle.score_cv_vs_mean(405, plot=False, max_expr_avg=100)
