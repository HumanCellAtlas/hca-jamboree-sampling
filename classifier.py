## Some functions for classification

from sklearn.ensemble import RandomForestClassifier

def RForest_top_vars(adata, label_key="included", n_top=36):
    """Using random forest to rank features and 
    return the feature index of the N features with top importance
    """
    RF_class = RandomForestClassifier(n_estimators=100, n_jobs=-1)
    RF_class.fit(adata.X, adata.obs[label_key])
    var_idx = np.argsort(RF_class.feature_importances_)[::-1][:n_top]
    
    return var_idx


def gate_prob(adata, train_set, test_set, var_ids, label_key, model):
    """ Function to determine the number of test cells occupying each
    reference sphere

    :param pd.Index test_set: Pandas index of test set observation names
    """
    if any(~train_set.isin(self.adata.obs_names)):
        raise ValueError(
            'Some of the cells in the test set are not in the AnnData object. '
            'Ensure that all the test cells are in the AnnData object'
        )
    if any(~test_set.isin(self.adata.obs_names)):
        raise ValueError(
            'Some of the cells in the test set are not in the AnnData object. '
            'Ensure that all the test cells are in the AnnData object'
        )
        
    var_idx = p.where(self.adata.var_names.isin(var_ids))[0]

    # Test and train data
    X_test = self.adata[_find_cell_indices(self.adata, test_set),:][:, var_idx].X
    X_train = self.adata[_find_cell_indices(self.adata, train_set),:][:, var_idx].X
    Y_train = self.adata[_find_cell_indices(self.adata, train_set),:].obs[label_key]
    
    model.fit(X_train, Y_train)
    pred_prob = model.predict(X_test)
    
    return pred_prob, pred_label, model


def DTree_gate(adata, label_key="included"):
    """
    
    """
    gates = None
    return gates
    