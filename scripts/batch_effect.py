# -*- coding: utf-8 -*-
"""
Corrección de batch effect con ComBat / NeuroCombat
Autor/es: Brandon Díaz Martín, Ana Hernandez Garcia, Martín Huqiao Astoreca Martín
"""

# ------------------------
# Librerías
# ------------------------
import numpy as np
import pandas as pd
import random
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.spatial.distance import pdist, squareform
import seaborn as sns

from combat.pycombat import pycombat   # pip install combat
from neuroCombat import neuroCombat    # pip install neurocombat

# ------------------------
# Semilla global
# ------------------------
SEED = 123
np.random.seed(SEED)
random.seed(SEED)


# ------------------------
# Funciones auxiliares
# ------------------------
def add_database_column(df):
    """Agrega columna Data_base según el prefijo del patient_id."""
    df["Data_base"] = df["patient_id"].apply(
        lambda x: "GTEX" if str(x).startswith("GTEX-") 
        else ("TCGA" if str(x).startswith("TCGA-") else "GTEX")
    )
    # Reordenar columna en posición 2
    col = df.pop("Data_base")
    df.insert(2, "Data_base", col)
    return df


def apply_combat(train_df, test_df, batch_col="Grupo"):
    """Aplica ComBat a train y test, devuelve DataFrames corregidos."""
    train_expr = train_df.drop(columns=["patient_id", batch_col, "Data_base"]).T
    test_expr  = test_df.drop(columns=["patient_id", batch_col]).T

    batch_train = train_df[batch_col].values
    batch_test  = test_df[batch_col].values

    train_corr = pycombat(train_expr, batch_train).T
    test_corr  = pycombat(test_expr, batch_test).T

    # Reconstruir DataFrames con índices jerárquicos
    train_corr.index = pd.MultiIndex.from_arrays(
        [train_df["patient_id"], train_df[batch_col]], names=["patient_id", batch_col]
    )
    train_corr.columns = train_expr.index

    test_corr.index = pd.MultiIndex.from_arrays(
        [test_df["patient_id"], test_df[batch_col]], names=["patient_id", batch_col]
    )
    test_corr.columns = test_expr.index

    return train_corr, test_corr


def plot_pca(before, after, batch, title="Train", cmap="Dark2"):
    """Genera gráfico PCA antes y después de ComBat."""
    pca = PCA(n_components=2, random_state=SEED)
    X_pca_before = pca.fit_transform(before)
    X_pca_after  = pca.fit_transform(after)

    categories = pd.Categorical(batch).categories
    colors = plt.cm.get_cmap(cmap, len(categories))
    legend_patches = [mpatches.Patch(color=colors(i), label=cat) for i, cat in enumerate(categories)]

    fig, ax = plt.subplots(1, 2, figsize=(12,5))

    ax[0].scatter(X_pca_before[:,0], X_pca_before[:,1], 
                  c=pd.Categorical(batch).codes, cmap=cmap)
    ax[0].set_title(f"{title} - Antes de ComBat")
    ax[0].set_xlabel("PC1"); ax[0].set_ylabel("PC2")
    ax[0].legend(handles=legend_patches, title="Batch")

    ax[1].scatter(X_pca_after[:,0], X_pca_after[:,1], 
                  c=pd.Categorical(batch).codes, cmap=cmap)
    ax[1].set_title(f"{title} - Después de ComBat")
    ax[1].set_xlabel("PC1"); ax[1].set_ylabel("PC2")

    plt.show()


def batch_r2(X, batch):
    """Calcula R² de batch en cada componente."""
    r2_scores = []
    for i in range(X.shape[1]):
        model = LinearRegression().fit(pd.get_dummies(batch), X[:, i])
        r2_scores.append(model.score(pd.get_dummies(batch), X[:, i]))
    return np.array(r2_scores)


def batch_distance_matrix(X, batch, metric="euclidean"):
    """Calcula matriz de distancias entre centroides de batch."""
    centroids = X.groupby(batch).mean()
    return pd.DataFrame(
        squareform(pdist(centroids, metric=metric)),
        index=centroids.index,
        columns=centroids.index
    )


# ------------------------
# Scrip Principal con pyCombat
# ------------------------
if __name__ == "__main__":
    # Cargar datos
    train = pd.read_csv(r"H:\TFE-TFM (Trabajo Fin de Estudios)/train_rna_ts_scaled.csv")
    test  = pd.read_csv(r"H:\TFE-TFM (Trabajo Fin de Estudios)/test_rna_scaled_export.csv")

    # Añadir columna de DB
    train = add_database_column(train)

    # Guardar copias antes
    RNA_train_before  = train.drop(columns=["Data_base"]).set_index(["patient_id", "Grupo"])
    RNA_test_before = test.set_index(["patient_id", "Grupo"])

    # Aplicar ComBat
    train_corr, test_corr = apply_combat(train, test)

    # Graficar PCA
    plot_pca(RNA_train_before.values, train_corr.values, train["Grupo"], title="Train", cmap="Dark2")
    plot_pca(RNA_test_before.values, test_corr.values, test["Grupo"], title="Test", cmap="autumn")

    # R²
    r2_train_before = batch_r2(RNA_train_before.values, train["Grupo"])
    r2_train_after  = batch_r2(train_corr.values, train["Grupo"])
    #---------------------------------------------------------------
    r2_before_test = batch_r2(RNA_test_before.values, test["Grupo"])
    r2_after_test  = batch_r2(test_corr.values, test["Grupo"])
    
    print("R² antes:", batch_r2(RNA_train_before.values, train["Grupo"]))
    print("R² después:", batch_r2(train_corr.values, train["Grupo"]))

    # Distancias entre centroides
    dist_train_before = batch_distance_matrix(RNA_train_before, train["Grupo"].values)
    dist_train_after  = batch_distance_matrix(train_corr, train["Grupo"].values)
    dist_train_before.to_csv("Distancias_train_centroides_preComBat.csv")
    dist_train_after.to_csv("Distancias_train_centroides_postComBat.csv")
    #-------------------------------------------------------------------
    dist_test_before = batch_distance_matrix(RNA_test_before, test["Grupo"].values)
    dist_test_after  = batch_distance_matrix(train_corr, train["Grupo"].values)
    dist_test_before.to_csv("Distancias_test_centroides_preComBat.csv")
    dist_test_after.to_csv("Distancias_test_centroides_postComBat.csv")

    # Guardar matrices corregidas
    train_corr.to_csv("RNA_train_batchCombat_corrected.csv")
    test_corr.to_csv("RNA_test_batchCombat_corrected.csv")
    print("Matrices corregidas guardadas.")

    # ------------------------
    # NeuroCombat (opcional)
    # ------------------------
    covars = pd.DataFrame({"batch": train["Grupo"]})
    RNA_expr = train.drop(columns=["patient_id", "Grupo", "Data_base"]).T
    RNA_corrected, _ = neuroCombat(dat=RNA_expr, covars=covars, batch_col="batch")
    RNA_corrected = RNA_corrected.T
