import os
from irods.session import iRODSSession
from irods.models import Collection, CollectionMeta
from irods.column import Criterion
import irods.client_configuration
import h5py
import anndata
import scipy.sparse
import getpass
import json
import pandas as pd

# required asof v1.1.9 for auto-close data objects if reading
irods.client_configuration.data_objects.auto_close = True


class iRanger(object):

    def __init__(self, irods_environment=None, password=None, password_file=None):
        if irods_environment is None:
            self.irods_environment = os.path.expanduser("~/.irods/irods_environment.json")
        if password is None and password_file is None:
            self.password = self._get_password()
        elif password is not None:
            self.password = password
        elif password_file is not None:
            password_file = os.path.expanduser(password_file)
            with open(password_file, mode="rt") as pf:
                self.password = pf.read().strip()
        self._check_connection()

    def _check_connection(self):
        with iRODSSession(irods_env_file=self.irods_environment, password=self.password) as session:
            print(f"Connected to host '{session.host}:{session.port}' server version {'.'.join(map(str,session.server_version))} as user '{session.username}'")

    def _get_password(self):
        return getpass.getpass(prompt="Sanger password: ")
    
    def _get_session(self):
        return iRODSSession(irods_env_file=self.irods_environment, password=self.password)

    def get_cellranger(self, collection_path):
        with self._get_session() as session:
            filtered_feature_bc_matrix = session.data_objects.get(
                os.path.join(collection_path, "filtered_feature_bc_matrix.h5")
            )
            print(f"Retrieving {filtered_feature_bc_matrix}")
            with filtered_feature_bc_matrix.open("r") as h5:
                with h5py.File(h5, "r") as f:
                    print("Reading matrix elements")
                    mtx_group = f["/matrix"]
                    # extract data from group
                    data = mtx_group["data"][:]
                    indices = mtx_group["indices"][:]
                    indptr = mtx_group["indptr"][:]
                    barcodes = mtx_group["barcodes"][:].astype(str)
                    features = mtx_group["features"]
                    gene_names = features["name"][:].astype(str)

                    mtx = scipy.sparse.csr_matrix((data, indices, indptr), shape=(len(barcodes), len(gene_names)))
                    
                    print("Reading var information")
                    var = {
                        "var_names": gene_names,
                        "feature_types": features["feature_type"][:].astype(str),
                    }
                    if "gene_id" not in features:
                        print("Read metadata specific to a feature-barcode matrix")
                        var["gene_ids"] = features["id"][:].astype(str)
                    else:
                        print("Read metadata specific to a probe-barcode matrix")
                        var.update(
                            {
                                "gene_ids": features["gene_id"].astype(str),
                                "probe_ids": features["id"].astype(str),
                            }
                        )
                    
                    print("Creating AnnData file")
                    adata = anndata.AnnData(X=mtx, obs={"obs_names": barcodes}, var=var)

                    print("Add irods metadata ('_irods') to unstructured observations")
                    adata.uns["_irods"] = {
                        "path": collection_path,
                        "metadata": [
                            {"name": meta_item.name, "value": meta_item.value}
                            for meta_item in session.collections.get(collection_path).metadata.items()
                        ],
                    }
                    print("Add h5 metadata ('_h5_metadata') to unstructured observations")
                    adata.uns["_h5_metadata"] = dict(f.attrs)
                    return adata

    def get_spaceranger(self, collection_path):
        adata = self.get_cellranger(collection_path)
        adata.uns["spatial"] = dict()
        library_id = str(adata.uns["_h5_metadata"].pop("library_ids")[0], "utf-8")
        adata.uns["spatial"][library_id] = dict()

        with self._get_session() as session:
            tissue_possitions = None
            scale_factors = None
            hires_image = None
            lowres_image = None
            if session.data_objects.exists(
                os.path.join(collection_path, "spatial", "tissue_positions.csv")
            ):
                tissue_possitions = session.data_objects.get(
                    os.path.join(collection_path, "spatial", "tissue_positions.csv")
                )
            else:
                tissue_possitions = session.data_objects.get(
                    os.path.join(collection_path, "spatial", "tissue_positions_list.csv")
                )
            scale_factors = session.data_objects.get(os.path.join(collection_path, "spatial", "scalefactors_json.json"))
            if session.data_objects.exists(os.path.join(collection_path, "spatial", "tissue_hires_image.png")):
                hires_image = session.data_objects.get(os.path.join(collection_path, "spatial", "tissue_hires_image.png"))
            if session.data_objects.exists(os.path.join(collection_path, "spatial", "tissue_lowres_image.png")):
                lowres_image = session.data_objects.get(os.path.join(collection_path, "spatial", "tissue_lowres_image.png"))
            
            adata.uns["spatial"][library_id]["images"] = dict()
            # with hires_image.open("r") as hires:
            #   adata.uns["spatial"][library_id]["images"]["hires"] = imread(hires)
            # with lowres_image.open("r") as lowres:
            with scale_factors.open("r") as sf:
                adata.uns["spatial"][library_id]["scalefactors"] = json.loads(sf.read())
            # read coordinates
            with tissue_possitions.open("r") as tp:
                positions = pd.read_csv(tp, header=0, index_col=0)
                positions.columns = [
                    "in_tissue",
                    "array_row",
                    "array_col",
                    "pxl_col_in_fullres",
                    "pxl_row_in_fullres",
                ]
            adata.obs = adata.obs.join(positions, how="left")
            adata.obsm["spatial"] = adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy()
            adata.obs.drop(columns=["pxl_row_in_fullres", "pxl_col_in_fullres"], inplace=True)
            adata.uns["spatial"][library_id]["metadata"] = {
                k: (
                    str(adata.uns["_h5_metadata"][k], "utf-8")
                    if isinstance(adata.uns["_h5_metadata"][k], bytes)
                    else adata.uns["_h5_metadata"][k]
                )
                for k in ("chemistry_description", "software_version")
                if k in adata.uns["_h5_metadata"]
            }
            return adata
