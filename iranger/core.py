import getpass
import json
import os
from importlib.metadata import version

import anndata
import h5py
import irods.client_configuration
import matplotlib.image
import numpy as np
import pandas as pd
import scipy.sparse
from irods.column import Criterion
from irods.models import Collection, CollectionMeta
from irods.session import iRODSSession

# required asof v1.1.9 for auto-close data objects if reading
irods.client_configuration.data_objects.auto_close = True


class iRanger(object):
    """
    A class for interacting with iRODS data collections, specifically designed
    to load and 10x tools outputs into AnnData objects.
    """

    def __init__(self, irods_environment=None, password=None, password_file=None, verbose=False):
        """
        Initializes an iRanger instance.

        Args:
            irods_environment (str): Path to the iRODS environment JSON file.
            password (str): Password string for iRODS login.
            password_file (str): Path to a file containing the iRODS password.
            verbose (bool): Flag to enable verbose logging.
        """
        self.verbose = verbose

        # read irods_environment
        if irods_environment is None:
            default_irods_environment = os.path.expanduser("~/.irods/irods_environment.json")
            self.log(f"No environment provided, using default irods environment file {default_irods_environment}", True)
            self.irods_environment = default_irods_environment
        else:
            custom_irods_environment = os.path.expanduser(irods_environment)
            self.log(f"Using custom irods environment file {custom_irods_environment}", True)
            self.irods_environment = custom_irods_environment
            self.log(f"Using {irods_environment}", True)

        if not os.path.exists(self.irods_environment):
            raise FileNotFoundError(f"iRODS environment file not found: {self.irods_environment}")

        # read password or password_file
        if password is None and password_file is None:
            self.log("No password provided, manual input required", True)
            self.password = self._get_password()
        elif password is not None:
            self.log("Using provided password", True)
            self.password = password
        elif password_file is not None:
            self.log("Using provided password_file", True)
            password_file = os.path.expanduser(password_file)
            with open(password_file, mode="rt") as pf:
                self.password = pf.read().strip()

        # smoke test to make sure connection works
        self.check_connection()

    def log(self, message, debug_message=False):
        """Prints log messages if verbose."""
        if self.verbose:
            print(message)
        elif not debug_message:
            print(message)

    def check_connection(self):
        """
        Performs a connection test to the iRODS server using current credentials to make sure they are working.
        """
        self.log("Checking iRODS connecting", True)
        with self._get_session() as session:
            self.log(
                f"Connected to host '{session.host}:{session.port}' server version {'.'.join(map(str, session.server_version))} as user '{session.username}'"
            )

    def _get_password(self):
        """
        Prompts the user for their password if none is provided.

        Returns:
            str: The password entered by the user.
        """
        return getpass.getpass(prompt="Sanger password: ")

    def _get_session(self):
        """
        Creates a new iRODS session with the provided environment and credentials.

        Returns:
            iRODSSession: The iRODS session object.
        """
        return iRODSSession(irods_env_file=self.irods_environment, password=self.password)

    def _read_feature_matrix_as_anndata(self, collection_path,count_file="filtered_feature_bc_matrix.h5"):
        """
        Reads a filtered_feature_bc_matrix.h5 file from a given 10x output in an iRODS collection.

        Args:
            collection_path (str): The iRODS collection path.

        Returns:
            anndata.AnnData: AnnData object with the parsed single-cell data.
        """
        self.log(f"Redaing feature matrix as AnnData", True)
        with self._get_session() as session:
            filtered_matrix_h5 = os.path.join(collection_path, count_file)
            self.log(f"Redaing {filtered_matrix_h5}", True)
            if not session.data_objects.exists(filtered_matrix_h5):
                raise FileNotFoundError(f"Missing '{filtered_matrix_h5}'")
            filtered_feature_bc_matrix = session.data_objects.get(filtered_matrix_h5)
            self.log(f"Retrieving {filtered_feature_bc_matrix}")
            with filtered_feature_bc_matrix.open("r") as h5:
                with h5py.File(h5, "r") as f:
                    self.log("Reading matrix elements from h5", True)
                    mtx_group = f["/matrix"]
                    # extract data from group
                    data = mtx_group["data"][:]
                    indices = mtx_group["indices"][:]
                    indptr = mtx_group["indptr"][:]
                    barcodes = mtx_group["barcodes"][:].astype(str)
                    features = mtx_group["features"]
                    gene_names = features["name"][:].astype(str)
                    self.log("Building CSR matrix", True)
                    mtx = scipy.sparse.csr_matrix((data, indices, indptr), shape=(len(barcodes), len(gene_names)))
                    self.log("Reading var elements", True)
                    var = dict()
                    var["var_names"] = gene_names
                    var["feature_types"] = features["feature_type"][:].astype(str)
                    if len(set(var["feature_types"])) != 1:
                        self.log(
                            f"Multiple feature_types. You may want to filter them to only have `var.feature_types == 'Gene Expression'`."
                        )

                    if "gene_id" not in features:
                        self.log("'gene_id' not in features, using 'id' as 'gene_id'", True)
                        var["gene_ids"] = features["id"][:].astype(str)
                    else:
                        self.log("using 'gene_id' as 'gene_id'", True)
                        var["gene_ids"] = (features["gene_id"].astype(str),)
                        self.log("using 'id' as 'probe_ids'", True)
                        var["probe_ids"] = features["id"].astype(str)

                    if "interval" in features:
                        self.log("found 'interval' in features (likely ATAC data)", True)
                        var["interval"] = np.array(features["interval"]).astype(str)

                    self.log("Creating AnnData file")
                    adata = anndata.AnnData(X=mtx, obs={"obs_names": barcodes}, var=var)

                    self.log("Adding irods metadata to uns", True)
                    adata.uns["_irods"] = {"path": collection_path, "metadata": []}
                    for item in session.collections.get(collection_path).metadata.items():
                        adata.uns["_irods"]["metadata"].append({item.name: item.value})

                    self.log("Adding package version information to uns", True)
                    adata.uns["_iranger_versions"] = {
                        "iranger": version("iranger"),
                        "anndata": version("anndata"),
                        "h5py": version("h5py"),
                        "pandas": version("pandas"),
                    }

                    self.log("Adding h5 metadata to uns", True)
                    adata.uns["_h5_metadata"] = dict(f.attrs)
                    return adata

    def read(self, collection_path,count_file="filtered_feature_bc_matrix.h5"):
        """
        Generic reader that detects the type of 10x analysis and calls the appropriate object reader.

        Args:
            collection_path (str): The iRODS collection path.

        Returns:
            anndata.AnnData: AnnData object with the data read from the specified collection.
        """
        self.log(f"Reading {collection_path}", True)
        # check that collection and feature matrix exist in the give path
        with self._get_session() as session:
            if not session.collections.exists(collection_path):
                raise FileNotFoundError(f"Collection '{collection_path}' not found")
            # sniff 'analysis_type' metadata
            self.log(f"Looking for analysis_type metadata to decide what format to read", True)
            collection = session.collections.get(collection_path)
            analysis_type = [m.value for m in collection.metadata.items() if m.name == "analysis_type"]
            if not analysis_type:
                raise KeyError(
                    "No 'analysis_type' metadata in iRODS collection. Use specific function to read the collection. For example: read_spaceranger('/seq/path/to/spaceranger')"
                )

        self.log(f"Detected analysis_type={analysis_type}")
        if "cellranger count" in analysis_type:
            return self.read_cellranger(collection_path,count_file)
        elif "spaceranger count" in analysis_type:
            return self.read_spaceranger(collection_path,count_file)
        elif "cellranger-arc count" in analysis_type:
            return self.read_cellranger_arc(collection_path,count_file)
        else:
            raise TypeError(
                f"Couldn't detect 10x analysis '{analysis_type}'. Use specific function to read the collection. For example: read_spaceranger(...)"
            )

    def read_cellranger(self, collection_path,count_file="filtered_feature_bc_matrix.h5"):
        """
        Reads a cellranger output.

        Args:
            collection_path (str): The iRODS collection path.

        Returns:
            anndata.AnnData: AnnData object with gene expression data.
        """
        adata = self._read_feature_matrix_as_anndata(collection_path,count_file)
        return adata

    def read_spaceranger(self, collection_path,count_file="filtered_feature_bc_matrix.h5"):
        """
        Reads a spaceranger output including spatial metadata and images.

        Args:
            collection_path (str): The iRODS collection path.

        Returns:
            anndata.AnnData: AnnData object with spatial transcriptomics data.
        """
        adata = self._read_feature_matrix_as_anndata(collection_path,count_file)
        adata.uns["spatial"] = dict()
        library_id = str(adata.uns["_h5_metadata"].pop("library_ids")[0], "utf-8")
        adata.uns["spatial"][library_id] = dict()

        with self._get_session() as session:
            tissue_possitions = None
            scale_factors = None
            hires_image = None
            lowres_image = None

            adata.uns["spatial"][library_id]["images"] = dict()

            # read images
            hires_image_path = os.path.join(collection_path, "spatial", "tissue_hires_image.png")
            if session.data_objects.exists(hires_image_path):
                hires_image = session.data_objects.get(hires_image_path)
                self.log(f"Reading {hires_image}", True)
                with hires_image.open("r") as hires:
                    adata.uns["spatial"][library_id]["images"]["hires"] = matplotlib.image.imread(hires)

            lowres_image_path = os.path.join(collection_path, "spatial", "tissue_lowres_image.png")
            if session.data_objects.exists(lowres_image_path):
                lowres_image = session.data_objects.get(lowres_image_path)
                self.log(f"Reading {lowres_image}", True)
                with lowres_image.open("r") as lowres:
                    adata.uns["spatial"][library_id]["images"]["lowres"] = matplotlib.image.imread(lowres)

            # read scale factors
            scale_factors = session.data_objects.get(os.path.join(collection_path, "spatial", "scalefactors_json.json"))
            self.log(f"Reading {scale_factors}", True)
            with scale_factors.open("r") as sf:
                adata.uns["spatial"][library_id]["scalefactors"] = json.loads(sf.read())

            # read coordinates
            if session.data_objects.exists(os.path.join(collection_path, "spatial", "tissue_positions.csv")):
                self.log(f"Found tissue_positions.csv", True)
                tpos = os.path.join(collection_path, "spatial", "tissue_positions.csv")
                tissue_possitions = session.data_objects.get(tpos)
            else:
                self.log(f"Found tissue_positions_list.csv", True)
                tpos = os.path.join(collection_path, "spatial", "tissue_positions_list.csv")
                tissue_possitions = session.data_objects.get(tpos)

            self.log(f"Reading {tissue_possitions}", True)
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

    def read_cellranger_arc(self, collection_path,count_file="filtered_feature_bc_matrix.h5"):
        """
        Reads cellranger-arc output including peak annotation from iRODS into an AnnData object.

        Args:
            collection_path (str): The iRODS collection path.

        Returns:
            anndata.AnnData: AnnData object with ATAC annotations.
        """
        adata = self._read_feature_matrix_as_anndata(collection_path,count_file)
        if "atac" not in adata.uns:
            adata.uns["atac"] = dict()

        peak_annotation_tsv = os.path.join(collection_path, "atac_peak_annotation.tsv")

        with self._get_session() as session:
            if session.data_objects.exists(peak_annotation_tsv):
                self.log("Parsing peak annotation file")
                peak_annotation = session.data_objects.get(peak_annotation_tsv)
                with peak_annotation.open("r") as pa:
                    adata.uns["atac"]["peak_annotation"] = pd.read_csv(pa, sep="\t")

            return adata

    def find(self, sample=None):
        """
        Searches iRODS collections by the 'sample' metadata field.

        Args:
            sample (str): Sample name to search for.

        Returns:
            list: A list of dictionaries containing paths and metadata of matched hits.
        """
        self.log(f"Searching: imeta qu -z seq -C sample = '{sample}'")
        with self._get_session() as session:
            results = []
            query = (
                session.query(Collection, CollectionMeta)
                .filter(Criterion("=", CollectionMeta.name, "sample"))
                .filter(Criterion("=", CollectionMeta.value, sample))
                .add_keyword("zone", "seq")
            )
            metadata = []
            for collection in query:
                collection_path = collection[Collection.name]
                for name, value in session.collections.get(collection_path).metadata.items():
                    metadata.append({name: value})
                results.append(
                    {
                        "path": collection_path,
                        "metadata": metadata,
                    }
                )
            return results
