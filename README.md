# iRODS ranger

## Install
```bash
pip install git+https://github.com/cellgeni/iranger.git
```

## Usage

```python
import iranger

# You can use your password directly
ir = iranger.setup(password='Your_iRODS_Password_Goes_Here')
# or alternatively put your password in a file and use that instead
ir = iranger.setup(password_file='~/my_password_file')

# 10x cellranger output
adata = ir.get_cellranger("/seq/26280/cellranger/cellranger302_count_26280_FCAImmP7555847_GRCh38-1_2_0")

# 10x spaceranger output
adata = ir.get_spaceranger("/seq/illumina/spaceranger/spaceranger130_count__WSSKNKCLsp12887269_GRCh38-2020-A")
print(adata)

```

### Embeded iRODS metadata

AnnData files will have an `_irods` key inside unstructured observations with the irods path and all the medata for the collection.

```python
>>> adata.uns['_irods']
{
  'path': '/seq/26280/cellranger/cellranger302_count_26280_FCAImmP7555847_GRCh38-1_2_0',
  'metadata': [{'name': 'library_type', 'value': 'Chromium single cell'}, {'name': 'study', 'value': 'FCA_ImmunoP'}, {'name': 'study_accession_number', 'value': 'EGAS00001002715'}, {'name': 'study_id', 'value': '5061'}, {'name': 'study_title', 'value': 'FCA_ImmunoP'}, {'name': 'id_run', 'value': '26280'}, {'name': 'sample', 'value': 'FCAImmP7555847'}, {'name': 'sample_id', 'value': '3775200'}, {'name': 'analysis_type', 'value': 'cellranger count'}, {'name': '10x:reference', 'value': '/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-GRCh38-1.2.0'}, {'name': '10x:pipeline', 'value': '/software/sciops/external/cellranger/3.0.2/cellranger'}, {'name': 'sample_lims', 'value': 'SQSCP'}, {'name': 'sample_uuid', 'value': '8b7789d6-7abf-11e8-8cef-68b599768938'}]
}
```

### iRODS environment file

The package requires an `irods_environment.json` file present in the machine. The standard location is:
```bash
~/.irods/irods_environment.json
```

Should you whish to have it in a different location you'll have to be explicit about it, for example:

```python
ir = iranger.setup(
    irods_environment='/path/to/my/irods_environment.json',
    password_file='/safe/path/to/my_password_file'
)
```
