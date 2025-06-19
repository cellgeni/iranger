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