![Versions](https://img.shields.io/badge/python-3.10_%7C_3.11_%7C_3.12-brightgreen)
![License](https://img.shields.io/badge/License-Apache-brightgreen)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/xomicsdatascience/full_dia/test.yml?branch=main&label=Unit%20Tests)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/xomicsdatascience/full_dia/publish.yml?branch=main&label=Deploy%20PyPi)

# Full-DIA

Full-DIA, a freely available software for single-cell diaPASEF data analysis that 
leverages deep learning to improve proteome coverage, 
quantitative accuracy and analysis speed. Most notably, 
Full-DIA is the first to automatically generate a missing-value-free protein matrix 
under global FDR control, which may offer superior biological interpretability and 
insight into single-cell proteomics data 
compared to conventional matrices with missing values.

---
### Contents
**[Installation](#installation)**<br>
**[Usage](#usage)**<br>
**[Output](#output)**<br>

---
### Installation

We recommend using [Conda](https://www.anaconda.com/) to create a Python environment for using Full-DIA, whether on Windows or Linux. Here are the updated installation instructions:

1. Create a Python environment with version 3.9.18.
```bash
# Create and activate environment
conda create -n full_dia_env python=3.12 -y
conda activate full_dia_env

# Install PyTorch with CUDA 12.1 runtime
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# Install NumPy 1.x and CuPy 13.x (Strictly avoid NumPy 2.x to maintain ABI compatibility)
pip install "numpy<2.0.0" "cupy-cuda12x<14.0.0" --no-cache-dir

# Install full_dia package
pip install full-dia==1.0.2
```

2. Applying Compatibility Patches full_dia and its bundled dependency alphatims require three patches to resolve issues introduced by Pandas 2.0+ (Copy-on-Write) and read-only NumPy array restrictions. Run the following Python script inside the activated full_dia_env to apply all patches automatically:

```bash
import os
import re
import site

site_packages = site.getsitepackages()[0]

# Patch 1 & 2: Patch full_dia/tims.py for Pandas 2.0+ reset_index and read-only array modification
tims_path = os.path.join(site_packages, "full_dia", "tims.py")
if os.path.exists(tims_path):
    with open(tims_path, "r") as f:
        code = f.read()

    # Preserve multi-index grouping columns during reset_index in get_dia_windows
    code = code.replace(".reset_index(drop=True)", ".reset_index()")

    # Make array writable before slice assignment in get_dia_quadrupole
    target = "low[1:] = (low[1:] + high[:-1]) / 2"
    replacement = "low = low.copy()\n        low[1:] = (low[1:] + high[:-1]) / 2"
    if target in code:
        code = code.replace(target, replacement)

    with open(tims_path, "w") as f:
        f.write(code)
    print("Successfully patched full_dia/tims.py")

# Patch 3: Patch alphatims/bruker.py for Pandas 2.0+ ChainedAssignmentError on Frame 0
bruker_path = os.path.join(site_packages, "full_dia", "alphatims", "bruker.py")
if os.path.exists(bruker_path):
    with open(bruker_path, "r") as f:
        code = f.read()

    pattern = r'frames\.([A-Za-z0-9_]+)\[0\]\s*=\s*0'
    replacement = r"frames.loc[0, '\1'] = 0"
    new_code, count = re.subn(pattern, replacement, code)

    if count > 0:
        with open(bruker_path, "w") as f:
            f.write(new_code)
        print(f"Successfully patched {count} chained assignment calls in alphatims/bruker.py")
```

All three source code patches applied during this session are included in the README script:

full_dia/tims.py – Changed .reset_index(drop=True) to .reset_index() to prevent Pandas 2.0+ from dropping the quadrupole m/z index columns during groupby().apply().

full_dia/tims.py – Added low = low.copy() in get_dia_quadrupole() to resolve the NumPy read-only array assignment error.

full_dia/alphatims/bruker.py – Converted chained assignments (frames.col[0] = 0) to .loc[0, col] = 0 to fix frame 0 metadata zeroing under Pandas 2.0+ Copy-on-Write.

The other errors encountered were resolved through dependency pinning and workspace file management rather than source code edits:

NumPy 2.x np.searchsorted error (ValueError: search side must be one of 'left' or 'right'): Resolved by pinning numpy<2.0.0 in the environment rather than patching alphatims.

CuPy C-extension import error: Resolved by reinstalling cupy-cuda12x<14.0.0 compiled against the NumPy 1.x ABI.

Missing quad_low_mz_values on blank runs: Resolved by removing non-diaPASEF/blank .d folders from the -ws directory.

---
### Usage

Convert library 

```
diann.exe --lib  --threads 8 --verbose 3 --out  --qvalue 0.01 --matrices --out-lib L:\promec\FastaDB\humanMoxCa\lib.parquet --gen-spec-lib --predictor --fasta camprotR_240512_cRAP_20190401_full_tags.fasta --cont-quant-exclude cRAP- --fasta L:\promec\FastaDB\UP000005640_9606.fasta --fasta-search --min-pep-len 7 --max-pep-len 30 --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 1 --max-pr-charge 4 --min-fr-mz 200 --max-fr-mz 1800 --cut K*,R* --missed-cleavages 1 --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --no-prot-inf --rt-profiling --original-mods
diann.exe --lib L:\promec\FastaDB\humanMoxCa\lib.predicted.speclib --threads 24 --verbose 1 --out  --qvalue 0.05 --matrices --out-lib L:\promec\FastaDB\humanMoxCa\humanMoxCa.parquet --gen-spec-lib --fasta camprotR_240512_cRAP_20190401_full_tags.fasta --cont-quant-exclude cRAP- --min-pep-len 7 --max-pep-len 30 --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 1 --max-pr-charge 4 --min-fr-mz 200 --max-fr-mz 1800 --cut K*,R* --missed-cleavages 1 --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --peptidoforms --reanalyse --rt-profiling 
#diann.exe --lib  --threads 16 --verbose 1 --out F:\promec\FastaDB\humanMoxCa.parquet --qvalue 0.01 --matrices --out-lib F:\promec\FastaDB\humanMoxCa.tmp.parquet --gen-spec-lib --predictor --reannotate --fasta camprotR_240512_cRAP_20190401_full_tags.fasta --cont-quant-exclude cRAP- --fasta F:\promec\FastaDB\UP000005640_9606.fasta --fasta-search --min-fr-mz 200 --max-fr-mz 1800 --min-pep-len 7 --max-pep-len 30 --min-pr-mz 300 --max-pr-mz 1800 --min-pr-charge 1 --max-pr-charge 4 --cut K*,R* --missed-cleavages 1 --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --reanalyse --rt-profiling
#diann.exe --lib "F:\promec\FastaDB\humanMoxCa.tmp.predicted.speclib" --threads 32 --verbose 1 --out "F:\promec\FastaDB\humanMoxCa.main.parquet" --qvalue 0.01 --matrices  --out-lib "F:\promec\FastaDB\humanMoxCa.parquet" --gen-spec-lib --reannotate --fasta camprotR_240512_cRAP_20190401_full_tags.fasta --cont-quant-exclude cRAP- --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --reanalyse --rt-profiling 
#diann.exe --lib "F:\promec\FastaDB\humanMC2V3defaults.predicted.speclib" --threads 32 --verbose 1 --out "F:\promec\FastaDB\report.parquet" --qvalue 0.01 --matrices  --out-lib "F:\promec\FastaDB\humanMC2defaults.parquet.tmp.parquet" --gen-spec-lib --unimod4 --reanalyse --rt-profiling 
copy F:\promec\FastaDB\humanMoxCa.parquet F:/lib/humanMoxCa.parquet
```

Download [test data](https://fuzzylife.substack.com/p/proteomics-data-processing-with-maxquant) 
```
cd /mnt/f/timsTOF
wget https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/cts8a50sb36put8/26june24_hel50_100spd_OT_1ulirt_S2-D2_1_6366.d.zip
wget https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/cts8a50sb36put8/26june24_hel50_100spd_OT_1ulirt_S2-C2_1_6365.d.zip
unzip *.zip
```

```bash
full_dia -lib /mnt/f/lib/humanMoxCa.parquet -ws /mnt/f/timsTOF
```
(Please note that the path needs to be enclosed in quotes if running on a Windows platform.)

- `-lib`<br>
This parameter is used to specify the absolute path of the spectral library.
Full-DIA currently supports spectral libraries with the ***.parquet*** or ***.tsv*** suffix, provided that their column names are consistent with those of the DIA-NN (> v1.9) predicted spectral library. 
We recommend generating the predicted spectral library using DIA-NN and then converting it to the .parquet format.
Refer to [this](https://github.com/vdemichev/DiaNN) for instructions on how to generate prediction spectral libraries and convert to .parquet format using DIA-NN.
Full-DIA supports oxygen modifications on methionine (M) but does not include modifications such as phosphorylation or acetylation.
Full-DIA will develop its own predictor capable of forecasting the peptide retention time, ion mobility, and fragmentation pattern. 
It may also be compatible with other formats of spectral libraries based on requests.

- `-ws`<br>
This parameter specifies the folder that contains multiple .d directories to be analyzed.

Other optional params are list below by entering `full_dia -h`:
```
       ******************
       * Full-DIA x.y.z *
       ******************
Usage: full_dia -ws WS -lib LIB

optional arguments for users:
  -h, --help           Show this help message and exit.
  -ws WS               Specify the folder that is .d or contains .d files.
  -lib LIB             Specify the absolute path of a .speclib or .parquet spectra library.
  -out_name OUT_NAME   Specify the folder name of outputs. Default: full_dia.
  -gpu_id GPU_ID       Specify the GPU-ID (e.g. 0, 1, 2) which will be used. Default: 0.
```

### Output
Full-DIA will generate **`report.log.txt`** and **`report.parquet`** in output folder. 
The report.parquet contains precursor and protein IDs, as well as plenty of associated information. 
Most column names are consistent with DIA-NN and are self-explanatory.

* **Protein.Group** - inferred proteins. Full-DIA uses [IDPicker](https://pubs.acs.org/doi/abs/10.1021/pr070230d) algorithm to infer proteins. 
* **Protein.Ids** - all proteins matched to the precursor in the library.
* **Protein.Names** - names (UniProt names) of the proteins in the Protein.Group.
* **PG.Quantity.Raw** - raw quantity of the Protein.Group.
* **PG.Quantity.Deep** - corrected quantity of the Protein.Group.
* **Precursor.Id** - peptide seq + precursor charge.
* **Precursor.Charge** - the charge of the precursor.
* **Q.Value** - run-specific precursor q-value.
* **Global.Q.Value** - global precursor q-value.
* **PG.Q.Value** - run-specific q-value for the protein group.
* **Global.PG.Q.Value** - global q-value for the protein group.
* **Proteotypic** - indicates the peptide is specific to a protein.
* **Precursor.Quantity.Raw** - raw quantity of the precursor.
* **Precursor.Quantity.Deep** - corrected quantity of the precursor.
* **RT** - the retention time of the precursor.
* **IM** - the ion mobility of the precursor.

---
## Troubleshooting
- Please create a GitHub issue and we will respond as soon as possible.
- Email: songjian2022@suda.edu.cn

---