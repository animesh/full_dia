# Library
Full-DIA will generate **`report.log.txt`** and **`report.parquet`** in output folder. 
The report.parquet contains precursor and protein IDs, as well as plenty of associated information. 
Most column names are consistent with DIA-NN and are self-explanatory.

* **Protein.Group** - inferred proteins. Full-DIA uses [IDPicker](https://pubs.acs.org/doi/abs/10.1021/pr070230d) algorithm to infer proteins. 
* **Protein.Ids** - all proteins matched to the precursor in the library.
* **Protein.Names** names (UniProt names) of the proteins in the Protein.Group.
* **PG.Quantity.Raw** raw quantity of the Protein.Group.
* **PG.Quantity.Deep** corrected quantity of the Protein.Group.
* **Precursor.Id** peptide seq + precursor charge.
* **Precursor.Charge** the charge of precursor.
* **Q.Value** run-specific precursor q-value.
* **Global.Q.Value** global precursor q-value.
* **PG.Q.Value** run-specific q-value for the protein group.
* **Global.PG.Q.Value** global q-value for the protein group.
* **Proteotypic** indicates the peptide is specific to a protein.
* **Precursor.Quantity.Raw** raw quantity of the precursor.
* **Precursor.Quantity.Deep** corrected quantity of the precursor.
* **RT** the retention time of the precursor.
* **IM** the ion mobility of the precursor.
