This is the Electronic Supplements for the paper:

> ## A Deep Learning Framework for Epileptic Seizure Detection based on Neonatal EEG Signals ##

by *Artur Gramacki & Jarosław Gramacki*

e-mails:  a.gramacki@issi.uz.zgora.pl, j.gramacki@ck.uz.zgora.pl

Paper details: [waiting for decision]

A complete repository consists of:
1. R and Python source files and complete output results obtained by the authors available at https://drive.google.com/file/d/1FwgR8GjZLwE3z8d36vL7XDAqloVtHgJa/view?usp=sharing (size about 770MB). 
2. Raw EDF files and CSV annotation files available at https://zenodo.org/record/4940267 (size about 4GB)

Download the repository and save all the files on your local or cloud-based disk. Remember not to change the directory structure as this is hard-coded in the R and Python codes. The HDF5 files generated in Step 5 (see ‘Replication of the results’ Section) are very large, as much as about 17GB. Therefore these files have been made available in separate zip archives. As it can be cumbersome to download a file that is 17GB in size, 15 smaller files were prepared. They are available for download at the links given in `working/inputs/__read_me__.txt`. Download these files and load them to the `working/inputs/` directory.

Then start by:
1. Reading a short tutorial in `R/Tutorial_on_how_to_run_R_codes.html`
2. Become familiar with the contents of the `Python/EEG_neonatal.ipynb` Jupyter notebook. 

Notes:
- To familiarize yourself with this notebook, we suggest that you first run it with the `COMPLETE_CALCULATIONS` global variable set to `False`. Leave the other global variables unchanged. 

- If you are working in [Google Colab](https://colab.research.google.com) set `GOOGLE_COLAB` global variable to `True`, otherwise set it to `False`. If `GOOGLE_COLAB = True` set also a correct value for `MOUNT_POINT` global variable.

- Two `HDF5` files (`working\inputs\expert_A_1sec_1chunk_64Hz.hdf5` and `working\inputs\expert_CC_10sec_20chunk_64Hz.hdf5`) are required and the files are already uploaded at the repository mentioned above. 

- Then you can safely run all the code cells (`Ctrl+F9` in [Google Colab](https://colab.research.google.com)).

In case of problems, the authors declare the necessary help for potential researchers.
