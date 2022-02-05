# A Deep Learning Framework for Epileptic Seizure Detection based on Neonatal EEG Signals

*Artur Gramacki & Jarosław Gramacki*

A version that only calculates one selected case `expert_A_1sec_1chunk_64Hz_fold_0`. 
Can be used to test codes quickly. Just execute `EEG_neonatal.ipynb` Jupyter notebook. 
In the [working/inputs](working/inputs) dir `expert_A_1sec_1chunk_64Hz.hdf5` file has been uploaded, so
Jupyter notebook can work properly.

If you want to calculated other cases, HDF5 files must be generated and uploaded to
the `working/inputs` dir. This task has been implemented in R. Follow the instructions 
in the paper (*Replication of the results* section). Read also the `R/Tutorial_on_how_to_run_R_codes.html` tutorial.

In the `working/inputs/illustrative_files` you can find txt versions of the 
`expert_A_1sec_1chunk_64Hz.hdf5` file (just for your information, txt files are not 
used during calculations).

The HDF5 files generated in Step 5 (see ‘Replication of the results’ Section) are very large, 
as much as 16.6GB. Therefore these files have been made available in separate zip archives.
See `working_inputs_dir.txt` file.

