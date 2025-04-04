annotators folder: Contains latest code for developping ECG annotators

ludb_set.RData: .RData file containing a list of all 200 LUDB samples. Contains 12 lead of ECG signal and accompanying 12 annotation leads

create_ludb_set.R: script to create the ludb_set.RData. Largely for reference. 

annotator_prep_functions.R: Contains functions for creating training and testing signal and annotation matrices, ready for ML input
    Input: 
        preferred lead (integer of value 1 thru 12)
        annotator style (see file, may leave blank)
        split (% samples in training set, may leave blank)
    Output:
        list of:
            training/testing ECG signal / annotations, in correct format for ML input
            training/testing samples: LUDB samples used in each set. Largely for reference, usually not needed
