# Deep learning identifies antigenic determinants of severe SARS-CoV-2 infection within T-cell repertoires

## Prepare Data
Since the repertoire data is large, it is not included in this repository. To replicate the anlayses found here, follow the following instructions.
 
 - Download the full repertoire data [here](https://clients.adaptivebiotech.com/pub/covid-2020). There is a file (INCOV067-AC-3_TCRB.tsv) that is corrupted from the download link (as of the date this repository was created) that needs to be downloaded separately from the ImmuneAccess portal and used to replace the corrupted one. 
 - All downloaded repertoire data should be put into a directory labeled "repertoires_raw" underneath Data/ImmuneCODE directory.
 - Run reformat_files.py under scripts/preprocessing to reformat the repertoire data to inputs that DeepTCR can use. Doing this should create another directory underneath Data/ImmuneCODE called repertoires.
 - Run organize_files.py under scripts/preprocessing to sort the files into their respective cohorts. Doing this should create another directory underneath Data/ImmuneCODE called repertoires_org.
 
 This steps should prepare and format the data in a way that the other analysis scripts in the repository can use.
 
 

