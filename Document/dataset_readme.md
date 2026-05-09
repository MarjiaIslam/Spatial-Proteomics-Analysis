# Additional data for Dayao et al. 'Deriving spatial features from in situ proteomics imaging to enhance cancer survival analysis'

See below for a description of each file

## `dataset_info.tar.gz` 

Extracting this tar file creates the `dataset_info` directory. This directory contains the following files:

- `cell_locations_and_labels.csv` CSV file where each row gives the sample id (`acquisition_id`), cell id, x-y coordinates, cluster id, and cluster label for each cell in the dataset.
- `labeled_arcsinh_norm_data.pkl` pickle file (readable by `pandas`) which yields a dataframe giving the normalized biomarker expression for each cell. This dataframe also includes columns for sample id, cell id, cluster id, and cluster label.
- `marker_names.csv` A list of protein biomarkers used in the paper.
- `qc_acq_ids_labeled.csv` A list of sample ids (acquisition ids) that have passed quality control.
- `sample_metadata.csv` CSV file where each row gives the metadata for each sample (e.g. acquisition id, patient id, tissue type, HPV status, survival length, survival status). 

## `all_k_neighborhood_mats.tar.gz` 

Extracting this tar file creates 20 directories: from `k1` to `k20`. These directories contain the neighborhood matrix features.  

Each directory contains a numpy array for each sample (identified by its acquisition id). This array is the neighborhood matrix using the value of `k` specified by the corresponding directory name. 

## `biomarker_expr_summary.tar.gz` 

Extracting this tar file creates the `expression_biomarkers` directory. This directory contains a CSV file for each sample with the name `{acq_id}_cell_info.csv`. These files are used for the average biomarker across image features.  

Each CSV file has two rows: `summed` and `average`. The `summed` row gives the summed pixel intensity across the entire sample for each biomarker. The `average` row gives the average (sum / number of pixels) pixel intensity across the entire sample for each biomarker.

## `patwa_comparisons.tar.gz`

Extracting this tar file creates the `comparisons` directory. This directory contains the features computed for Patwa et al. comparison method. It contains the following files:

- `binary_biomarker_expr_cells_qc_labeled.csv` CSV file indicating the biomarker positivities for each cell in the dataset.
- `biomarker_frac_positivity_qc_labeled.csv` CSV file giving the fraction of cells in each sample that are positive for each biomarker.
- `interaction_biomarker_features.csv` CSV file giving the number of pairwise biomarker interactions in each sample.

## `denvar.tar.gz` 

Extracting this tar file creates the `denvar` directory. This directory contains the features computed for the DenVar comparison method. For each biomarker, there are two files in this directory:

- `{biomarker}_DenVar_clusters.csv` CSV file listing the assigned DenVar cluster for each sample for the specified biomarker.
- `{biomarker}_JSD.RDS` RDS file (R object) that gives the Jenson-Shannon distance matrix between all samples for the specified biomarker.

## `rsf_risk_scores.tar.gz` 

Extracting this tar file creates the `rsf_risk_scores` directory. This directory contains the following files:

- `sample_risk_info.csv` CSV file giving the survival length, survival status, and output RSF risk score for each sample.
- `patient_risk_info.csv` CSV file giving the survival length, survival status, output RSF risk score, and assigned risk cohort (`low` or `high`) for each patient.

## `k_fns_norm_by_uw_qc_labeled.csv` 

This CSV files contains the K function features for each biomarker and sample. Each row corresponds to a normalized K function (as used in the paper), with each column indicating the function value at that value of `r` (in pixels).  

