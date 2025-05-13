#!/usr/bin/env nextflow


process spatial_preprocess {

    publishDir '/Users/sarah/Documents/ICB/Panpipes/15.nextflow/preprocess_spatial/preprocessed.data', mode: 'copy'

    input:
        val sample
        val norm_hvg_flavour
        val n_top_genes
        val filter_by_hvg
        val hvg_batch_key
        val squidpy_hvg_flavour
        val min_mean
        val max_mean
        val min_disp
        val theta
        val clip 
        val n_pcs

    output:
        path "$sample-preprocessed.zarr"

    script:
    """
    python run_preprocess_spatial.py --input_spatialdata /Users/sarah/Documents/ICB/Panpipes/15.nextflow/preprocess_spatial/preprocessed.data/$sample-filtered.zarr \
            --output_spatialdata $sample-preprocessed.zarr \
            --figdir ./figures \
            --norm_hvg_flavour $norm_hvg_flavour --filter_by_hvg $filter_by_hvg --hvg_batch_key $hvg_batch_key --squidpy_hvg_flavour $squidpy_hvg_flavour \
            --min_mean $min_mean --max_mean $max_mean --min_disp $min_disp --theta $theta --clip $clip \
            --n_pcs $n_pcs

    """
}

process concatenate {

    publishDir '/Users/sarah/Documents/ICB/Panpipes/15.nextflow/preprocess_spatial/concatenated.data', mode: 'copy'

    input:
        path samples

    output:
        path "concatenated.zarr"

    script:
    """
    python concatenation_spatial.py 
    """
}




/*
 * Pipeline parameters
 */
params.sample= ["V1_Human_Lymph_Node", "V1_Human_Lymph_Node2"] /*change this*/

params.norm_hvg_flavour = "squidpy"
params.n_top_genes = 2000
params.filter_by_hvg= "False"
params.hvg_batch_key= "None" 
params.squidpy_hvg_flavour="seurat"
params.min_mean= 0.05
params.max_mean=1.5
params.min_disp=0.5
params.theta=100
params.clip="None" 
params.n_pcs=50

params.concat="True"


workflow {

    /* Run Preprocessing */
    samples = Channel.of(params.sample)
                         .flatten()

    spatial_preprocess(samples, params.norm_hvg_flavour, params.n_top_genes, params.filter_by_hvg,
                        params.hvg_batch_key, params.squidpy_hvg_flavour, params.min_mean,
                        params.max_mean, params.min_disp, params.theta, params.clip,
                        params.n_pcs)

    /* Run Concatenation */
    if (params.concat == 'True') {
        concatenate(spatial_preprocess.out.collect())
    }
}
