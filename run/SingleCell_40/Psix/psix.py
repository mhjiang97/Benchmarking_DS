import psix
psix_object = psix.Psix()
psix_object.junctions2psi(
        sj_dir='~/scRNA_AS/SJout/v2/',
        intron_file='psix_annotation.tab',
        tpm_file='~/gene_tpm',
        save_files_in='psix_output/'
        )
psix_object.run_psix(latent='~/gene_pca1_2', 
                     n_random_exons=2000, 
                     n_neighbors=10,
                     n_jobs=10,
                     )
psix_object.save_psix_object('psix_output/')