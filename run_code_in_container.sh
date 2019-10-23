docker run --rm \
-v $PWD/data:/temp/data \
-v $PWD/job_1234:/job_1234 \
-w /job_1234 \
-it genepattern/conos:1.1 \
Rscript /module/run_conos.R \
--file_list '/temp/data/file_list.txt' --k 40 --perplexity 30 --pagoda_odgenes 6000 --projection_method 'PCA' --ncomps 100 --conos_odgenes 6000
