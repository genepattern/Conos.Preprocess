docker run --rm \
-w /job_1234 \
-v $PWD/job_1234:/job_1234 \
-v $PWD/data:/temp/data \
-it genepattern/conos:1.1 \
Rscript /module/run_conos.R \
--file_list '/temp/data/file_list.txt' \
--k 40 \
--perplexity 30 \
--pagoda_odgenes 6000 \
--projection_method 'PCA' \
--ncomps 100 \
--conos_odgenes 6000
