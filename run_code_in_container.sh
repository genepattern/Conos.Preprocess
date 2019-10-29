docker build -t genepattern/conos:1.2 .

docker run --rm \
-v $PWD/data:/temp/data \
-v $PWD/job_1234:/job_1234 \
-w /job_1234 \
-it genepattern/conos:1.2 \
time Rscript /module/run_conos_preprocess.R \
--file_list '/temp/data/file_list.txt' --k 40 --perplexity 30 --pagoda_odgenes 6000 --projection_method 'PCA' --ncomps 100 --conos_odgenes 6000

docker run --rm \
-v $PWD/data:/temp/data \
-v $PWD/job_1234:/job_1234 \
-w /job_1234 \
-it genepattern/conos:1.2 \
time Rscript /module/run_conos_cluster.R \
--conos_object conos_object.rds --runleiden True --leiden_resolution 1.0 --runwalktrap True --walktrap_steps 10 --umap_distance 0.05 --umap_spread 5.0
