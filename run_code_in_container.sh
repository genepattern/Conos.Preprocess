docker run --rm \
-w /job_1234 \
-v $PWD/job_1234:/job_1234 \
-v $PWD/data:/temp/data \
-v $PWD/src:/module \
-it genepattern/conos:1.0 \
Rscript /module/run_conos.R --file_list '/temp/data/file_list.txt'
