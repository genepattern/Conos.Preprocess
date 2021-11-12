FROM rocker/tidyverse:4.0.5

# Dockerfile adapted from kharchenkolab/conos/blob/master/docker/Dockerfile

LABEL authors="Viktor Petukhov <viktor.s.petuhov@ya.ru>, Evan Biederstedt <evan.biederstedt@gmail.com>" \
    version.image="1.4.4" \
    version.pagoda2="1.0.6" \
    description="tidyverse image R 4.0.5 to run pagoda2 with Rstudio"

RUN apt-get update --yes && apt-get install --no-install-recommends --yes \
  build-essential \
  cmake \
  git \
  less \
  libcurl4-openssl-dev \
  libssl-dev \
  libgsl0-dev \
  libeigen3-dev \
  libcairo2-dev \
  libxt-dev \
  libgtk2.0-dev \
  xvfb  \
  xauth \
  xfonts-base \
  libz-dev \
  libhdf5-dev \
  time

RUN R -e 'chooseCRANmirror(ind=1); install.packages("BiocManager")'

RUN R -e 'install.packages("optparse",repos = "http://cran.us.r-project.org")'

RUN R -e 'BiocManager::install(c("AnnotationDbi", "BiocGenerics", "GO.db", "pcaMethods", "org.Dr.eg.db", "org.Hs.eg.db", "org.Mm.eg.db", "scde", "BiocParallel"))'

RUN R -e "install.packages('Seurat',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "install.packages('entropy',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "install.packages('p2data',dependencies=TRUE, repos='https://kharchenkolab.github.io/drat/', type='source')"

RUN R -e "install.packages('pagoda2',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "install.packages('conosPanel',dependencies=TRUE, repos='https://kharchenkolab.github.io/drat/', type='source')"

RUN R -e 'devtools::install_github("kharchenkolab/leidenAlg")'

RUN R -e 'devtools::install_github("kharchenkolab/conos@v1.4.4")'

RUN mkdir /module
RUN mkdir /job_1234

# We will add some code to /module later

#ADD companion_script.py /module/companion_script.py
#ADD call_discover.py /module/call_discover.py
#ADD null_module.sh /module/null_module.sh

COPY src /module

# build using this:
# docker build -t genepattern/conos:2.0 .

# default command
CMD ["Rscript", "--version"]
