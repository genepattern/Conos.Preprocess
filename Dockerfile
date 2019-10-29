from vpetukhov/conos:version-1.1.1
# using https://hub.docker.com/r/vpetukhov/conos/tags on 2019-10-16

MAINTAINER Edwin Juarez <ejuarez@ucsd.edu>

ENV LANG=C LC_ALL=C

USER root

RUN apt-get install time

RUN mkdir /module
RUN mkdir /job_1234

# We will add some code to /module later

#ADD companion_script.py /module/companion_script.py
#ADD call_discover.py /module/call_discover.py
#ADD null_module.sh /module/null_module.sh

RUN R -e 'install.packages("optparse",repos = "http://cran.us.r-project.org")'

COPY src /module

# build using this:
# docker build -t genepattern/conos:1.2 .
