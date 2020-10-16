from continuumio/miniconda3:4.8.2

WORKDIR /tmp


RUN conda install pycryptosat \
 && conda config --set sat_solver pycryptosat \
 && conda config --set channel_priority strict

# make that two steps (base and then tools)
ADD workflow/envs/environment.yml environment.yml
RUN conda env create -f environment.yml
ADD parser parser
RUN pip install ./parser
RUN mkdir -p /app/workdir
ADD workflow/scripts/  /app/scripts
VOLUME /app/workdir
WORKDIR /app/workdir

