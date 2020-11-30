from continuumio/miniconda3:4.8.2

WORKDIR /tmp


RUN conda install pycryptosat \
 && conda config --set sat_solver pycryptosat \
 && conda config --set channel_priority strict

ENV NCOV ncov-qc
# make that two steps (base and then tools)
ADD workflow/envs/environment.yml environment.yml
ADD parser parser
RUN conda env create --name $NCOV -f environment.yml \
    && rm -rf /opt/conda/pkgs/
RUN /opt/conda/envs/$NCOV/bin/pip install ./parser
ENV PYTHONNOUSERSITE=1
RUN mkdir -p /app/workdir
ADD workflow /app
VOLUME /app/workdir
WORKDIR /app/workdir
# To be set with "singularity {shell,exec}"
RUN  mkdir -p /.singularity.d/env \
  && echo ". /etc/profile.d/conda.sh" >> /etc/bash.bashrc \
  && echo ". /etc/profile.d/conda.sh" >>  /.singularity.d/env/999-conda.sh \
  && echo "conda activate $NCOV" >>  /etc/bash.bashrc \
  && echo "conda activate $NCOV" >> /.singularity.d/env/999-conda.sh \
  && rm /root/.bashrc && rm /bin/sh && ln -s /bin/bash /bin/sh
#docker build --tag ncov-tools .  &&  singularity build ../.snakemake/singularity/7f0a6cf10df7320052c57c5468fcf292.simg docker-daemon://ncov-tools:latest
