FROM rapidsai/rapidsai:cuda11.4-base-ubuntu20.04-py3.8
# rapidsai/rapidsai:22.02-cuda11.0-base-ubuntu18.04-py3.8

COPY ./requirements.txt /tmp/requirements.txt

RUN source activate rapids && \
    pip install -r /tmp/requirements.txt

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/rapids/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name rapids > rapids_container.yml