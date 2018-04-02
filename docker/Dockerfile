FROM ubuntu:latest
MAINTAINER Adrian Schneider <post@adrian-schneider.de>

# Install Miniconda3
RUN apt-get -qq update && apt-get -qq -y install curl bzip2 \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3 \
    && conda update conda \
    && apt-get -qq -y remove curl bzip2 \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && conda clean --all --yes

ENV PATH /opt/conda/bin:$PATH

# Install standard Python stuff
RUN conda update --all && conda install -y numpy scipy jupyter matplotlib nose

# Install PY-ADOL-C
# RUN apt-get update && apt-get install -y wget git build-essential autotools-dev libtool autoconf
# RUN cd /usr/local && wget https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.bz2 && \
#     tar --bzip2 -xf boost_1_66_0.tar.bz2 && cd boost_1_66_0 && ./bootstrap.sh && ./b2 install
# RUN mkdir install && cd install && git clone https://github.com/b45ch1/pyadolc.git && cd pyadolc && pip install -e .

# Algopy
RUN pip install algopy
# Py-ADOLC

VOLUME /forschungspraktikum/ /certificates/
EXPOSE 8888
ADD ../jupyter_notebook_config.py /etc/jupyter/
WORKDIR /forschungspraktikum/
ENTRYPOINT pip install -e /forschungspraktikum/ && cd demo && jupyter notebook --allow-root
