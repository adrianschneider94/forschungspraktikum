FROM ubuntu:latest
MAINTAINER Adrian Schneider <post@adrian-schneider.de>

# Install basic stuff
RUN apt-get update && apt-get install -y build-essential g++ python-dev autotools-dev libicu-dev build-essential \
    libbz2-dev cmake git wget curl

# Install Boost
RUN cd && wget https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.gz \
    && tar xzvf boost_1_66_0.tar.gz \
    && cd boost_1_66_0 \
    && ./bootstrap.sh --prefix=/usr/local \
    && ./b2 \
    && ./b2 install

# Install CppAD
RUN cd && git clone https://github.com/coin-or/CppAD.git \
    && cd CppAD \
    && cmake -D cppad_prefix=/usr/local . \
    && make install

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


# Install Python packages
RUN conda install -y numpy scipy jupyter matplotlib nose

# Install PyCppAD

ENTRYPOINT /bin/bash
