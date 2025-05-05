FROM python:3.12-slim

ARG project=vcfdelta

RUN mkdir /${project} /apps

WORKDIR /${project}

RUN apt-get update \
&& apt-get install -y build-essential autoconf automake libtool zlib1g-dev libbz2-dev liblzma-dev libssl-dev wget bzip2

RUN cd /apps \
&& wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 \
&& tar -jxvf bcftools-1.21.tar.bz2 \
&& cd bcftools-1.21 \
&& ./configure --prefix=/apps/bcftools \
&& make \
&& make install \
&& cd -\
&& rm -r bcftools-1.21 bcftools-1.21.tar.bz2

ENV PATH="/apps/bcftools/bin:$PATH"

COPY requirements.txt /${project}

RUN pip install --upgrade pip \
&& pip install --no-cache-dir -r /${project}/requirements.txt

COPY src /${project}/src

ENTRYPOINT [ "python", "src/main.py" ]