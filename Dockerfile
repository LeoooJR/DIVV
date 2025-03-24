FROM python:3.12-slim

RUN mkdir /vcfdelta /apps

WORKDIR /vcfdelta

COPY requirements.txt /vcfdelta

RUN apt-get update \
&& apt-get install -y build-essential autoconf automake libtool zlib1g-dev libbz2-dev liblzma-dev libssl-dev wget bzip2

RUN pip install --upgrade pip \
&& pip install --no-cache-dir -r requirements.txt

RUN cd /apps \
&& wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 \
&& tar -jxvf bcftools-1.21.tar.bz2 \
&& cd bcftools-1.21 \
&& ./configure --prefix=/apps/bcftools \
&& make \
&& make install \
&& cd - \
&& rm -r bcftools-1.21 bcftools-1.21.tar.bz2

ENV PATH="/apps/bcftools/bin:$PATH"

COPY src /vcfdelta