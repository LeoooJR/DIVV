FROM python:3.12-slim

RUN mkdir /vcfdelta

WORKDIR /vcfdelta

COPY src /vcfdelta/
COPY requirements.txt /vcfdelta/

RUN apt-get update \
&& apt-get install -y autoconf automake libtool zlib1g-dev libbz2-dev liblzma-dev libssl-dev

RUN pip install --upgrade pip \
&& pip install --no-cache-dir -r requirements.txt