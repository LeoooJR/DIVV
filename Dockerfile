ARG python_version=3.12

FROM python:${python_version}-slim AS builder

ARG project=divv
ARG branch=main
ARG binary=divv

RUN apt-get update \
&& apt-get install -y git build-essential autoconf automake libtool zlib1g-dev libbz2-dev liblzma-dev libssl-dev bzip2 \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/*

WORKDIR /${project}

RUN git clone --depth 1 --branch ${branch} https://github.com/LeoooJR/DIVV.git /${project} \
&& ./install.sh \
&& pip install --upgrade pip \
&& pip install --no-cache-dir -r /${project}/requirements.txt \
&& pip install --no-cache-dir pyinstaller \
&& pyinstaller --clean --onefile --paths src --add-data src/templates:templates --add-data src/htslib/bin/:htslib/bin/ --workpath build --distpath dist -n ${binary} src/main.py

FROM python:${python_version}-slim

ARG project=divv
ARG branch=main
ARG binary=divv

WORKDIR /${project}

ENV PATH="$PATH:/${project}/bin"

COPY --from=builder /${project}/dist/${binary} /${project}/bin 

ENTRYPOINT [ "${binary}" ]