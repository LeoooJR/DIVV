ARG python_version=3.12

FROM python:${python_version}-slim AS builder

ARG project=divv
ARG branch=main

RUN apt-get update \
&& apt-get install -y git build-essential autoconf automake libtool zlib1g-dev libbz2-dev liblzma-dev libssl-dev bzip2 \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/*

WORKDIR /${project}

RUN git clone --depth 1 --branch ${branch} https://github.com/LeoooJR/DIVV.git /${project} \
&& ./install.sh /${project}/src/ \
&& pip install --upgrade pip \
&& pip install --no-cache-dir -r /${project}/requirements.txt \
&& pip install --no-cache-dir pyinstaller \
&& pyinstaller --clean --onefile --paths src --add-data src/templates:templates --add-data src/htslib/bin/:htslib/bin/ --hidden-import=scipy._cyutility --workpath build --distpath dist -n divv src/main.py

FROM python:${python_version}-slim

ARG project=divv
ARG branch=main

WORKDIR /${project}

ENV PATH="$PATH:/${project}/bin"

RUN mkdir bin

COPY --from=builder /${project}/dist/divv /${project}/bin 

ENTRYPOINT [ "divv" ]