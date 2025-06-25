ARG project=vcfdelta
ARG branch=main
ARG binary=vcfdelta

FROM python:3.12-slim as builder

RUN apt-get update \
&& apt-get install -y git build-essential autoconf automake libtool zlib1g-dev libbz2-dev liblzma-dev libssl-dev bzip2 \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/*

WORKDIR /${project}

RUN git clone --depth 1 --branch ${branch} https://github.com/LeoooJR/VCFDelta.git /${project} \
&& ./install.sh \
&& pip install --upgrade pip \
&& pip install --no-cache-dir -r /${project}/requirements.txt \
&& pip install --no-chache-dir pyinstaller \
&& pyinstaller --clean --onefile --paths src --add-data src/templates:templates --add-data src/htslib/bin/:htslib/bin/ --workpath build --distpath dist -n ${binary} src/main.py

FROM python:3.12-slim

WORKDIR /${project}

ENV PATH="$PATH:/${project}/bin"

COPY --from=builder /${project}/dist/${binary} /${project}/bin 

ENTRYPOINT [ "${binary}" ]