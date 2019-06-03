FROM alpine:3.9.2
RUN apk update upgrade --no-cache && apk add bash musl-dev libgcc libgomp libstdc++  boost-filesystem
COPY w384/ /384/
COPY w96/ /96/
ADD bwa /usr/local/bin/bwa
ADD scripts/multibwa.sh /usr/local/bin/multibwa.sh
ADD scripts/start.sh /usr/local/bin/start.sh
ENV NWELLS 96
ENTRYPOINT ["start.sh"]
