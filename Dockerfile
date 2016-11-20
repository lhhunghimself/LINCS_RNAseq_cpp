FROM alpine:edge
RUN apk upgrade && apk add --no-cache musl-dev libgcc libgomp libstdc++  boost-filesystem
COPY w384/ /384/
COPY w96/ /96/
ADD scripts/start.sh start.sh
ENV NWELLS 96
ENTRYPOINT ["/start.sh"]
