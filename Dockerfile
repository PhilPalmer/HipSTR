FROM alpine:3.8
MAINTAINER Pablo Prieto <pablo@lifebit.ai>

LABEL \
    description="HipSTR container from Lifebit"

RUN apk update && apk add bash && apk add --no-cache --virtual .build-deps \
	git build-base zlib-dev bzip2-dev xz-dev \
	&& git clone https://github.com/HipSTR-Tool/HipSTR \
	&& cd HipSTR && make && mv HipSTR /usr/local/bin/

ENTRYPOINT ["HipSTR"]
