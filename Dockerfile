FROM lpryszcz/redundans:v0.13c

# metadata
LABEL base.image="lpryszcz/redundans"
LABEL version="1"
LABEL software="Karyon"
LABEL software.version="0.10"
LABEL description=""
LABEL website="https://github.com/Gabaldonlab/karyon"
LABEL license="GNU General Public License"

MAINTAINER Manu Molina (CRG)

RUN apt-get update
RUN apt-get install -y git
RUN git clone https://github.com/Gabaldonlab/karyon.git /root/src/karyon

#Set working directory
WORKDIR /root/src/karyon
## ENTRYPOINT ["sh", "./scripts/install.sh"]