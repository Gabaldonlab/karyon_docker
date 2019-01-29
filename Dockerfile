FROM lpryszcz/redundans:v0.13c
MAINTAINER Manu Molina (CRG)

ARG GATK_VERSION=4.0.12.0
ARG SPAdes_VERSION=3.9.0


COPY ./src/ /root/src/karyon
RUN mkdir -p /root/src/karyon/shared
RUN echo 'alias karyon="python /root/src/karyon/scripts/karyon.py"' >> ~/.bashrc

WORKDIR /root/src/karyon/dependencies/
# --------------------------------------- 
# 	RUN echo "Installing Git + wget + nano..."
# 	RUN apt-get update
# 	RUN apt-get install -y software-properties-common
# 	RUN add-apt-repository -y ppa:webupd8team/java
# 	RUN apt-get update
# 	RUN apt-get install -y git && \
# 		apt install -y wget && \
# 		apt-get install -y unzip && \
# 		apt-get install -y nano
# 	RUN apt-get install -y build-essential
# 	RUN apt-get install -y oracle-java8-installer
# 	RUN apt install -y python-pip
# 	RUN echo "Git + wget + nano installed"
# ---------------------------------------
# 	RUN echo "Installing GATK"
# 	RUN wget https://github.com/broadinstitute/gatk/releases/download/$GATK_VERSION/gatk-$GATK_VERSION.zip
# 	RUN unzip gatk-$GATK_VERSION.zip
# 	RUN echo "GATK installed"
# 	# --------------------------------------- 
# 	RUN echo "Installing Samtools"
# 	RUN sh install_samtools.sh
# 	RUN echo "Samtools installed"
# 	# ---------------------------------------
# 	# SOAPdeNovo (folder)
# 	RUN echo "Installing SOAPdeNovo"
# 	RUN echo "SOAPdeNovo installed"
# 	# --------------------------------------- 
# 	# SPAdes (.tar.gz)
# 	RUN echo "Installing SPAdes"
# 	RUN tar -xvf SPAdes-$SPAdes_VERSION-Linux.tar.gz
# 	RUN echo "SPAdes installed"
# 	# --------------------------------------- 
# 	# Trimmomatic (folder)
# 	RUN echo "Installing Trimmomatic"
# 	RUN echo "Trimmomatic installed"
# 	# --------------------------------------- 
# 	# Picard-tools (folder)
# 	RUN echo "Installing Picard-tools"
# 	RUN echo "Picard-tools installed"
# 	# ---------------------------------------
# 	# nQuire (folder)
# 	RUN echo "Installing nQuire"
# 	RUN echo "nQuire installed"
# 	# ---------------------------------------
# 	# anaconda_ete (folder)
# 	RUN echo "Installing anaconda_ete"
# 	RUN echo "anaconda_ete installed"
# 	# ---------------------------------------
# 	# Python libraries
# 	RUN pip install --upgrade pip
# 	RUN pip install numpy
# 	RUN pip install biopython
# 	RUN pip install psutil
# 	RUN pip install pysam
# 	RUN python -m pip install --user scipy matplotlib ipython jupyter pandas sympy nose seaborn
# 	# ---------------------------------------
RUN echo "Installing Samtools"
RUN apt-get update
RUN apt-get install -y gcc
RUN apt-get install -y make
RUN apt-get install -y libbz2-dev
RUN apt-get install -y bzip2
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y libncurses5-dev 
RUN apt-get install -y libncursesw5-dev
RUN apt-get install -y liblzma-dev

RUN sh install_samtools.sh
RUN echo "Samtools installed"

WORKDIR "/root/src/karyon"

# Clean cache
RUN apt-get clean
RUN set -x; rm -rf /var/lib/apt/lists/*