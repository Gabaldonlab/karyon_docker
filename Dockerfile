FROM lpryszcz/redundans
MAINTAINER Manu Molina (CRG)

COPY ./src/ /root/src/karyon
RUN mkdir -p /root/src/karyon/shared
RUN echo 'alias karyon="python /root/src/karyon/src/scripts/karyon.py"' >> ~/.bashrc

WORKDIR /root/src/karyon/dependencies/
# --------------------------------------- 
RUN echo "Installing Git + wget + nano..."
RUN apt-get update && \
	apt-get install -y git && \
	apt install -y wget && \
	apt-get install -y unzip && \
	apt-get install -y nano
RUN apt-get install -y build-essential
RUN apt install -y python-pip
RUN echo "Git + wget + nano installed"
# ---------------------------------------
RUN echo "Installing GATK"
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.0.12.0/gatk-4.0.12.0.zip
RUN unzip gatk-4.0.12.0.zip
RUN echo "GATK installed"
# --------------------------------------- 
RUN echo "Installing Samtools"
RUN sh install_samtools.sh
RUN echo "Samtools installed"
# ---------------------------------------
# SOAPdeNovo (folder)
RUN echo "Installing SOAPdeNovo"
RUN echo "SOAPdeNovo installed"
# --------------------------------------- 
# SPAdes (.tar.gz)
RUN echo "Installing SPAdes"
RUN tar -xvf SPAdes-3.9.0-Linux.tar.gz
RUN echo "SPAdes installed"
# --------------------------------------- 
# Trimmomatic (folder)
RUN echo "Installing Trimmomatic"
RUN echo "Trimmomatic installed"
# --------------------------------------- 
# Picard-tools (folder)
RUN echo "Installing Picard-tools"
RUN echo "Picard-tools installed"
# ---------------------------------------
# nQuire (folder)
RUN echo "Installing nQuire"
RUN echo "nQuire installed"
# ---------------------------------------
# anaconda_ete (folder)
RUN echo "Installing anaconda_ete"
RUN echo "anaconda_ete installed"
# ---------------------------------------
# Python libraries
RUN pip install --upgrade pip
RUN pip install numpy
RUN pip install biopython
RUN pip install psutil
RUN pip install pysam
RUN python -m pip install --user scipy matplotlib ipython jupyter pandas sympy nose seaborn
# ---------------------------------------
WORKDIR "/root/src/karyon"