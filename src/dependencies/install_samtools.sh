###################################################
# Info: Shell script to install Samtools
# Author: Manu Molina
# Date: 23/01/2019
# version: 1.0
###################################################
apt-get update
apt-get install -y gcc
apt-get install -y make
apt-get install -y libbz2-dev
apt-get install -y zlib1g-dev
apt-get install -y libncurses5-dev 
apt-get install -y libncursesw5-dev
apt-get install -y liblzma-dev

echo "Installing htslib-1.9"
tar -vxjf htslib-1.9.tar.bz2
cd htslib-1.9
make

echo "Installing samtools-1.9"
cd ..
tar -vxjf samtools-1.9.tar.bz2
cd samtools-1.9
make

echo "Installing bcftools-1.9"
cd ..
tar -vxjf bcftools-1.9.tar.bz2
cd bcftools-1.9
make

mv htslib-1.9 /usr/bin/htslib
mv samtools-1.9 /usr/bin/samtools
mv bcftools-1.9 /usr/bin/bcftools

echo "Export PATH"
export PATH="$PATH:/usr/bin/htslib"
export PATH="$PATH:/usr/bin/samtools"
export PATH="$PATH:/usr/bin/bcftools"