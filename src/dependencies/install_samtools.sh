###################################################
# Info: Shell script to install Samtools
# Author: Manu Molina
# Date: 23/01/2019
# version: 1.0
###################################################
docker_workdir=/root/src/karyon/dependencies


htslib_version=1.9
samtools_version=1.9
bcftools_version=1.9

echo "Installing htslib-$htslib_version"
tar -vxjf $docker_workdir/htslib-$htslib_version.tar.bz2
cd htslib-$htslib_version
$docker_workdir/htslib-$htslib_version/make
cd ..

# echo "Installing samtools-$samtools_version"
# tar -vxjf $docker_workdir/samtools-$samtools_version.tar.bz2
# cd samtools-$samtools_version
# make
# cd ..

# echo "Installing bcftools-$bcftools_version"
# tar -vxjf $docker_workdir/bcftools-$bcftools_version.tar.bz2
# cd bcftools-$bcftools_version
# make
# cd ..

pwd
ls -la

mv $docker_workdir/htslib-$htslib_version /usr/bin/htslib
# mv $docker_workdir/samtools-$samtools_version /usr/bin/samtools
# mv $docker_workdir/bcftools-$bcftools_version /usr/bin/bcftools

echo "Export PATH"
export PATH="$PATH:/usr/bin/htslib"
export PATH="$PATH:/usr/bin/samtools"
# export PATH="$PATH:/usr/bin/bcftools"