#!/bin/bash

apt-get update
apt-get -y upgrade
apt-get install -y htop wget nano vim less zip
apt-get install -y perl gcc default-jre default-jdk apt-utils
apt-get install -y make python-pip
apt-get install -y python-tk

apt-get install -y git
pip install --upgrade pip

apt-get clean

##########################
# INSTALL SRATOOLKIT NCBI
##########################
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1-3/sratoolkit.2.8.1-3-ubuntu64.tar.gz
tar -xzf sratoolkit.2.8.1-3-ubuntu64.tar.gz
mv sratoolkit.2.8.1-3-ubuntu64/ /etc/
ln -s /etc/sratoolkit.2.8.1-3-ubuntu64/bin/* /bin/
rm sratoolkit.2.8.1-3-ubuntu64.tar.gz

##########################
# INSTALL SEQPREP
##########################
apt-get install -y seqprep

##########################
# INSTALL CONDA
##########################
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
chmod +x Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh -b -p $HOME/miniconda2
export PATH="$HOME/miniconda2/bin:$PATH"
conda update --all --yes
rm Miniconda2-latest-Linux-x86_64.sh

cpan install GD::Graph
conda install -c conda-forge libgd -y
conda install -c bioconda perl-gdtextutil -y

##########################
# INSTALL BOWTIE2
##########################
apt-get install bowtie2

##########################
# INSTALL FASTQC
##########################
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip
chmod 755 FastQC/fastqc
mv FastQC/ /etc/
ln -s /etc/FastQC/fastqc /bin/fastqc
rm fastqc_v0.11.5.zip

##########################
# INSTALLl FASTQ-SCREEM
##########################
wget https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_v0.11.4.tar.gz
tar -xzf fastq_screen_v*
mv fastq_screen_v*/ /etc/
ln -s /etc/fastq_screen_v*/fastq_screen /bin/fastq_screen
rm fastq_screen_v*

##########################
# INSTALL CUTADAPT
##########################
pip install --user --upgrade cutadapt
ln -s ~/.local/bin/cutadapt /bin/cutadapt

##########################
# INSTALL TrimGalore
##########################
wget https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.zip -O trim_galore.zip
unzip trim_galore.zip
mv TrimGalore*/ /etc/
ln -s /etc/TrimGalore-*/trim_galore /bin/trim_galore
rm trim_galore.zip

##########################
# INSTALL BIOPYTHON
##########################
apt-get -y install python-biopython
apt-get -y install python-biopython-doc

##########################
# INSTALL HMMER
##########################
wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
tar zxf hmmer-3.1b2-linux-intel-x86_64.tar.gz
mv hmmer-3.1b2-linux-intel-x86_64/ /etc/
cd /etc/hmmer-3.1b2-linux-intel-x86_64/
./configure
make
make check
make install

##########################
# INSTALL QUIIME 1
##########################
# ** Não está sendo instalado que execulto o scritp shell "sudo bash Meta_install.sh"
conda create -y -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda

##########################
# INSTALL QIIME 2
##########################
# ** Não está sendo instalado que execulto o scritp shell "sudo bash Meta_install.sh"
conda create -y -n qiime2-2017.7 --file https://data.qiime2.org/distro/core/qiime2-2017.7-conda-linux-64.txt

##########################
# INSTALL Parallel-META
##########################
wget http://bioinfo.single-cell.cn/Released_Software/parallel-meta/3.4.1/parallel-meta-3.4.1-src-64.tar.gz -O parallel-meta.tar.gz
tar -xzvf parallel-meta.tar.gz
mv parallel-meta /etc/
cd /etc/parallel-meta
make
export ParallelMETA=/etc/parallel-meta
export PATH="$PATH:$ParallelMETA/bin"
ln -s /etc/parallel-meta/bin/* /bin/
source ~/.bashrc
Rscript $ParallelMETA/Rscript/PM_Config.R
rm parallel-meta.tar.gz

##########################
# INSTALL FragGeneScan
##########################
wget https://sourceforge.net/projects/fraggenescan/files/FragGeneScan1.30.tar.gz
tar zxf FragGeneScan1.30.tar.gz
mv FragGeneScan1.30/ /etc/
cd /etc/FragGeneScan1.30/
make clean
make fgs
ln -s /etc/FragGeneScan1.30/FragGeneScan /bin/FragGeneScan

##########################
# INSTALL InterProScan
##########################
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.25-64.0/interproscan-5.25-64.0-64-bit.tar.gz -O interproscan.tar.gz
tar -xvzf interproscan.tar.gz
mv interproscan-5.25-64.0/ /etc/
ln -s /etc/interproscan-5.25-64.0/interproscan.sh /bin/interproscan.sh

##########################
# INSTALL Kaiju
##########################
git clone https://github.com/bioinformatics-centre/kaiju.git
cd kaiju/src
make
export PATH="$PATH:$kaiju/bin"
cd ../
mkdir kaijudb
cd kaijudb
makeDB.sh -r -v
makeDB.sh -p 
makeDB.sh -n 
makeDB.sh -e 

