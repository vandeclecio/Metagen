#!/bin/bash

sudo su -s
apt-get update
apt-get -y upgrade
apt-get install -y htop wget nano vim less zip
apt-get install -y perl gcc default-jre default-jdk apt-utils
apt-get install -y make python-pip
apt-get install -y python-tk
apt-get install -y r-base
pip install --upgrade pip

apt-get install -y apache2 openssh-server cadvisor
# Pacote que contem o programa para gerar o passwd
apt-get install -y whois
# Criando usuario e dando permissão
useradd -m Heisenberg
mkpasswd Heisenberg > senha.txt
usermod -a -G sudo Heisenberg
chsh -s /bin/bash Heisenberg

mkdir /home/Heisenberg/public_html
chown Heisenberg /home/Heisenberg/public_html
chgrp Heisenberg /home/Heisenberg/public_html
cd /etc/apache2/mods-enabled/

ln -s ../mods-available/userdir.conf
ln -s ../mods-available/userdir.load
/etc/init.d/apache2 start > /home/Haisenberg/ip.txt
mkdir /home/Heisenberg/public_html/Genoma
mv /home/senha.txt /home/Heisenberg/

cd /root

apt-get clean

###########################
cpan install GD::Graph
conda install -c conda-forge libgd
conda install -c bioconda perl-gdtextutil
###########################

##########################
# INSTALL SRATOOLKIT NCBI
##########################
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1-3/sratoolkit.2.8.1-3-ubuntu64.tar.gz
tar -xzf sratoolkit.2.8.1-3-ubuntu64.tar.gz
mv sratoolkit.2.8.1-3-ubuntu64/ /etc/
ln -s /etc/sratoolkit.2.8.1-3-ubuntu64/bin/* /bin/

##########################
# INSTALL SEQPREP
##########################
apt-get install -y seqprep

##########################
# INSTALL CONDA
##########################
wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
chmod +x Anaconda3-5.0.1-Linux-x86_64.sh
bash Anaconda3-5.0.1-Linux-x86_64.sh



wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
chmod +x Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh -b -p $HOME/miniconda2
export PATH="$HOME/miniconda2/bin:$PATH"
conda update --all --yes

##########################
# INSTALL BOWTIE2
##########################
conda install -c bioconda bowtie2
or 
apt-get install bowtie2

##########################
# INSTALL FASTQC
##########################
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip
chmod 755 FastQC/fastqc
mv FastQC/ /etc/
ln -s /etc/FastQC/fastqc /bin/fastqc

##########################
# INSTALL CUTADAPT
##########################
pip install --user --upgrade cutadapt
ln -s ~/.local/bin/cutadapt /bin/cutadapt

##########################
# INSTALL TrimGalore
##########################
wget https://github.com/FelixKrueger/TrimGalore/archive/0.4.3.tar.gz -O trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
mv TrimGalore-0.4.3/ /etc/
ln -s /etc/TrimGalore-0.4.3/trim_galore /bin/trim_galore

##########################
# INSTALL BIOPYTHON
##########################
RUN apt-get -y install python-biopython
RUN apt-get -y install python-biopython-doc

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

cd /root

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

cd /root

##########################
# INSTALL QUIIME 1
##########################
conda create -y -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda
# pip install qiime #<--depois que entrar no ambiente do qiime esse comando deve ser executado
# ou conda install -c bioconda qiime

##########################
# INSTALL QIIME 2
##########################
conda create -y -n qiime2-2017.7 --file https://data.qiime2.org/distro/core/qiime2-2017.7-conda-linux-64.txt

##########################
# INSTALL Parallel-META
##########################
wget http://bioinfo.single-cell.cn/Released_Software/parallel-meta/3.4.1/parallel-meta-3.4.1-src-64.tar.gz -O parallel-meta.tar.gz
tar -xzvf parallel-meta.tar.gz
mv parallel-meta /etc/
cd /etc/parallel-meta
make
# export ParallelMETA=/etc/parallel-meta
# export PATH="$PATH:$ParallelMETA/bin"
# ln -s /etc/parallel-meta/bin/* /bin/
source ~/.bashrc
Rscript $ParallelMETA/Rscript/PM_Config.R


##########################
# INSTALL InterProScan
##########################
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.25-64.0/interproscan-5.25-64.0-64-bit.tar.gz -O interproscan.tar.gz
tar -xvzf interproscan.tar.gz
mv interproscan-5.25-64.0/ /etc/
ln -s /etc/interproscan-5.25-64.0/interproscan.sh /bin/interproscan.sh


cd /root
rm *.zip
rm *.gz
rm *.sh

