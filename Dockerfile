FROM ubuntu:16.04
MAINTAINER bioinfo@imd.ufrn.br
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get install -y htop wget nano vim less zip
RUN apt-get install -y perl gcc default-jre default-jdk apt-utils
RUN apt-get install -y make python-pip
RUN apt-get install -y python-tk
RUN apt-get install -y r-base
RUN pip install --upgrade pip

#######################
# INSTALL APACH
#######################
RUN apt-get install -y apache2
# Pacote que contem o programa para gerar o passwd
RUN apt-get install -y whois
# Criando usuario e dando permissão
RUN useradd -m Heisenberg
RUN mkpasswd Heisenberg > /home/Heisenberg/senha.txt
RUN usermod -a -G sudo Heisenberg
RUN chsh -s /bin/bash Heisenberg

# Criando pasta public_html
RUN mkdir /home/Heisenberg/public_html
RUN chown Heisenberg /home/Heisenberg/public_html
RUN chgrp Heisenberg /home/Heisenberg/public_html
WORKDIR /etc/apache2/mods-enabled/

RUN ln -s ../mods-available/userdir.conf
RUN ln -s ../mods-available/userdir.load
# Ver configuraçẽos de IP para o apache
RUN /etc/init.d/apache2 start > /home/Heisenberg/ip.txt
RUN mkdir /home/Heisenberg/public_html/Genoma
RUN mv /home/senha.txt /home/Heisenberg/

WORKDIR /root

RUN apt-get clean

######################
# INSTALL SRATOOLKIT NCBI
######################
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1-3/sratoolkit.2.8.1-3-ubuntu64.tar.gz
RUN tar -xzf sratoolkit.2.8.1-3-ubuntu64.tar.gz
RUN mv sratoolkit.2.8.1-3-ubuntu64/ /etc/
RUN ln -s /etc/sratoolkit.2.8.1-3-ubuntu64/bin/* /bin/

######################
# INSTALL SEQPREP
######################
RUN apt-get install -y seqprep

####################
# INSTALL CONDA
####################
RUN wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
RUN chmod +x Miniconda2-latest-Linux-x86_64.sh
RUN bash Miniconda2-latest-Linux-x86_64.sh -b -p $HOME/miniconda2
RUN export PATH="$HOME/miniconda2/bin:$PATH"
RUN conda update --all --yes

#####################
# INSTALL FASTQC
#####################
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
RUN unzip fastqc_v0.11.5.zip
RUN chmod 755 FastQC/fastqc
RUN mv FastQC/ /etc/
RUN ln -s /etc/FastQC/fastqc /bin/fastqc

######################
# INSTALL cutadapt
######################
RUN pip install --user --upgrade cutadapt
RUN ln -s ~/.local/bin/cutadapt /bin/cutadapt

######################
# INSTALL TrimGalore
######################
wget https://github.com/FelixKrueger/TrimGalore/archive/0.4.3.tar.gz -O trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
mv TrimGalore-0.4.3/ /etc/
ln -s /etc/TrimGalore-0.4.3/trim_galore /bin/trim_galore

###################
# INSTALL BIOPYTHON
###################
RUN apt-get -y install python-biopython
RUN apt-get -y install python-biopython-doc

###################
# INSTALL HMMER
###################
RUN wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
RUN tar zxf hmmer-3.1b2-linux-intel-x86_64.tar.gz
RUN mv hmmer-3.1b2-linux-intel-x86_64/ /etc/
WORKDIR /etc/hmmer-3.1b2-linux-intel-x86_64/
RUN ./configure
RUN make
RUN make check
RUN make install

WORKDIR /root

###################
# INSTALL FragGeneScan
###################
RUN wget https://sourceforge.net/projects/fraggenescan/files/FragGeneScan1.30.tar.gz
RUN tar zxf FragGeneScan1.30.tar.gz
RUN mv FragGeneScan1.30/ /etc/
WORKDIR  "/etc/FragGeneScan1.30/"
RUN make clean
RUN make fgs
RUN ln -s /etc/FragGeneScan1.30/FragGeneScan /bin/FragGeneScan

WORKDIR /root

###################
# INSTALL QUIIME 1
###################
RUN conda create -y -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda
#pip install qiime #<--depois que entrar no ambiente do qiime esse comando deve ser executado
# ou conda install -c bioconda qiime

##################
# INSTALL QIIME 2
##################
RUN conda create -y -n qiime2-2017.7 --file https://data.qiime2.org/distro/core/qiime2-2017.7-conda-linux-64.txt

##########################
# INSTALL Parallel-META
##########################
RUN wget http://bioinfo.single-cell.cn/Released_Software/parallel-meta/3.4.1/parallel-meta-3.4.1-src-64.tar.gz -O parallel-meta.tar.gz
RUN tar -xzvf parallel-meta.tar.gz
RUN mkdir /opt/tools
RUN mv parallel-meta /opt/tools/
WORKDIR /opt/tools/parallel-meta
RUN make
RUN export ParallelMETA=/opt/tools/parallel-meta
RUN export PATH="$PATH:$ParallelMETA/bin"
RUN source ~/.bashrc
RUN Rscript $ParallelMETA/Rscript/PM_Config.R

WORKDIR /root

##########################
# INSTALL InterProScan
##########################
RUN wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.25-64.0/interproscan-5.25-64.0-64-bit.tar.gz -O interproscan.tar.gz
RUN tar -xvzf interproscan.tar.gz
RUN interproscan-5.25-64.0/ /etc/
RUN ln -s /etc/interproscan-5.25-64.0/interproscan.sh /bin/interproscan.sh

WORKDIR /root
RUN rm *.zip
RUN rm *.gz
RUN rm *.sh


