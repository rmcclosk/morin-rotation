How to set up a developer workstation
=====================================

This guide assumes you are working on a system where you do not have
administrator privelages. If you want to install software for all users,
change the prefix argument to configure/make accordingly.


Basic settings
--------------

Copy the file "settings\_template.R" to "settings.R" and change the values
according to your environment. You won't know the values for some of the
variables until _after_ completing all the installation steps below.


Setting up exomeCNV
-------------------

1. Install R.

        wget http://cran.stat.sfu.ca/src/base/R-3/R-3.1.1.tar.gz
        tar xf R-3.1.1.tar.gz
        cd R-3.1.1
        ./configure --prefix=$HOME && make && make install

2. Install samtools (requires ncurses development libraries).

        wget -O samtools-1.0.tar.bz2 http://sourceforge.net/projects/samtools/files/latest/download?source=files 
        tar xf samtools-1.0.tar.bz2
        cd samtools-1.0
        make && make prefix=$HOME install

3. Register for an account with the Broad Institute (in order to download
   GATK). Go to this website.

        https://www.broadinstitute.org/gatk/download

   It will prompt you to log in or register. Click "register", follow the
   prompts, and wait for a confirmation email. 
   
4. Install GATK. Once you have an account, go back to the download page, click
   "GATK" (you don't need Queue). Navigate to where the archive is, then do

        tar xf GenomeAnalysisTK-3.2-2.tar.bz2

   This will produce a jar file. To use GATK more easily, copy the jar file to
   $HOME/bin, then create an executable script called "gatk" in $HOME/bin with
   the following contents.

        #!/bin/sh
        java -jar $HOME/bin/GenomeAnalysis.jar "$@"

5. Download the human genome reference sequence. This is quite a large file (~850 MB).

        wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
        gunzip human_g1k_v37.fasta.gz

6. Retrieve the exome list from genesis. Its location is

        /projects/rmorin/analysis/aligned/DLBCL_EXOME/copy_number/SureSelect_regions.list

Setting up HMMcopy
------------------

You should already have R installed from the previous step.

1. Install HMMcopy.

        wget http://compbio.bccrc.ca/files/2013/12/HMMcopy.zip
        unzip HMMcopy.zip
        rm -rf __MACOSX
        cd HMMcopy
        cmake -DCMAKE_INSTALL_PREFIX:PATH=$HOME .
        make

2. Copy the binaries to $HOME/bin.

        mkdir -p $HOME/bin
        cp bin/* $HOME/bin

3. Add the following line to your ~/.bash_profile.

        export PATH=$PATH:$HOME/bin

4. Install the HMMcopy R package. In the R console:

        source("http://bioconductor.org/biocLite.R")
        biocLite("HMMcopy")
