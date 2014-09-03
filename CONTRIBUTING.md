How to set up a developer workstation
=====================================

This guide assumes you are working on a system where you do not have
administrator privelages. If you want to install software for all users,
change the --prefix argument to configure accordingly.

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
   prompts, and wait for a confirmation email. Once you get one, go back to the
   download page, click "GATK" (you don't need Queue). Navigate to where the
   archive is, then do

        tar xf GenomeAnalysisTK-3.2-2.tar.bz2

   This will produce a jar file. To use GATK more easily, copy the jar file to
   $HOME/bin, then create an executable script called "gatk" in $HOME/bin with
   the following contents.

        #!/bin/sh
        java -jar $HOME/bin/GenomeAnalysis.jar "$@"

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
