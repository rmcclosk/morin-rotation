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

5. Install the DNAcopy R package. From the R console:

        source("http://bioconductor.org/biocLite.R")
        biocLite("DNAcopy")

7. Install ExomeCNV. First get the file.

        wget http://cran.r-project.org/src/contrib/Archive/ExomeCNV/ExomeCNV_1.4.tar.gz

    Then from the R console:

        install.packages("ExomeCNV_1.4.tar.gz", repos=NULL, type="source")

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

Setting up Lumpy
----------------

1. Install yaha.

        wget http://faculty.virginia.edu/irahall/support/yaha/YAHA.0.1.82.tar.gz
        tar xf YAHA.0.1.82.tar.gz
        mv yaha $HOME/bin

2. Install Lumpy.

        wget -O lumpy-sv-0.2.7.tar.gz https://github.com/arq5x/lumpy-sv/releases/download/0.2.7/lumpy-sv-0.2.7.tar.gz
        tar xf lumpy-sv-0.2.7.tar.gz
        cd lumpy-sv-0.2.7
        make
        cp scripts/* $HOME/bin

Setting up TITAN
----------------

1. Install the TitanCNA R package. From R:

        source("http://bioconductor.org/biocLite.R")
        biocLite("BiocUpgrade")
        biocLite("TitanCNA")

   Note that you must be running BioConductor >= 2.14, which is only available
   for R 3.1. The instructons on the TITAN website state that you need R 3.0.2,
   but this is incorrect.

2. Install bcftools.

        wget -O bcftools-1.1.tar.bz2 http://sourceforge.net/projects/samtools/files/samtools/1.1/bcftools-1.1.tar.bz2/download
        tar xf bcftools-1.1.tar.bz2
        cd bcftools-1.1
        make
        cp bcftools $HOME/bin

3. Download TITANRunner.

        wget http://compbio.bccrc.ca/files/2013/07/TITANRunner-0.1.1.zip
        unzip TITANRunner-0.1.1.zip
        cd TITANRunner

4. Obtain a local copy of dbSNP in VCF format.

        wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/common_all.vcf.gz

5. Copy ``config_default.cfg`` to ``config.cfg`` and change the values to suit
   your environment.

6. Follow the installation instructions
   (here)[http://compbio.bccrc.ca/software/titan/titan-installation/].

Other tools
-----------

1. Download Picard tools.

        wget -O picard-tools-1.119.zip http://sourceforge.net/projects/picard/files/latest/download?source=files
        unzip picard-tools-1.119.zip
        cd picard-tools-1.119
        cp * $HOME/bin
