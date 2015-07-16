HiCdat: Hi-C data analysis tool
=========================

User guide including a tutorial for data pre-processing:

[user guide](https://github.com/MWSchmid/HiCdat/blob/master/other/userGuide.pdf?raw=true)

Binaries (note that the MacOSX binary was built on 10.10.1 and in addition tested on 10.8.5):

[windows 7 64bit binary](https://github.com/MWSchmid/HiCdat/blob/master/other/windows_64bit.zip?raw=true)

[MacOSX 64bit binary](https://github.com/MWSchmid/HiCdat/blob/master/other/mac_64bit.zip?raw=true)

[linux 64bit binary](https://github.com/MWSchmid/HiCdat/blob/master/other/linux_64bit.zip?raw=true)

Data for pre-processing tutorial:

[pre-processing data (full data set, ~20 Gb)](http://www.botinst.uzh.ch/static/HiCat/At_pre-process_tutorial.zip)

[pre-processing data (reduced data set, ~5 Gb)](http://www.botinst.uzh.ch/static/HiCat/At_pre-process_tutorial_small.zip)



For the tutorial in R, download the package and the archives below, unpack the two archives, and open "HiCdat-tutorial-arabidopsis.R" in a text editor and follow the instructions. Note that HiCdatR requires the R libraries "gplots", "randomizeBE", "MASS", and "HiCseg". You can install them with install.packages("gplots"), install.packages("randomizeBE"), install.packages("MASS"), and install.packages("HiCdat"). The HiCdat package HiCdatR can be installed with install.packages("/path/to/HiCdatR_0.99.0.tar.gz", repos=NULL, type = "source").

[R-package](https://github.com/MWSchmid/HiCdat/blob/master/other/HiCdatR_0.99.0.tar.gz?raw=true)

[R-Scripts (including the R-tutorial script)](https://github.com/MWSchmid/HiCdat/blob/master/other/Rscripts.zip?raw=true)

[files required for the tutorial in R](https://github.com/MWSchmid/HiCdat/blob/master/other/At_tutorial_files.zip?raw=true)

If you encounter problems, please contact me.

NOTE: if you encounter problems installing one of the R packages (other than HiCdatR), try to install it via [Bioconductor](http://http://www.bioconductor.org/):

```R
source("http://bioconductor.org/biocLite.R")
biocLite("insertNameOfPackageHere")
```

NOTE: On linux, GLIBC needs to be at least version 2.14. Biolinux6 has a lower version.

=========================
About HiCdat:

HiCdat is an implementation of all data analysis approaches employed in <a class="reference external" href="http://www.sciencedirect.com/science/article/pii/S1097276514006029">Grob et al. 2014</a>.
Importantly, HiCdat is focussed on analysis of larger structural features of chromosomes and on comparative studies. 



<figure>
  <img src="https://raw.githubusercontent.com/MWSchmid/HiCdat/master/other/figure1.png" alt="Schematic HiCdat workflow">
  <figcaption>
  <strong>Schematic HiCdat workflow.</strong>
(<strong>A-B</strong>) After sequencing and initial quality checks have been performed, the read-ends (f: forward, r: reverse) are aligned separately to a reference genome. (<strong>C-D</strong>) After pairing the separately aligned read-ends, each end is mapped to genomic fragments, which are either genomic bins with a fixed size or restriction fragments with variable size. (<strong>E</strong>) Genomic fragments can be associated with various data types to test for correlation and enrichment of Hi-C data with genomic and epigenetic features. (<strong>F</strong>) Finally, the data can be conveniently analyzed in R using HiCdatR.
  </figcaption>
</figure>

