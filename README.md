HiCdat: Hi-C data analysis tool
=========================

User guide including a tutorial for data pre-processing:

[user guide](https://github.com/MWSchmid/HiCdat/blob/master/other/userGuide.pdf?raw=true)

Binaries:

[windows 7 64bit binary](https://github.com/MWSchmid/HiCdat/blob/master/other/windows_64bit.zip?raw=true)

[MacOSX (10.9) 64bit binary](https://github.com/MWSchmid/HiCdat/blob/master/other/mac64_bit.zip?raw=true)

[linux 64bit binary](https://github.com/MWSchmid/HiCdat/blob/master/other/linux_64bit.zip?raw=true)

Data for pre-processing tutorial:

[pre-processing data (large!)](http://www.botinst.uzh.ch/static/HiCdat/At_pre-process_tutorial.zip)

For the tutorial in R, download the package and the archives below, unpack the two archives, and open "HiCdat-tutorial-arabidopsis.R" in a text editor and follow the instructions. Note that HiCdatR requires the R libraries "gplots", "randomizeBE", and "MASS". You can install them with install.packages("gplots"), install.packages("randomizeBE"), and install.packages("MASS"). The HiCdat package HiCdatR can be installed with install.packages("/path/to/HiCdatR_0.99.0.tar.gz", repos=NULL, type = "source").

[R-package](https://github.com/MWSchmid/HiCdat/blob/master/other/HiCdatR_0.99.0.tar.gz?raw=true)

[R-Scripts (including the R-tutorial script)](https://github.com/MWSchmid/HiCdat/blob/master/other/Rscripts.zip?raw=true)

[files required for the tutorial in R](https://github.com/MWSchmid/HiCdat/blob/master/other/At_tutorial_files.zip?raw=true)

If you encounter problems, please contact me.

NOTE: On linux, GLIBC needs to be at least version 2.14. Biolinux6 has a lower version.

=========================
About HiCdat:

HiCdat is an implementation of all data analysis approaches employed in <a class="reference external" href="http://www.sciencedirect.com/science/article/pii/S1097276514006029">Grob et al. 2014</a>.
Importantly, HiCdat is focussed on analysis of larger structural features of chromosomes and on comparative studies. 



<figure>
  <img src="https://raw.githubusercontent.com/MWSchmid/HiCdat/master/other/figure1.png" alt="Schematic HiCdat workflow">
  <figcaption>
  <strong>Schematic HiCdat workflow.</strong>
(<strong>A-B</strong>) After sequencing and initial quality checks have been performed, the read-ends (f: forward, r: reverse) are aligned separately to a reference genome. (<strong>C-D</strong>) After  merging the separated read-ends, each end is mapped to genomic fragments, which are either genomic bins with a fixed size or restriction fragments with variable size. (<strong>E</strong>) Genomic fragments can be associated with various data types to test for correlation and enrichment of Hi-C data with genomic and epigenetic features. (<strong>F</strong>) Finally, the data can be conveniently analyzed in R using HiCdatR.
  </figcaption>
</figure>

