HiCat: Hi-C data analysis tool
=========================

User guide including a tutorial for data pre-processing:

<a class="reference external" href="https://github.com/MWSchmid/HiCat/blob/master/userGuide.pdf?raw=true">user guide</a>

Binaries:

<a class="reference external" href="https://github.com/MWSchmid/HiCat/blob/master/windows_64bit.zip?raw=true">windows 7 64bit binary</a>

<a class="reference external" href="https://github.com/MWSchmid/HiCat/blob/master/mac64bit.zip?raw=true">MacOSX (10.9) 64bit binary</a>

<a class="reference external" href="https://github.com/MWSchmid/HiCat/blob/master/linux_64bit.zip?raw=true">linux 64bit binary</a>

Data for pre-processing tutorial:

<a class="reference external" href="http://www.botinst.uzh.ch/static/HiCat/At_pre-process_tutorial.zip">pre-processing data (large!)</a>

For the tutorial in R, download the archives below, unpack both of them, and open "HiCat-tutorial-arabidopsis.R" in a text editor and follow the instructions. Note that HiCat requires the R libraries "gplots", "randomizeBE", and "MASS". You can install them with install.packages("gplots"), install.packages("randomizeBE"), and install.packages("MASS"). The HiCat package HiCatR can be installed with install.packages("/path/to/HiCatR_0.99.0.tar.gz", repos=NULL, type = "source").

<a class="reference external" href="https://github.com/MWSchmid/HiCat/blob/master/HiCatR_0.99.0.tar.gz?raw=true">R-package</a>

<a class="reference external" href="https://github.com/MWSchmid/HiCat/blob/master/Rscripts.zip?raw=true">R-Scripts (including the R-tutorial script)</a>

<a class="reference external" href="https://github.com/MWSchmid/HiCat/blob/master/At_tutorial_files.zip?raw=true">files required for the tutorial in R</a>

If you encounter problems, please contact me.

NOTE: On linux, GLIBC needs to be at least version 2.14. Biolinux6 has a lower version.

=========================
About HiCat:

HiCat is an implementation of all data analysis approaches employed in <a class="reference external" href="http://www.sciencedirect.com/science/article/pii/S1097276514006029">Grob et al. 2014</a>.
Importantly, HiCat is focussed on analysis of larger structural features of chromosomes and on comparative studies. 



<figure>
  <img src="https://raw.githubusercontent.com/MWSchmid/HiCat/master/figure1.png" alt="Schematic HiCat workflow">
  <figcaption>
  <strong>Schematic HiCat workflow.</strong>
(<strong>A-B</strong>) After sequencing and initial quality checks have been performed, the read-ends (f: forward, r: reverse) are aligned separately to a reference genome. (<strong>C-D</strong>) After  merging the separated read-ends, each end is mapped to genomic fragments, which are either genomic bins with a fixed size or restriction fragments with variable size. (<strong>E</strong>) Genomic fragments can be associated with various data types to test for correlation and enrichment of Hi-C data with genomic and epigenetic features. (<strong>F</strong>) Finally, the data can be conveniently analyzed in R using the provided script.
  </figcaption>
</figure>

