HiCat
=====

a Hi-C data analysis tool

This page is under construction (today is the 25th of September 2014).

It will soon host a fast and easy-to-use GUI tool for Hi-C data pre-processing and functions and examples for Hi-C data analysis in R.

HiCat is an implementation of all data analysis approaches employed in <a class="reference external" href="http://www.sciencedirect.com/science/article/pii/S1097276514006029">Grob et al. 2014</a>.
Importantly, HiCat is focussed on analysis of larger structural features of chromosomes and on comparative studies. 



<figure>
  <img src="https://raw.githubusercontent.com/MWSchmid/HiCat/master/figure1.png" alt="Schematic HiCat workflow">
  <figcaption>
  <strong>Schematic HiCat workflow.</strong>
(<strong>A-B</strong>) After sequencing and initial quality checks have been performed, the read-ends (f: forward, r: reverse) are aligned separately to a reference genome. (<strong>C-D</strong>) After  merging the separated read-ends, each end is mapped to genomic fragments, which are either genomic bins with a fixed size or restriction fragments with variable size. (<strong>E</strong>) Genomic fragments can be associated with various data types to test for correlation and enrichment of Hi-C data with genomic and epigenetic features. (<strong>F</strong>) Finally, the data can be conveniently analyzed in R using the provided script.
  </figcaption>
</figure>

