# Typing

StaphSCAN includes several typing methods, each of which can be run as a stand-alone module.

For clarity, all these methods and their respective parameters are presented together on this page

## MLST

``-m assembly`` 

All genomes identified as *Staphylococcus aureus* are subject to MLST using the seven-locus typing scheme described [here](https://pubmlst.org/organisms/staphylococcus-aureus) 

A copy of the MLST alleles and ST definitions is stored in the **/data** directory of this module.

### Parameters

Hits are filtered by min. identity (``default: 90.0``) and min. coverage (``default: 90.0``). 

### Output

In the detailed report, both the ST and the allelic profile are reported. 

In the simple report, only the ST is reported. 

StaphSCAN reports also imperfect matches, following this criteria:

* Imperfect hits (identity or coverage < 100%) are indicated with a  ``*``
* Imprecise Sequence Types are annotated by reporting the closest matching ST, followed by the number of differing loci (n-locus variants, up to 2): ``ST1-1LV`` (closest match is ST1, with 1 different allele).

