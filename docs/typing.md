# Typing

StaphSCAN includes several typing methods, each of which can be run as a stand-alone module.

For clarity, all these methods and their respective parameters are presented together on this page

## Multi Locus Sequence Typing

``-m mlst`` 

All genomes identified as *Staphylococcus aureus* are subject to MLST using the seven-locus typing scheme described [here](https://pubmlst.org/organisms/staphylococcus-aureus) 

A copy of the MLST alleles and ST definitions is stored in the **/data** directory of this module.

### Parameters

Hits are filtered by min. identity (``default: 90.0``) and min. coverage (``default: 90.0``). 

### Output

Both the ST and the allelic profile are reported. 

StaphSCAN reports also imperfect matches, following this criteria:

* Imperfect hits (identity or coverage < 100%) are indicated with a  ``*``
* Imprecise Sequence Types are annotated by reporting the closest matching ST, followed by the number of differing loci (n-locus variants, up to 2): ``ST1-1LV`` (closest match is ST1, with 1 different allele).

## spa typing

``-m spa``

The spa-typing is a method based on the characterization of the repeat regions of Staphylococcus protein A gene (spa).This method is widely used for rapid typing of MRSA, particularly in hospital and surveillance settings.

For more information visit [here](https://spa.ridom.de/index.shtml). 

A local copy of [Ridom database](https://spa.ridom.de/spatypes.shtml) is distributed with this module and stored the modules's **/data** directory.

Genome assemblies are screened for the presence of the *spa* gene X-region by simulating PCR amplification with multiple published primer sets. Each primer set is tested against all contigs, and both forward and reverse-complement orientations are evaluated. 

When multiple primer sets yield a valid amplicon, the first detected product is used for downstream analysis. Assemblies in which no valid amplicon is detected are reported as spa-negative.

The amplified X-region is scanned to identify spa repeat units using a curated database of known repeat sequences. Detected repeats are recorded sequentially to generate a repeat pattern, which is then compared against the reference spa type database.

### Output

* If the identified repeat patterns match a known spa type, the corresponding type is reported.

* Patterns not present in the reference database were classified as novel and reported together with their repeat composition.

* Assemblies in which an X-region was amplified but no known repeat units were detected were reported as “Unknown".

### Limitation

Spa typing is dependent on genome assembly quality, and fragmentation or sequencing errors within the spa X-region may result in spa-negative or Unknown calls. Novel or divergent repeat patterns not present in the reference database are reported as Novel. As with all in silico typing approaches, results may differ from laboratory-based spa typing in cases of mixed populations or incomplete assemblies.

