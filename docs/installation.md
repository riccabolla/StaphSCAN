# Installation

## Dependencies
StaphSCAN requires the following software to be installed and available in your path:

* [Python](https://www.python.org/) v3.10 or later
* [Mash](https://github.com/marbl/Mash) v2.0 or later
* [Blast](https://blast.ncbi.nlm.nih.gov/) v2.17.0

And the following libraries:

* [Biopython](https://biopython.org/) v1.75 or later
* [Pandas]()

### Install StaphSCAN

Install with conda (available soon):

```bash
conda create -n staphscan -c bioconda staphscan -y
conda activate staphscan
```

CLone repo:

```bash
conda create -n staphscan -c bioconda python=3.10 mash blast biopython -y
conda activate staphscan
git clone https://github.com/riccabolla/StaphSCAN.git
cd StaphSCAN/
python main.py -h
```


