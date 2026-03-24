# Installation

## Dependencies
StaphSCAN requires the following software to be installed and available in your path:

* [Python](https://www.python.org/) v3.10 or later
* [Mash](https://github.com/marbl/Mash) v2.0 or later
* [Blast](https://blast.ncbi.nlm.nih.gov/) v2.17.0

And the following libraries:

* [Biopython](https://biopython.org/) v1.83 or later
* [Pandas](https://pandas.pydata.org/)
* [Rauth](https://rauth.readthedocs.io/en/latest/)
* [Requests](https://pypi.org/project/requests/)

### Install StaphSCAN

Install with conda:

```bash
conda create -n staphscan -c bioconda staphscan -y

conda activate staphscan
```

CLone repo:

```bash
conda create -n staphscan -c bioconda python=3.10 mash pandas blast biopython -y

conda activate staphscan

git clone https://github.com/riccabolla/StaphSCAN.git

cd StaphSCAN/

python main.py -h
```

### MLST database update

Once installed, you have to update the MLST database if you want to run it with the most recent alleles and profiles. <br>

To do that, you have first to:

1) Create an account on https://pubmlst.org/ 
2) Subscribe to the *Staphylococcus aureus* database  (you can do it by My Account -> Database registrations -> search for Staphylococcus aureus typing database -> Register)
3) Generate your API keys (My Account -> API keys -> Submit)

Once you have done this, you can run `staphscan --mlst_update`

The first time, it will require you to enter your API keys.

```bash
#Example

Data will be downloaded to: /path/to/staphscan/module/mlst/data/

--- PubMLST Authentication Required ---
1. Register/log in at https://pubmlst.org/site-accounts
2. Go to your profile and generate a Client Key and Client Secret.

Enter your PubMLST Client Key: #you client key
Enter your PubMLST Client Secret (hidden): #your client secret code

Requesting temporary token...

Please visit this URL to authorize StaphScan:
https://pubmlst.org/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_saureus_seqdef&page=authorizeClient&oauth_token=#randomtoken

Enter the verification code from the website: #your verification code

Requesting access token...
Authentication successful! Credentials saved.

Fetching S. aureus scheme metadata...
Downloading profiles.tsv...
Downloading arcC.fas...
Downloading aroE.fas...
Downloading glpF.fas...
Downloading gmk.fas...
Downloading pta.fas...
Downloading tpi.fas...
Downloading yqiL.fas...

Successfully updated S. aureus MLST database!

```