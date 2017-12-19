# omsTools: A Tool Box for NGS analysis


omsTools is a collection of scripts for NGS (mainly for RNAseq) analysis.

Currently it includes following scripts:

Name | Function |
--------------------| ----------------------------------------|
oms_lncRNA_classify |Classify lncRNA according to its relative location to mRNA|

## Installation

A number of scripts will be added to omsTools not long. My suggestion is to install omsTools in a virtual environment, so that you could keep up with the updating.

Create a venv folder and activate the environment.

```bash
cd /path/to/your/environment/
mkdir om004_venv
virtualenv om004_venv
. om004_venv/bin/activate
```


Download source code.
```bash
cd /path/to/your/download/
git clone https://github.com/bioShaun/omsTools.git
```

Install omsTools.
```bash
cd /path/to/your/download/omsTools
pip install -e .
```

Update omsTools.
```bash
# remember to activate your virtual environment first.
cd /path/to/your/download/omsTools
git pull origin master
pip install -e .
```

## Usage

After installation, input following command will output the usage of the script you need.

```bash
name_of_script --help
```
 
