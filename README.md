# OM004script: A Tool Box for NGS analysis


OM004script is a collection of scripts for NGS (mainly for RNAseq) analysis.

Currently it includes following scripts:

Name | Function |
--------------------| ----------------------------------------|
oms_lncRNA_classify |Classify lncRNA according to its relative location to mRNA|

## Installation

A number of scripts will be added to OM004script not long. My suggestion is to install OM004script in a virtual environment, so that you could keep up with the updating.

Create a venv folder and activate the environment.

```
cd /path/to/your/environment/
mkdir om004_venv
virtualenv om004_venv
. om004_venv/bin/activate
```


Download source code.
```
cd /path/to/your/download/
git clone https://github.com/bioShaun/OM004script.git
```

Install OM004script.
```
cd /path/to/your/download/OM004script
pip install -e .
```

Update OM004script.
```
# remember to activate your virtual environment first.
cd /path/to/your/download/OM004script
git pull origin master
pip install -e .
```

## Usage

After installation, input following command will output the usage of the script you need.

```
name_of_script --help
```

## 
