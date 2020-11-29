# PDBTOOL

## Table of contents
* [General info](#general-info)
* [Technologies](#technologies)
* [Setup](#setup)

## General info
This project started as simple PDB (Protein Data Bank) file parser that could be used as a bridge to some more complicated analysis of protein structures.
	
## Technologies
Project is created with:
* Python version: 3.5
* Matplotlib version: 3.0.3
* Scipy version: 1.4.1
	
## Setup
Probably the most predictable way to use pdbtool is to install it in virtual environment.  You can also probably configure it within your Linux environment, but some version issues with scipy and matplotlib may get tricky, so do it only if you are comfortable with debugging and/or enjoy it.  

```
$ cd wherever/your/sandbox/is
$ python3 -m venv pdbtool
$ cd pdbtool
$ . ./bin/activate
$ git clone https://github.com/pozharski/pdbtool.git
$ pip3 install --upgrade pip
$ pip3 install scipy matplotlib
```
Good to go.  For additional convenience, run this brief bash script

```
#!/bin/bash
cd bin
for name in $(ls ../pdbtool/scripts/*.py)
do 
fname=$(basename $name)
fname=${fname%.*}
ln -s $name $fname
done
```
This generates symbolic links so that whenever you activate the environment, all the scripts you can run in scripts subfolder are immediately available without the need to provide full path to them. 
