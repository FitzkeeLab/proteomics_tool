# proteomics_tool
A Simple Tool for Analyzying Proteomics Data
Nicholas Fitzkee, 7/28/2025

## Introduction

These are two programs designed to help analyze proteomics data on nanoparticle surfaces. Generally, a proteomics facility will provide you with a list of UniProt IDs and abundances. This software is designed to simplify the following

1.	Catalogging properties like pI and MW
2.	Obtaining general functional features of the proteins
3.	Creating histograms that summarize what proteins were identified on the nanoparticle

## Usage

To use this program, you will need to have the Bio package (import Bio) installed on a recent Python 3 interpreter. To pull data from UniProt, you will need an active internet connection. If you received these files via email, you may need to rename files ending in “.py.txt” to just “.py” since some email programs don’t allow “.py” files to be shared.

The software is run in two parts. The first part, “01_uniprot_fetch.py” takes one command-line argument: a file name containing a list of UniProt IDs. In this file (and in other files used by these programs), lines starting with “#” are ignored, and anything after a “#” symbol is also ignored. One UniProt ID per line is expected. An example input file is provided as “01_input_uniprot_ids.txt.” To run thus, execute:

```
python 01_uniprot_fetch.py 01_input_uniprot_ids.txt
```

This will produce two files, “go_functions.txt” and “go_processes.txt.” These files contain

a.	The Uniprot ID
b.	The MW in Da
c.	The pI
d.	A list of all GO functions and biological functions (depending on the file), where the process/function is delimited with parentheses

In a semi-interactive approach, your responsibility is to review the files (generally only processes or functions, but both are included to help you) and come up with a list of simplified terms used to describe the proteins identified. This culled file should be saved; it has the exact same format, but only one GO term is expected (multiple words are allowed, and the parentheses should still be there to group words). An example culled process file is provided as “02_input_go_processes.txt.”

Once you have reviewed the process/function data and compiled a culled list, you can run the second program. It takes the culled list along with a filename containing UniProt/abundance pairs. For example:

```
python 02_uniprot_crunch.py 02_input_go_processes.txt 02_input_test.txt
```

As before, the input files can contain comments, designated by “#.” The Uniprot IDs found in the second file (“02_input_test.txt”) must all exist in the first file (“02_input_go_processes.txt”), though some IDs from the first file may be missing from the second.

The abundances are assumed to be log2 abundances by default, but they can also be raw, scaled, or normalized abundances. There is a flag in the “crunch” program that allows you to set whether the abundances are transformed before using (e.g., 2x). Whether you use this flag will depend on how your Proteomics Facility has processed the data; please consult them for the correct approach. For us, we do want to transform the data in this way (you may not).

You can also look at the crunch file to customize the histogram of pI values and MW.

## Output

The program will output three files. The name of the files will be based on the name of the second file you provided as input.

1.	“…-pi_hist.out.txt” This file contains the protein abundance organized by pI. The more proteins within the pI range, the higher that bin will be. By default, it is scaled so that all bins add up to 100.
2.	“…-mw_hist.out.txt” This file contains the protein abundance sorted by molecular weight.
3.	“…-proc_hist.out.txt” This isn’t really a histogram, but it contains the abundance of each protein in the classes defined by the first input file. If “other” is defined as a process and appears in the list, it is placed at the end of the file.
These files can be imported into a graphing program, such as Prism, and displayed as either a histogram or a stacked bar chart.

## Citation

If you use this software, please cite the following article:

> Shaikh T., Amarasekara D.L., Hulugalla K., Torgall V., Garrigues R.J., Mayatt R., Werfel T.A.,
> Zeczycki T.N., Fitzkee N.C. (2025) “Precision control of nanoparticle behavior with engineered
> biomimetic protein coronas.” In Preparation.


