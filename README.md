![image](https://github.com/dnazerg666/OrthoVenn2/raw/master/image.jpg)

# OrthoVenn2

OrthoVenn2 is a tool for comparison and annotation of orthologous gene clusters among multiple species.

Genome wide analysis of orthologous clusters is an important component of comparative genomics studies. Identifying the overlap among orthologous clusters can enable us to elucidate the function and evolution of proteins across multiple species.

web server: http://www.bioinfogenome.net/OrthoVenn

# Dependencies

There is a number of additional dependencies not provided by OrthoVenn2 authors. Additional programs include:

1. diamond: a sequence aligner for protein and translated DNA searches (https://github.com/bbuchfink/diamond).<br/>
2. orthAgogue: an agile tool for the rapid prediction of orthology relations (https://github.com/samyeaman/orthagogue).<br/>
3. mcl: a cluster algorithm for graphs (https://www.micans.org/mcl).<br/>

# Installation

1. Download the above Dependencies, and copy excutable files to a bin folder.<br/>
2. Download OrthoVenn2 and upzip.<br/>
3. Enter the directory after extracting.<br/>
4. type "python ./OrthoVenn2".<br/>

# bin

The project folder contains python files of OrthoVenn2 and the bin folder contains the executable files of the dependencies, all of them has tested in Ubuntu 18.04 (64-bit).

# Usage

Usage:<br/>
-i the path of protein fasta files folder<br/>
-w the path of output folder<br/>
-b the path of bin folder<br/>

example: python ./OrthoVenn2.py -i ./example_data -w ./work_dir -b ./bin

Note:
The protein sequence file must be fasta format (*.fasta) or compressed fasta (*.fasta.tar.gz, *.fasta.zip).
The description line is distinguished from the sequence data by a greater-than (">") symbol in the first column. The first word following the ">" symbol is the identifier of the sequence, and the rest of the line is the description (both are optional). There should be no space between the ">" and the first letter of the identifier. 

