#!/usr/bin/env python

# Very simple script (no multithreading) to download
# the complete non redundant

import ftplib
import re


ftp = ftplib.FTP("ftp.ncbi.nih.gov")
ftp.login("anonymous")
ftp.cwd("/refseq/release/complete")

# it would probably be better to use startsWith and check the last
# letters to make sure it ends with faa.gz, instead of using regular 
# expression
nr_re = re.compile("complete.nonredundant_protein.(\d)*.protein.faa.gz")
nr_filelist = (i for i in ftp.nlst() if not re.match(nr_re, i) is None)

for f in nr_filelist:
    output_file = open(f, "w")
    ftp.retrbinary("RETR "+f, output_file.write)
    output_file.close()
