#!/usr/bin/python

# Map peptides on protein sequence and return coordinates

import sys
import getopt
from Bio import SeqIO
import pandas as pd
import re

if __name__ == "__main__":

    # reading arguments
        
    unixOptions = "f:d:"  
    gnuOptions = ["file=","db="] # two arguments - fasta file with peptide sequences (--file) and fasta of protein database (db)
    fullCmdArguments = sys.argv
    argumentList = fullCmdArguments[1:]

    try:  
        arguments, values = getopt.getopt(argumentList, unixOptions, gnuOptions)
    except getopt.error as err:  
        print (str(err))
        sys.exit(2) 

    
    db_input = ''
    file_input = ''
    for currentArgument, currentValue in arguments:  
        if currentArgument in ("-f", "--file"):
            file_input = currentValue                     
        if currentArgument in ("-d", "--db"):
            db_input = currentValue

    #reading protein database
    dict_db1 = SeqIO.index(db_input, "fasta")

    #reading peptide fasta file
    results = []

    for record in SeqIO.parse(file_input, 'fasta'):
        results.append([record.id, record.seq])

    # mapping peptides to protein sequence
    print('Start mapping...')
    err_cnt = 0
    pep_list = []
    for record in results:
        try:
            seq=dict_db1[record[0]].seq # protein seq
            res = re.search(pattern=str(record[1]),string=str(seq))
            #print(res)
            pep_start = res.start()
            pep_stop = res.end()
            pep_list.append([record[0],str(record[1]),pep_start,pep_stop, len(seq)]) # add protein ID, peptide and ccordinates in a list
            
        except:
            err_cnt+=1
            continue
    
    print(f"{err_cnt} proteins were missed")

    # writing a table
    df = pd.DataFrame(columns=['protein','peptide','start','stop','prot_length'], data=pep_list)

    # out table
    df.to_csv(f"mapped_pep_lst_from_{file_input.split('/')[-1]}.csv", index=False)
