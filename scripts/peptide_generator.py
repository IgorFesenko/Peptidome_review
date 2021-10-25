#!/usr/bin/python

# Generate a list of random peptides from protein database based on peptide fasta

import sys
import getopt
from Bio import SeqIO
import pandas as pd
import numpy as np

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
        results.append([record.id, len(record.seq)])
    

    # Generating peptide list based on existing
    print('Start peptide generating...')
    err_cnt = 0
    pep_list = []
    for _ in range(1,11): # 10 repeats
        for record in results:
            try:
                seq=dict_db1[record[0]].seq # protein seq
                length = record[1] # peptide length
                pep_start = np.random.randint(1, len(seq)-length) # random peptide starts
                new_pep = seq[pep_start:pep_start+length] # new peptide sequence
                pep_list.append([record[0],new_pep,pep_start,pep_start+length,len(seq)]) # add protein ID and generated peptide in a list
                
            except:
                err_cnt+=1
                continue
    
    print(f'{len(pep_list)} peptides were generated')
    print(f"{err_cnt} proteins were missed")
    print()

    # writing a table
    df = pd.DataFrame(columns=['protein','peptide','start','stop','prot_length'], data=pep_list)

    # out table
    df.to_csv(f"generated_pep_lst_from{file_input.split('/')[-1]}.csv", index=False)


