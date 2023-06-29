#! /usr/bin/python3
import fileinput
import sys

# written by Dominik Bujna, modified by Brona Brejova
# original name cmscanToBAM.py

target_number = {}
#outputfile = open(sys.argv[2], "w+")
for line in fileinput.input():
    columns = line.split()
    #check if the line contains an entry
    if len(columns) > 16:
        if columns[5].isdigit():
            #take out the relevant data
            target_name = columns[0]
            if target_name not in target_number.keys():
                target_number[target_name] = 1
            else:
                target_number[target_name] += 1
            unique_name = target_name + "-" + str(target_number[target_name])
            accession = columns[1]
            query_name = columns[2]
            strand = columns[9]
            truncated = columns[10]
            e_value = columns[15]
            if strand == "+":
                start = int(columns[7]) - 1
                end = int(columns[8])
            else:
                start = int(columns[8]) - 1
                end = int(columns[7])

            length = end - start
            score = int(float(columns[14]) * 10)
            description = "truncated: " + truncated + ", E-value: " + e_value
            #format the output
            output = "\t".join([query_name, str(start), str(end),
                                unique_name, str(score), strand,
                                str(start), str(start), "0,0,0",
                                "1", str(end-start), "0", accession,
                                description])                                
            print(output)
