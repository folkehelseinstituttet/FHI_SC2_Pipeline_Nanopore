#!/usr/bin/env python

'''Script for converting a Nextclade output to Karoline Bragstad-format'''

import sys, csv

# Custom function for handling output
def handleoutput(aasubs, aadels):
    '''Takes a row's values for aa mutations and aa deletions and returns in dict form'''
    outputdic = {k: [] for k in ["ORF1a", "ORF1b", "S", "S2", "S3", "S4", "S5", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF14", "ORF10"]}
    subs = aasubs.split(",")
    dels = aadels.split(",") if not aadels.split(",") == [''] else []
    
    for elem in subs + dels:
        if elem == '':
            continue
        prot, pos = elem.split(":")[0], elem.split(":")[1]
        outputdic[prot] += [pos]
    return outputdic

# Output headers
output_headers = ["name", "ORF1a", "ORF1b", "S", "S2", "S3", "S4", "S5", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF14", "ORF10"]

# Processing each row in the CSV file
final_results = []
input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

with open(input_file_path) as infile:
    clades = csv.reader(infile, delimiter=";")
    header = next(clades)  # Read header

    for row in clades:
        name = row[1]  # Using 'seqName' as the 'name' column
        muts = row[29]  # 'aaSubstitutions' column
        dels = row[30]  # 'aaDeletions' column
        allchanges = handleoutput(muts, dels)

        # Applying corrected segmentation logic for S parts
        initial_S = allchanges['S']
        allchanges['S'] = initial_S[:22]
        allchanges['S2'] = initial_S[22:44]
        allchanges['S3'] = initial_S[44:66]
        allchanges['S4'] = initial_S[66:88]
        allchanges['S5'] = initial_S[88:]

        # Constructing the output for each row
        outstring = name
        for prot in output_headers[1:]:  # Skip "name"
            outstring += "\t" + ";".join(allchanges[prot])
        final_results.append(outstring)

# Write the header and final results to the output file
with open(output_file_path, 'w') as outfile:
    # Writing header
    outfile.write("\t".join(output_headers) + "\n")
    
    # Writing each processed row to the file
    for result in final_results:
        outfile.write(result + "\n")
