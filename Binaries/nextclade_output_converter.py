#!/usr/bin/env python

'''Script for converting a Nextclade output to Karoline Bragstad-format'''

import sys, csv, re

def handleoutput(aasubs, aadels):
    '''Takes a row's values for aa mutations and aa deletions and returns in dic form'''
    outputdic = {k: [] for k in ["ORF1a", "ORF1b", "S", "S2", "S3", "S4", "S5", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF14", "ORF10", "ORF1a_1"]}
    subs = aasubs.split(",")
    dels = aadels.split(",") if not aadels.split(",") == [''] else []

    for elem in subs + dels:
        if elem == '':
            continue
        prot, pos = elem.split(":")[0], elem.split(":")[1]
        outputdic[prot] += [pos]
    return outputdic

with open(sys.argv[1]) as infile:
    clades = csv.reader(infile, delimiter=";")
    header = next(clades)
    namecol = 1
    try:
        aasubcol = header.index("aaSubstitutions")
        aadelcol = header.index("aaDeletions")
    except ValueError:
        aasubcol = 17
        aadelcol = 19

    output_headers = ["name", "ORF1a", "ORF1a_1", "ORF1b", "S", "S2", "S3", "S4", "S5", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF14", "ORF10"]
    print("\t".join(output_headers))

    for row in clades:
        name = row[namecol]
        muts = row[aasubcol]
        dels = row[aadelcol]
        allchanges = handleoutput(muts, dels)

        if len(allchanges['S']) <= 22:
            allchanges['S2'] = []
            allchanges['S3'] = []
            allchanges['S4'] = []
            allchanges['S5'] = []
        elif len(allchanges['S']) <= 44:
            allchanges['S2'] = allchanges['S'][22:]
            allchanges['S'] = allchanges['S'][:22]
            allchanges['S3'] = []
            allchanges['S4'] = []
            allchanges['S5'] = []
        elif len(allchanges['S']) <= 66:
            allchanges['S2'] = allchanges['S'][22:44]
            allchanges['S'] = allchanges['S'][:22]
            allchanges['S3'] = allchanges['S'][44:]
            allchanges['S4'] = []
            allchanges['S5'] = []
        elif len(allchanges['S']) <= 88:
            allchanges['S2'] = allchanges['S'][22:44]
            allchanges['S'] = allchanges['S'][:22]
            allchanges['S3'] = allchanges['S'][44:66]
            allchanges['S4'] = allchanges['S'][66:]
            allchanges['S5'] = []
        else:
            allchanges['S2'] = allchanges['S'][22:44]
            allchanges['S'] = allchanges['S'][:22]
            allchanges['S3'] = allchanges['S'][44:66]
            allchanges['S4'] = allchanges['S'][66:88]
            allchanges['S5'] = allchanges['S'][88:]

        # Splitting 'ORF1a' into 'ORF1a' and 'ORF1a_1'
        if len(allchanges['ORF1a']) <= 22:
            allchanges['ORF1a_1'] = []
        else:
            allchanges['ORF1a_1'] = allchanges['ORF1a'][22:]
            allchanges['ORF1a'] = allchanges['ORF1a'][:22]

        outstring = name
        for prot in output_headers[1:]:  # Skip "name"
            outstring += "\t" + ";".join(allchanges[prot])
        print(outstring)
