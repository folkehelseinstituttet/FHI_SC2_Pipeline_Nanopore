#!/usr/bin/env python

'''Script for converting a Nextclade output to Karoline Bragstad-format'''

import sys, csv, re

def handleoutput(aasubs, aadels):
    '''Takes a row's values for aa mutations and aa deletions and returns in dic form'''
    outputdic = {k: [] for k in ["ORF1a","ORF1a_1", "ORF1b", "S", "S2", "S3", "S4", "S5", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF14", "ORF10"]}
    subs = aasubs.split(",")
    dels = aadels.split(",") if not aadels.split(",") == [''] else []
    # if dels == ['']:
    #     dels = []
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

    output = ["name", "ORF1a", "ORF1a_1", "ORF1b", "S", "S2", "S3", "S4", "S5", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF14", "ORF10"]
    outstring = "\t".join(output)
    print(outstring)

    for row in clades:
        name = row[namecol]
        muts = row[aasubcol]
        dels = row[aadelcol]
        allchanges = handleoutput(muts, dels)

        # Use the original list for all splits
        original_S = allchanges['S']

        if len(original_S) <= 22:
            allchanges['S'] = original_S
            allchanges['S2'] = []
            allchanges['S3'] = []
            allchanges['S4'] = []
            allchanges['S5'] = []
        elif len(original_S) <= 44:
            allchanges['S'] = original_S[:22]
            allchanges['S2'] = original_S[22:]
            allchanges['S3'] = []
            allchanges['S4'] = []
            allchanges['S5'] = []
        elif len(original_S) <= 66:
            allchanges['S'] = original_S[:22]
            allchanges['S2'] = original_S[22:44]
            allchanges['S3'] = original_S[44:]
            allchanges['S4'] = []
            allchanges['S5'] = []
        elif len(original_S) <= 88:
            allchanges['S'] = original_S[:22]
            allchanges['S2'] = original_S[22:44]
            allchanges['S3'] = original_S[44:66]
            allchanges['S4'] = original_S[66:]
            allchanges['S5'] = []
        else:
            allchanges['S'] = original_S[:22]
            allchanges['S2'] = original_S[22:44]
            allchanges['S3'] = original_S[44:66]
            allchanges['S4'] = original_S[66:88]
            allchanges['S5'] = original_S[88:]

        # Splitting 'ORF1a' into 'ORF1a' and 'ORF1a_1'
        if 'ORF1a' in allchanges:
            if len(allchanges['ORF1a']) <= 22:
                allchanges['ORF1a_1'] = []
            else:
                allchanges['ORF1a_1'] = allchanges['ORF1a'][22:]
                allchanges['ORF1a'] = allchanges['ORF1a'][:22]


        outstring = ''
        for prot in ["ORF1a", "ORF1a_1", "ORF1b", "S", "S2", "S3", "S4", "S5", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF14", "ORF10"]:
            outstring += ";".join(allchanges[prot])
            outstring += "\t"
        outstring = outstring.rstrip()
        print(name + "\t" + outstring)
