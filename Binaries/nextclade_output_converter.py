#!/usr/bin/env python

'''Script for converting a Nextclade output to Karoline Bragstad-format'''

import sys, csv, re

def handleoutput(aasubs, aadels):
    '''Takes a row's values for aa mutations and aa deletions and returns in dic form'''
    outputdic = {k: [] for k in ["ORF1a", "ORF1b", "S", "S2", "S3", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF14", "ORF10"]}
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

    output = ["name", "ORF1a", "ORF1b", "S", "S2", "S3", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF14", "ORF10"]
    outstring = "\t".join(output)
    print(outstring)

    for row in clades:
        name = row[namecol]
        muts = row[aasubcol]
        dels = row[aadelcol]
        allchanges = handleoutput(muts, dels)

        if len(allchanges['S']) <= 22:
            allchanges['S2'] = []
            allchanges['S3'] = []
        elif len(allchanges['S']) <= 44:
            allchanges['S2'] = allchanges['S'][22:]
            allchanges['S'] = allchanges['S'][:22]
            allchanges['S3'] = []
        else:
            allchanges['S3'] = allchanges['S'][44:]
            allchanges['S2'] = allchanges['S'][22:44]
            allchanges['S'] = allchanges['S'][:22]

        outstring = ''
        for prot in ["ORF1a", "ORF1b", "S", "S2", "S3", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF14", "ORF10"]:
            outstring += ";".join(allchanges[prot])
            outstring += "\t"
        outstring = outstring.rstrip()
        print(name + "\t" + outstring)
