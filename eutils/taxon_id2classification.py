#!/usr/bin/env python

if __name__ == '__main__':
    import argparse
    import sequence_id2scientific_classification
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",'--taxon_id', type=str, help="taxon ncbi id", nargs='+')

    args = parser.parse_args()


    classif = sequence_id2scientific_classification.taxon_id2scientific_classification(args.taxon_id)
    #print classif

    for accession in classif:
        data = classif[accession]
        try:
            superkingdom = data["superkingdom"]
        except:
            superkingdom = "-"
        try:
            phylum = data["phylum"]
        except:
            phylum = "-"
        try:
            class_ = data["class"]
        except:
            class_ = "-"
        try:
            order = data["order"]
        except:
            order = "-"
        try:
            family = data["family"]
        except:
            family = "-"
        try:
            genus = data["genus"]
        except:
            genus = "-"
        try:
            species = data["species"]
        except:
            species = "-"

        print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (superkingdom, phylum, class_, order, family, genus, species)