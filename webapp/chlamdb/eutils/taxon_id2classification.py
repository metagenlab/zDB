#!/usr/bin/env python

if __name__ == '__main__':
    import argparse
    import sequence_id2scientific_classification
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--taxon_id', type=str,
                        help="taxon ncbi id", nargs='+')
    parser.add_argument("-r", '--rank', action="store_true",
                        help="get taxon rank ")

    args = parser.parse_args()

    classif = sequence_id2scientific_classification.taxon_id2scientific_classification(
        args.taxon_id)
    # print classif
    if not args.rank:

        for accession in classif:
            data = classif[accession]
            try:
                superkingdom = data["superkingdom"]
            except Exception:
                superkingdom = "-"
            try:
                phylum = data["phylum"]
            except Exception:
                phylum = "-"
            try:
                class_ = data["class"]
            except Exception:
                class_ = "-"
            try:
                order = data["order"]
            except Exception:
                order = "-"
            try:
                family = data["family"]
            except Exception:
                family = "-"
            try:
                genus = data["genus"]
            except Exception:
                genus = "-"
            try:
                species = data["species"]
            except Exception:
                species = "-"

            print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (superkingdom, phylum, class_, order, family, genus, species)
    else:
        rank_dico = sequence_id2scientific_classification.taxon_id2scientific_classification(
            args.taxon_id, True)
        for taxon_id in rank_dico:
            print '%s\t%s\t%s' % (taxon_id, rank_dico[taxon_id][0], rank_dico[taxon_id][1])
