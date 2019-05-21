#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# convert embl file to gbk
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: mars 1015
# ---------------------------------------------------------------------------

import parse_blastoutfmt6
import gbk2circos

def blast2circos(gbk_file, accession2n_eukaryotes, accession2n_hits, accession2bacteria, accession2n_chlamydiales, accession2n_non_chlamydiales, draft_coordinates=False):
    from Bio import SeqIO
    with open(gbk_file, "r") as gbk:

        circos_file_n_euk = open("circos_n_blast_hits_euk.txt", "w")
        circos_file_n_hits = open("circos_n_blast_hits.txt", "w")
        circos_file_n_bact = open("circos_n_blast_hits_bact.txt", "w")
        circos_file_n_chlamydiales = open("circos_n_chlamydiales.txt", "w")
        circos_file_n_non_chlamydiales = open("circos_n_non_chlamydiales.txt", "w")


        record = SeqIO.read(gbk, "genbank")
        draft_data = gbk2circos.circos_fasta_draft_misc_features(record)

        for feature in record.features:
            if feature.type == "CDS":
                for i in draft_data:
                    # determine to which contig the feature belong

                    if feature.location.start >= i[1] and feature.location.end <= i[2]:
                        if draft_coordinates:
                            contig = i[0]
                            start = feature.location.start - i[1]
                            end = feature.location.end - i[1]
                        else:
                            contig = i[0]
                            start = feature.location.start
                            end = feature.location.end

                        if 'pseudo' in feature.qualifiers:
                            continue
                        try:

                            circos_file_n_euk.write('%s %s %s %s id=%s\n' % (contig,
                                                                            start,
                                                                            end,
                                                                            accession2n_eukaryotes[feature.qualifiers["locus_tag"][0]],
                                                                            feature.qualifiers["locus_tag"][0]))
                            circos_file_n_hits.write('%s %s %s %s id=%s\n' % (contig,
                                                                            start,
                                                                            end,
                                                                            accession2n_hits[feature.qualifiers["locus_tag"][0]],
                                                                            feature.qualifiers["locus_tag"][0]))
                            circos_file_n_bact.write('%s %s %s %s id=%s\n' % (contig,
                                                                            start,
                                                                            end,
                                                                            accession2bacteria[feature.qualifiers["locus_tag"][0]],
                                                                            feature.qualifiers["locus_tag"][0]))
                            circos_file_n_chlamydiales.write('%s %s %s %s id=%s\n' % (contig,
                                                                            start,
                                                                            end,
                                                                            accession2n_chlamydiales[feature.qualifiers["locus_tag"][0]],
                                                                            feature.qualifiers["locus_tag"][0]))
                            circos_file_n_non_chlamydiales.write('%s %s %s %s id=%s\n' % (contig,
                                                                            start,
                                                                            end,
                                                                            accession2n_non_chlamydiales[feature.qualifiers["locus_tag"][0]],
                                                                            feature.qualifiers["locus_tag"][0]))


                        except:

                            circos_file_n_euk.write('%s %s %s %s id=\n' % (contig,
                                                                            start,
                                                                            end,
                                                                            0))
                            circos_file_n_hits.write('%s %s %s %s id=no_hit\n' % (contig,
                                                                            start,
                                                                            end,
                                                                            0))
                            circos_file_n_bact.write('%s %s %s %s id=no_hit\n' % (contig,
                                                                            start,
                                                                            end,
                                                                            0))

                            circos_file_n_chlamydiales.write('%s %s %s %s id=no_hit\n' % (contig,
                                                                            start,
                                                                            end,
                                                                            0))
                            circos_file_n_non_chlamydiales.write('%s %s %s %s id=no_hit\n' % (contig,
                                                                            start,
                                                                            end,
                                                                            0))




if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_tab', type=str, help="input fasta tab files", nargs='+')
    parser.add_argument("-g", '--gbk_file', type=str, help="input gbk files")
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)
    parser.add_argument("-t", '--taxonomy_dico', type=str, help="taxonomy dico", default=False)

    args = parser.parse_args()

    blast_results = parse_blastoutfmt6.parse_blast_outfmt6(*args.input_tab)

    accession2n_hits = parse_blastoutfmt6.blast_dico2n_blast_hits(blast_results)
    accession2n_eukaryote_hits = parse_blastoutfmt6.blast_dico2n_eukaryote_hits(blast_results)
    accession2n_bacterial_hits = parse_blastoutfmt6.blast_dico2n_bacterial_hits(blast_results)

    with open(args.taxonomy_dico) as f:
        taxid2classification=json.load(f)

    classif, sum_order_count, chlamydiales_count, non_chlamydiales_count = parse_blastoutfmt6.investigate_classification(blast_results, taxid2classification)

    blast2circos(args.gbk_file, accession2n_eukaryote_hits, accession2n_hits, accession2n_bacterial_hits, chlamydiales_count, non_chlamydiales_count)