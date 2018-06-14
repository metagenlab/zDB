

def circos_orthology(all_record_list, ref_record_and_location, target_record_and_location_list, location = "assets/circos"):
    import gbk2circos
    import os
    import re
    import shell_command
    #print location

    with open(os.path.join(location, "circos.kar"), "w") as contig_file:
        i = 0
        for record in all_record_list:
            if i%2 == 0:
                col = 1
            else:
                col = 2
            i+=1
            # band chr_name band_name band_label start end color

            description = record.description
            description = re.sub(", complete genome\.", "", description)
            description = re.sub(", complete genome", "", description)
            description = re.sub(", complete sequence\.", "", description)
            description = re.sub("strain ", "", description)
            description = re.sub("str\. ", "", description)
            description = re.sub(" complete genome sequence\.", "", description)
            description = re.sub(" complete genome\.", "", description)
            description = re.sub(" chromosome", "", description)
            description = re.sub(" DNA", "S.", description)
            description = re.sub("Merged record from ", "", description)
            description = re.sub(", wgs", "", description)
            description = re.sub("Candidatus ", "", description)
            description = re.sub(".contig.0_1, whole genome shotgun sequence.", "", description)
            description = re.sub(" ", "-", description)
            description = re.sub("Chlamydia_", "C_", description)
            description = re.sub("Chlamydophila_", "C_", description)
            description = re.sub("Simkania_", "S_", description)
            description = re.sub("Parachlamydia_", "P_", description)
            description = re.sub("Neochlamydia_", "N_", description)
            description = re.sub("Protochlamydia_", "P_", description)
            description = re.sub("Waddlia_", "W_", description)
            description = re.sub("Estrella_", "E_", description)
            description = re.sub("Methylacidiphilum_", "M_", description)
            description = re.sub("Criblamydia_", "C_", description)
            description = re.sub("Methylacidiphilum_", "M_", description)


            accession = record.id.split(".")[0]

            contig_file.write('chr - %s %s %s %s spectral-5-div-%s\n' % (accession, description, 0, len(record), col))
            #chr - Rhab_1 Rhab_1 0 125191 spectral-5-div-4


    with open(os.path.join(location, "circos.link"), "w") as link_file:
        for link in target_record_and_location_list:
            line = '%s %s %s %s %s %s\n' % (ref_record_and_location[0], ref_record_and_location[1],ref_record_and_location[2], link[0], link[1], link[2])
            link_file.write(line)

    circos_conf = gbk2circos.Circos_config("circos.kar", show_ticks="no", show_tick_labels="no", ideogram_spacing=100, label_radius=0.01, radius=0.45)

    circos_conf.add_link("circos.link", thickness=3)

    with open(os.path.join(location, "circos.config"), "w") as f:
        f.write(circos_conf.get_file())

    cmd = "circos -outputfile %s -outputdir %s -conf %s" % ("circos_ortho", location, os.path.join(location, "circos.config"))
    #print cmd
    (stdout, stderr, return_code) = shell_command.shell_command(cmd)
    #print stdout
    #print stderr



