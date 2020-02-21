#!/usr/bin/python



def locus2hydrophobicity_plot(biodb, locus):

    import pylab
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    # Kyte & Doolittle index of hydrophobicity
    kd = { 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
           'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
           'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
           'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }

    sql = 'select locus_tag, translation from orthology_detail where locus_tag="%s"' % (locus)
    data = server.adaptor.execute_and_fetchall(sql,)[0]
    locus = data[0]
    seq = data[1]
    num_residues = len(seq)

    #print matplotlib.get_cachedir()

    values = []
    for residue in seq:
        values.append(kd[residue])

    # I've computed the values for each residue.  Now compute the moving
    # average for a window size of 19.  Transmembrane helixes are often
    # about 20 residues in length so peaks in this plot suggest
    # transmembrane regions.

    window_size = 19
    half_window = int(round((window_size-1)/2,0))

    # Precompute the offsets list for better performance.
    offsets = range(-half_window, half_window+1)

    y_data = []
    for i in range(half_window, num_residues-half_window):
        average_value = 0.0
        for offset in offsets:
            average_value += values[i+offset]
        y_data.append(average_value / window_size)

    # Exclude the first and last residues.
    # The final list is in biologist coordinates (the first residue is
    # index 1)
    x_data = range(half_window, half_window+len(y_data))
    fig, ax = pylab.subplots(nrows=1, ncols=1, figsize=(10,5))
    ax.plot(x_data, y_data, linewidth=1.0)


    # Draw a reasonable cutoff for membrane prediction
    # Value of 1.6 taken from
    #   http://arbl.cvmbs.colostate.edu/molkit/hydropathy/
    #pylab.axhline(y=1.6)

    # Draw the known helix and ribbon ranges

    sql2 = 'select start, stop from interpro where locus_tag="%s" and signature_accession="TRANSMEMBRANE"' % (locus)
    try:
        transmembrane_data = server.adaptor.execute_and_fetchall(sql2,)
    except:
        transmembrane_data = []
    for transmembrane in transmembrane_data:
        pylab.axvspan(transmembrane[0], transmembrane[1], facecolor="yellow", alpha=0.4) # helix

    # Show exactly the length of the sequence
    pylab.axis(xmin = 1, xmax = num_residues)

    pylab.xlabel("residue number")
    pylab.ylabel("hydrophobicity (moving average over %d values)" % window_size)

    record_id = '%s' % locus
    pylab.title("K&D hydrophobicity for " + record_id)
    return fig   # save the figure to file
    #pylab.close(fig)

#locus2hydrophobicity_plot("chlamydia_04_16","WCW_RS07525")
