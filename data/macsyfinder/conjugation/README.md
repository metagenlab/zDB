# Conjugation (CONJScan)

This folder contains the definitions for the 8 types of conjugative systems (MPF B, C, F, FA, FATA, G, I, T). It's better to call them one at a time, otherwise in case of tandem or integration of one system into another, macsyfinder might rejects both.

The model `CONJ` is a model that is called by other models, but one can use it to detect systems that are degraded or for which other accessory proteins are absent.

The model `MOB` detects any MOB proteins without `virb4` nearby (`forbidden`) and possibly `T4CP` (`accessory`)

The basic command is:

    macsyfinder mpf_type \
                --db-type ordered_replicon \
                -d Macsyfinder_models/Data/Conjugation/DEF \
                -p Macsyfinder_models/Data/Conjugation/HMM \
                --profile-suffix .HMM \
                --sequence-db my_protein_file \

# How to delimit ICEs ?

In the file `Tutorial_ICE.ipynb`, one can find how to detect conjugative systems in different chromosomes (or ICEs), and how to delimit them using comparative genomics.

# Download

To download only this module: [click here](https://minhaskamal.github.io/DownGit/#/home?url=https://github.com/gem-pasteur/Macsyfinder_models/tree/master/Data/Conjugation)

# References

- Abby Sophie S., Néron Bertrand, Ménager Hervé, Touchon Marie, Rocha Eduardo P. C. (2014). MacSyFinder: A Program to Mine Genomes for Molecular Systems with an Application to CRISPR-Cas Systems. In PLoS ONE, 9 (10), pp. e110726. [doi:10.1371/journal.pone.0110726](http://dx.doi.org/10.1371/journal.pone.0110726)

- Abby Sophie S., Cury Jean, Guglielmini Julien, Néron Bertrand, Touchon Marie, Rocha Eduardo P. C. (2016). Identification of protein secretion systems in bacterial genomes. In Scientific Reports, 6, pp. 23080. [doi:10.1038/srep23080](http://dx.doi.org/10.1038/srep23080)

If you use the tutorial:

- Jean Cury, Marie Touchon, Eduardo P. C. Rocha, Integrative And Conjugative Elements And Their Hosts: Composition, Distribution, And Organization. Nucleic Acids Res., [doi:10.1093/nar/gkx607](https://doi.org/10.1093/nar/gkx607).
