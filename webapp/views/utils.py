def safe_replace(string, search_string, replace_string):
    if string:
        return string.replace(search_string, replace_string)
    return string


title2page = {
    'COG Ortholog': ['fam_cog'],
    'Comparisons: Antimicrobial Resistance': [
        'amr_comparison', 'index_comp_amr', 'entry_list_amr', 'extract_amr'],
    'Comparisons: Clusters of Orthologous groups (COGs)': [
        'cog_barchart', 'index_comp_cog', 'cog_phylo_heatmap',
        'cog_comparison', 'entry_list_cog', 'extract_cog', 'heatmap_COG',
        'pan_genome_cog', 'plot_heatmap_cog', 'venn_cog'],
    'Comparisons: Kegg Orthologs (KO)': [
        'entry_list_ko', 'extract_ko', 'heatmap_ko', 'index_comp_ko',
        'ko_comparison', 'module_barchart', 'pan_genome_ko', 'plot_heatmap_ko',
        'venn_ko'],
    'Comparisons: PFAM domains': [
        'entry_list_pfam', 'extract_pfam', 'heatmap_pfam', 'pan_genome_pfam',
        'pfam_comparison', 'venn_pfam', 'index_comp_pfam', 'plot_heatmap_pfam'],
    'Comparisons: orthologous groups': [
        'extract_orthogroup', 'heatmap_orthogroup', 'index_comp_orthogroup',
        'orthogroup_comparison', 'pan_genome_orthogroup', 'plot_heatmap_orthogroup',
        'venn_orthogroup'],
    'Genome alignments: Circos plot': ['circos'],
    'Genome alignments: Plot region': ['plot_region'],
    'Genomes: table of contents': ['extract_contigs', 'genomes_intro'],
    'Homology search: Blast': ['blast'],
    'Kegg Ortholog': ['fam_ko'],
    'Kegg metabolic pathways': ['kegg_genomes'],
    'Kegg module': ['KEGG_module_map'],
    'Metabolism: kegg based': [
        'KEGG_mapp_ko', 'kegg', 'kegg_genomes_modules', 'kegg_module',
        'kegg_module_subcat', 'module_comparison'],
    'Pfam domain': ['fam_pfam'],
    'Phylogeny': ['phylogeny_intro'],
    'Table of content: Genomes': ['genomes'],
}

page2title = {}
for value, keys in title2page.items():
    page2title.update({key: value for key in keys})
