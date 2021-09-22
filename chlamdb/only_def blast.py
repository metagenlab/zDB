def blast(request):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)
    #server = manipulate_biosqldb.load_db(db)
    blast_form_class = make_blast_form(db) 

    if request.method == 'POST': 

        form = blast_form_class(request.POST)

        if form.is_valid():  
            from Bio.Blast.Applications import NcbiblastpCommandline
            from Bio.Blast.Applications import NcbiblastnCommandline
            from Bio.Blast.Applications import NcbitblastnCommandline
            from Bio.Blast.Applications import NcbiblastxCommandline
            from Bio.Seq import reverse_complement, translate
            from tempfile import NamedTemporaryFile
            from io import StringIO
            from Bio.Blast import NCBIXML
          
            
            import os
            from chlamdb.biosqldb import shell_command
            import re


            input_sequence = form.cleaned_data['blast_input']

            target_accession = form.cleaned_data['target'] #form.cleaned_data returns a dictionary of validated form input fields and their values, where string primary keys are returned as objects 
            blast_type = form.cleaned_data['blast']

            target_accession_pl=''
            if target_accession.endswith(' plasmid'):
                target_accession_pl=target_accession #keep for plasmid filtering
                target_accession= int (target_accession.replace(' plasmid', ''))
                

            unknown_format = False

            print('target_accession' , target_accession)

            if '>' in input_sequence:
                my_record = [i for i in SeqIO.parse(StringIO(input_sequence), 'fasta')]
                seq = my_record[0]
                                          
            else:
                input_sequence == input_sequence.rstrip(os.linesep)
                seq = Seq(input_sequence)
                                
            dna = set("ATGCN") #add ambiguous N
            prot = set('ACDEFGHIKLMNPQRSTVWY') #X
            check_seq_DNA = set(seq) - dna
            check_seq_prot = set(seq) - prot
            
            if  not check_seq_DNA:
                try: my_record[0].description = "DNA"
                except: my_record = [SeqRecord(seq, id="INPUT", description="DNA")]
                seq_type= my_record[0].description
                    
            elif  not check_seq_prot:
                try: my_record[0].description = "Protein"
                except: my_record = [SeqRecord(seq, id="INPUT", description="Protein")]
                seq_type= my_record[0].description
                print('my_record_changed', my_record)
            else:
                unknown_format = True
                                                   
                
            if not unknown_format:
                if seq_type == 'DNA' and blast_type in ["blastp", "tblastn"]:
                    wrong_format = True
                elif seq_type =='Protein' and blast_type in ["blastn_ffn", "blastn_fna", "blastx"]:
                    wrong_format = True

                else:
                    query_file = NamedTemporaryFile(mode='w')
                    SeqIO.write(my_record, query_file, "fasta")
                    query_file.flush()
                    
                    
                    if target_accession =='all':
                        key_dict = 'merged'
                    else:
                        dictionary_acc_names=db.get_taxon_id_to_filenames() #here db replaces 'db_utils.DB' that is used to create the link because it has already been done
                        key_dict=dictionary_acc_names[int(target_accession)]               
                        print('dictionary_key', key_dict)  
                        print('dictionary_acc_names', dictionary_acc_names)            
                    if blast_type=='blastn_ffn':
                        blastType = 'locus'
                        blastdb = settings.BLAST_PATH + "ffn/%s_seq" % (key_dict)
                        blast_cline = NcbiblastnCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
                        print('BLASTDB', blastdb)
                    if blast_type=='blastn_fna':
                        blastType = 'genome'
                        blastdb = settings.BLAST_PATH  + "fna/%s_seq" % (key_dict)
                        print(blastdb)
                        blast_cline = NcbiblastnCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
                    if blast_type=='blastp':
                        blastType = 'locus'
                        blastdb = settings.BLAST_PATH  + "faa/%s_seq" % (key_dict)
                        blast_cline = NcbiblastpCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
                    if blast_type=='tblastn':
                        blastType = 'genome'
                        blastdb = settings.BLAST_PATH  + "fna/%s_seq" % (key_dict)
                        blast_cline = NcbitblastnCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
                        blast_cline2 = NcbitblastnCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=5)
                    if blast_type=='blastx':
                        blastType = 'locus'
                        blastdb = settings.BLAST_PATH  + "faa/%s_seq" % (key_dict)
                        blast_cline = NcbiblastxCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
                    
                    
                    blast_stdout, blast_stderr = blast_cline()
                    print('FINISH FIRST PART')
                    print('blast_stdout', blast_stdout)
             
             
             
                     #part below is for tblastn      

                    if blast_type=='tblastn':
                        from Bio.SeqUtils import six_frame_translations

                        blast_stdout2, blast_stderr2 = blast_cline2()

                        blast_records = NCBIXML.parse(StringIO(blast_stdout2))
                        all_data = []
                        best_hit_list = []
                        for record in blast_records:
                            for n, alignment in enumerate(record.alignments):
                                print('ALIGNMENTS', alignment)
                                accession = alignment.title.split(' ')[0] #before it was 1 but it went to the beginning of the name of the sample, with 0 we take the contig name
                                print('ACCESSION', accession)
                                sql = 'select description from bioentry where accession="%s" ' % accession

                                #description = server.adaptor.execute_and_fetchall(sql,)[0][0]

                    #            for n2, hsp in enumerate(alignment.hsps):
                    #                if n == 0 and n2 == 0:
                    #                    best_hit_list.append([record.query, hsp.sbjct_start, hsp.sbjct_end])
                                
                    #                start = hsp.sbjct_start
                    #                end = hsp.sbjct_end
                    #                if start > end:
                    #                    start = hsp.sbjct_end
                    #                    end = hsp.sbjct_start
                    #                length = end-start
                    #                leng = end-start
                    #                seq = manipulate_biosqldb.location2sequence(db, db.db_name, accession, db, start, end, leng) #replaced biodb wtih db
                                    
                    #                anti = reverse_complement(seq)
                    #                comp = anti[::-1]
                    #                length = len(seq)
                    #                frames = {}
                    #                for i in range(0, 3):
                    #                    fragment_length = 3 * ((length-i) // 3)
                    #                    tem1 = translate(seq[i:i+fragment_length], 1)
                    #                    frames[i+1] = '<span style="color: #181407;">%s</span><span style="color: #bb60d5;">%s</span><span style="color: #181407;">%s</span>' % (tem1[0:100], tem1[100:len(tem1)-99], tem1[len(tem1)-99:])
                    #                    tmp2 = translate(anti[i:i+fragment_length], 1)[::-1]
                    #                    frames[-(i+1)] = tmp2
                    #                all_data.append([accession, start, end, length, frames[1], frames[2], frames[3], frames[-1], frames[-2], frames[-3], description, seq])   #all_data is required in the hyml file to display the second graph
                    #    if len(best_hit_list) > 0:
                    #        fig_list = []
                    #        for best_hit in best_hit_list:
                    #            accession = best_hit[0]
                    #            print['ACCESSION', accession]
                    #            best_hit_start = best_hit[1]
                    #            best_hit_end = best_hit[2]
                    #            temp_location = os.path.join(settings.BASE_DIR, "assets/temp/")
                    #            temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
                    #            name = 'temp/' + os.path.basename(temp_file.name)
                    #            fig_list.append([accession, name])
                    #            orthogroup_list = mysqldb_plot_genomic_feature.location2plot(db, #UNDERSTAND WHAT WAS THE DIFFERENCE BETWEEN db AND biodb IN THE PREVIOUS TOP CODE #I removed biodb
                    #                                                                        target_accession,
                    #                                                                        temp_file.name,
                    #                                                                        best_hit_start-15000,
                    #                                                                        best_hit_end+15000,
                    #                                                                        cache,
                    #                                                                        color_locus_list = [],
                    #                                                                        region_highlight=[best_hit_start, best_hit_end])

                             #part above is for tblastn


                    no_match = re.compile('.* No hits found .*', re.DOTALL) #I still do not know how it knows that there are no hits bit it works properly
                    

                    if no_match.match(blast_stdout):
                        print ("no blast hit")
                        blast_no_hits = blast_stdout
                    elif len(blast_stderr) != 0:
                        print ("blast error")
                        blast_err = blast_stderr #linked to the html file
                    else:
                    #PART FOR PLASMIDS
                        if target_accession_pl.endswith('plasmid'):
                            from Bio import SearchIO
                            blast_record = list(SearchIO.parse(StringIO(blast_stdout), 'blast-text'))#'blast-text'
                            all_locus_tag = []
                            for query in blast_record:
                                print('QUERY', query)
                                for hit in query:
                                    locus_tag = hit.id
                                    all_locus_tag.append(locus_tag)

                            locus_filter = '"' + '","'.join(all_locus_tag) + '"'
                            print('locus_filter', locus_filter)
                            print('all_locus_tag', all_locus_tag)
                            print('blast_record', blast_record)
                            print('I WANT PLASMID')

                            dictionary_pl_contig_bioentry= db.get_bioentry_id_to_plasmid_contigs() #dictionary with bioentry_id of locus tags linked to plasmids and the locus_tag names as values
                            print('contig_pl_dict', dictionary_pl_contig_bioentry)  #CHECK THE FUNCTION IN db_utils
                            #key_dict=dictionary_acc_names[int(target_accession)] 
                            set_all_locus_tag= set(all_locus_tag)
                            print(set_all_locus_tag, 'set_all_locus_tag')
                            dict_values= list(dictionary_pl_contig_bioentry.values())
                            print('dict_values', dict_values)

                            lst_plasmid_output=(list(set(all_locus_tag).intersection(list(dictionary_pl_contig_bioentry.values())))) #list of locus_tags present in the blast output and classified as plasmid-linked
                            print('plasmid_output', lst_plasmid_output)

                            if  len(lst_plasmid_output) == 0:
                                print ('No hits on plasmids')
                                no_plasmid_hits = True #linked to the html file
                            else: 
                                for query in blast_record:
                                    for hit in query:
                                        print('HIT', hit)
                                        print(dir(hit))
                                        print('hit_map')    
                                        print('hsps', hit.hsps)    
                                        print('fragments',hit.fragments)    
                                        print('index', hit.index)    
                                        print('iterhits',hit.iterhits)   

                                        for HSPFragment in query:
                                            print(dir(HSPFragment.hsps))
                                            print(HSPFragment.index)
                                           # for hit in query:
                                            #    locus_tag = hit.id
                                             #   if locus_tag in lst_plasmid_output:
                                            blast_file_pl = settings.BASE_DIR + '/assets/temp/blast_pl.xml'
                                            filtered_out_blast=open(blast_file_pl, 'w')
                                            print("Foobar")
                                            filtered_out_blast.write(HSPFragment.hsps)

                                            print('RE-OUTPUT', filtered_out_blast)
                       
                        #sql = 'select locus_tag, product from orthology_detail_%s where locus_tag in (%s)' % (db, locus_filter) #REDO WITH THE NEW SQL DATABASE STRUCTURE
                        #locus2product = manipulate_biosqldb.to_dict(db.server.adaptor.execute_and_fetchall(sql,))

                        rand_id = id_generator(6) #it generates a random id of 6 character
                        
                        blast_file_l = settings.BASE_DIR + '/assets/temp/%s.xml' % rand_id #this file changes every time you run it and it contians the blast output
                        f = open(blast_file_l, 'w')
                        f.write(blast_stdout)
                        #print('stout', blast_stdout)
                        #print('f', f)
                        f.close()
                        

                        asset_blast_path = '/assets/temp/%s.xml' % rand_id
                        js_out = True #it could be for java script

            envoi= True #here the data are passed when it is done correctly

    else:  
        form = blast_form_class()

    return render(request, 'chlamdb/blast.html', my_locals(locals()))