
class CircosJs():
    '''
    Generate javascript code for circos plots

    - self.add_heatmap: add a new heatmap
    - self.get_js_code: return complete js code as string
    '''
    
    
    
    def __init__(self,
                 figure_div_id="heatmapChart"):
        
        self.contigs_data = []
        self.heatmap_tracks = '{'
        self.highlight_tracks = '{'
        self.lineplot_tracks = '{'
        self.histogram_tracks = '{'
        self.last_inner_radius = None
        self.last_outer_radius = 1.01

        self.js_template = ''' 
        
        var width = 900;

        var configuration = {
                  innerRadius: width / 2 - 80,
                  outerRadius: width / 2 - 70,
                  ticks: {display: true},
                  gap: 0,
                  tooltipContent: function (d) {
                  return `<h5> TEST: ${d.label}</h5>`
                  },
                  labels: {
                    position: 'center',
                    display: false,
                    size: 20,
                    color: 'red',
                    radialOffset: 25
                  },
              ticks: {
              display: false,
              color: 'grey',
              spacing: 50000,
              labels: true,
              labelSpacing: 1,
              labelSuffix: 'Kb',
              labelDenominator: 1000,
              labelDisplay0: true,
              labelSize: '10px',
              labelColor: '#000000',
              labelFont: 'default',
              majorSpacing: 5,
              size: {
                minor: 0,
                major: 0,
              }
              },
              events: {'click.alert': function (datum, index, nodes, event) {
                    window.alert(datum.label)
              },
              'mouseover': function (datum, index, nodes, event) {
                    return 'test'
              },
              }
        }

        var contig_data = %s ;


        var myCircos = new Circos({
            container: '#heatmapChart',
            width: width,
            height: width,
            defaultTrackWidth: 10
        });

        
        myCircos.layout(contig_data, configuration);

        %s

        %s
        
        %s
        
        %s

        myCircos.render();
  
        '''
        
        self.template_heatmap = '''
        var heatmap_data_lst = %s ;

        var i = 0;
        for (var key in heatmap_data_lst) { 
            i = i+1;
            var data = heatmap_data_lst[key][0];
            console.log(heatmap_data_lst[key][1]);
            var conf = heatmap_data_lst[key][1];
            conf["tooltipContent"] = function (d) { if (d.value == "false"){return None} else if (d.value == -1) {return `<h5>${d.locus_tag}</h5>`} else {return `<h5>${d.locus_tag}: ${d.value} %% identity</h5>` }}
            conf["events"] = {"click.alert": function (datum, index, nodes, event) { if (datum.locus_tag){window.location="/locusx/" + datum.locus_tag} else {return None};}}

            myCircos.heatmap(key, data, conf);
        }
        '''
        
        self.template_highlight = '''
        var highlight_data_lst = %s ;

        var i = 0;
        for (var key in highlight_data_lst) { 
            i = i+1;
            var data = highlight_data_lst[key][0];
            console.log(highlight_data_lst[key][1]);
            var conf = highlight_data_lst[key][1];
            conf["tooltipContent"] = function (d) {return `<h5>Locus: ${d.locus_tag} <br>Gene: ${d.gene} <br>Product: ${d.product}</h5>`}
            conf["events"] = {"click.alert": function (datum, index, nodes, event) { window.location="/locusx/" + datum.locus_tag;}}
            myCircos.highlight(key, data, conf);
        }
        
        '''
        
        self.template_lineplot = '''
        var lineplot_data_lst = %s ;

        var i = 0;
        for (var key in lineplot_data_lst) { 
            i = i+1;
            var data = lineplot_data_lst[key][0];
            console.log(lineplot_data_lst[key][1]);
            var conf = lineplot_data_lst[key][1];
            myCircos.line(key, data, conf);
        }
        
        '''

        self.template_histogram = '''
        var histogram_data_lst = %s ;

        var i = 0;
        for (var key in histogram_data_lst) { 
            i = i+1;
            var data = histogram_data_lst[key][0];
            console.log(histogram_data_lst[key][1]);
            var conf = histogram_data_lst[key][1];
            myCircos.histogram(key, data, conf);
        }
        
        '''

        self.heatmap_configuration_template = '{outerRadius: %s, innerRadius: %s, min: 0, max: 100, color: %s, logScale: false}'
        
        self.histogram_configuration_template = '{outerRadius: %s, innerRadius: %s, color: %s, logScale: false}'
        
        self.lineplot_configuration_template = '{outerRadius: %s, innerRadius: %s, color: "%s", logScale: false, fillColor: "%s", fill: true, direction: "out"}'
    

    def add_contigs_data(self, df_bioentry):
      '''
      Input df with the following columns:
      - length
      - accession
      
      The index is used as id
      '''     
      i = 0
      for bioentry, row in df_bioentry.iterrows():
          if i % 2 == 0:
              self.contigs_data.append({"len": row.length, "color": "#8dd3c7", "label": row.accession, "id": bioentry})
          else:
              self.contigs_data.append({"len": row.length, "color": "#fb8072", "label": row.accession, "id": bioentry})
          i+=1
      
    
    def add_gene_track(self, df, label, highlight_list = ["OJCDPLGI_00091", "OJCDPLGI_00096"], sep=0, radius_diff=0.04):
      '''
      Input df columns:
      - bioentry_id
      - start_pos
      - end_pos
      - locus_ref
      '''
    
      highlight_data = [{"block_id": str(row.bioentry_id), "start":row.start_pos, "end": row.end_pos, "color": row.color, "locus_tag": row.locus_ref, "gene": f"{row.gene}", "product": f'{row.gene_product}'} if row.locus_ref not in highlight_list else  {"block_id": str(row.bioentry_id), "start":row.start_pos, "end": row.end_pos, "color": "magenta", "locus_tag": row.locus_ref,  "gene": f"{row.gene}", "product": f'{row.gene_product}'} for n, row in df.iterrows()]
      
      if not self.last_inner_radius:
        outer_radius = 0.98
        inner_radius = 0.98 - radius_diff
        #conf["color"] = color
      else:
        outer_radius = self.last_inner_radius - sep
        inner_radius  = self.last_inner_radius - sep - radius_diff
        #conf["color"] = color
      
      self.last_inner_radius = inner_radius

      conf = self.heatmap_configuration_template % (outer_radius, inner_radius, 'function (d) {return d.color}') # 

      self.highlight_tracks += f'"{label}": [%s,%s],' % (highlight_data, conf)
    
    
    def add_heatmap_track(self, df, label, color="comp", sep=0, radius_diff=0.06):
        '''
        Input df columns:
        - bioentry_id
        - start_pos
        - end_pos
        - identity
        - locus_tag
        
        Locus without homologs in the target genome should have a locus_tag with a value of None
        '''
        heatmap_data = [{"block_id": row.bioentry_id, "start":row.start_pos, "end": row.end_pos, "value": row.identity, "locus_tag": row.locus_tag} if row.locus_tag is not None else {"block_id": row.bioentry_id, "start":row.start_pos, "end": row.end_pos, "value": "false"} for n, row in df.iterrows()]
          
        if not self.last_inner_radius:
            outer_radius = 0.98
            inner_radius = 0.98 - radius_diff
            #conf["color"] = color
        else:
            outer_radius = self.last_inner_radius - sep
            inner_radius  = self.last_inner_radius - sep - radius_diff
            #conf["color"] = color
          
        self.last_inner_radius = inner_radius

        if color == 'comp':
              col = '''function(datum, index) {var color_scale = chroma.scale('YlOrRd').domain([30,100]); if (datum.value == "false") {return "lightblue"} else {var col = color_scale((datum.value)).hex();return col}}'''
        else:
              col = color

        conf = self.heatmap_configuration_template % (outer_radius, inner_radius, col)
              
        self.heatmap_tracks += f'"{label}": [%s,%s],' % (heatmap_data, conf)
    

    def add_line_track(self, bioentry_df, label, fillcolor="red", sep=0, radius_diff=0.06, windows=500):
        '''
        Minimal Input df columns:
        - bioentry_id
        - seq
        '''
        from Bio.SeqUtils import GC
        
        ordered_seqs = [bioentry.seq for n,bioentry in bioentry_df.iterrows()]
        concat_seq = ''.join(ordered_seqs)
        average_gc = GC(concat_seq)
        
        linedata_data = []
        for index, bientry in bioentry_df.iterrows():
            for i in range(0, len(bientry.seq), windows):
                start = i
                stop = i + windows
                gc = GC(bientry.seq[start:stop]) #- average_gc
                if stop > len(bientry.seq):
                    stop = len(bientry.seq)
                if stop - start < 500:
                    break
                #print("gc-diff", gc)
                linedata_data.append({"block_id": str(index), "position": start, "value": gc}) 
                linedata_data.append({"block_id": str(index), "position": stop, "value": gc}) 
          
        if not self.last_inner_radius:
            outer_radius = 0.98
            inner_radius = 0.98 - radius_diff
            #conf["color"] = color
        else:
            outer_radius = self.last_inner_radius - sep
            inner_radius  = self.last_inner_radius - sep - radius_diff
            #conf["color"] = color
          
        self.last_inner_radius = inner_radius

        conf = self.lineplot_configuration_template % (outer_radius, inner_radius, fillcolor, fillcolor)
              
        self.lineplot_tracks += f'"{label}": [%s,%s],' % (linedata_data, conf)
 
 
    def add_histogram_track(self, df, label, fillcolor="blue", sep=0, radius_diff=0.06, outer=False):
        '''
        Minimal Input df columns:
        - bioentry_id
        - start_pos 
        - end_pos 
        - value
        '''
        histogram_data = [{"block_id": f'{row.bioentry_id}', "start":row.start_pos, "end": row.end_pos, "value": row.value} for n, row in df.iterrows()]
        
        if not outer:
            if not self.last_inner_radius:
                outer_radius = 0.98
                inner_radius = 0.98 - radius_diff
                #conf["color"] = color
            else:
                outer_radius = self.last_inner_radius - sep
                inner_radius  = self.last_inner_radius - sep - radius_diff
                #conf["color"] = color
            self.last_inner_radius = inner_radius
        else:
            inner_radius = self.last_outer_radius + sep
            outer_radius  = self.last_outer_radius + sep + radius_diff
            self.last_outer_radius = outer_radius
          
        

        col = '''function(datum, index) {console.log(datum); if (datum.value < 6) {return "red"} else {return "lightblue"}}'''
        print(col)

        conf = self.histogram_configuration_template % (outer_radius, inner_radius, col)
              
        self.histogram_tracks += f'"{label}": [%s,%s],' % (histogram_data, conf)
    
    
    def get_js_code(self,):
        
        if len(self.heatmap_tracks)> 1:
            heat =  self.heatmap_tracks[0:-1] + '}'
            heatmap_code = self.template_heatmap % heat
        else:
            heatmap_code = ''
            
        if len(self.highlight_tracks) > 1:           
            high = self.highlight_tracks[0:-1] + '}'
            highlight_code = self.template_highlight % high
        else:
            highlight_code = ''
        
        if len(self.lineplot_tracks) > 1:
            line = self.lineplot_tracks[0:-1] + '}'
            lineplot_code = self.template_lineplot % line
        else:
            lineplot_code = ''
            
        if len(self.histogram_tracks) > 1:
            histogram = self.histogram_tracks[0:-1] + '}'
            histogram_code = self.template_histogram % histogram
        else:
            histogram_code = ''
        
        return self.js_template % (self.contigs_data, heatmap_code, highlight_code, lineplot_code, histogram_code)
