#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2016
# ---------------------------------------------------------------------------

def successive_cutof_search(biodb,
                            dist, 
                            orthogroup,
                            cutof1, 
                            cutof2, 
                            cutof3,
                            cutof4):
    too_much_hits = False
    all_groups_profile = find_profile_links_recusrsive(biodb, [orthogroup], cutof1, 0, dist)
    cutoff = cutof1
    if all_groups_profile == False:
        all_groups_profile = find_profile_links_recusrsive(biodb, [orthogroup], cutof2, 0, dist)
        cutoff = cutof2
        if all_groups_profile == False:
            #print 'cotoff 0 #######################'
            all_groups_profile = find_profile_links_recusrsive(biodb, [orthogroup], cutof3, 0, dist)
            cutoff = cutof3
            #print all_groups_profile
            if all_groups_profile == False:
                #print 'cotoff 0 #######################'
                all_groups_profile = find_profile_links_recusrsive(biodb, [orthogroup], cutof4, 0, dist)
                cutoff = cutof4
                #print all_groups_profile
                if all_groups_profile == False:
                    all_groups_profile = []
                    too_much_hits = True
    return all_groups_profile, cutoff, too_much_hits



def find_links_recusrsive(biodb, all_connected_seqfeatures, ratio_cutoff=0.5, n_comp_cutoff=1):

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    if biodb == 'chlamydia_04_16':
        filter = '"' + '","'.join(all_connected_seqfeatures) + '"'

        sql = 'select locus_1, locus_2 from interactions_colocalization_table_locus where (locus_1 in (%s) or locus_2 in (%s)) and ' \
              ' (ratio >= %s and n_comparisons >= %s)' % (filter, filter, ratio_cutoff, n_comp_cutoff)
    else:
        all_connected_seqfeatures = [str(i) for i in all_connected_seqfeatures]
        filter = ','.join(all_connected_seqfeatures)

        sql = 'select locus_1, locus_2 from interactions_colocalization_table_locus where (locus_1 in (%s) or locus_2 in (%s)) and ' \
              ' (ratio >= %s and n_comparisons >= %s)' % (filter, filter, ratio_cutoff, n_comp_cutoff)

        #print sql
    data = server.adaptor.execute_and_fetchall(sql,)

    all_groups = []
    for i in data:
        if i[0] not in all_groups:
            all_groups.append(i[0])
        if i[1] not in all_groups:
            all_groups.append(i[1])
    if len(all_groups) > len(all_connected_seqfeatures):
        return find_links_recusrsive(biodb, all_groups, ratio_cutoff, n_comp_cutoff)
    else:
        return all_groups


def find_profile_links_recusrsive(biodb, 
                                  all_connected_groups, 
                                  max_dist=1, 
                                  max_iterations=0,
                                  distance="eucl"):
    '''
    all_connected_groups: orthogroup list under consideration
    max_dist : distance cutoff
    max_iterations : max number of iteration search to identify connections
    distance : either eucl (euclidean) or jac (jaccard)
    '''
    
    if max_iterations >= 1:
        return all_connected_groups
    else:
        max_iterations+=1
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    filter = '"' + '","'.join(all_connected_groups) + '"'
    sql = 'select * from interactions_phylo_profiles_%s_dist where (group_1 in (%s) or group_2 in (%s)) and ' \
          ' (euclidian_dist <= %s)' % (distance,
                                       filter, 
                                       filter, 
                                       max_dist)


    data = server.adaptor.execute_and_fetchall(sql,)
    #print 'n hits', len(data)
    if len(data)>30:
        #print 'too much hits, return False'
        return False
    all_groups = []
    for i in data:
        if i[0] not in all_groups:
            all_groups.append(i[0])
        if i[1] not in all_groups:
            all_groups.append(i[1])
    if len(all_groups) > len(all_connected_groups):
        return find_profile_links_recusrsive(biodb, 
                                             all_groups, 
                                             max_dist, 
                                             max_iterations, 
                                             distance=distance)
    else:
        return all_groups


class orthogroup2network:

    def __init__(self, orthogroup, threshold=995):
        import MySQLdb
        import os
        sqlpsw = os.environ['SQLPSW']

        mysql_host = 'localhost'
        mysql_user = 'root'
        mysql_pwd = sqlpsw
        mysql_db = 'string'
        conn = MySQLdb.connect(host=mysql_host, # your host, usually localhost
                               user=mysql_user, # your username
                               passwd=mysql_pwd, # your password
                               db=mysql_db) # name of the data base
        self.cursor = conn.cursor()
        self.reference_link_data = self.get_reference_link_data(orthogroup, threshold)

        #print self.link_data

    def get_reference_link_data(self, orthogroup, threshold):
        sql = 'select * from detailed_cog_links_v10 where group1="%s" and combined_score >%s' % (orthogroup, threshold)
        self.cursor.execute(sql)
        link_data = self.cursor.fetchall()
        return link_data

    def get_group2_link_data(self,threshold):
        group2 = [i[1] for i in self.reference_link_data]
        sql_filter = '(%)' % ''.join(group2)
        self.group2_link_data = []
        for group in group2:
            sql = 'select * from detailed_cog_links_v10 where group1="%s" and group2 in %s and combined_score >%s' % (group, sql_filter, threshold)
            self.cursor.execute(sql)
            self.group2_link_data.append(sql)

        self.main_cotoscape_js_template = '''

            $(function(){ // on dom ready

            var cy = cytoscape({
              layout: {
                name: 'cose',
                padding: 150
              },
              container: document.getElementById('cy'),
              boxSelectionEnabled: false,
              autounselectify: true,

              style: cytoscape.stylesheet()
                .selector('node')
                  .css({
                    'shape': 'data(faveShape)',
                    'width': 'mapData(weight, 100, 200, 100, 10)',
                    'content': 'data(name)',
                    'text-valign': 'center',
                    'text-outline-width': 0,
                    'text-outline-color': 'data(faveColor)',
                    'background-color': 'data(faveColor)',
                    'color': '#000000'
                  })
                .selector(':selected')
                  .css({
                    'border-width': 3,
                    'border-color': '#333'
                  })
                .selector('edge')
                  .css({
                    'opacity': 0.666,
                    'width': 'mapData(strength, 7, 7, 4, 7)',
                    /*'target-arrow-shape': 'triangle',*/
                    /*'source-arrow-shape': 'circle',*/
                    'line-color': 'data(faveColor)',
                    'source-arrow-color': 'data(faveColor)',
                    'target-arrow-color': 'data(faveColor)'
                  })
                .selector('edge.questionable')
                  .css({
                    'line-style': 'dotted',
                    'target-arrow-shape': 'diamond'
                  })
                .selector('.faded')
                  .css({
                    'opacity': 0.25,
                    'text-opacity': 0
                  }),

              elements: {
                nodes: [
               %s
              ],
              edges: [
              %s
              ]
            },


            /*ready: function(){
              window.cy = this;

              // giddy up
            }*/
          });

          cy.$('node').on('click', function(e){
            var ele = e.cyTarget;
            console.log('clicked ' + ele.id());
          });

            cy.nodes().forEach(function(n){
              var g = n.data('name');

              n.qtip({
                content: [
                  {
                    name: 'NCBI COGS',
                    url: 'http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=' + g
                  },
                  {
                    name: 'eggnog',
                    url: 'http://eggnogdb.embl.de/#/app/results?target_nogs='+ g
                  },
                  {
                    name: 'STRING',
                    url: 'http://string.embl.de/newstring_cgi/show_network_section.pl?all_channels_on=1&interactive=yes&network_flavor=evidence&targetmode=cogs&identifier=' + g
                  }
                ].map(function( link ){
                  return '<a target="_blank" href="' + link.url + '">' + link.name + '</a>';
                }).join('<br />\\n'),
                position: {
                  my: 'bottom left',
                  at: 'top right'
                },
                style: {

                  classes: 'qtip-bootstrap',
                  width: 180,
                  height: 130,
                  lineHeight: 5.5,
                  tip: {
                    width: 22,
                    height: 22,
                  }
                }
              });
            });

          }); // on dom ready
        '''

        self.node_template = '''
        { data: { id: '%s', name: '%s', weight: 10, faveColor: '#6FB1FC', faveShape: 'ellipse' } }
        ''' % (cog, cog)
        self.edge_template = '''
        { data: { source: '%s', target: '%s', faveColor: '#6FB1FC', strength: %s } }

        ''' % (group1, group2, combined_score)


        def get_node_template():
            # reference node
            node_list = [self.node_template % (self.link_data[0][0])]
            for link in self.link_data:
                node_list.append(self.node_template % (self.link_data[0][1]))
            return ','.join(node_list)



def generate_network(biodb,
                     seqfeature_id_list,
                     target_locus_list,
                     ratio_limit,
                     scale_link=True,
                     width=160,
                     height=70,
                     interpro=False,
                     annot=False):
    from chlamdb.biosqldb import manipulate_biosqldb
    import re


    server, db = manipulate_biosqldb.load_db(biodb)
    filter =  ','.join([str(i) for i in seqfeature_id_list])
    sql = 'select * from interactions_colocalization_table_locus where (locus_1 in (%s) or locus_2 in (%s)) and ratio>=%s;' % (filter, filter, ratio_limit)

    data = server.adaptor.execute_and_fetchall(sql,)

    # filtering singletons
    locus_keep = []
    for i in data:
        if i[0] not in locus_keep:
            locus_keep.append(i[0])
        if i[1] not in locus_keep:
            locus_keep.append(i[1])
    # group_1    | group_2    | n_links | n_comparisons | ratio
    #print 'number of locus keep:', len(locus_keep)
    sources = '''
    <meta charset=utf-8 />
    <title>Visual style</title>
    <meta name="viewport" content="user-scalable=no, initial-scale=1.0, minimum-scale=1.0, maximum-scale=1.0, minimal-ui">

    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap-theme.min.css">
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-slider/4.10.3/css/bootstrap-slider.min.css">
    <script src="http://marvl.infotech.monash.edu/webcola/cola.v3.min.js"></script>
    <script src="cytoscape-cola.js"></script>
    <link href="style.css" rel="stylesheet" />

    '''

    body = '''
<body>
  <div id="cy"></div>
</body>
    '''


    # layout cola
    template = '''

<script>

$(function(){ // on dom ready

var cy = cytoscape({
layout: {
	name: 'cose',
    maxSimulationTime: 6000,
    randomize: true,
    nodeSpacing: 30,
    animate: true
  },
  container: document.getElementById('cy'),
  boxSelectionEnabled: false,
  autounselectify: true,

  style: cytoscape.stylesheet()
    .selector('node')
      .css({
        'shape': 'data(faveShape)',
        'width': %s,
	    'height': %s,
        'content': 'data(name)',
        'text-valign': 'center',
        'text-outline-width': 0,
        'text-outline-color': 'data(faveColor)',
        'background-color': 'data(faveColor)',
        'color': '#000000',
	    'font-size': 18

      })
    .selector(':selected')
      .css({
        'border-width': 3,
        'border-color': '#333'
      })
    .selector('edge')
      .css({
        'opacity': 0.666,
        'width': 'data(strength)',
        /*'target-arrow-shape': 'triangle',*/
        /*'source-arrow-shape': 'circle',*/
        'line-color': 'data(faveColor)',
        'source-arrow-color': 'data(faveColor)',
        'target-arrow-color': 'data(faveColor)',
	'curve-style': 'bezier',
        'control-point-distance': '20px',
	'control-point-weight': '0.7', // '0': curve towards source node, '1': towards target node.

      })
    .selector('edge.questionable')
      .css({
        'line-style': 'dotted',
        'target-arrow-shape': 'diamond'
      })
    .selector('.faded')
      .css({
        'opacity': 0.25,
        'text-opacity': 0
      }),

  elements: {
    nodes: [
      %s
    ],
    edges: [
      %s
    ]
  },


  /*ready: function(){
    window.cy = this;

    // giddy up
  }*/
});



  cy.$('node').forEach(function(n){
    var g = n.data('annot');
    var l = n.data('name');

    n.qtip({
      content: [
        {
          name: g,
          url: '/chlamdb/locusx/%s/' + l + '/T'
        }
      ].map(function( link ){
        return '<a target="_blank" href="' + link.url + '">' + link.name + '</a>';
      }).join('<br />\\n'),
      position: {
        my: 'bottom left',
        at: 'top right'
      },
      style: {

        classes: 'qtip-bootstrap',
        width: 180,
        height: 130,
        lineHeight: 1,
        tip: {
          width: 12,
          height: 12,
        }
      }
    });
  });
  cy.edges().forEach(function(n){
    var g = n.data('n_comp');
    var h = n.data('n_links');
    n.qtip({
      content: [
        {
          name: 'Num:' + g + " / " + h
        },

      ].map(function( link ){
        return '<a target="_blank" href="' + link.url + '">' + link.name + '</a>';
      }).join('<br />\\n'),
      position: {
        my: 'bottom left',
        at: 'top right'
      }
    });

  });

    $("#cy_png").click(function(e) {
      var pngData = cy.png({maxWidth:3000,maxHeight:6000, full: true});

    // put the png data in an img tag
    $('#cy_png').attr('href', pngData);

    });


}); // on dom ready



</script>


    '''

    node_template1 = """{ data: { id: '%s', name: '%s', weight: 1, faveColor: '#%s', faveShape: 'ellipse' } }""" # 6FB1FC
    node_template2 = """{ data: { id: '%s', name: '%s', weight: 1, faveColor: '#%s', faveShape: 'ellipse', annot: '%s' } }""" # 6FB1FC
    edge_template = """{ data: { source: '%s', target: '%s', faveColor: '#%s', strength: %s , n_comp:%s, n_links: %s} }"""
    #print interpro,interpro
    node_list = []
    for node in locus_keep:#locus_tag_list:
        if node in target_locus_list:
            color = "FA58F4"
        else:
            color="6FB1FC"
        if annot:
            #print interpro[node][0],
            try:
                node_annot = '%s<br>%s: %s' % (annot[node][1], interpro[node][0], interpro[node][1])
            except:
                node_annot = '%s<br>' % (annot[node][1])
            node_annot = re.sub('\'', '', node_annot)
            if 'hyp' in node_annot and not 'IPR' in node_annot:
                color="F7BE81"
            node_list.append(node_template2 % (node, node, color, node_annot))

        else:
            node_list.append(node_template1 % (node, node, "6FB1FC"))

    edge_list = []
    # group_1    | group_2    | n_links | n_comparisons | ratio
    for edge in data:
        if float(edge[2])< 5:
            if scale_link:
                edge_list.append(edge_template % (edge[0], edge[1], "47de47", float(edge[2]), edge[2], edge[3]))
            else:
                edge_list.append(edge_template % (edge[0], edge[1], "47de47", 2, edge[2], edge[3]))
        else:
            if scale_link:
                edge_list.append(edge_template % (edge[0], edge[1], "ff0404", (float(edge[2]/6)), edge[2], edge[3]))
            else:
                edge_list.append(edge_template % (edge[0], edge[1], "47de47", 10, edge[2], edge[3]))

    return template % (width, height, ',\n'.join(node_list), ',\n'.join(edge_list), biodb)


def generate_network_string(biodb,
                     locus_tag_list,
                     target_locus_list,
                     taxon_id,
                     ratio_limit,
                     width=160,
                     height=70,
                     interpro=False,
                     annot=False):
    from chlamdb.biosqldb import manipulate_biosqldb
    import re


    server, db = manipulate_biosqldb.load_db(biodb)
    filter = '"' + '","'.join(locus_tag_list) + '"'
    sql = 'select t2.locus_tag,t3.locus_tag,t1.label_1,t1.label_2,t1.global_score ' \
          ' from string_interactions t1 inner join custom_tables_locus2seqfeature_id t2 ' \
          ' on t1.seqfeature_id_1=t2.seqfeature_id inner join custom_tables_locus2seqfeature_id t3 ' \
          ' on t1.seqfeature_id_2=t3.seqfeature_id where t1.taxon_id=%s and global_score>=%s and t2.locus_tag in (%s)' \
          ' and t3.locus_tag in (%s) ;' % (taxon_id,
                                           ratio_limit,
                                           filter,
                                           filter,)

    data = server.adaptor.execute_and_fetchall(sql,)

    # filtering singletons
    locus_keep = []
    for i in data:
        if i[0] not in locus_keep:
            locus_keep.append(i[0])
        if i[1] not in locus_keep:
            locus_keep.append(i[1])
    # 0 locus_tag
    # 1 locus_tag
    # 2 label_1
    # 3 label_2
    # 4 global_score
    #print 'number of locus keep:', len(locus_keep)
    sources = '''
    <meta charset=utf-8 />
    <title>Visual style</title>
    <meta name="viewport" content="user-scalable=no, initial-scale=1.0, minimum-scale=1.0, maximum-scale=1.0, minimal-ui">

    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap-theme.min.css">
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-slider/4.10.3/css/bootstrap-slider.min.css">
    <script src="http://marvl.infotech.monash.edu/webcola/cola.v3.min.js"></script>
    <script src="cytoscape-cola.js"></script>
    <link href="style.css" rel="stylesheet" />

    '''

    body = '''
<body>
  <div id="cy"></div>
</body>
    '''


    # layout cola
    template = '''

<script>

$(function(){ // on dom ready

var cy = cytoscape({
layout: {
	name: 'cose',
    maxSimulationTime: 6000,
    randomize: true,
    nodeSpacing: 30,
    animate: true
  },
  container: document.getElementById('cy'),
  boxSelectionEnabled: false,
  autounselectify: true,

  style: cytoscape.stylesheet()
    .selector('node')
      .css({
        'shape': 'data(faveShape)',
        'width': %s,
	    'height': %s,
        'content': 'data(name)',
        'text-valign': 'center',
        'text-outline-width': 0,
        'text-outline-color': 'data(faveColor)',
        'background-color': 'data(faveColor)',
        'color': '#000000',
	    'font-size': 18

      })
    .selector(':selected')
      .css({
        'border-width': 3,
        'border-color': '#333'
      })
    .selector('edge')
      .css({
        'opacity': 0.666,
        'width': 'data(strength)',
        /*'target-arrow-shape': 'triangle',*/
        /*'source-arrow-shape': 'circle',*/
        'line-color': 'data(faveColor)',
        'source-arrow-color': 'data(faveColor)',
        'target-arrow-color': 'data(faveColor)',
	'curve-style': 'bezier',
        'control-point-distance': '20px',
	'control-point-weight': '0.7', // '0': curve towards source node, '1': towards target node.

      })
    .selector('edge.questionable')
      .css({
        'line-style': 'dotted',
        'target-arrow-shape': 'diamond'
      })
    .selector('.faded')
      .css({
        'opacity': 0.25,
        'text-opacity': 0
      }),

  elements: {
    nodes: [
      %s
    ],
    edges: [
      %s
    ]
  },


  /*ready: function(){
    window.cy = this;

    // giddy up
  }*/
});



  cy.$('node').forEach(function(n){
    var g = n.data('annot');
    var l = n.data('name');

    n.qtip({
      content: [
        {
          name: g,
          url: '/chlamdb/locusx/%s/' + l + '/T'
        }
      ].map(function( link ){
        return '<a target="_blank" href="' + link.url + '">' + link.name + '</a>';
      }).join('<br />\\n'),
      position: {
        my: 'bottom left',
        at: 'top right'
      },
      style: {

        classes: 'qtip-bootstrap',
        width: 180,
        height: 130,
        lineHeight: 1,
        tip: {
          width: 12,
          height: 12,
        }
      }
    });
  });
  cy.edges().forEach(function(n){
    var g = n.data('n_comp');
    var h = n.data('n_links');
    n.qtip({
      content: [
        {
          name: 'Num:' + g + " / " + h
        },

      ].map(function( link ){
        return '<a target="_blank" href="' + link.url + '">' + link.name + '</a>';
      }).join('<br />\\n'),
      position: {
        my: 'bottom left',
        at: 'top right'
      }
    });

  });

    $("#cy_png").click(function(e) {
      var pngData = cy.png({maxWidth:3000,maxHeight:6000, full: true});

    // put the png data in an img tag
    $('#cy_png').attr('href', pngData);

    });


}); // on dom ready



</script>


    '''

    node_template1 = """{ data: { id: '%s', name: '%s', weight: 1, faveColor: '#%s', faveShape: 'ellipse' } }""" # 6FB1FC
    node_template2 = """{ data: { id: '%s', name: '%s', weight: 1, faveColor: '#%s', faveShape: 'ellipse', annot: '%s' } }""" # 6FB1FC
    edge_template = """{ data: { source: '%s', target: '%s', faveColor: '#%s', strength: %s , global_score:%s} }"""
    #print interpro,interpro
    node_list = []
    for node in locus_keep:#locus_tag_list:
        if node in target_locus_list:
            color = "FA58F4"
        else:
            color="6FB1FC"
        if annot:
            #print interpro[node][0],
            try:
                node_annot = '%s<br>%s: %s' % (annot[node][1], interpro[node][0], interpro[node][1])
            except:
                node_annot = '%s<br>' % (annot[node][1])
            node_annot = re.sub('\'', '', node_annot)
            if 'hyp' in node_annot and not 'IPR' in node_annot:
                color="F7BE81"
            node_list.append(node_template2 % (node, node, color, node_annot))

        else:
            node_list.append(node_template1 % (node, node, "6FB1FC"))

    edge_list = []
    # 0 locus_tag
    # 1 locus_tag
    # 2 label_1
    # 3 label_2
    # 4 global_score
    # """{ data: { source: '%s', target: '%s', faveColor: '#%s', strength: %s , global_score:%s} }"""
    for edge in data:
        #if float(edge[2])< 1:
        edge_list.append(edge_template % (edge[0], edge[1], "47de47", (10-float(edge[4]))*2, edge[4]))

        #else:
        #   edge_list.append(edge_template % (edge[0], edge[1], "47de47", 10, edge[2], edge[3]))

    return template % (width, height, ',\n'.join(node_list), ',\n'.join(edge_list), biodb)




def generate_network_profile(biodb,
                             group_list,
                             reference_group_list,
                             euclidian_distance_limit=1.5,
                             scale_link=True,
                             interpro=False,
                             annot=False,
                             distance="eucl"):
    
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    filter = '"' + '","'.join(group_list) + '"'
    # group_1     | group_2     | euclidian_dist
    sql = 'select * from interactions_phylo_profiles_%s_dist where (group_1 in (%s) and group_2 in (%s)) and euclidian_dist<=%s;' % (distance,
                                                                                                                                      filter, 
                                                                                                                                      filter, 
                                                                                                                                      euclidian_distance_limit)

    data = server.adaptor.execute_and_fetchall(sql,)

    #print '##############'
    #print data

    # group_1    | group_2    | n_links | n_comparisons | ratio

    sources = '''
    <meta charset=utf-8 />
    <title>Visual style</title>
    <meta name="viewport" content="user-scalable=no, initial-scale=1.0, minimum-scale=1.0, maximum-scale=1.0, minimal-ui">
    <script src="http://ajax.googleapis.com/ajax/libs/jquery/1/jquery.min.js"></script>
    <script src="http://cytoscape.github.io/cytoscape.js/api/cytoscape.js-latest/cytoscape.min.js"></script>

<script src="http://cdnjs.cloudflare.com/ajax/libs/qtip2/2.2.0/jquery.qtip.min.js"></script>
<link href="http://cdnjs.cloudflare.com/ajax/libs/qtip2/2.2.0/jquery.qtip.min.css" rel="stylesheet" type="text/css" />
<script src="https://cdn.rawgit.com/cytoscape/cytoscape.js-qtip/2.7.0/cytoscape-qtip.js"></script>

    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap-theme.min.css">
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/bootstrap-slider/4.10.3/css/bootstrap-slider.min.css">

    <script src="cytoscape-cola.js"></script>
    <link href="style.css" rel="stylesheet" />

    '''

    body = '''
<body>
  <div id="cy"></div>
</body>
    '''


    # breadthfirst
    #maxSimulationTime: 6000,
    #randomize: true,
    #nodeSpacing: 50,
    #animate: true,
    template = '''

<script>

$(function(){ // on dom ready

var cy = cytoscape({
layout: {
	name: 'cose',
    idealEdgeLength: 100,
    nodeOverlap: 200,
    nodeSpacing: 50,
    randomize: true
  },
  container: document.getElementById('cy'),
  boxSelectionEnabled: false,
  autounselectify: true,

  style: cytoscape.stylesheet()
    .selector('node')
      .css({
        'shape': 'data(faveShape)',
        'width': 80,
	    'height': 50,
        'content': 'data(name)',
        'text-valign': 'center',
        'text-outline-width': 0,
        'text-outline-color': 'data(faveColor)',
        'background-color': 'data(faveColor)',
        'color': '#000000',
	    'font-size': 12

      })
    .selector(':selected')
      .css({
        'border-width': 3,
        'border-color': '#333'
      })
    .selector('edge')
      .css({
        'opacity': 0.666,
        'width': 'data(strength)',
        /*'target-arrow-shape': 'triangle',*/
        /*'source-arrow-shape': 'circle',*/
        'line-color': 'data(faveColor)',
        'source-arrow-color': 'data(faveColor)',
        'target-arrow-color': 'data(faveColor)',
	'curve-style': 'bezier',
        'control-point-distance': '20px',
	'control-point-weight': '0.7', // '0': curve towards source node, '1': towards target node.

      })
    .selector('edge.questionable')
      .css({
        'line-style': 'dotted',
        'target-arrow-shape': 'diamond'
      })
    .selector('.faded')
      .css({
        'opacity': 0.25,
        'text-opacity': 0
      }),

  elements: {
    nodes: [
      %s
    ],
    edges: [
      %s
    ]
  },


  /*ready: function(){
    window.cy = this;

    // giddy up
  }*/
});


  cy.nodes().forEach(function(n){
    var g = n.data('annot');

    n.qtip({
      content: [
        {
          name: g
        },
        {
          name: 'eggnog',
          url: 'http://eggnogdb.embl.de/#/app/results?target_nogs='+ g
        },
        {
          name: 'STRING',
          url: 'http://string.embl.de/newstring_cgi/show_network_section.pl?all_channels_on=1&interactive=yes&network_flavor=evidence&targetmode=cogs&identifier=' + g
        }
      ].map(function( link ){
        return '<a target="_blank" href="' + link.url + '">' + link.name + '</a>';
      }).join('<br />\\n'),
      position: {
        my: 'bottom left',
        at: 'top right'
      },
      style: {
        classes: 'qtip-bootstrap',
        width: 180,
        height: 130,
        lineHeight: 5.5,
        tip: {
          width: 22,
          height: 22,
        }
      }
    });
  });
  cy.edges().forEach(function(n){
    var g = n.data('n_comp');
    var h = n.data('n_links');
    n.qtip({
      content: [
        {
          name: 'Num:' + g + " / " + h
        },

      ].map(function( link ){
        return '<a target="_blank" href="' + link.url + '">' + link.name + '</a>';
      }).join('<br />\\n'),
      position: {
        my: 'bottom left',
        at: 'top right'
      },
      style: {
        style: {classes: "MyQtip"},
      }
    });
  });




}); // on dom ready



</script>


    '''
    node_template = """{ data: { id: '%s', name: '%s', weight: 1, faveColor: '#%s', faveShape: 'ellipse'} }""" # 6FB1FC
    edge_template = """{ data: { source: '%s', target: '%s', faveColor: '#%s', strength: %s , n_comp:%s, n_links: %s} }"""

    node_list = []
    for node in group_list:
        if node in reference_group_list:

            node_list.append(node_template % (node, node, "fd5713"))
        else:
            if annot:
                node_annot = '%s\t%s' % (annot[node][0],annot[node][1])
                node_list.append(node_template % (node, node, "6FB1FC", ""))
            else:
                node_list.append(node_template % (node, node, "6FB1FC"))

    edge_list = []
    # group_1    | group_2    | n_links | n_comparisons | ratio
    for edge in data:
        if edge[0] == edge[1]:
            # self link
            continue
        if scale_link:
            if euclidian_distance_limit == 0:
                fscale = 1
            else:
                fscale = euclidian_distance_limit
            edge_list.append(edge_template % (edge[0], edge[1], "ff0404", (fscale-float(edge[2]) + 0.1)*4, edge[2], euclidian_distance_limit))
        else:
            edge_list.append(edge_template % (edge[0], edge[1], "ff0404", 10, edge[2], euclidian_distance_limit))
    return template % (',\n'.join(node_list), ',\n'.join(edge_list))


def get_subgraph(biodb, seqfeature_id_list, ratio_limit, target_locus):

    import networkx as nx
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)
    myfilter =  ','.join([str(i) for i in seqfeature_id_list])
    sql = 'select * from interactions_colocalization_table_locus where (locus_1 in (%s) or locus_2 in (%s)) and ratio>=%s;' % (myfilter, myfilter, ratio_limit)

    data = server.adaptor.execute_and_fetchall(sql,)
    #print 'ex', data[0]

    all_verticles = []
    for i in data:
        all_verticles.append(i[0])
        all_verticles.append(i[1])
    all_verticles = list(set(all_verticles))

    G = nx.Graph()
    G.add_nodes_from(all_verticles)

    label = {}
    for n, orthogroup in enumerate(all_verticles):
        label[n] = orthogroup # orthology_chlamydia_04_16
    #print label[10], label[2388]

    for i in data:
        index_1 = i[0]
        index_2 = i[1]
        weight = i[2]
        G.add_edge(index_1, index_2, weight=weight)

    graph_list = [i for i in nx.connected_component_subgraphs(G)]
    #print 'n subgraphs', len(graph_list)
    locus_keep = []
    for i in range(0, len(graph_list)):

        n_list = graph_list[i].node.keys()
        common = set(n_list).intersection(target_locus)
        if len(common) > 0:
            locus_keep+=n_list
    return locus_keep


if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--biodb', type=str, help="biodb")


    args = parser.parse_args()
    #print 'bonjour'
    #n = orthogroup2network("COG0593")

    #all_groups = find_links_recusrsive('chlamydia_07_16_austr', ["group_1808"], 0.5)
    #print generate_network('chlamydia_07_16_austr', all_groups, "group_1808")
    #print all_groups

    #all_groups = find_profile_links_recusrsive('chlamydia_07_16_austr', ["group_1808"])


    #print all_groups
