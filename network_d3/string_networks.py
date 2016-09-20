#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-  


# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2016
# ---------------------------------------------------------------------------


def find_links_recusrsive(biodb, all_connected_groups, ratio_cutoff=0.5, n_comp_cutoff=1):

    import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    filter = '"' + '","'.join(all_connected_groups) + '"'
    sql = 'select group_1, group_2 from interactions.colocalization_table_%s where (group_1 in (%s) or group_2 in (%s)) and ' \
          ' (ratio >= %s and n_comparisons >= %s)' % (biodb, filter, filter, ratio_cutoff, n_comp_cutoff)

    print sql

    data = server.adaptor.execute_and_fetchall(sql,)

    all_groups = []
    for i in data:
        if i[0] not in all_groups:
            all_groups.append(i[0])
        if i[1] not in all_groups:
            all_groups.append(i[1])
    if len(all_groups) > len(all_connected_groups):
        return find_links_recusrsive(biodb, all_groups, ratio_cutoff, n_comp_cutoff)
    else:
        return all_groups


def find_profile_links_recusrsive(biodb, all_connected_groups, max_euclidian_dist=1):
    print 'find prifiles'
    import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    filter = '"' + '","'.join(all_connected_groups) + '"'
    sql = 'select * from comparative_tables.phylogenetic_profiles_euclidian_distance_%s where (group_1 in (%s) or group_2 in (%s)) and ' \
          ' (euclidian_dist <= %s)' % (biodb, filter, filter, max_euclidian_dist)

    print sql

    data = server.adaptor.execute_and_fetchall(sql,)
    print 'n hits', len(data)
    if len(data)>50:
        return False
    all_groups = []
    for i in data:
        if i[0] not in all_groups:
            all_groups.append(i[0])
        if i[1] not in all_groups:
            all_groups.append(i[1])
    if len(all_groups) > len(all_connected_groups):
        return find_profile_links_recusrsive(biodb, all_groups, max_euclidian_dist)
    else:
        return all_groups


class orthogroup2network:

    def __init__(self, orthogroup, threshold=995):
        import MySQLdb
        mysql_host = 'localhost'
        mysql_user = 'root'
        mysql_pwd = 'estrella3'
        mysql_db = 'string'
        conn = MySQLdb.connect(host=mysql_host, # your host, usually localhost
                               user=mysql_user, # your username
                               passwd=mysql_pwd, # your password
                               db=mysql_db) # name of the data base
        self.cursor = conn.cursor()
        self.reference_link_data = self.get_reference_link_data(orthogroup, threshold)
        
        print self.link_data
        
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


    
def generate_network(biodb, group_list, reference_group, ratio_limit):
    import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    filter = '"' + '","'.join(group_list) + '"'
    sql = 'select * from interactions.colocalization_table_%s where (group_1 in (%s) or group_2 in (%s)) and ratio>=%s;' % (biodb, filter, filter, ratio_limit)

    data = server.adaptor.execute_and_fetchall(sql,)



    # group_1    | group_2    | n_links | n_comparisons | ratio

    sources = '''
    <meta charset=utf-8 />
    <title>Visual style</title>
    <meta name="viewport" content="user-scalable=no, initial-scale=1.0, minimum-scale=1.0, maximum-scale=1.0, minimal-ui">
    <script src="http://ajax.googleapis.com/ajax/libs/jquery/1/jquery.min.js"></script>
    <script src="http://cytoscape.github.io/cytoscape.js/api/cytoscape.js-latest/cytoscape.min.js"></script>

    <script src="http://cdnjs.cloudflare.com/ajax/libs/qtip2/2.2.0/jquery.qtip.min.js"></script>
    <link href="http://cdnjs.cloudflare.com/ajax/libs/qtip2/2.2.0/jquery.qtip.min.css" rel="stylesheet" type="text/css" />
    <script src="https://cdn.rawgit.com/cytoscape/cytoscape.js-qtip/2.2.5/cytoscape-qtip.js"></script>

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



    template = '''

<script>

$(function(){ // on dom ready

var cy = cytoscape({
layout: {
	name: 'cola',
    maxSimulationTime: 6000,
    randomize: true,
    nodeSpacing: 50,
    animate: true
  },
  container: document.getElementById('cy'),
  boxSelectionEnabled: false,
  autounselectify: true,

  style: cytoscape.stylesheet()
    .selector('node')
      .css({
        'shape': 'data(faveShape)',
        'width': 160,
	    'height': 70,
        'content': 'data(name)',
        'text-valign': 'center',
        'text-outline-width': 0,
        'text-outline-color': 'data(faveColor)',
        'background-color': 'data(faveColor)',
        'color': '#000000',
	    'font-size': 22

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
    node_template = """{ data: { id: '%s', name: '%s', weight: 1, faveColor: '#%s', faveShape: 'ellipse' } }""" # 6FB1FC
    edge_template = """{ data: { source: '%s', target: '%s', faveColor: '#%s', strength: %s , n_comp:%s, n_links: %s} }"""

    node_list = []
    for node in group_list:
        if node == reference_group:
            node_list.append(node_template % (node, node, "fd5713"))
        else:
            node_list.append(node_template % (node, node, "6FB1FC"))

    edge_list = []
    # group_1    | group_2    | n_links | n_comparisons | ratio
    for edge in data:
        if float(edge[2])< 20:
            edge_list.append(edge_template % (edge[0], edge[1], "47de47", 10, edge[2], edge[3]))
        else:
            edge_list.append(edge_template % (edge[0], edge[1], "ff0404", float(edge[2])/10, edge[2], edge[3]))

    return template % (',\n'.join(node_list), ',\n'.join(edge_list))


def generate_network_profile(biodb, group_list, reference_group, euclidian_distance_limit=1.5):
    import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    filter = '"' + '","'.join(group_list) + '"'
    # group_1     | group_2     | euclidian_dist
    sql = 'select * from comparative_tables.phylogenetic_profiles_euclidian_distance_%s where (group_1 in (%s) or group_2 in (%s)) and euclidian_dist<=%s;' % (biodb, filter, filter, euclidian_distance_limit)

    data = server.adaptor.execute_and_fetchall(sql,)



    # group_1    | group_2    | n_links | n_comparisons | ratio

    sources = '''
    <meta charset=utf-8 />
    <title>Visual style</title>
    <meta name="viewport" content="user-scalable=no, initial-scale=1.0, minimum-scale=1.0, maximum-scale=1.0, minimal-ui">
    <script src="http://ajax.googleapis.com/ajax/libs/jquery/1/jquery.min.js"></script>
    <script src="http://cytoscape.github.io/cytoscape.js/api/cytoscape.js-latest/cytoscape.min.js"></script>

    <script src="http://cdnjs.cloudflare.com/ajax/libs/qtip2/2.2.0/jquery.qtip.min.js"></script>
    <link href="http://cdnjs.cloudflare.com/ajax/libs/qtip2/2.2.0/jquery.qtip.min.css" rel="stylesheet" type="text/css" />
    <script src="https://cdn.rawgit.com/cytoscape/cytoscape.js-qtip/2.2.5/cytoscape-qtip.js"></script>

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



    template = '''

<script>

$(function(){ // on dom ready

var cy = cytoscape({
layout: {
	name: 'cola',
    maxSimulationTime: 6000,
    randomize: true,
    nodeSpacing: 50,
    animate: true
  },
  container: document.getElementById('cy'),
  boxSelectionEnabled: false,
  autounselectify: true,

  style: cytoscape.stylesheet()
    .selector('node')
      .css({
        'shape': 'data(faveShape)',
        'width': 160,
	    'height': 70,
        'content': 'data(name)',
        'text-valign': 'center',
        'text-outline-width': 0,
        'text-outline-color': 'data(faveColor)',
        'background-color': 'data(faveColor)',
        'color': '#000000',
	    'font-size': 22

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
    node_template = """{ data: { id: '%s', name: '%s', weight: 1, faveColor: '#%s', faveShape: 'ellipse' } }""" # 6FB1FC
    edge_template = """{ data: { source: '%s', target: '%s', faveColor: '#%s', strength: %s , n_comp:%s, n_links: %s} }"""

    node_list = []
    for node in group_list:
        if node == reference_group:
            node_list.append(node_template % (node, node, "fd5713"))
        else:
            node_list.append(node_template % (node, node, "6FB1FC"))

    edge_list = []
    # group_1    | group_2    | n_links | n_comparisons | ratio
    for edge in data:
        edge_list.append(edge_template % (edge[0], edge[1], "ff0404", float(edge[2]), edge[2], euclidian_distance_limit))

    return template % (',\n'.join(node_list), ',\n'.join(edge_list))



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