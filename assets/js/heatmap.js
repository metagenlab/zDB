var classesNumber = 10,
    cellSize = 14;

//#########################################################
function heatmap_display(url, heatmapId, paletteName) {


    //##########################################################################
    // Patrick.Brockmann@lsce.ipsl.fr
    //##########################################################################
    
    //==================================================
    // References
    // http://bl.ocks.org/Soylent/bbff6cc507dca2f48792
    // http://bost.ocks.org/mike/selection/
    // http://bost.ocks.org/mike/join/
    // http://stackoverflow.com/questions/9481497/understanding-how-d3-js-binds-data-to-nodes
    // http://bost.ocks.org/mike/miserables/
    // http://bl.ocks.org/ianyfchang/8119685

    //==================================================
    var tooltip = d3.select(heatmapId)
        .append("div")
        .style("position", "absolute")
        .style("visibility", "hidden");

    //==================================================
    // http://bl.ocks.org/mbostock/3680958
    function zoom() {
    	svg.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
    }

    // define the zoomListener which calls the zoom function on the "zoom" event constrained within the scaleExtents
    var zoomListener = d3.behavior.zoom().scaleExtent([0.1, 3]).on("zoom", zoom);

    //==================================================
    var viewerWidth = $(document).width();
    var viewerHeight = $(document).height();
    var viewerPosTop = 150;
    var viewerPosLeft = 400;

    var legendElementWidth = cellSize * 2;

    // http://bl.ocks.org/mbostock/5577023
    var colors = colorbrewer[paletteName][classesNumber];

    // http://bl.ocks.org/mbostock/3680999
    var svg;

    //==================================================
    d3.json(url, function(error, data) {

        //console.log(data);
        var arr = data.data;
        var row_number = arr.length;
        var col_number = arr[0].length;
        //console.log(col_number, row_number);

        var colorScale = d3.scale.quantize()
            .domain([0.0, 1.0])
            .range(colors);

        svg = d3.select(heatmapId).append("svg")
            .attr("width", viewerWidth)
            .attr("height", viewerHeight)
	    .call(zoomListener)
            .append("g")
            .attr("transform", "translate(" + viewerPosLeft + "," + viewerPosTop + ")");

        svg.append('defs')
            .append('pattern')
            .attr('id', 'diagonalHatch')
            .attr('patternUnits', 'userSpaceOnUse')
            .attr('width', 4)
            .attr('height', 4)
            .append('path')
            .attr('d', 'M-1,1 l2,-2 M0,4 l4,-4 M3,5 l2,-2')
            .attr('stroke', '#000000')
            .attr('stroke-width', 1);

        var rowSortOrder = false;
        var colSortOrder = false;

        var rowLabels = svg.append("g")
            .attr("class", "rowLabels")
            .selectAll(".rowLabel")
            .data(data.index)
            .enter().append("text")
            .text(function(d) {
                return d.count > 1 ? d.join("/") : d;
            })
            .attr("x", 0)
            .attr("y", function(d, i) {
                return (i * cellSize);
            })
            .style("text-anchor", "end")
            .attr("transform", function(d, i) {
                return "translate(-3," + cellSize / 1.5 + ")";
            })
            .attr("class", "rowLabel mono")
            .attr("id", function(d, i) {
                return "rowLabel_" + i;
            })
            .on('mouseover', function(d, i) {
                d3.select('#rowLabel_' + i).classed("hover", true);
            })
            .on('mouseout', function(d, i) {
                d3.select('#rowLabel_' + i).classed("hover", false);
            })
            .on("click", function(d, i) {
                rowSortOrder = !rowSortOrder;
                sortByValues("r", i, rowSortOrder);
                d3.select("#order").property("selectedIndex", 0);
                //$("#order").jqxComboBox({selectedIndex: 0});
            });

        var colLabels = svg.append("g")
            .attr("class", "colLabels")
            .selectAll(".colLabel")
            .data(data.columns)
            .enter().append("text")
            .text(function(d) {
                //d.shift();
                return d;
            })
            .attr("x", 0)
            .attr("y", function(d, i) {
                return (i * cellSize);
            })
            .style("text-anchor", "left")
            .attr("transform", function(d, i) {
                return "translate(" + cellSize / 2 + ", -3) rotate(-90) rotate(45, 0, " + (i * cellSize) + ")";
            })
            .attr("class", "colLabel mono")
            .attr("id", function(d, i) {
                return "colLabel_" + i;
            })
            .on('mouseover', function(d, i) {
                d3.select('#colLabel_' + i).classed("hover", true);
            })
            .on('mouseout', function(d, i) {
                d3.select('#colLabel_' + i).classed("hover", false);
            })
            .on("click", function(d, i) {
                colSortOrder = !colSortOrder;
                sortByValues("c", i, colSortOrder);
                d3.select("#order").property("selectedIndex", 0);
            });

        var row = svg.selectAll(".row")
            .data(data.data)
            .enter().append("g")
            .attr("id", function(d) {
                return d.idx;
            })
            .attr("class", "row");

        var j = 0;
        var heatMap = row.selectAll(".cell")
            .data(function(d) {
                j++;
                return d;
            })
            .enter().append("svg:rect")
            .attr("x", function(d, i) {
                return i * cellSize;
            })
            .attr("y", function(d, i, j) {
                return j * cellSize;
            })
            .attr("rx", 4)
            .attr("ry", 4)
            .attr("class", function(d, i, j) {
                return "cell bordered cr" + j + " cc" + i;
            })
            .attr("row", function(d, i, j) {
                return j;
            })
            .attr("col", function(d, i, j) {
                return i;
            })
            .attr("width", cellSize)
            .attr("height", cellSize)
            .style("fill", function(d) {
                if (d != null) return colorScale(d);
                else return "url(#diagonalHatch)";
            })
            .on('mouseover', function(d, i, j) {
                d3.select('#colLabel_' + i).classed("hover", true);
                d3.select('#rowLabel_' + j).classed("hover", true);
                if (d != null) {
                    tooltip.html('<div class="heatmap_tooltip">' + d.toFixed(3) + '</div>');
                    tooltip.style("visibility", "visible");
                } else
                    tooltip.style("visibility", "hidden");
            })
            .on('mouseout', function(d, i, j) {
                d3.select('#colLabel_' + i).classed("hover", false);
                d3.select('#rowLabel_' + j).classed("hover", false);
                tooltip.style("visibility", "hidden");
            })
            .on("mousemove", function(d, i) {
                tooltip.style("top", (d3.event.pageY - 55) + "px").style("left", (d3.event.pageX - 60) + "px");
            })
            .on('click', function() {
                //console.log(d3.select(this));
            });

        var legend = svg.append("g")
            .attr("class", "legend")
            .attr("transform", "translate(0,-300)")
            .selectAll(".legendElement")
            .data([0.0, 0.1, 0.2,0.3,0.4,0.5, 0.6, 0.7, 0.8, 0.9])
            .enter().append("g")
            .attr("class", "legendElement");

        legend.append("svg:rect")
            .attr("x", function(d, i) {
                return legendElementWidth * i;
            })
            .attr("y", viewerPosTop)
            .attr("class", "cellLegend bordered")
            .attr("width", legendElementWidth)
            .attr("height", cellSize / 2)
            .style("fill", function(d, i) {
                return colors[i];
            });

        legend.append("text")
            .attr("class", "mono legendElement")
            .text(function(d) {
                return "≥" + Math.round(d * 100) / 100;
            })
            .attr("x", function(d, i) {
                return legendElementWidth * i;
            })
            .attr("y", viewerPosTop + cellSize);

        //==================================================
        // Change ordering of cells
        function sortByValues(rORc, i, sortOrder) {
            var t = svg.transition().duration(1000);
            var values = [];
            var sorted;
            d3.selectAll(".c" + rORc + i)
                .filter(function(d) {
                    if (d != null) values.push(d);
                    else values.push(-999); // to handle NaN
                });
            //console.log(values);		
            if (rORc == "r") { // sort on cols
                sorted = d3.range(col_number).sort(function(a, b) {
                    if (sortOrder) {
                        return values[b] - values[a];
                    } else {
                        return values[a] - values[b];
                    }
                });
                t.selectAll(".cell")
                    .attr("x", function(d) {
                        var col = parseInt(d3.select(this).attr("col"));
                        return sorted.indexOf(col) * cellSize;
                    });
                t.selectAll(".colLabel")
                    .attr("y", function(d, i) {
                        return sorted.indexOf(i) * cellSize;
                    })
                    .attr("transform", function(d, i) {
                        return "translate(" + cellSize / 2 + ", -3) rotate(-90) rotate(45, 0, " + (sorted.indexOf(i) * cellSize) + ")";
                    });
            } else { // sort on rows
                sorted = d3.range(row_number).sort(function(a, b) {
                    if (sortOrder) {
                        return values[b] - values[a];
                    } else {
                        return values[a] - values[b];
                    }
                });
                t.selectAll(".cell")
                    .attr("y", function(d) {
                        var row = parseInt(d3.select(this).attr("row"));
                        return sorted.indexOf(row) * cellSize;
                    });
                t.selectAll(".rowLabel")
                    .attr("y", function(d, i) {
                        return sorted.indexOf(i) * cellSize;
                    })
                    .attr("transform", function(d, i) {
                        return "translate(-3," + cellSize / 1.5 + ")";
                    });
            }
        }

        //==================================================
        d3.select("#order").on("change", function() {
	    var newOrder = d3.select("#order").property("value");	
            changeOrder(newOrder, heatmapId);
        });

        //==================================================
        d3.select("#palette")
            .on("keyup", function() {
		var newPalette = d3.select("#palette").property("value");
		if (newPalette != null)						// when interfaced with jQwidget, the ComboBox handles keyup event but value is then not available ?
                	changePalette(newPalette, heatmapId);
            })
            .on("change", function() {
		var newPalette = d3.select("#palette").property("value");
                changePalette(newPalette, heatmapId);
            });
    });

    //==================================================
}

//#########################################################
function changeOrder(newOrder, heatmapId) {
    var svg = d3.select(heatmapId);
    var t = svg.transition().duration(1000);
    if (newOrder == "sortinit_col") { // initial sort on cols (alphabetically if produced like this)
        t.selectAll(".cell")
            .attr("x", function(d) {
                var col = parseInt(d3.select(this).attr("col"));
                return col * cellSize;
            });
        t.selectAll(".colLabel")
            .attr("y", function(d, i) {
                return i * cellSize;
            })
            .attr("transform", function(d, i) {
                return "translate(" + cellSize / 2 + ", -3) rotate(-90) rotate(45, 0, " + (i * cellSize) + ")";
            });
    } else if (newOrder == "sortinit_row") { // initial sort on rows (alphabetically if produced like this)
        t.selectAll(".cell")
            .attr("y", function(d) {
                var row = parseInt(d3.select(this).attr("row"));
                return row * cellSize;
            });
        t.selectAll(".rowLabel")
            .attr("y", function(d, i) {
                return i * cellSize;
            })
            .attr("transform", function(d, i) {
                return "translate(-3," + cellSize / 1.5 + ")";
            });
    } else if (newOrder == "sortinit_col_row") { // initial sort on rows and cols (alphabetically if produced like this)
        t.selectAll(".cell")
            .attr("x", function(d) {
                var col = parseInt(d3.select(this).attr("col"));
                return col * cellSize;
            })
            .attr("y", function(d) {
                var row = parseInt(d3.select(this).attr("row"));
                return row * cellSize;
            });
        t.selectAll(".colLabel")
            .attr("y", function(d, i) {
                return i * cellSize;
            })
            .attr("transform", function(d, i) {
                return "translate(" + cellSize / 2 + ", -3) rotate(-90) rotate(45, 0, " + (i * cellSize) + ")";
            });
        t.selectAll(".rowLabel")
            .attr("y", function(d, i) {
                return i * cellSize;
            })
            .attr("transform", function(d, i) {
                return "translate(-3," + cellSize / 1.5 + ")";
            });
    }
}

//#########################################################
function changePalette(paletteName, heatmapId) {
    var colors = colorbrewer[paletteName][classesNumber];
    var colorScale = d3.scale.quantize()
        .domain([0.0, 1.0])
        .range(colors);
    var svg = d3.select(heatmapId);
    var t = svg.transition().duration(500);
    t.selectAll(".cell")
        .style("fill", function(d) {
                if (d != null) return colorScale(d);
                else return "url(#diagonalHatch)";
        })
    t.selectAll(".cellLegend")
        .style("fill", function(d, i) {
            return colors[i];
        });
}
