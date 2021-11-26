
// Note: the regions should be a dict containing the following elements:
//   start: start of the region
//   end: end of the region
//   features: list of the different features, with each feature having the following elements:
//   	type:  either CDS, tmRNA, tRNA, rRNA, pseudo
//   	gene: name of the gene, or empty string if none
//   	start: start of the coding region
//   	end: end of the coding region
//   	strand: either +1 or -1
//   	locus_tag
//   highlight: a list of locus tag to highlight in red
function createGenomicRegion(div, regions, connections, highlight, window_size) {
	const text_field_size = 40;
	const margin = { top:5, right: 5, bottom:5, left:5 };
	const max_arrow_size = 350;
	const arrow_tube_size = 12;
	const default_width = 600;
	const base_line_width = 2;
	const regions_vertical_interval = 20;
	const arrow_height = 30;
	const diagram_vertical_size = 2*arrow_height + base_line_width;
	const region_height = text_field_size+diagram_vertical_size;

	var total_height = region_height*regions.length+regions_vertical_interval*(regions.length-1);

	var svg = div
		.append("svg")
		.attr("width",  default_width+margin.right+margin.left)
		.attr("height", total_height+margin.top+margin.bottom)
		.append("g")
		.attr("transform", "translate(" + margin.left + "," + margin.top + ")");

	var Tooltip = div
		.append("div")
		.style("opacity", 0)
		.attr("class", "tooltip")
		.style("background-color", "white")
		.style("border", "solid")
		.style("border-width", "2px")
		.style("border-radius", "5px")
		.style("padding", "5px");

	function get_feature_points(feature, region_start, region_end) {
		let start = feature.start;
		let stop = feature.end;
		let top_y = 0;
		let bot_y = arrow_height;

		if(stop-start <= max_arrow_size) {
			top_left = {x: start , y:top_y};
			bottom_left = {x:start, y:bot_y};
			spike = {x:stop, y:(top_y+bot_y)/2};
			points = [top_left, bottom_left, spike];
		} else {
			let diff = (arrow_height-arrow_tube_size)/2;
			arrow_top = {x: stop-max_arrow_size, y:top_y};
			arrow_bot = {x: stop-max_arrow_size, y:bot_y};
			arrow_spike = {x: stop, y:(top_y+bot_y)/2};
			tube_beg_top = {x:start, y: diff};
			tube_beg_bot = {x:start, y: diff+arrow_tube_size};
			tube_end_top = {x: stop-max_arrow_size, y:diff}
			tube_end_bot = {x: stop-max_arrow_size, y:diff+arrow_tube_size}
			points = [tube_beg_top, tube_end_top, arrow_top, arrow_spike, arrow_bot, tube_end_bot, tube_beg_bot];
		}

		if(feature.strand == -1) {
			let midpoint = (feature.start+feature.end)/2;
			for(point in points) {
				let x_coord = points[point].x;
				points[point].y += arrow_height+base_line_width;
				points[point].x -= 2*(x_coord-midpoint);
				if(points[point].x > region_end) {
					points[point].x = region_end;
				}
			}
		}

		for(var i=0; i<points.length;i++) {
			cur_point = points[i];
			if(cur_point.x < region_start) {
				cur_point.x = region_start;
			} else if (cur_point.x > region_end) {
				cur_point.x = region_end;
			}
		}
		return points;
	}

	// maybe not needed as a function anymore
	function draw_genes_arrow(svg, region, x_scale, y_scale) {
		// if a gene is partly in the region, only draw those
		// whose arrow can be completely drawn, discard the other ones.
		let start=region.start;
		let end=region.end;
		filtered_features = region.features.filter(function(d) {
			if(d.strand==1) {
				return d.end <= end-max_arrow_size;
			} else {
				return d.end >= start+max_arrow_size;
			}
		});

		svg.append("g")
		.selectAll("polygon")
		.data(filtered_features)
		.enter()
		.append("polygon")
		.attr("points", function(d) {
			let points  = get_feature_points(d, region.start, region.end)
				.map(function(x) {
					return [x_scale(x.x), y_scale(x.y)].join(" ");
			});
			return points.join(",");
		})
		.style("stroke-width", 2)
		.style("opacity", .9)
		.style("fill", function(d) {
			if(highlight!=null && highlight.find(el => el===d.locus_tag)!==undefined) {
				return "red";
			} else if(d.type=="CDS") {
				return "green";
			} else if(d.type=="tmRNA" || d.type=="tRNA" || d.type=="rRNA") {
				return "orange";
			} else if(d.type=="pseudo") {
				return "black";
			}
		})
		.on("mouseover", mouseover_feature)
		.on("mouseleave", mouseleave_feature)
		.on("click", click_feature);
	}


	function add_genes_name(svg, region, x_scale, y_text_pos) {
		let filtered_features = region.features.filter(function(d) {
			let start = region.start;
			let stop = region.end;
			return (d.start>=start) & (d.end<=stop);
		});
		svg.selectAll("text")
			.data(filtered_features)
			.enter()
			.append("text")
			.attr("x", d => x_scale((d.start+d.end)/2))
			.attr("y", y_text_pos)
			.attr("text-anchor", "start")
			.attr("transform", function(d ) {
				return "rotate(315,"+x_scale((d.start+d.end)/2)+","+y_text_pos+")";
			})
			.attr("id", d => d.locus_tag)
			.text(d => d.gene);
	}


	function load_axis(svg, start, end, x_scale, y_scale) {
		svg.append("line")
			.style("stroke", "black")
			.style("stroke-width", base_line_width)
			.attr("x1", x_scale(start))
			.attr("y1", y_scale(arrow_height+base_line_width/2))
			.attr("x2", x_scale(end))
			.attr("y2", y_scale(arrow_height+base_line_width/2));

		// border lines
		svg.append("line")
			.style("stroke", "black")
			.style("stroke-width", 2)
			.attr("x1", x_scale(start))
			.attr("y1", y_scale(0))
			.attr("x2", x_scale(start))
			.attr("y2", y_scale(2*arrow_height+base_line_width));

		// border lines
		svg.append("line")
			.style("stroke", "black")
			.style("stroke-width", 2)
			.attr("x1", x_scale(end))
			.attr("y1", y_scale(0))
			.attr("x2", x_scale(end))
			.attr("y2", y_scale(2*arrow_height+base_line_width));

		// background rectangle
		svg.append("rect")
			.attr("x", x_scale(start))
			.attr("y", y_scale(0))
			.attr("width", x_scale(end)-x_scale(start))
			.attr("height", y_scale(2*arrow_height+base_line_width)-y_scale(0))
			.attr("opacity", .1)
			.attr("fill", "gray");
	}


	function mouseover_feature(d) {
		let pos = d3.mouse(this);
		Tooltip.style("opacity", 1)
			.html(d.product + " " + d.locus_tag)
			.style("left", pos[0] + "px")
			.style("top", pos[1] + "px");
		d3.select(this)
			.style("stroke", "red");
		d3.select("#"+d.locus_tag)
			.style("stroke", "red");
	}


	function click_feature(d) {
		window.open("/locusx/"+d.locus_tag);
	}


	function mouseleave_feature(d) {
		Tooltip.style("opacity", 0);
		d3.select(this)
			.style("stroke", "none");
		d3.select("#"+d.locus_tag)
			.style("stroke", "none");
	}


	let max_region_size = d3.max(regions.map(d => d.end-d.start));
	for(var i=0; i<regions.length; i++) {
		let current_region = regions[i];
		let region_size = current_region.end-current_region.start;
		let y_base_pos = region_height*i + i*regions_vertical_interval;
		let ratio = max_region_size/window_size;
		let this_region_width = ratio*(region_size/max_region_size)*default_width;
		let blank_space = default_width-this_region_width;
		let x_scale = d3.scale.linear().
			domain([current_region.start, current_region.end]).
			range([blank_space/2, default_width-(blank_space/2)]);
		let y_scale = d3.scale.linear().
			domain([0, diagram_vertical_size]).
			range([y_base_pos+text_field_size, y_base_pos+diagram_vertical_size+text_field_size]);

		load_axis(svg, current_region.start, current_region.end, x_scale, y_scale);
		draw_genes_arrow(svg, current_region, x_scale, y_scale);
		add_genes_name(svg, current_region, x_scale, y_base_pos+text_field_size);
	}
}
