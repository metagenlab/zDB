
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
function createGenomicRegion(div, regions, connections, highlight, window_size, ident_range) {
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

	function draw_genes_arrow(svg, region, x_scale, y_scale) {
		// if a gene is partly in the region, only draw those
		// whose arrow can be completely drawn, discard the other ones.
		let start=region.start;
		let end=region.end;
		let filtered_features = region.features.filter(function(d) {
			if(d.strand==1) {
				return d.end <= end-max_arrow_size;
			} else {
				return d.end >= start+max_arrow_size;
			}
		});
		let locus_to_position = {};

		svg.append("g")
		.selectAll("polygon")
		.data(filtered_features)
		.enter()
		.append("polygon")
		.attr("points", function(d) {
			let points = get_feature_points(d, region.start, region.end);

			// ugly
			locus_to_position[d.locus_tag] = [x_scale(d.start), x_scale(d.end)];
			return points.map(x => x_scale(x.x)+" "+y_scale(x.y)).join(",");
		})
		.attr("id", d => d.locus_tag)
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

		return locus_to_position;
	}

	function mouseover_link(d) {
		d3.select(this)
			.style("stroke", "blue")
			.style("stroke-width", 2);
		let top_locus_tag = d3.select("#"+d.top_feature);
		let bot_locus_tag = d3.select("#"+d.bottom_feature);
		top_locus_tag.style("stroke", "blue");
		bot_locus_tag.style("stroke", "blue");
		let pos = d3.mouse(this);
		Tooltip.style("opacity", 1)
			.html(d.ident + "% identity")
			.style("left", pos[0] + "px")
			.style("top", pos[1] + "px");
	}

	function mouseleave_link(d) {
		Tooltip.style("opacity", 0).
			style("left", "-1px").
			style("top", "-1px");
		d3.select(this)
			.style("stroke", "none");
		let top_locus_tag = d3.select("#"+d.top_feature);
		let bot_locus_tag = d3.select("#"+d.bottom_feature);
		top_locus_tag.style("stroke", "none");
		bot_locus_tag.style("stroke", "none");
		Tooltip.style("opacity", 0);
	}

	function mouseclick_link(d) {
		window.open("/orthogroup/group_"+d.group);
	}

	function add_genes_name(svg, region, x_scale, y_text_pos) {
		let filtered_features = region.features.filter(function(d) {
			let start = region.start;
			let stop = region.end;
			return (d.start>=start) & (d.end<=stop);
		});
		svg.append("g")
			.selectAll("text")
			.data(filtered_features)
			.enter()
			.append("text")
			.attr("x", d => x_scale((d.start+d.end)/2))
			.attr("y", y_text_pos)
			.attr("text-anchor", "start")
			.attr("transform", function(d ) {
				return "rotate(315,"+x_scale((d.start+d.end)/2)+","+y_text_pos+")";
			})
			.attr("id", d => d.locus_tag+"_gene_name")
			.text(d => d.gene);
	}


	function load_axis(svg, region, x_scale, y_scale) {
		start = region.start;
		end = region.end;
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

		if("name" in region) {
			svg.append("text").
				attr("x", x_scale(start)+3).
				attr("y", y_scale(2*arrow_height+base_line_width)).
				attr("fill", "gray").
				attr("opacity", 1).
				style("font-size", "10px").
				text(region.name);
		}
	}


	function mouseover_feature(d) {
		let pos = d3.mouse(this);
		Tooltip.style("opacity", 1)
			.html(d.product)
			.style("left", pos[0] + "px")
			.style("top", pos[1] + "px");
		d3.select(this)
			.style("stroke", "red");
		d3.select("#"+d.locus_tag+"_gene_name")
			.style("stroke", "red");
	}


	function click_feature(d) {
		window.open("/locusx/"+d.locus_tag);
	}


	function mouseleave_feature(d) {
		Tooltip.style("opacity", 0).
			style("left", "-1px").
			style("top", "-1px");
		d3.select(this)
			.style("stroke", "none");
		d3.select("#"+d.locus_tag+"_gene_name")
			.style("stroke", "none");
	}


	let max_region_size = d3.max(regions.map(d => d.end-d.start));
	let prev_gene_pos = null;
	let prev_region = null;
	let prev_xscale = null;
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

		load_axis(svg, current_region, x_scale, y_scale);
		let curr_gene_pos = draw_genes_arrow(svg, current_region, x_scale, y_scale);
		let this_g = svg.append("g");
		if(i>=1 && connections != null) {
			// messy code... I certainly hope to not have to debug it...

			let ident_scale = d3.scale.linear().
				domain(ident_range).
				range([.2, .9]);
			let curr_connections = connections[i-1];
			let prev_start = prev_xscale(prev_region.start);
			let this_start = x_scale(current_region.start);
			let prev_end = prev_xscale(prev_region.end);
			let this_end = x_scale(current_region.end);
			for(var locus_tag in curr_gene_pos) {
				if(!(locus_tag in curr_connections)) {
					continue;
				}
				let connection_data = curr_connections[locus_tag];
				let connects_to = connection_data[0];
				let orthogroup = connection_data[1];
				let identity = connection_data[2];
				if(!(connects_to in prev_gene_pos)) {
					continue;
				}
				let top_pos = prev_gene_pos[connects_to];
				let curr_pos = curr_gene_pos[locus_tag];
				let top_pos0 = top_pos[0];
				let top_pos1 = top_pos[1];
				let curr_pos1 = curr_pos[1];
				let curr_pos0 = curr_pos[0];
				let all_points = [];
				if(top_pos0<prev_start && curr_pos0<this_start) {
					all_points.push([this_start, y_base_pos+text_field_size]);
					all_points.push([prev_start, y_base_pos-regions_vertical_interval]);
				} else if(top_pos0<prev_start) {
					slope = (regions_vertical_interval+text_field_size)/(curr_pos0-top_pos0);
					all_points.push([curr_pos0, y_base_pos+text_field_size]);
					all_points.push([prev_start,
						y_base_pos-(slope*(curr_pos0-prev_start))]);
					all_points.push([prev_start,
						y_base_pos-regions_vertical_interval]);
				} else if(curr_pos0<this_start) {
					slope = (regions_vertical_interval+text_field_size)/(curr_pos0-top_pos0);
					all_points.push([this_start, y_base_pos+text_field_size]);
					all_points.push([this_start,
						y_base_pos-regions_vertical_interval+slope*(this_start-top_pos0)]);
					all_points.push([top_pos0, y_base_pos-regions_vertical_interval]);
				} else {
					all_points.push([curr_pos0, y_base_pos+text_field_size]);
					all_points.push([top_pos0, y_base_pos-regions_vertical_interval]);
				}

				if(top_pos1>prev_end && curr_pos1>this_end) {
					all_points.push([prev_end, y_base_pos-regions_vertical_interval]);
					all_points.push([this_end, y_base_pos+text_field_size]);
				}else if (curr_pos1>this_end) {
					slope = (regions_vertical_interval+text_field_size)/(curr_pos1-top_pos1); 
					all_points.push([top_pos1, y_base_pos-regions_vertical_interval]);
					all_points.push([this_end, y_base_pos-regions_vertical_interval
						+slope*(this_end-top_pos1)]);
					all_points.push([this_end, y_base_pos+text_field_size]);
				} else if(top_pos1>prev_end) {
					slope = (regions_vertical_interval+text_field_size)/(curr_pos1-top_pos1);
					all_points.push([prev_end, y_base_pos-regions_vertical_interval]);
					all_points.push([prev_end, y_base_pos+text_field_size
						-slope*(curr_pos1-prev_end)]);
					all_points.push([curr_pos1, y_base_pos+text_field_size]);
				} else {
					all_points.push([top_pos1, y_base_pos-regions_vertical_interval]);
					all_points.push([curr_pos1, y_base_pos+text_field_size]);
				}

				this_g.append("polygon")
					.datum({group: orthogroup, top_feature: connects_to,
						bottom_feature: locus_tag, ident: identity})
					.attr("points",
						all_points.map(d => d[0]+" "+d[1]).join(",")
					)
					.style("fill-opacity", ident_scale(identity))
					.style("fill", "gray")
					.style("stroke-opacity", .5)
					.on("mouseover", mouseover_link)
					.on("mouseleave", mouseleave_link)
					.on("click", mouseclick_link);
			}
		}
		add_genes_name(svg, current_region, x_scale, y_base_pos+text_field_size);
		prev_gene_pos = curr_gene_pos;
		prev_region = current_region;
		prev_xscale = x_scale;
	}
}
