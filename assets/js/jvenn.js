/***************************************************************
*  Copyright notice
*
*  (c) 2014 PF bioinformatique de Toulouse
*  All rights reserved
* 
*
*  This script is an adaptation of the venny script developed by
*  Juan Carlos Oliveros, BioinfoGP, CNB-CSIC:
*  Oliveros, J.C. (2007) VENNY. An interactive tool for comparing 
*  lists with Venn Diagrams.
*  http://bioinfogp.cnb.csic.es/tools/venny/index.html.
*  It is distributed under the terms of the GNU General Public 
*  License as published by the Free Software Foundation; either 
*  version 2 of the License, or (at your option) any later version.
*
*  The GNU General Public License can be found at
*  http://www.gnu.org/copyleft/gpl.html.
*
*  This script is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  This copyright notice MUST APPEAR in all copies of the script!
***************************************************************/
(function($) {
	$.fn.jvenn = function(options) {
        var defaults = {
            series: [{
            	name: 'Actors',
            	data: ["Marilyn Monroe", "Arnold Schwarzenegger", "Jack Nicholson", "Barbra Streisand", "Robert de Niro", "Dean Martin", "Harrison Ford"]
            }, {
            	name: 'Singers',
            	data: ["Freddy Mercury", "Barbra Streisand", "Dean Martin", "Ricky Martin", "Celine Dion", "Marilyn Monroe"]
            }],
//            series: [{
//            	name: 'sample1',
//            	data: ["Otu1", "Otu2", "Otu3", "Otu4", "Otu5", "Otu6", "Otu7"],
//            	values: [5, 15, 250, 20, 23, 58, 89]
//            }, {
//            	name: 'sample2',
//            	data: ["Otu1", "Otu2", "Otu5", "Otu7", "Otu8", "Otu9"],
//            	values: [90, 300, 10, 2, 45, 9]
//            }],
//            series: [{
//            	name: {A: 'Actors',
//                     B: 'Singers'},
//            	data: {A: ["Arnold Schwarzenegger", "Jack Nicholson", "Robert de Niro", "Harrison Ford"], B: ["Freddy Mercury", "Ricky Martin", "Celine Dion"], AB: ["Marilyn Monroe", "Barbra Streisand", "Dean Martin"]},
//            	values: {A: 4, B: 3, AB: 3}
//            }],
            fnClickCallback: function() {
            	var value = "";
            	if (this.listnames.length == 1) {
            		value += "Elements only in ";
            	} else {
            		value += "Common elements in ";
            	}
            	for (name in this.listnames) {
            		value += this.listnames[name] + " ";
            	}
            	value += ":\n";
            	for (val in this.list) {
            		value += this.list[val] + "\n";
            	}
            	alert(value);
            },
            disableClick: false,
            useValues: false,
            exporting: true,
            displayMode: 'classic',
            displayStat: false,
            displaySwitch: true,
            shortNumber: true,
            searchInput: null,
            searchStatus: null,
            searchMinSize:1,
            // Colors must be RGB
            //       green         , blue            , red            , yellow          , orange         , brown
            colors: ["rgb(0,102,0)","rgb(90,155,212)","rgb(241,90,96)","rgb(250,220,91)","rgb(255,117,0)","rgb(192,152,83)"],
            fontSize: "12px",
            fontFamily: "Arial"
        };  
		var opts = $.extend(defaults, options),
			__canvasWidth = 500,
			__canvasHeight = 415,
			__context = null;
		
		function drawEllipse(x, y, r, w, h, a, fillcolor) {
			__context.beginPath();
			__context.save();
			__context.translate(x, y);
			__context.rotate(a*Math.PI/180);
			__context.scale(w, h);
			__context.arc(0, 0, r, 0, Math.PI * 2);
			__context.fillStyle = fillcolor;
			__context.fill();
			__context.restore();
		};
		
		function drawTriangle(x1, y1, x2, y2, x3, y3, fillcolor) {
			__context.beginPath();
			__context.save();
			__context.moveTo(x1,y1);
		    __context.lineTo(x2,y2);
		    __context.lineTo(x3,y3);
			__context.fillStyle = fillcolor;
			__context.fill();
			__context.restore();
		};
		
		function drawLine(x1, y1, x2, y2, strokecolor) {
		    __context.lineWidth = 1;
		    __context.beginPath();
		    __context.moveTo(x1, y1);
		    __context.lineTo(x2, y2);
			__context.strokeStyle = strokecolor;
		    __context.stroke();
		};
		
		function drawCircle(x, y, r, w, h, a, strokecolor) {
			__context.beginPath();
			__context.save();
			__context.translate(x, y);
			__context.rotate(a*Math.PI/180);
			__context.scale(w, h);
			__context.arc(0, 0, r, 0, Math.PI * 2);
		    __context.lineWidth = 3;
			__context.strokeStyle = strokecolor;
		    if(getOpacity(strokecolor) >= 0.6) {  // strokecolor must be rgba
			    __context.fillStyle = changeOpacity(strokecolor, 0.1);
				__context.fill();
		    }
		    __context.stroke();
			__context.restore();
		};
		
		function drawDoubleCircle(strokecolor) {
			var counterClockwise = false;
			__context.strokeStyle = strokecolor;
			__context.lineWidth = 3;
			// Circles
			__context.beginPath();
			__context.arc(120, 205, 90, 0.25 * Math.PI, 1.75 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(370, 205, 90, 1.25 * Math.PI, 0.75 * Math.PI, counterClockwise);
		   	__context.stroke();
			// Arcs
			__context.beginPath();
			__context.arc(245, 360, 110, 1.3 * Math.PI,  1.7 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(245, 50, 110, 0.3 * Math.PI,  0.7 * Math.PI, counterClockwise);
			__context.stroke();
			// Fill if hover
			if(getOpacity(strokecolor) >= 0.6) {  // strokecolor must be rgba
	    		__context.save();
	    		__context.beginPath();
	    		__context.moveTo(0,0);
	    		__context.lineTo(500,0);
	    		__context.lineTo(500,415);
	    		__context.lineTo(0,415);
	    		__context.closePath();
	    		__context.clip();
	    		__context.strokeStyle = 'rgba(0,0,0,0)';
	    		__context.lineCap = 'butt';
	    		__context.lineJoin = 'miter';
	    		__context.miterLimit = 4;
	    		__context.save();
	    		__context.restore();
	    		__context.save();
	    		__context.restore();
	    		__context.save();
	    		__context.translate(-11,-633);
	    		__context.save();
	    		__context.strokeStyle = 'rgba(0,0,0,0)';
	    		__context.translate(0,636);
	    		__context.beginPath();
	    		__context.moveTo(116.68996,288.55847);
	    		__context.bezierCurveTo(114.18984,288.06339,109.64415,287.21717,106.58843999999999,286.67797);
	    		__context.bezierCurveTo(91.905866,284.08711,74.873379,272.58587,62.202665,256.70632);
	    		__context.bezierCurveTo(50.11906,241.56257,44.228048,226.01748,43.257799,206.71498);
	    		__context.bezierCurveTo(41.560634,172.95091,57.502272,144.04252,87.68464900000001,126.1518);
	    		__context.bezierCurveTo(113.85468,110.63942,151.67555000000002,111.15554999999999,176.27348,127.36075);
	    		__context.bezierCurveTo(180.25579000000002,129.9843,185.18402,133.59167,187.22511,135.37712);
	    		__context.bezierCurveTo(194.61144,141.83832999999998,205.94366,148.05088999999998,218.69964,152.63209999999998);
	    		__context.bezierCurveTo(245.95179,162.41950999999997,276.56952,161.09978999999998,302.53531,149.01852999999997);
	    		__context.bezierCurveTo(312.66749999999996,144.30425999999997,319.21466,140.23029999999997,324.17127999999997,135.55560999999997);
	    		__context.bezierCurveTo(333.89739,126.38270999999997,351.79472999999996,118.14201999999997,367.20779999999996,115.73977999999997);
	    		__context.bezierCurveTo(378.37728999999996,113.99892999999997,395.02293999999995,115.10837999999997,406.09866999999997,118.33190999999997);
	    		__context.bezierCurveTo(427.39104999999995,124.52893999999996,447.81205,140.89770999999996,458.49307999999996,160.32946999999996);
	    		__context.bezierCurveTo(466.21541999999994,174.37854999999996,468.14135,182.79986999999997,468.18706999999995,202.71742999999995);
	    		__context.bezierCurveTo(468.21646999999996,215.53506999999996,467.84246999999993,219.76843999999994,466.13687999999996,225.92328999999995);
	    		__context.bezierCurveTo(460.34646999999995,246.81870999999995,443.97394999999995,267.47014999999993,423.77633999999995,279.35468999999995);
	    		__context.bezierCurveTo(411.54527999999993,286.55159999999995,399.98404999999997,289.09283999999997,379.83469999999994,289.01338999999996);
	    		__context.bezierCurveTo(354.90229999999997,288.91508999999996,338.08241999999996,282.66270999999995,322.01485999999994,267.52031999999997);
	    		__context.bezierCurveTo(312.0940299999999,258.17073,292.45028999999994,249.89035999999996,272.32759999999996,246.57576999999998);
	    		__context.bezierCurveTo(254.94397999999995,243.71236999999996,234.25916999999995,245.61216,218.58554999999996,251.51169);
	    		__context.bezierCurveTo(206.28085999999996,256.14317,191.57338999999996,264.45267,188.72360999999995,268.38324);
	    		__context.bezierCurveTo(186.89117999999996,270.91063,172.96095999999994,280.26608,166.96462999999994,282.99643);
	    		__context.bezierCurveTo(164.03704999999994,284.32946999999996,158.23248999999996,286.16623,154.06560999999994,287.07811999999996);
	    		__context.bezierCurveTo(145.28290999999993,289.00015999999994,123.31454999999994,289.87026999999995,116.68995999999993,288.55846999999994);
	    		__context.lineTo(116.68995999999993,288.55846999999994);
	    		__context.closePath();
			    __context.fillStyle = changeOpacity(strokecolor, 0.1);  // strokecolor must be rgba
	    		__context.fill();
	    		__context.stroke();
	    		__context.restore();
	    		__context.restore();
	    		__context.restore();
		    }
		};
		
		function drawCross(strokecolor) {
			var counterClockwise = false;
			__context.strokeStyle = strokecolor;
			__context.lineWidth = 3;
			__context.beginPath();
			__context.arc(145, 205, 40, 0.35 * Math.PI, 1.65 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(345, 205, 40, 1.35 * Math.PI, 0.65 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(245, 100, 40, 0.8 * Math.PI, 2.2 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(245, 310, 40, 1.82 * Math.PI, 1.18 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(182, 140, 35, 1.84 * Math.PI, 0.7 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(308, 140, 35, 0.3 * Math.PI, 1.18 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(182, 270, 35, 1.3 * Math.PI, 0.2 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(308, 270, 35, 0.8 * Math.PI, 1.7 * Math.PI, counterClockwise);
			__context.stroke();
			// Fill if hover
			if(getOpacity(strokecolor) >= 0.6) {  // strokecolor must be rgba
				__context.save();
				__context.beginPath();
				__context.moveTo(0,0);
				__context.lineTo(500,0);
				__context.lineTo(500,415);
				__context.lineTo(0,415);
				__context.closePath();
				__context.clip();
				__context.strokeStyle = 'rgba(0,0,0,0)';
				__context.lineCap = 'butt';
				__context.lineJoin = 'miter';
				__context.miterLimit = 4;
				__context.save();
				__context.restore();
				__context.save();
				__context.restore();
				__context.save();
				__context.translate(-140,-212);
				__context.save();
				__context.strokeStyle = 'rgba(0,0,0,0)';
				__context.translate(139,217);
				__context.beginPath();
				__context.moveTo(234.95394,341.16933);
				__context.bezierCurveTo(231.80731999999998,340.135,227.61570999999998,338.28931,225.63925999999998,337.0678);
				__context.bezierCurveTo(221.28248999999997,334.37516999999997,214.37806999999998,326.70736999999997,211.83772,321.74032);
				__context.bezierCurveTo(207.9414,314.12196,208.45042999999998,293.37118,212.67719,287.51943);
				__context.bezierCurveTo(215.87649,283.09018,218.95801,274.65101,219.59244,268.58113);
				__context.bezierCurveTo(220.45250000000001,260.35249,218.0029,251.63114999999996,212.48275,243.26837999999998);
				__context.bezierCurveTo(206.02273000000002,233.48175999999998,198.89664000000002,229.43794999999997,185.95006,228.21201);
				__context.bezierCurveTo(178.17816000000002,227.47606,172.06373000000002,228.9391,162.89332000000002,233.72894);
				__context.bezierCurveTo(156.28128,237.18251999999998,155.23208000000002,237.42228,146.73088,237.42244);
				__context.bezierCurveTo(135.51808000000003,237.42266,128.41302000000002,234.90269,121.10778000000002,228.33465999999999);
				__context.bezierCurveTo(114.30638000000002,222.21962,111.33466000000001,216.76216,109.78745000000002,207.54525999999998);
				__context.bezierCurveTo(108.21951000000003,198.20489999999998,109.69855000000003,187.35654,113.42449000000002,180.86914);
				__context.bezierCurveTo(116.70310000000002,175.16054,123.50535000000002,169.17129,130.61112000000003,165.73660999999998);
				__context.bezierCurveTo(135.84135000000003,163.2085,137.81975000000003,162.82216,145.72073000000003,162.78605);
				__context.bezierCurveTo(154.43018000000004,162.74624999999997,155.21524000000002,162.9348,164.40855000000002,167.27446999999998);
				__context.bezierCurveTo(173.55745000000002,171.59319999999997,174.43144,171.80518999999998,183.14876,171.82016);
				__context.bezierCurveTo(191.71748000000002,171.83496,192.78155,171.59418,200.06856000000002,167.99274);
				__context.bezierCurveTo(206.71571000000003,164.70754,208.40559000000002,163.31437,211.70876,158.39629);
				__context.bezierCurveTo(216.86766,150.71523,219.00004,144.97677,219.64970000000002,137.02646);
				__context.bezierCurveTo(220.25703000000001,129.59413999999998,217.81806000000003,120.60798999999999,213.36866000000003,113.88449999999999);
				__context.bezierCurveTo(208.31047000000004,106.24103999999998,207.60430000000002,86.789031,212.06743000000003,78.04059099999999);
				__context.bezierCurveTo(215.36468000000002,71.57742099999999,222.48250000000002,64.79066099999999,229.83172000000002,61.10251999999999);
				__context.bezierCurveTo(236.05208000000002,57.98088699999999,237.41276000000002,57.68862899999999,245.72583000000003,57.68862899999999);
				__context.bezierCurveTo(253.41441000000003,57.68862899999999,255.67446000000004,58.09343499999999,260.37304,60.31215099999999);
				__context.bezierCurveTo(277.08479,68.203604,286.51407,86.23497699999999,282.57463,102.76767);
				__context.bezierCurveTo(281.86319000000003,105.75339,279.29147,112.06009,276.85971,116.78255);
				__context.bezierCurveTo(273.65874,122.99881,272.27383000000003,127.0461,271.84239,131.44528);
				__context.bezierCurveTo(270.33639000000005,146.80139,279.61874,162.9296,293.20300000000003,168.5596);
				__context.bezierCurveTo(306.55802000000006,174.09458999999998,314.0903,173.65240999999997,329.06341000000003,166.45448);
				__context.bezierCurveTo(336.00568000000004,163.11717,337.40021,162.8096,345.73093000000006,162.77846);
				__context.bezierCurveTo(353.96010000000007,162.74776,355.45355000000006,163.05339,361.47832000000005,166.00153);
				__context.bezierCurveTo(369.53210000000007,169.94254,376.64791,177.20452,380.3403400000001,185.25103000000001);
				__context.bezierCurveTo(382.5644300000001,190.09773,382.9971500000001,192.49377,383.0256800000001,200.12014000000002);
				__context.bezierCurveTo(383.0562800000001,208.31117000000003,382.74275000000006,209.86176000000003,379.8593200000001,215.77750000000003);
				__context.bezierCurveTo(376.0651900000001,223.56169000000003,370.1637800000001,229.56756000000004,362.0334900000001,233.91890000000004);
				__context.bezierCurveTo(351.4550500000001,239.58048000000002,339.7398500000001,239.17006000000003,326.6582200000001,232.67959000000005);
				__context.bezierCurveTo(310.6073100000001,224.71590000000003,294.4869300000001,227.24494000000004,281.98739000000006,239.68774000000005);
				__context.bezierCurveTo(270.03863000000007,251.58229000000006,268.4080000000001,270.18789000000004,277.92910000000006,285.99336000000005);
				__context.bezierCurveTo(284.3299400000001,296.61906000000005,285.13914000000005,307.05507000000006,280.51249000000007,319.31069);
				__context.bezierCurveTo(277.5959900000001,327.03626,268.9701400000001,335.93369,260.6715500000001,339.77631);
				__context.bezierCurveTo(253.04107000000008,343.30956000000003,243.1026400000001,343.84789,234.95394000000007,341.16933);
				__context.lineTo(234.95394000000007,341.16933);
				__context.closePath();
			    __context.fillStyle = changeOpacity(strokecolor, 0.1);  // strokecolor must be rgba
				__context.fill();
				__context.stroke();
				__context.restore();
				__context.restore();
				__context.restore();
		    }
		};
		
		function drawCross2(strokecolor) {
			var counterClockwise = false;
			__context.strokeStyle = strokecolor;
			__context.lineWidth = 3;
			__context.beginPath();
			__context.arc(151, 205, 18, 0.4 * Math.PI, 1.6 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(167, 238, 18, 1.35 * Math.PI, 2.42 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(180, 272, 18, 0.18 * Math.PI, 1.33 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(212, 289, 18, 1.15 * Math.PI, 2.1 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(245, 302, 18, 1.85 * Math.PI, 1.15 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(279, 289, 18, 0.87 * Math.PI, 1.85 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(311, 272, 18, 1.65 * Math.PI, 0.85 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(325, 238, 18, 0.54 * Math.PI, 1.68 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(340, 205, 18, 1.37 * Math.PI, 0.59 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(325, 172, 18, 0.35 * Math.PI, 1.4 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(311, 138, 18, 1.15 * Math.PI, 2.35 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(279, 122, 18, 0.17 * Math.PI, 1.1 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(245, 108, 18, 0.8 * Math.PI, 2.15 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(212, 122, 18, 1.9 * Math.PI, 0.85 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(180, 138, 18, 0.65 * Math.PI, 1.85 * Math.PI, counterClockwise);
			__context.stroke();
			__context.beginPath();
			__context.arc(167, 172, 18, 1.6 * Math.PI, 0.7 * Math.PI, counterClockwise);
			__context.stroke();
			// Fill if hover
		    if(getOpacity(strokecolor) >= 0.6) {  // strokecolor must be rgba
				__context.save();
				__context.beginPath();
				__context.moveTo(0,0);
				__context.lineTo(500,0);
				__context.lineTo(500,415);
				__context.lineTo(0,415);
				__context.closePath();
				__context.clip();
				__context.strokeStyle = 'rgba(0,0,0,0)';
				__context.lineCap = 'butt';
				__context.lineJoin = 'miter';
				__context.miterLimit = 4;
				__context.save();
				__context.restore();
				__context.save();
				__context.restore();
				__context.save();
				__context.translate(-4,-637);
				__context.save();
				__context.strokeStyle = 'rgba(0,0,0,0)';
				__context.translate(0,637);
				__context.beginPath();
				__context.moveTo(242.80929,316.76895);
				__context.bezierCurveTo(239.2558,315.3947,234.83257,308.84895,234.18254000000002,304.00262000000004);
				__context.bezierCurveTo(233.85998,301.59777,234.27994,298.62012000000004,235.27891000000002,296.22928);
				__context.bezierCurveTo(237.39948,291.15406,237.29592000000002,285.99823000000004,234.97058,280.87787000000003);
				__context.bezierCurveTo(231.58546,273.42388000000005,224.01393000000002,268.42150000000004,216.05280000000002,268.37923);
				__context.bezierCurveTo(210.33377000000002,268.34883,202.51030000000003,272.78013,198.27794000000003,278.44704);
				__context.bezierCurveTo(196.41545000000002,280.94083,193.36619000000002,284.00917000000004,191.50181000000003,285.26557);
				__context.bezierCurveTo(180.81154000000004,292.46974,165.85368000000003,281.37461,168.71196000000003,268.36097);
				__context.bezierCurveTo(169.66342000000003,264.02899,174.42957000000004,258.84685,179.89636000000004,256.20042);
				__context.bezierCurveTo(188.48844000000005,252.04107000000002,193.00795000000005,242.34352,190.58184000000006,233.27254000000002);
				__context.bezierCurveTo(189.87870000000007,230.64358000000001,188.21908000000005,227.07100000000003,186.89381000000006,225.33348);
				__context.bezierCurveTo(182.51221000000007,219.58889000000002,172.48125000000005,216.46373,166.00572000000005,218.82574);
				__context.bezierCurveTo(160.88401000000005,220.69395,153.47334000000006,221.17446999999999,149.58323000000004,219.89061);
				__context.bezierCurveTo(147.65941000000004,219.25570000000002,144.69894000000005,217.15717,143.00442000000004,215.22721);
				__context.bezierCurveTo(140.25313000000003,212.09367,139.92347000000004,211.02317000000002,139.92347000000004,205.22251);
				__context.bezierCurveTo(139.92347000000004,199.42185,140.25313000000003,198.35135,143.00442000000004,195.21781);
				__context.bezierCurveTo(148.08184000000003,189.43495,155.59116000000003,188.31217999999998,166.73325000000003,191.66995999999997);
				__context.bezierCurveTo(171.89873000000003,193.22662999999997,176.76581000000002,192.32402999999996,182.19737000000003,188.80213999999998);
				__context.bezierCurveTo(186.99980000000002,185.68819,191.44125000000003,178.0236,191.44125000000003,172.85002999999998);
				__context.bezierCurveTo(191.44125000000003,165.91124999999997,187.02605000000003,157.91411999999997,181.72540000000004,155.25199999999998);
				__context.bezierCurveTo(168.68097000000003,148.70074999999997,164.62205000000003,137.69753999999998,171.78539000000004,128.30592);
				__context.bezierCurveTo(178.71368000000004,119.22245999999998,190.67567000000003,120.01398999999998,197.23888000000002,129.99017999999998);
				__context.bezierCurveTo(204.51058000000003,141.04331,214.32951000000003,144.72351999999998,224.37555000000003,140.16123);
				__context.bezierCurveTo(235.23365000000004,135.23013999999998,238.86734000000004,125.44275999999999,234.78543000000002,112.12207);
				__context.bezierCurveTo(233.16446000000002,106.83225999999999,234.54797000000002,101.56208,238.76789000000002,96.95178);
				__context.bezierCurveTo(248.33151000000004,86.503446,267.39820000000003,96.55735,264.708,110.63005);
				__context.bezierCurveTo(262.24834000000004,123.49681,263.43732,130.24138,269.17704000000003,135.9811);
				__context.bezierCurveTo(273.84920000000005,140.65326,279.80664,142.46487,287.54524000000004,141.56671);
				__context.bezierCurveTo(290.80029,141.18892,293.25246000000004,140.53831,292.99449000000004,140.12091);
				__context.bezierCurveTo(292.73652000000004,139.70352,293.68403000000006,138.74196,295.10007,137.98412000000002);
				__context.bezierCurveTo(296.51611,137.22628000000003,298.81193,134.66566000000003,300.2019,132.29385000000002);
				__context.bezierCurveTo(303.47878000000003,126.70229000000002,307.77937000000003,123.14413000000002,312.28166000000004,122.29950000000002);
				__context.bezierCurveTo(323.19238000000007,120.25263000000002,332.5219000000001,130.29284,330.36964000000006,141.76535);
				__context.bezierCurveTo(329.6006300000001,145.86456,323.8285400000001,152.29884,317.5065000000001,156.10419000000002);
				__context.bezierCurveTo(312.2969800000001,159.23989,308.17124000000007,167.72107000000003,308.8642700000001,173.86973);
				__context.bezierCurveTo(309.6787300000001,181.09556,313.7704700000001,186.95133,320.2047800000001,190.09934);
				__context.bezierCurveTo(326.1904100000001,193.02784000000003,330.20042000000007,193.13284000000002,338.9627100000001,190.59052000000003);
				__context.bezierCurveTo(344.9917000000001,188.84125000000003,350.1988400000001,189.89377000000002,354.59001000000006,193.74927000000002);
				__context.bezierCurveTo(368.2933100000001,205.78092000000004,353.68808000000007,226.52374000000003,336.55395000000004,219.36465);
				__context.bezierCurveTo(323.93852000000004,214.09358,308.93266000000006,224.16784,308.71798000000007,238.05247);
				__context.bezierCurveTo(308.6646800000001,241.50290999999999,311.2253800000001,248.07176,313.86683000000005,251.26035);
				__context.bezierCurveTo(314.59179000000006,252.13548,316.54865000000007,253.7222,318.21540000000005,254.7864);
				__context.bezierCurveTo(324.03371000000004,258.50135,328.85284000000007,263.67991,329.97401,267.42204);
				__context.bezierCurveTo(332.95992,277.3881,325.11943,287.55831,314.45977000000005,287.54616999999996);
				__context.bezierCurveTo(309.88386,287.54116999999997,305.3854400000001,284.53792999999996,300.84666000000004,278.45842999999996);
				__context.bezierCurveTo(295.84636000000006,271.76072,288.36647000000005,267.87215,281.88796,268.60236999999995);
				__context.bezierCurveTo(274.64614,269.41862999999995,268.80722000000003,273.50679999999994,265.64137,279.97756999999996);
				__context.bezierCurveTo(263.21128,284.94451999999995,262.99769,286.21504999999996,263.64589,291.84819999999996);
				__context.bezierCurveTo(265.39685000000003,307.06499999999994,264.64547,310.09781999999996,257.87161000000003,315.15523999999994);
				__context.bezierCurveTo(254.81776000000002,317.43526999999995,246.76347000000004,318.2981699999999,242.80929000000003,316.76894999999996);
				__context.lineTo(242.80929000000003,316.76894999999996);
				__context.closePath();
				__context.fillStyle = changeOpacity(strokecolor, 0.1); // strokecolor must be rgba
				__context.fill();
				__context.stroke();
				__context.restore();
				__context.restore();
				__context.restore();
		    }
		};
		
		function drawRoundRect(x, y, width, height, strokecolor) {
		    var radius = 10;
		    __context.beginPath();
		    __context.moveTo(x + radius, y);
		    __context.lineTo(x + width - radius, y);
		    __context.quadraticCurveTo(x + width, y, x + width, y + radius);
		    __context.lineTo(x + width, y + height - radius);
		    __context.quadraticCurveTo(x + width, y + height, x + width - radius, y + height);
		    __context.lineTo(x + radius, y + height);
		    __context.quadraticCurveTo(x, y + height, x, y + height - radius);
		    __context.lineTo(x, y + radius);
		    __context.quadraticCurveTo(x, y, x + radius, y);
		    __context.closePath();
		    __context.lineWidth = 3;
		    __context.strokeStyle = strokecolor;
		    if(getOpacity(strokecolor) >= 0.6) {  // strokecolor must be rgba
			    __context.fillStyle = changeOpacity(strokecolor, 0.1); // strokecolor must be rgba
				__context.fill();
		    }
		    __context.stroke();
		};
		
		function drawAxis(ctx, startx, starty, endx, endy) {
			ctx.beginPath();
			ctx.moveTo(startx, starty);
			ctx.lineTo(endx, endy);
			ctx.closePath();
			ctx.stroke();
		}
		
		function drawRectangle(ctx, x, y, w, h, fillcolor, strokecolor) {
			ctx.beginPath();
			ctx.rect(x, y, w, h);
			ctx.closePath();
			//ctx.save();
			ctx.lineWidth = 0.75;
			ctx.shadowColor = "rgba(0,0,0, 0.4)";
			ctx.shadowBlur = 7;
			ctx.shadowOffsetX = 2;
			ctx.shadowOffsetY = -2;
			ctx.strokeStyle = strokecolor;
			ctx.stroke();
			ctx.fillStyle = fillcolor;
			ctx.fill();
			//ctx.restore();
		}
		
		function changeOpacity(rgba, opacity) {
			var	colorStr = rgba.slice(rgba.indexOf('(') + 1, rgba.indexOf(')'));
			var	colorArr = colorStr.split(',');
			return "rgba("+colorArr[0] + "," + colorArr[1] + "," + colorArr[2] + "," + opacity + ")";
		}
		
		function getOpacity(rgba) {
			var	colorStr = rgba.slice(rgba.indexOf('(') + 1, rgba.indexOf(')'));
			var	colorArr = colorStr.split(',');
			return colorArr[3];
		}
		
		function clearCanvas() {
			__context.clearRect(0, 0, __canvasWidth, __canvasHeight);
		}

		function placeStat(vennSize) {
			var axiscolor = "rgba(0,0,0, 0.7)";
			/* 
			 * Bar chart
			*/
			var	h = 120,
				xmargin = 70,
				ymargin = 415,
				xspacer = 20,
				barwidth = (370-(vennSize*xspacer))/vennSize,
				ytext = 265;
			
			// Data
			var	data = sizeOfClass(),
				dataplot = new Array();
			var max = 0;
			for (var i=0; i<vennSize; i++) {
				max = Math.max(max, data[i]);
			}
			for (var i=0; i<vennSize; i++) {
				dataplot.push(data[i]/max * (h-50));
			}
		
			// Draw the bar chart
			__context.fillStyle = "#000";
			__context.font = 'italic 10pt Arial';
			__context.textAlign = 'center';
		    __context.fillText('Size of each list', 250, ymargin+27);
		    __context.font = 'normal 8pt Arial';
		    __context.textAlign = 'right';
		    __context.fillText(0, 45, ymargin + h + 2);
		    __context.fillText(max/2, 45, ymargin + ((h+55)/2));
		    __context.fillText(max, 45, ymargin + 54);
		    __context.textAlign = 'left';
			for (var i=0; i<vennSize; i++) {
				drawRectangle(__context,
							  xmargin + i*barwidth + i*xspacer,
							  ymargin + h-dataplot[i],
							  barwidth,
							  dataplot[i],
							  changeOpacity(opts.colors[i], 0.5),
							  opts.colors[i]);
				if(h-dataplot[i]+15 <= h-10) {
					__context.textAlign = 'right';
					__context.fillStyle = 'white';
					__context.fillText(data[i], barwidth + 65 + i*barwidth + i*xspacer, ymargin + h-dataplot[i]+15);
				}
				else {
					__context.textAlign = 'right';
					__context.fillText(data[i], barwidth + 65 + i*barwidth + i*xspacer, ymargin + h-dataplot[i]-10);
				}
				if(i%2 && vennSize>2) {	ytext = ymargin + h + 25; }
				else { ytext = ymargin + h + 15; }
				__context.fillStyle = "#000";
				__context.textAlign = 'center';
				__context.fillText($("#"+cpref+"-label"+(i+1)).html(), (xmargin + i*barwidth + i*xspacer) + barwidth/2, ytext, 200);
				__context.strokeStyle = axiscolor;
				__context.lineWidth = 0.4;
				drawAxis(__context, (xmargin + i*barwidth + i*xspacer) + barwidth/2, ymargin + h, (xmargin + i*barwidth + i*xspacer) + barwidth/2, ymargin + h+5);
			}
			
			// Draw the x and y axes
			__context.lineWidth = 1;
			__context.strokeStyle = axiscolor;
			drawAxis(__context, 50, ymargin + h, 50, ymargin + 35); 
			drawAxis(__context, 50, ymargin + h, 450, ymargin + h);
			__context.lineWidth = 0.4;
			drawAxis(__context, 47, ymargin + (h+48)/2, 53, ymargin + (h+48)/2);
			drawAxis(__context, 47, ymargin + 50, 53, ymargin + 50);
			drawTriangle(50, ymargin+25,  46,ymargin+35,   54,ymargin+35, axiscolor);
			drawTriangle(460,ymargin+h, 450,ymargin+h-4, 450,ymargin+h+4, axiscolor);
			
			/* 
			 * Stacked bar chart
			*/
			// Data
			var	data2 = countByNbClass(),
				data2plot = new Array();
			xspacer = 2,
			xmargin = 60;
			ymargin += h + 50;
			var	maxwidth = 390 + xspacer,
				sum = 0;
			for (var i=0; i<vennSize; i++) {
				sum += data2[i];
			}
			for (var i=0; i<vennSize; i++) {
				data2plot.push(data2[i]/sum * maxwidth);
			}

			// Draw the bar chart
			__context.font = 'italic 10pt Arial';
			__context.textAlign = 'center';
		    __context.fillText('Number of elements: specific (1) or shared by 2, 3, ... lists', 250, ymargin);
		    __context.font = 'normal 8pt Arial';
			var	xprev = 0,
				ylegend = 0;
			
			for (var i=vennSize-1; i>=0; i--) {
				if(data2plot[i] == 0) { continue };
			    drawRectangle(__context,
			    			  xmargin + xprev,
							  ymargin + 15,
							  data2plot[i] - xspacer,
							  20,
							  changeOpacity("rgb(156,106,156)", (1/(i+1.5))),
							  "rgba(0,0,0,0.5)");
				__context.textAlign = 'center';
				__context.fillStyle = 'white';
				if((data2plot[i] - xspacer) < 25) {
					__context.fillStyle = '#000';
					__context.fillText(i+1, (data2plot[i] - xspacer)/2  + xmargin + xprev, ymargin + 57 + ylegend);
					__context.textAlign = 'left';
					__context.fillText("(" + data2[i] + ")", (data2plot[i] - xspacer)/2  + xmargin + xprev + 6, ymargin + 57 + ylegend);
					if(ylegend != 48) { ylegend += 12; }
					else { ylegend = 0; }
				}
				else {
					if(ylegend >= 24) { ylegend = 0; }
					__context.fillText(data2[i], (data2plot[i] - xspacer)/2  + xmargin + xprev, ymargin + 29);
					__context.fillStyle = '#000';
					__context.fillText(i+1, (data2plot[i] - xspacer)/2  + xmargin + xprev, ymargin + 57 + ylegend);
					ylegend = 0;
				}
				__context.strokeStyle = axiscolor;
				__context.lineWidth = 0.4;
				drawAxis(__context, (data2plot[i] - xspacer)/2  + xmargin + xprev, ymargin+40, (data2plot[i] - xspacer)/2  + xmargin + xprev, ymargin+45);
				xprev += data2plot[i];
			}

		    // Draw the x axis
		    __context.lineWidth = 1;
			__context.strokeStyle = axiscolor;
			drawAxis(__context, 50, ymargin+40, 460, ymargin+40);
		}
		
		function placeNumber(div, left, top, space) {
			// if shortNumber option && lenght>space
            if(opts.shortNumber  &&  div.text().length > space) {
				div.html("<span title=" + div.text() + ">?</span>");
			}
            var l = left - ((div.text().length-1)/2 * 4);
			div.css("left", l).css("top", top);
		}
		
		function transpose (a) {
			// Calculate the width and height of the Array
			var	w = a.length ? a.length : 0,
				h = a[0] instanceof Array ? a[0].length : 0;
			for(i=0; i<w; i++) {
				if(a[i] instanceof Array  &&  a[i].length > h) {
					h = a[i].length;
				}
			}
			// In case it is a zero matrix, no transpose routine needed.
			if(h === 0 || w === 0) { return []; }
			var i, j, t = [];
			// Loop through every item in the outer array (height)
			for(i=0; i<h; i++) {
				// Insert a new row (array)
				t[i] = [];
				// Loop through every item per item in outer array (width)
				for(j=0; j<w; j++) {
					// Save transposed data.
					t[i][j] = a[j][i];
				}
			}
			return t;
		}
		
		function placeClassicVenn(vennSize) {
			var	grey = "rgba(0,0,0,0.1)";
						
			if (vennSize == 6) {
				drawTriangle(0,11,    254,160, 174,235, changeOpacity(opts.colors[0], $("#"+cpref+"-label1").css('opacity')));
				drawTriangle(188,0,   134,242, 236,202, changeOpacity(opts.colors[1], $("#"+cpref+"-label2").css('opacity')));
				drawTriangle(338,52,  135,123, 191,242, changeOpacity(opts.colors[2], $("#"+cpref+"-label3").css('opacity')));
				drawTriangle(500,260, 163,117, 134,219, changeOpacity(opts.colors[3], $("#"+cpref+"-label4").css('opacity')));
				drawTriangle(250,415, 133,150, 203,67,  changeOpacity(opts.colors[4], $("#"+cpref+"-label5").css('opacity')));
				drawTriangle(11,307,  263,81,  214,220, changeOpacity(opts.colors[5], $("#"+cpref+"-label6").css('opacity')));
				
				$("#"+cpref+"-label1").css("left",  35).css("top",  10).css("color", opts.colors[0]);
				$("#"+cpref+"-label2").css("left", 200).css("top",   5).css("color", opts.colors[1]);
				$("#"+cpref+"-label3").css("left", 335).css("top",  60).css("color", opts.colors[2]);
				$("#"+cpref+"-label4").css("left", 410).css("top", 200).css("color", opts.colors[3]);
				$("#"+cpref+"-label5").css("left", 255).css("top", 385).css("color", opts.colors[4]);
				$("#"+cpref+"-label6").css("left",  30).css("top", 300).css("color", opts.colors[5]);
				placeNumber($("#"+cpref+"-resultC100000"),  98,  90, 6);
				placeNumber($("#"+cpref+"-resultC010000"), 187,  50, 4);
				placeNumber($("#"+cpref+"-resultC001000"), 280,  80, 6);
				placeNumber($("#"+cpref+"-resultC000100"), 320, 205, 6);
				placeNumber($("#"+cpref+"-resultC000010"), 212, 272, 5);
				placeNumber($("#"+cpref+"-resultC000001"), 100, 240, 6);
				placeNumber($("#"+cpref+"-resultC111111"), 185, 170, 6);
				drawLine(140, 80, 166,110, grey);
				placeNumber($("#"+cpref+"-resultC110000"), 132,  60, 6);
				placeNumber($("#"+cpref+"-resultC101000"), 142, 116, 2);
				drawLine( 75,180, 145,185, grey);
				drawLine( 75,180,  65,175, grey);
				placeNumber($("#"+cpref+"-resultC100100"),  55, 157, 6);
				placeNumber($("#"+cpref+"-resultC100010"), 140, 145, 2);
				drawLine( 75,200, 142,190, grey);
				drawLine( 75,200,  65,195, grey);
				placeNumber($("#"+cpref+"-resultC100001"),  55, 177, 6);
				drawLine(230, 80, 212,115, grey);
				placeNumber($("#"+cpref+"-resultC011000"), 230,  60, 6);
				placeNumber($("#"+cpref+"-resultC010100"), 225, 190, 1);
				placeNumber($("#"+cpref+"-resultC010010"), 193,  82, 3);
				placeNumber($("#"+cpref+"-resultC010001"), 145, 218, 2);
				drawLine(295,145, 235,180, grey);
				placeNumber($("#"+cpref+"-resultC001100"), 302, 132, 6);
				drawLine(275,270, 193,233, grey);
				placeNumber($("#"+cpref+"-resultC001010"), 275, 268, 6);
				placeNumber($("#"+cpref+"-resultC001001"), 232, 113, 4);
				placeNumber($("#"+cpref+"-resultC000110"), 216, 211, 1);
				drawLine( 75,220, 140,205, grey);
				drawLine( 75,220,  65,215, grey);
				placeNumber($("#"+cpref+"-resultC000101"),  55, 197, 6);
				drawLine(150,270, 183,230, grey);
				placeNumber($("#"+cpref+"-resultC000011"), 145, 268, 6);
				
				placeNumber($("#"+cpref+"-resultC111110"), 170, 130, 6);
				$("#"+cpref+"-resultC111110").hide();
				placeNumber($("#"+cpref+"-resultC111101"), 218, 155, 6);
				$("#"+cpref+"-resultC111101").hide();
				placeNumber($("#"+cpref+"-resultC111011"), 205, 125, 6);
				$("#"+cpref+"-resultC111011").hide();
				placeNumber($("#"+cpref+"-resultC110111"), 160, 190, 6);
				$("#"+cpref+"-resultC110111").hide();
				placeNumber($("#"+cpref+"-resultC101111"), 180, 212, 6);
				$("#"+cpref+"-resultC101111").hide();
				placeNumber($("#"+cpref+"-resultC011111"), 208, 197, 6);
				$("#"+cpref+"-resultC011111").hide();
				
				placeNumber($("#"+cpref+"-resultC111100"), 225, 170, 6);
				$("#"+cpref+"-resultC111100").hide();
				placeNumber($("#"+cpref+"-resultC111010"), 182, 113, 6);
				$("#"+cpref+"-resultC111010").hide();
				placeNumber($("#"+cpref+"-resultC110110"), 149, 163, 6);
				$("#"+cpref+"-resultC110110").hide();
				placeNumber($("#"+cpref+"-resultC101110"), 152, 140, 6);
				$("#"+cpref+"-resultC101110").hide();
				placeNumber($("#"+cpref+"-resultC011110"), 215, 195, 6);
				$("#"+cpref+"-resultC011110").hide();
				placeNumber($("#"+cpref+"-resultC111001"), 212, 128, 6);
				$("#"+cpref+"-resultC111001").hide();
				placeNumber($("#"+cpref+"-resultC110101"), 150, 194, 6);
				$("#"+cpref+"-resultC110101").hide();
				placeNumber($("#"+cpref+"-resultC110011"), 168, 216, 6);
				$("#"+cpref+"-resultC110011").hide();
				placeNumber($("#"+cpref+"-resultC101101"), 226, 145, 6);
				$("#"+cpref+"-resultC101101").hide();
				placeNumber($("#"+cpref+"-resultC101011"), 178, 216, 6);
				$("#"+cpref+"-resultC101011").hide();
				placeNumber($("#"+cpref+"-resultC100111"), 178, 214, 6);
				$("#"+cpref+"-resultC100111").hide();
				placeNumber($("#"+cpref+"-resultC011101"), 218, 180, 6);
				$("#"+cpref+"-resultC011101").hide();
				placeNumber($("#"+cpref+"-resultC011011"), 206, 122, 6);
				$("#"+cpref+"-resultC011011").hide();
				placeNumber($("#"+cpref+"-resultC010111"), 215, 198, 6);
				$("#"+cpref+"-resultC010111").hide();
				placeNumber($("#"+cpref+"-resultC001111"), 200, 208, 6);
				$("#"+cpref+"-resultC001111").hide();
				
				placeNumber($("#"+cpref+"-resultC111000"), 156, 105, 6);
				$("#"+cpref+"-resultC111000").hide();
				placeNumber($("#"+cpref+"-resultC110100"), 145, 174, 6);
				$("#"+cpref+"-resultC110100").hide();
				placeNumber($("#"+cpref+"-resultC110010"), 163, 100, 6);
				$("#"+cpref+"-resultC110010").hide();
				placeNumber($("#"+cpref+"-resultC101100"), 238, 152, 6);
				$("#"+cpref+"-resultC101100").hide();
				placeNumber($("#"+cpref+"-resultC101010"), 148, 128, 6);
				$("#"+cpref+"-resultC101010").hide();
				placeNumber($("#"+cpref+"-resultC100110"), 145, 158, 6);
				$("#"+cpref+"-resultC100110").hide();
				placeNumber($("#"+cpref+"-resultC011100"), 226, 178, 6);
				$("#"+cpref+"-resultC011100").hide();
				placeNumber($("#"+cpref+"-resultC011010"), 195, 105, 6);
				$("#"+cpref+"-resultC011010").hide();
				placeNumber($("#"+cpref+"-resultC010110"), 218, 198, 6);
				$("#"+cpref+"-resultC010110").hide();
				placeNumber($("#"+cpref+"-resultC001110"), 202, 216, 6);
				$("#"+cpref+"-resultC001110").hide();
				placeNumber($("#"+cpref+"-resultC110001"), 160, 213, 6);
				$("#"+cpref+"-resultC110001").hide();
				placeNumber($("#"+cpref+"-resultC101001"), 220, 132, 6);
				$("#"+cpref+"-resultC101001").hide();
				placeNumber($("#"+cpref+"-resultC100101"), 139, 182, 6);
				$("#"+cpref+"-resultC100101").hide();
				placeNumber($("#"+cpref+"-resultC100011"), 170, 220, 6);
				$("#"+cpref+"-resultC100011").hide();
				placeNumber($("#"+cpref+"-resultC011001"), 212, 122, 6);
				$("#"+cpref+"-resultC011001").hide();
				placeNumber($("#"+cpref+"-resultC010101"), 146, 203, 6);
				$("#"+cpref+"-resultC010101").hide();
				placeNumber($("#"+cpref+"-resultC010011"), 164, 220, 6);
				$("#"+cpref+"-resultC010011").hide();
				placeNumber($("#"+cpref+"-resultC001101"), 233, 138, 6);
				$("#"+cpref+"-resultC001101").hide();
				placeNumber($("#"+cpref+"-resultC001011"), 185, 218, 6);
				$("#"+cpref+"-resultC001011").hide();
				placeNumber($("#"+cpref+"-resultC000111"), 209, 207, 6);
				$("#"+cpref+"-resultC000111").hide();
				
			} else if (vennSize == 5) {
				drawEllipse(214,230,10,18.6,9.5,25,  changeOpacity(opts.colors[0], $("#"+cpref+"-label1").css('opacity')));
				drawEllipse(232,187,10,18.6,9.5,98,  changeOpacity(opts.colors[1], $("#"+cpref+"-label2").css('opacity')));
				drawEllipse(273,196,10,18.6,9.5,170, changeOpacity(opts.colors[2], $("#"+cpref+"-label3").css('opacity')));
				drawEllipse(282,238,10,18.6,9.5,62,  changeOpacity(opts.colors[3], $("#"+cpref+"-label4").css('opacity')));
				drawEllipse(242,260,10,18.6,9.5,134, changeOpacity(opts.colors[4], $("#"+cpref+"-label5").css('opacity')));
				
				$("#"+cpref+"-label1").css("left",   0).css("top", 100).css("color",   opts.colors[0]);
				$("#"+cpref+"-label2").css("left", 310).css("top",  15).css("color",  opts.colors[1]);
				$("#"+cpref+"-label3").css("left", 450).css("top", 120).css("color", opts.colors[2]);
				$("#"+cpref+"-label4").css("left", 387).css("top", 400).css("color", opts.colors[3]);
				$("#"+cpref+"-label5").css("left",  40).css("top", 400).css("color",  opts.colors[4]);
				$("#"+cpref+"-label6").css("left", -1000).css("top", -2200);
				placeNumber($("#"+cpref+"-resultC100000"),  75, 150, 6);
				placeNumber($("#"+cpref+"-resultC010000"), 245,  30, 6);
				placeNumber($("#"+cpref+"-resultC001000"), 415, 162, 6);
				placeNumber($("#"+cpref+"-resultC000100"), 350, 370, 6);
				placeNumber($("#"+cpref+"-resultC000010"), 132, 370, 6);
				placeNumber($("#"+cpref+"-resultC110000"), 152, 118, 2);
				placeNumber($("#"+cpref+"-resultC101000"), 110, 200, 6);
				placeNumber($("#"+cpref+"-resultC100100"), 350, 295, 6);
				placeNumber($("#"+cpref+"-resultC100010"), 125, 273, 2);
				placeNumber($("#"+cpref+"-resultC011000"), 309,  94, 3);
				placeNumber($("#"+cpref+"-resultC010100"), 215,  82, 6);
				placeNumber($("#"+cpref+"-resultC010010"), 195, 340, 6);
				placeNumber($("#"+cpref+"-resultC001100"), 378, 232, 2);
				placeNumber($("#"+cpref+"-resultC001010"), 360, 140, 5);
				placeNumber($("#"+cpref+"-resultC000110"), 261, 347, 2);
				placeNumber($("#"+cpref+"-resultC111000"), 148, 180, 4);
				placeNumber($("#"+cpref+"-resultC110100"), 168, 113, 2);
				placeNumber($("#"+cpref+"-resultC110010"), 176, 295, 6);
				placeNumber($("#"+cpref+"-resultC101100"), 367, 248, 2);
				placeNumber($("#"+cpref+"-resultC101010"), 127, 256, 2);
				placeNumber($("#"+cpref+"-resultC100110"), 305, 300, 4);
				placeNumber($("#"+cpref+"-resultC011100"), 240, 110, 6);
				placeNumber($("#"+cpref+"-resultC011010"), 317, 108, 2);
				placeNumber($("#"+cpref+"-resultC010110"), 248, 338, 2);
				placeNumber($("#"+cpref+"-resultC001110"), 345, 180, 4);
				placeNumber($("#"+cpref+"-resultC111100"), 185, 140, 6);
				placeNumber($("#"+cpref+"-resultC111010"), 158, 245, 6);
				placeNumber($("#"+cpref+"-resultC110110"), 250, 310, 6);
				placeNumber($("#"+cpref+"-resultC101110"), 330, 240, 6);
				placeNumber($("#"+cpref+"-resultC011110"), 290, 133, 6);
				placeNumber($("#"+cpref+"-resultC111110"), 245, 210, 6);				
				$("#"+cpref+"-resultC000001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111111").css("left", -1000).css("top", -2200);

			} else if (vennSize == 4) {	
				drawEllipse(181,238,10,18.5,11.5,40,  changeOpacity(opts.colors[0], $("#"+cpref+"-label1").css('opacity')));
				drawEllipse(242,177,10,18.5,11.5,40,  changeOpacity(opts.colors[1], $("#"+cpref+"-label2").css('opacity')));
				drawEllipse(259,177,10,18.5,11.5,140, changeOpacity(opts.colors[2], $("#"+cpref+"-label3").css('opacity')));
				drawEllipse(320,238,10,18.5,11.5,140, changeOpacity(opts.colors[3], $("#"+cpref+"-label4").css('opacity')));
				
				$("#"+cpref+"-label1").css("left", 5).css("top", 70).css("color",   opts.colors[0]);
				$("#"+cpref+"-label2").css("left", 85).css("top", 5).css("color",   opts.colors[1]);
				$("#"+cpref+"-label3").css("left", 350).css("top", 5).css("color",  opts.colors[2]);
				$("#"+cpref+"-label4").css("left", 428).css("top", 70).css("color", opts.colors[3]);
				$("#"+cpref+"-label5").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label6").css("left", -1000).css("top", -2200);
				placeNumber($("#"+cpref+"-resultC100000"),  55, 190, 6);
				placeNumber($("#"+cpref+"-resultC010000"), 140,  60, 6);
				placeNumber($("#"+cpref+"-resultC001000"), 335,  60, 6);
				placeNumber($("#"+cpref+"-resultC000100"), 430, 190, 6);
				placeNumber($("#"+cpref+"-resultC110000"), 105, 120, 6);
				placeNumber($("#"+cpref+"-resultC101000"), 130, 260, 6);
				placeNumber($("#"+cpref+"-resultC100100"), 245, 340, 6);
				placeNumber($("#"+cpref+"-resultC011000"), 245,  90, 6);
				placeNumber($("#"+cpref+"-resultC010100"), 365, 260, 6);
				placeNumber($("#"+cpref+"-resultC001100"), 385, 120, 6);
				placeNumber($("#"+cpref+"-resultC111000"), 160, 170, 6);
				placeNumber($("#"+cpref+"-resultC110100"), 310, 290, 6);
				placeNumber($("#"+cpref+"-resultC101100"), 180, 290, 6);
				placeNumber($("#"+cpref+"-resultC011100"), 330, 170, 6);
				placeNumber($("#"+cpref+"-resultC111100"), 245, 220, 6);
				$("#"+cpref+"-resultC000010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111111").css("left", -1000).css("top", -2200);
			
			} else if (vennSize == 3) {
				drawEllipse(171,142,120,1,1,0, changeOpacity(opts.colors[0], $("#"+cpref+"-label1").css('opacity')));
				drawEllipse(327,142,120,1,1,0, changeOpacity(opts.colors[1], $("#"+cpref+"-label2").css('opacity')));
				drawEllipse(249,271,120,1,1,0, changeOpacity(opts.colors[2], $("#"+cpref+"-label3").css('opacity')));
				
				$("#"+cpref+"-label1").css("left", 55).css("top", 5).css("color",    opts.colors[0]);
				$("#"+cpref+"-label2").css("left", 380).css("top", 5).css("color",   opts.colors[1]);
				$("#"+cpref+"-label3").css("left", 220).css("top", 400).css("color", opts.colors[2]);
				$("#"+cpref+"-label4").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label5").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label6").css("left", -1000).css("top", -2200);
				placeNumber($("#"+cpref+"-resultC100000"), 120, 100, 8);
				placeNumber($("#"+cpref+"-resultC010000"), 360, 100, 8);
				placeNumber($("#"+cpref+"-resultC001000"), 245, 330, 8);
				placeNumber($("#"+cpref+"-resultC110000"), 245, 100, 8);
				placeNumber($("#"+cpref+"-resultC101000"), 170, 220, 8);
				placeNumber($("#"+cpref+"-resultC011000"), 320, 220, 8);
				placeNumber($("#"+cpref+"-resultC111000"), 245, 175, 8);
				$("#"+cpref+"-resultC000100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111111").css("left", -1000).css("top", -2200);
			
			} else if (vennSize == 2) {	
				drawEllipse(171,206,140,1,1,0, changeOpacity(opts.colors[0], $("#"+cpref+"-label1").css('opacity')));
				drawEllipse(327,206,140,1,1,0, changeOpacity(opts.colors[1], $("#"+cpref+"-label2").css('opacity')));
				
				$("#"+cpref+"-label1").css("left", 95).css("top", 40).css("color",  opts.colors[0]);
				$("#"+cpref+"-label2").css("left", 360).css("top", 40).css("color", opts.colors[1]);
				$("#"+cpref+"-label3").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label4").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label5").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label6").css("left", -1000).css("top", -2200);
				placeNumber($("#"+cpref+"-resultC100000"), 120, 195, 10);
				placeNumber($("#"+cpref+"-resultC010000"), 360, 195, 10);
				placeNumber($("#"+cpref+"-resultC110000"), 250, 195, 10);
				$("#"+cpref+"-resultC001000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111111").css("left", -1000).css("top", -2200);
			
			} else {
				drawEllipse(246,210,140,1,1,0, changeOpacity(opts.colors[0], $("#"+cpref+"-label1").css('opacity')));
				
				$("#"+cpref+"-label1").css("left", 225).css("top", 30).css("color", opts.colors[0]);
				$("#"+cpref+"-label2").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label3").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label4").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label5").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label6").css("left", -1000).css("top", -2200);
				placeNumber($("#"+cpref+"-resultC100000"), 245, 200, 10);
				$("#"+cpref+"-resultC010000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111111").css("left", -1000).css("top", -2200);
			}
		}
		
		function placeEdwardsVenn(vennSize) {
			if (vennSize == 6) {
				drawCircle(246,206,90,1,1,0, changeOpacity(opts.colors[0], $("#"+cpref+"-label1").css('opacity')));
				drawRoundRect(245,25,250,365,changeOpacity(opts.colors[1], $("#"+cpref+"-label2").css('opacity')));
				drawRoundRect(5,205,490,185,changeOpacity(opts.colors[2], $("#"+cpref+"-label3").css('opacity')));
				drawDoubleCircle(changeOpacity(opts.colors[3], $("#"+cpref+"-label4").css('opacity')));
				drawCross(changeOpacity(opts.colors[4], $("#"+cpref+"-label5").css('opacity')));
				drawCross2(changeOpacity(opts.colors[5], $("#"+cpref+"-label6").css('opacity')));
				
				$("#"+cpref+"-label1").css("left", 297).css("top", 90).css("color", opts.colors[0]);
				$("#"+cpref+"-label2").css("left", 267).css("top",  0).css("color", opts.colors[1]);
				$("#"+cpref+"-label3").css("left",  27).css("top",400).css("color", opts.colors[2]);
				$("#"+cpref+"-label4").css("left",  55).css("top", 90).css("color", opts.colors[3]);
				$("#"+cpref+"-label5").css("left", 282).css("top", 45).css("color", opts.colors[4]);
				$("#"+cpref+"-label6").css("left", 297).css("top",302).css("color", opts.colors[5]);
				placeNumber($("#"+cpref+"-resultC100000"), 204, 122, 1);
				placeNumber($("#"+cpref+"-resultC010000"), 390,  70, 6);
				placeNumber($("#"+cpref+"-resultC001000"),  85, 330, 6);
				placeNumber($("#"+cpref+"-resultC000100"),  85, 150, 6);
				placeNumber($("#"+cpref+"-resultC000010"), 225,  70, 3);
				placeNumber($("#"+cpref+"-resultC110000"), 280, 122, 1);
				placeNumber($("#"+cpref+"-resultC101000"), 204, 269, 1);
				placeNumber($("#"+cpref+"-resultC100100"), 172, 155, 1);
				placeNumber($("#"+cpref+"-resultC100010"), 218, 118, 1);
				placeNumber($("#"+cpref+"-resultC011000"), 390, 330, 6);
				placeNumber($("#"+cpref+"-resultC010100"), 390, 150, 6);
				placeNumber($("#"+cpref+"-resultC010010"), 258,  70, 3);
				placeNumber($("#"+cpref+"-resultC001100"),  85, 240, 6);
				placeNumber($("#"+cpref+"-resultC001010"), 225, 320, 3);
				placeNumber($("#"+cpref+"-resultC000110"), 120, 178, 2);
				placeNumber($("#"+cpref+"-resultC111000"), 280, 269, 1);
				placeNumber($("#"+cpref+"-resultC110100"), 314, 155, 1);
				placeNumber($("#"+cpref+"-resultC110010"), 266, 118, 1);
				placeNumber($("#"+cpref+"-resultC101100"), 172, 235, 1);
				placeNumber($("#"+cpref+"-resultC101010"), 219, 274, 1);
				placeNumber($("#"+cpref+"-resultC100110"), 167, 171, 1);
				placeNumber($("#"+cpref+"-resultC011100"), 390, 240, 6);
				placeNumber($("#"+cpref+"-resultC011010"), 258, 320, 3);
				placeNumber($("#"+cpref+"-resultC010110"), 365, 178, 2);
				placeNumber($("#"+cpref+"-resultC001110"), 120, 212, 2);
				placeNumber($("#"+cpref+"-resultC111100"), 315, 235, 1);
				placeNumber($("#"+cpref+"-resultC111010"), 266, 274, 1);
				placeNumber($("#"+cpref+"-resultC110110"), 320, 171, 1);
				placeNumber($("#"+cpref+"-resultC101110"), 166, 218, 1);
				placeNumber($("#"+cpref+"-resultC011110"), 365, 212, 2);
				placeNumber($("#"+cpref+"-resultC111110"), 320, 218, 1);
				placeNumber($("#"+cpref+"-resultC000001"), 180, 120, 1);
				placeNumber($("#"+cpref+"-resultC100001"), 199, 134, 1);
				placeNumber($("#"+cpref+"-resultC010001"), 304, 120, 1);
				placeNumber($("#"+cpref+"-resultC001001"), 179, 271, 1);
				placeNumber($("#"+cpref+"-resultC000101"), 169, 132, 1);
				placeNumber($("#"+cpref+"-resultC000011"), 234,  95, 1);
				placeNumber($("#"+cpref+"-resultC110001"), 286, 134, 1);
				placeNumber($("#"+cpref+"-resultC101001"), 199, 256, 1);
				placeNumber($("#"+cpref+"-resultC100101"), 190, 149, 1);
				placeNumber($("#"+cpref+"-resultC100011"), 228, 135, 3);
				placeNumber($("#"+cpref+"-resultC011001"), 306, 272, 1);
				placeNumber($("#"+cpref+"-resultC010101"), 317, 132, 1);
				placeNumber($("#"+cpref+"-resultC010011"), 250,  95, 1);
				placeNumber($("#"+cpref+"-resultC001101"), 167, 260, 1);
				placeNumber($("#"+cpref+"-resultC001011"), 233, 296, 1);
				placeNumber($("#"+cpref+"-resultC000111"), 143, 187, 1);
				placeNumber($("#"+cpref+"-resultC111001"), 285, 256, 1);
				placeNumber($("#"+cpref+"-resultC110101"), 295, 149, 1);
				placeNumber($("#"+cpref+"-resultC110011"), 255, 135, 3);
				placeNumber($("#"+cpref+"-resultC101101"), 190, 242, 1);
				placeNumber($("#"+cpref+"-resultC101011"), 227, 255, 3);
				placeNumber($("#"+cpref+"-resultC100111"), 213, 178, 6);
				placeNumber($("#"+cpref+"-resultC011101"), 318, 259, 1);
				placeNumber($("#"+cpref+"-resultC011011"), 250, 296, 1);
				placeNumber($("#"+cpref+"-resultC010111"), 342, 187, 1);
				placeNumber($("#"+cpref+"-resultC001111"), 143, 203, 1);
				placeNumber($("#"+cpref+"-resultC111101"), 295, 240, 1);
				placeNumber($("#"+cpref+"-resultC111011"), 256, 255, 3);
				placeNumber($("#"+cpref+"-resultC110111"), 270, 178, 6);
				placeNumber($("#"+cpref+"-resultC101111"), 213, 212, 6);
				placeNumber($("#"+cpref+"-resultC011111"), 342, 203, 1);
				placeNumber($("#"+cpref+"-resultC111111"), 270, 212, 6);
				
			} else if (vennSize == 5) {
				drawCircle(246,206,90,1,1,0, changeOpacity(opts.colors[0], $("#"+cpref+"-label1").css('opacity')));
				drawRoundRect(245,25,250,365,changeOpacity(opts.colors[1], $("#"+cpref+"-label2").css('opacity')));
				drawRoundRect(5,205,490,185, changeOpacity(opts.colors[2], $("#"+cpref+"-label3").css('opacity')));
				drawDoubleCircle(changeOpacity(opts.colors[3], $("#"+cpref+"-label4").css('opacity')));
				drawCross(changeOpacity(opts.colors[4], $("#"+cpref+"-label5").css('opacity')));
			
				$("#"+cpref+"-label1").css("left", 297).css("top", 90).css("color", opts.colors[0]);
				$("#"+cpref+"-label2").css("left", 267).css("top",  0).css("color", opts.colors[1]);
				$("#"+cpref+"-label3").css("left",  27).css("top",400).css("color", opts.colors[2]);
				$("#"+cpref+"-label4").css("left",  55).css("top", 90).css("color", opts.colors[3]);
				$("#"+cpref+"-label5").css("left", 282).css("top", 45).css("color", opts.colors[4]);
				$("#"+cpref+"-label6").css("left", -1000).css("top", -2200);
				placeNumber($("#"+cpref+"-resultC100000"), 200, 130, 3);
				placeNumber($("#"+cpref+"-resultC010000"), 390,  70, 6);
				placeNumber($("#"+cpref+"-resultC001000"),  85, 330, 6);
				placeNumber($("#"+cpref+"-resultC000100"),  85, 150, 6);
				placeNumber($("#"+cpref+"-resultC000010"), 225,  80, 4);
				placeNumber($("#"+cpref+"-resultC110000"), 284, 130, 3);
				placeNumber($("#"+cpref+"-resultC101000"), 200, 260, 3);
				placeNumber($("#"+cpref+"-resultC100100"), 185, 150, 4);
				placeNumber($("#"+cpref+"-resultC100010"), 228, 130, 3);
				placeNumber($("#"+cpref+"-resultC011000"), 390, 330, 6);
				placeNumber($("#"+cpref+"-resultC010100"), 390, 150, 6);
				placeNumber($("#"+cpref+"-resultC010010"), 258,  80, 4);
				placeNumber($("#"+cpref+"-resultC001100"),  85, 240, 6);
				placeNumber($("#"+cpref+"-resultC001010"), 225, 310, 4);
				placeNumber($("#"+cpref+"-resultC000110"), 133, 178, 6);
				placeNumber($("#"+cpref+"-resultC111000"), 284, 260, 3);
				placeNumber($("#"+cpref+"-resultC110100"), 300, 150, 4);
				placeNumber($("#"+cpref+"-resultC110010"), 255, 130, 3);
				placeNumber($("#"+cpref+"-resultC101100"), 185, 240, 4);
				placeNumber($("#"+cpref+"-resultC101010"), 228, 260, 3);
				placeNumber($("#"+cpref+"-resultC100110"), 213, 178, 6);
				placeNumber($("#"+cpref+"-resultC011100"), 390, 240, 6);
				placeNumber($("#"+cpref+"-resultC011010"), 258, 310, 4);
				placeNumber($("#"+cpref+"-resultC010110"), 355, 178, 6);
				placeNumber($("#"+cpref+"-resultC001110"), 133, 212, 6);
				placeNumber($("#"+cpref+"-resultC111100"), 300, 240, 4);
				placeNumber($("#"+cpref+"-resultC111010"), 256, 260, 3);
				placeNumber($("#"+cpref+"-resultC110110"), 270, 178, 6);
				placeNumber($("#"+cpref+"-resultC101110"), 213, 212, 6);
				placeNumber($("#"+cpref+"-resultC011110"), 355, 212, 6);
				placeNumber($("#"+cpref+"-resultC111110"), 270, 212, 6);
				$("#"+cpref+"-resultC000001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111111").css("left", -1000).css("top", -2200);
				
			} else if (vennSize == 4) {
				drawCircle(246,206,90,1,1,0,  changeOpacity(opts.colors[0], $("#"+cpref+"-label1").css('opacity')));
				drawRoundRect(245,25,250,365, changeOpacity(opts.colors[1], $("#"+cpref+"-label2").css('opacity')));
				drawRoundRect(5,205,490,185,  changeOpacity(opts.colors[2], $("#"+cpref+"-label3").css('opacity')));
				drawDoubleCircle(changeOpacity(opts.colors[3], $("#"+cpref+"-label4").css('opacity')));

				$("#"+cpref+"-label1").css("left", 267).css("top", 90).css("color", opts.colors[0]);
				$("#"+cpref+"-label2").css("left", 267).css("top",  0).css("color", opts.colors[1]);
				$("#"+cpref+"-label3").css("left",  27).css("top",400).css("color", opts.colors[2]);
				$("#"+cpref+"-label4").css("left",  27).css("top", 90).css("color", opts.colors[3]);
				$("#"+cpref+"-label5").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label6").css("left", -1000).css("top", -2200);
				placeNumber($("#"+cpref+"-resultC100000"), 215, 130, 6);
				placeNumber($("#"+cpref+"-resultC010000"), 390,  70, 6);
				placeNumber($("#"+cpref+"-resultC001000"),  85, 330, 6);
				placeNumber($("#"+cpref+"-resultC000100"),  85, 150, 6);
				placeNumber($("#"+cpref+"-resultC110000"), 265, 130, 6);
				placeNumber($("#"+cpref+"-resultC101000"), 215, 260, 6);
				placeNumber($("#"+cpref+"-resultC100100"), 200, 170, 6);
				placeNumber($("#"+cpref+"-resultC011000"), 390, 330, 6);
				placeNumber($("#"+cpref+"-resultC010100"), 390, 150, 6);
				placeNumber($("#"+cpref+"-resultC001100"),  85, 240, 6);
				placeNumber($("#"+cpref+"-resultC111000"), 265, 260, 6);
				placeNumber($("#"+cpref+"-resultC110100"), 280, 170, 6);
				placeNumber($("#"+cpref+"-resultC101100"), 200, 220, 6);
				placeNumber($("#"+cpref+"-resultC011100"), 390, 240, 6);
				placeNumber($("#"+cpref+"-resultC111100"), 280, 220, 6);
				$("#"+cpref+"-resultC000010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100010").css("left", -1000).css("top", -2200);			
				$("#"+cpref+"-resultC010010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111111").css("left", -1000).css("top", -2200);
			
			} else if (vennSize == 3) {
				drawCircle(246,206,110,1,1,0, changeOpacity(opts.colors[0], $("#"+cpref+"-label1").css('opacity')));
				drawRoundRect(245,25,250,365, changeOpacity(opts.colors[1], $("#"+cpref+"-label2").css('opacity')));
				drawRoundRect(5,205,490,185,  changeOpacity(opts.colors[2], $("#"+cpref+"-label3").css('opacity')));
				
				$("#"+cpref+"-label1").css("left", 95).css("top", 75).css("color",  opts.colors[0]);
				$("#"+cpref+"-label2").css("left", 267).css("top",  0).css("color", opts.colors[1]);
				$("#"+cpref+"-label3").css("left",  27).css("top",400).css("color", opts.colors[2]);
				$("#"+cpref+"-label4").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label5").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label6").css("left", -1000).css("top", -2200);
				placeNumber($("#"+cpref+"-resultC100000"), 190, 150, 8);
				placeNumber($("#"+cpref+"-resultC010000"), 390,  70, 8);
				placeNumber($("#"+cpref+"-resultC001000"),  85, 330, 8);
				placeNumber($("#"+cpref+"-resultC110000"), 290, 150, 8);
				placeNumber($("#"+cpref+"-resultC101000"), 190, 240, 8);
				placeNumber($("#"+cpref+"-resultC011000"), 390, 330, 8);
				placeNumber($("#"+cpref+"-resultC111000"), 290, 240, 8);
				$("#"+cpref+"-resultC000100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111111").css("left", -1000).css("top", -2200);
			
			} else if (vennSize == 2) {
				drawCircle(246,206,110,1,1,0, changeOpacity(opts.colors[0], $("#"+cpref+"-label1").css('opacity')));
				drawRoundRect(245,25,250,365, changeOpacity(opts.colors[1], $("#"+cpref+"-label2").css('opacity')));
				
				$("#"+cpref+"-label1").css("left", 95).css("top", 75).css("color",  opts.colors[0]);
				$("#"+cpref+"-label2").css("left", 267).css("top",  0).css("color", opts.colors[1]);
				$("#"+cpref+"-label3").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label4").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label5").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label6").css("left", -1000).css("top", -2200);
				placeNumber($("#"+cpref+"-resultC100000"), 180, 195, 10);
				placeNumber($("#"+cpref+"-resultC010000"), 420, 195, 10);
				placeNumber($("#"+cpref+"-resultC110000"), 290, 195, 10);
				$("#"+cpref+"-resultC001000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111111").css("left", -1000).css("top", -2200);
			
			} else {
				drawCircle(246,210,140,1,1,0, changeOpacity(opts.colors[0], $("#"+cpref+"-label1").css('opacity')));
				
				$("#"+cpref+"-label1").css("left", 225).css("top", 30).css("color", opts.colors[0]);
				$("#"+cpref+"-label2").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label3").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label4").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label5").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-label6").css("left", -1000).css("top", -2200);
				placeNumber($("#"+cpref+"-resultC100000"), 245, 200, 10);
				$("#"+cpref+"-resultC010000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111000").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111100").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111010").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111110").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC000111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111001").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC100111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC010111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC001111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111101").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111011").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC110111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC101111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC011111").css("left", -1000).css("top", -2200);
				$("#"+cpref+"-resultC111111").css("left", -1000).css("top", -2200);
			}
		}
		
		function addLegend(div, vennSize) {
			$t = div;
			var i=1;
			var hide = "";
			var div_legend = '<div class="module-legend">';
			$("*[id^="+cpref+"-label]").each(function(){
				div_legend += '<div id="'+cpref+'-item-'+i+'" name="' + $(this).text() + '" class="leg-items" style="opacity:0.5; background-color:' + $(this).css("color") + hide + '">';
				div_legend += '<span style="background-color:white; margin-left:0px; padding:0px 2px 0px 2px; transition: margin-left .3s ease-in-out;">off</span></div>';
				i += 1;
				if(i>vennSize) {
					hide = ";display:none";
				}
			});
			div_legend += '</div>';
			$t.append(div_legend);
			
			$("*[id^="+cpref+"-item]").hover(function(){
		        $(this).css('opacity', 0.75);
		        $(this).css('box-shadow',  '0px 0px 8px 1px lightgrey');
		    },function(){
		    	if($(this).children("span").text() === 'off') {
					$(this).css('opacity', 0.5);
		    	}
				$(this).css('box-shadow',  'none');
		    });
			$("*[id^="+cpref+"-item]").click(function(){
				var visible_id = null ;
				if ($(this).children("span").text() === 'off') {
					visible_id = select( $(this).attr('name') );
				}
				else {
					visible_id = unselect( $(this).attr('name') );
				}
				// Clear search
				if (opts.searchInput != null) {
					opts.searchInput.val("");
					if (opts.searchStatus != null) {
						opts.searchStatus.text("");
					}
				}
				// Draw
				clearCanvas();
            	if (opts.displayMode == 'edwards') {
                   	placeEdwardsVenn(vennSize);
                } else {
                   	placeClassicVenn(vennSize);
                }
            	if (opts.displayStat) {
            		placeStat(vennSize);
                }
            	if (visible_id != cpref+"-resultC000000") { 
            		$("#"+cpref+"-" + visible_id).show();
            	}
			});
		}
		
		function fillListVenn() {
			var classified=new Array();
			var actualList=new Array();
			actualList[0]=new Array();
			actualList[1]=new Array();
			actualList[2]=new Array();
			actualList[3]=new Array();
			actualList[4]=new Array();
			actualList[5]=new Array();
			
			for (m=0;m<opts.series.length;m++) {
				actualList[m]=new Array();
				var list = opts.series[m].data;
				for (t=0;t<list.length;t++) {
					if (list[t].length>0) {
						if (actualList[m][list[t]]) {
							actualList[m][list[t]]++;
						} else {
							actualList[m][list[t]]=1;
						}
						classified[list[t]]="C";
					}
				}
			}
			for (t=0;t<6;t++) {
				for (tt in actualList[t]) {
					if (classified[tt]) {
						classified[tt]=classified[tt]+"1";
					}
				}
				
				for (cl in classified) {
					if (classified[cl].length<t+2) {
						classified[cl]=classified[cl]+"0";
					}
				}
			}
			for (cl in classified) {
				var value = parseInt($("#"+cpref+"-result"+classified[cl]).html());
				if (opts.useValues) {
					for (var m=0;m<opts.series.length;m++) {
						// Is the value in the list
						var index_val = opts.series[m].data.indexOf(cl);
						if (index_val != -1) {
							value += parseInt(opts.series[m].values[index_val]);
						}
					}
				} else {
					value += 1;
				}
				$("#"+cpref+"-result"+classified[cl]).html(value);
				$("#"+cpref+"-result"+classified[cl]).removeClass("number-empty");
			}
			
			// Update the labels
			if (opts.series.length == 6) {
				$("#"+cpref+"-label1").html(opts.series[0].name);
				$("#"+cpref+"-label2").html(opts.series[1].name);
				$("#"+cpref+"-label3").html(opts.series[2].name);
				$("#"+cpref+"-label4").html(opts.series[3].name);
				$("#"+cpref+"-label5").html(opts.series[4].name);
				$("#"+cpref+"-label6").html(opts.series[5].name);
			} else if (opts.series.length == 5) {
				$("#"+cpref+"-label1").html(opts.series[0].name);
				$("#"+cpref+"-label2").html(opts.series[1].name);
				$("#"+cpref+"-label3").html(opts.series[2].name);
				$("#"+cpref+"-label4").html(opts.series[3].name);
				$("#"+cpref+"-label5").html(opts.series[4].name);
			} else if (opts.series.length == 4) {
				$("#"+cpref+"-label1").html(opts.series[0].name);
				$("#"+cpref+"-label2").html(opts.series[1].name);
				$("#"+cpref+"-label3").html(opts.series[2].name);
				$("#"+cpref+"-label4").html(opts.series[3].name);
			} else if (opts.series.length == 3) {
				$("#"+cpref+"-label1").html(opts.series[0].name);
				$("#"+cpref+"-label2").html(opts.series[1].name);
				$("#"+cpref+"-label3").html(opts.series[2].name);
			} else if (opts.series.length == 2) {
				$("#"+cpref+"-label1").html(opts.series[0].name);
				$("#"+cpref+"-label2").html(opts.series[1].name);
			} else if (opts.series.length == 1) {
				$("#"+cpref+"-label1").html(opts.series[0].name);
			}
			
			if (!opts.disableClick) {
				// Add some eventlistener
				$("*[id^="+cpref+"-resultC]").mouseover(function(){
	            	$(this).addClass("number-over");
	            });
	            $("*[id^="+cpref+"-resultC]").mouseout(function(){
	            	$(this).removeClass("number-over");
	            }); 
	            $("*[id^="+cpref+"-resultC]").click(opts.fnClickCallback);
			}
            
            // Add info to the number
            $("*[id^="+cpref+"-resultC]").each(function(){
            	this.listnames = new Array();
            	for (var i=cpref.length+7; i<$(this).attr("id").length; i++) {
            		if ($(this).attr("id").substring(i+1,i+2) == "1") {
            			try { this.listnames.push(opts.series[i-(cpref.length+7)].name); } catch(err) { }
            		}
            	}
            	this.empty = true;
				this.list = new Array();
				var cvalue = $(this).attr("id").substring(cpref.length+7,cpref.length+14);
				for (cl in classified) {
					if (classified[cl]==cvalue) {
						this.list.push(cl);
						this.empty = false;
					}
				}
            });
		}
		
		// Return an Array with number of common and specific element (x in 1 class, y in 2 class...) 
		function countByNbClass() {
			var data = new Array(0,0,0,0,0,0);
			$("*[id^="+cpref+"-resultC]").each(function(){
				var n = 0;
				for (var i=cpref.length+7; i<$(this).attr("id").length; i++) {
					n += $(this).attr("id").substring(i+1,i+2) == "1";
				}
				var val = $(this).text();
				if(val == '?') {
					val = $(this).children("span").attr('title');
				}
				data[n-1] += parseInt(val);
			});
			return data;
		}
		
		// Return an Array with size of each class 
		function sizeOfClass() {
			var data = new Array(0,0,0,0,0,0);
			$("*[id^="+cpref+"-resultC]").each(function(){
				for (var i=cpref.length+7; i<$(this).attr("id").length; i++) {
					if ($(this).attr("id").substring(i+1,i+2) == "1") {
						
						var val = $(this).text();
						if(val == '?') {
							val = $(this).children("span").attr('title');
						}
						data[i-(cpref.length+7)] += parseInt(val);
					}
				}
			});
			return data;
		}
		
		function setValuesforFillCountVenn (values) {
			if (values.A) { $("#"+cpref+"-resultC100000").html(values.A); }
			if (values.B) { $("#"+cpref+"-resultC010000").html(values.B); }
			if (values.C) { $("#"+cpref+"-resultC001000").html(values.C); }
			if (values.D) { $("#"+cpref+"-resultC000100").html(values.D); }
			if (values.E) { $("#"+cpref+"-resultC000010").html(values.E); }
			if (values.F) { $("#"+cpref+"-resultC000001").html(values.F); }
			if (values.AB) { $("#"+cpref+"-resultC110000").html(values.AB); }
			if (values.AC) { $("#"+cpref+"-resultC101000").html(values.AC); }
			if (values.AD) { $("#"+cpref+"-resultC100100").html(values.AD); }
			if (values.AE) { $("#"+cpref+"-resultC100010").html(values.AE); }
			if (values.AF) { $("#"+cpref+"-resultC100001").html(values.AF); }
			if (values.BC) { $("#"+cpref+"-resultC011000").html(values.BC); }
			if (values.BD) { $("#"+cpref+"-resultC010100").html(values.BD); }
			if (values.BE) { $("#"+cpref+"-resultC010010").html(values.BE); }
			if (values.BF) { $("#"+cpref+"-resultC010001").html(values.BF); }
			if (values.CD) { $("#"+cpref+"-resultC001100").html(values.CD); }
			if (values.CE) { $("#"+cpref+"-resultC001010").html(values.CE); }
			if (values.CF) { $("#"+cpref+"-resultC001001").html(values.CF); }
			if (values.DE) { $("#"+cpref+"-resultC000110").html(values.DE); }
			if (values.DF) { $("#"+cpref+"-resultC000101").html(values.DF); }
			if (values.EF) { $("#"+cpref+"-resultC000011").html(values.EF); }
			if (values.ABC) { $("#"+cpref+"-resultC111000").html(values.ABC); }
			if (values.ABD) { $("#"+cpref+"-resultC110100").html(values.ABD); }
			if (values.ABE) { $("#"+cpref+"-resultC110010").html(values.ABE); }
			if (values.ABF) { $("#"+cpref+"-resultC110001").html(values.ABF); }
			if (values.ACD) { $("#"+cpref+"-resultC101100").html(values.ACD); }
			if (values.ACE) { $("#"+cpref+"-resultC101010").html(values.ACE); }
			if (values.ACF) { $("#"+cpref+"-resultC101001").html(values.ACF); }
			if (values.ADE) { $("#"+cpref+"-resultC100110").html(values.ADE); }
			if (values.ADF) { $("#"+cpref+"-resultC100101").html(values.ADF); }
			if (values.AEF) { $("#"+cpref+"-resultC100011").html(values.AEF); }
			if (values.BCD) { $("#"+cpref+"-resultC011100").html(values.BCD); }
			if (values.BCE) { $("#"+cpref+"-resultC011010").html(values.BCE); }
			if (values.BCF) { $("#"+cpref+"-resultC011001").html(values.BCF); }
			if (values.BDE) { $("#"+cpref+"-resultC010110").html(values.BDE); }
			if (values.BDF) { $("#"+cpref+"-resultC010101").html(values.BDF); }
			if (values.BEF) { $("#"+cpref+"-resultC010011").html(values.BEF); }
			if (values.CDE) { $("#"+cpref+"-resultC001110").html(values.CDE); }
			if (values.CDF) { $("#"+cpref+"-resultC001101").html(values.CDF); }
			if (values.CEF) { $("#"+cpref+"-resultC001011").html(values.CEF); }
			if (values.DEF) { $("#"+cpref+"-resultC000111").html(values.DEF); }
			if (values.ABCD) { $("#"+cpref+"-resultC111100").html(values.ABCD); }
			if (values.ABCE) { $("#"+cpref+"-resultC111010").html(values.ABCE); }
			if (values.ABCF) { $("#"+cpref+"-resultC111001").html(values.ABCF); }
			if (values.ABDE) { $("#"+cpref+"-resultC110110").html(values.ABDE); }
			if (values.ABDF) { $("#"+cpref+"-resultC110101").html(values.ABDF); }
			if (values.ACDE) { $("#"+cpref+"-resultC101110").html(values.ACDE); }
			if (values.ACDF) { $("#"+cpref+"-resultC101101").html(values.ACDF); }
			if (values.BCDE) { $("#"+cpref+"-resultC011110").html(values.BCDE); }
			if (values.BCDF) { $("#"+cpref+"-resultC011101").html(values.BCDF); }
			if (values.CDEF) { $("#"+cpref+"-resultC001111").html(values.CDEF); }
			if (values.BDEF) { $("#"+cpref+"-resultC010111").html(values.BDEF); }
			if (values.BCEF) { $("#"+cpref+"-resultC011011").html(values.BCEF); }
			if (values.ADEF) { $("#"+cpref+"-resultC100111").html(values.ADEF); }
			if (values.ACEF) { $("#"+cpref+"-resultC101011").html(values.ACEF); }
			if (values.ABEF) { $("#"+cpref+"-resultC110011").html(values.ABEF); }
			if (values.ABCDE) { $("#"+cpref+"-resultC111110").html(values.ABCDE); }
			if (values.ABCDF) { $("#"+cpref+"-resultC111101").html(values.ABCDF); }
			if (values.ABCEF) { $("#"+cpref+"-resultC111011").html(values.ABCEF); }
			if (values.ABDEF) { $("#"+cpref+"-resultC110111").html(values.ABDEF); }
			if (values.ACDEF) { $("#"+cpref+"-resultC101111").html(values.ACDEF); }
			if (values.BCDEF) { $("#"+cpref+"-resultC011111").html(values.BCDEF); }
			if (values.ABCDEF) { $("#"+cpref+"-resultC111111").html(values.ABCDEF); }
		}

		function fillCountVenn() {
			
			// Update the labels
			if (opts.series[0].name.A) { $("#"+cpref+"-label1").html(opts.series[0].name.A); }
			if (opts.series[0].name.B) { $("#"+cpref+"-label2").html(opts.series[0].name.B); }
			if (opts.series[0].name.C) { $("#"+cpref+"-label3").html(opts.series[0].name.C); }
			if (opts.series[0].name.D) { $("#"+cpref+"-label4").html(opts.series[0].name.D); }
			if (opts.series[0].name.E) { $("#"+cpref+"-label5").html(opts.series[0].name.E); }
			if (opts.series[0].name.F) { $("#"+cpref+"-label6").html(opts.series[0].name.F); }
			
			// if values and data are given
			if (opts.series[0].hasOwnProperty("values") && opts.series[0].hasOwnProperty("data")) {
				setValuesforFillCountVenn(opts.series[0].values);
			// if only the data
			} else if (opts.series[0].hasOwnProperty("data")) {
				var values = {};
				for (var letter in opts.series[0].data) {
					values[letter] = opts.series[0].data[letter].length;
				}
				setValuesforFillCountVenn(values);
			// if only the values
			} else if(opts.series[0].hasOwnProperty("values")) {
				setValuesforFillCountVenn(opts.series[0].values);
			}

			// Add info to the number
            $("*[id^="+cpref+"-resultC]").each(function(){
            	this.listnames = new Array();
            	for (var i=cpref.length+7; i<$(this).attr("id").length; i++) {
            		if ($(this).attr("id").substring(i+1,i+2) == "1") {
            			try {
            				if(i-(cpref.length+7) == 0)		{ this.listnames.push(opts.series[0].name.A); }
            				else if(i-(cpref.length+7) == 1)	{ this.listnames.push(opts.series[0].name.B); } 
            				else if(i-(cpref.length+7) == 2)	{ this.listnames.push(opts.series[0].name.C); }
            				else if(i-(cpref.length+7) == 3)	{ this.listnames.push(opts.series[0].name.D); }
            				else if(i-(cpref.length+7) == 4)	{ this.listnames.push(opts.series[0].name.E); }
            				else if(i-(cpref.length+7) == 5)	{ this.listnames.push(opts.series[0].name.F); }
            				}
            			catch(err) { }
            		}
            	}
            	// if data are provided
            	if (opts.series[0].hasOwnProperty("data")) {
            		var cvalue = "";
    				for (var i=cpref.length+7; i<$(this).attr("id").length; i++) {
    					if ($(this).attr("id").substring(i+1,i+2) == "1") {
    						if (i == (cpref.length+7)) {cvalue += "A";}
    						else if (i == (cpref.length+8)) {cvalue += "B";}
    						else if (i == (cpref.length+9)) {cvalue += "C";}
    						else if (i == (cpref.length+10)) {cvalue += "D";}
    						else if (i == (cpref.length+11)) {cvalue += "E";}
    						else if (i == (cpref.length+12)) {cvalue += "F";}
    					}
    				}
    				this.list = opts.series[0].data[cvalue];
					if (!opts.disableClick) {
						// Add some eventlistener
						$(this).mouseover(function(){
			            	$(this).addClass("number-over");
			            });
			            $(this).mouseout(function(){
			            	$(this).removeClass("number-over");
			            }); 
			            $(this).click(opts.fnClickCallback);
					}
            	}
					
				if($(this).text() > 0) {
            		$(this).removeClass("number-empty");
            	}
            });
            
		}
		
		function getVennType() {
			// If more than 1 sample, it's a list venn
			if (opts.series.length > 1) {
				return (new Array("list", opts.series.length));
			} else {
				if (opts.series[0] != undefined && opts.series[0].name.A) {
					var count = 0;
					for (i in opts.series[0].name){ count++; }
					return (new Array("count", count));
				} else {
					return (new Array("list", opts.series.length));
				}
			}
		}
		
		function addSVGText(text, x, y, color, opacity, ffamily, fsize, fweight, isLabel) {
			if (isLabel) {
				x = parseInt(x.replace("px", "")) + 20;
				y = parseInt(y.replace("px", "")) + 14;
			} else {
				x = parseInt(x.replace("px", "")) + 5;
				y = parseInt(y.replace("px", "")) + 14;
			}
			return '<text fill="'+color+'" fill-opacity="'+opacity+'" stroke="none" font-family="'+ffamily+'" font-size="'+fsize+'"pt font-style="normal" font-weight="'+fweight+'" text-decoration="none" x="'+x+'"pt y="'+y+'"pt text-anchor="middle" dominant-baseline="alphabetic">'+text+'</text>';
		}
		
		function addExportModule(div, extraheight, type){
			$t = div;
			
			var div_export = '<div id="'+cpref+'-module-export" style="position: relative; left:475px; top: -'+(418+extraheight)+'px; width: 25px; height: 20px;">';
			div_export += '<a id="'+cpref+'-data-export" href="" download="" style="display:none;"></a>';
			div_export += '<canvas id="'+cpref+'-canvasExport" style="border:1px solid white" width="25" height="20"></canvas>';
        	div_export += '<div id="'+cpref+'-export-menu" style="position: relative;width:150px; height:30px; display:none; right:123px; top:-4px;">';
        	div_export += '<div style="box-shadow: 3px 3px 10px rgb(136, 136, 136); border: 1px solid rgb(160, 160, 160); background: none repeat scroll 0% 0% rgb(255, 255, 255);padding: 5px 0px;">';
        	div_export += '<div id="'+cpref+'-format-png" style="padding: 0px 10px; background: none repeat scroll 0% 0% transparent; color: rgb(48, 48, 48); font-size: 12px;">Download PNG image</div>';
        	div_export += '<div id="'+cpref+'-format-svg" style="padding: 0px 10px; background: none repeat scroll 0% 0% transparent; color: rgb(48, 48, 48); font-size: 12px;">Download SVG image</div>';
        	if (type == "list") {
        		div_export += '<div id="'+cpref+'-format-csv" style="padding: 0px 10px; background: none repeat scroll 0% 0% transparent; color: rgb(48, 48, 48); font-size: 12px;">Download CSV lists</div>';
        	}
        	div_export += '</div>';
			div_export += '</div>';
			div_export += '</div>';
			$t.append(div_export);
			
			//draw canvas button
			var ecanvas = $("#"+cpref+"-canvasExport")[0],
				econtext = ecanvas.getContext("2d");
			for (i=0;i<3;i++){
				econtext.lineWidth = 3;
				econtext.beginPath();
				econtext.lineCap="round";
				econtext.moveTo(5,5+i*5.2);
				econtext.lineTo(20,5+i*5.2);
				econtext.strokeStyle = "#666";
				econtext.stroke();
			}
			
			var select_form = $("#"+cpref+"-export-menu");
			var ceColorOri = $("#"+cpref+"-canvasExport").css('background');
			$("#"+cpref+"-format-png").click(function(event) {
				$("#"+cpref+"-canvasExport").css('background', ceColorOri);
				$("#"+cpref+"-canvasExport").css('border-color', "white");
				$("#"+cpref+"-canvasExport").hover(function() {
					$(this).css('background', 'linear-gradient(to bottom, white, #AECEFF) repeat scroll 0 0 transparent');
					$(this).css('border', '1px solid #6688AA');
					$(this).css('border-radius', '3px');
				}, function() {
					$(this).css('background', ceColorOri);
					$(this).css('border-color', "white");
				});
				select_form.hide();
				var export_canvas = document.createElement( "canvas" );
				export_ctx = export_canvas.getContext( "2d" );
				export_canvas.width = __canvasWidth;
				export_canvas.height = __canvasHeight;
				var img = document.createElement( "img" ),
					$expDiv = $("<div></div>");
				$expDiv.append(__context.getSerializedSvg(true));
				
				$("[id^="+cpref+"-resultC]").each(function(){					
					if($(this).css("display") != "none") {
						if ($(this).html().indexOf("<span") === 0) {
							value = "?";
						} else {
							value = $(this).html();
						}
						$expDiv.children("svg").children("g").append(addSVGText(value, 
							$(this).css("left"),  $(this).css("top"),
							$(this).css("color"), $(this).css("opacity"),
							$(this).css("font-family"), $(this).css("font-size"),
							$(this).css("font-weight"), false));
					}
				});
				$("[id^="+cpref+"-label]").each(function(){
					$expDiv.children("svg").children("g").append(addSVGText($(this).html(), 
							$(this).css("left"),  $(this).css("top"), 
							$(this).css("color"), $(this).css("opacity"),
							$(this).css("font-family"), $(this).css("font-size"),
							$(this).css("font-weight"), true));
				});

				img.setAttribute( "src", "data:image/svg+xml;base64," + window.btoa($expDiv.html()));
				img.setAttribute('height', '700px');
				img.setAttribute('width', '500px');
				img.onload = function() {
					export_ctx.drawImage( img, 0, 0 );
					var export_img = export_canvas.toDataURL("image/png");
					$('#'+cpref+'-data-export').attr("href", export_img);
					$('#'+cpref+'-data-export').attr("download", "jVenn_chart.png");
					$('#'+cpref+'-data-export')[0].click();
				};
			});
			$("#"+cpref+"-format-svg").click(function(event) {
				$("#"+cpref+"-canvasExport").css('background', ceColorOri);
				$("#"+cpref+"-canvasExport").css('border-color', "white");
				$("#"+cpref+"-canvasExport").hover(function() {
					$(this).css('background', 'linear-gradient(to bottom, white, #AECEFF) repeat scroll 0 0 transparent');
					$(this).css('border', '1px solid #6688AA');
					$(this).css('border-radius', '3px');
				}, function() {
					$(this).css('background', ceColorOri);
					$(this).css('border-color', "white");
				});
				select_form.hide();
				
				var $expDiv = $("<div></div>");
				$expDiv.append(__context.getSerializedSvg(true));
				
				$("[id^="+cpref+"-resultC]").each(function(){					
					if($(this).css("display") != "none") {
						if ($(this).html().indexOf("<span") === 0) {
							value = "?";
						} else {
							value = $(this).html();
						}
						$expDiv.children("svg").children("g").append(addSVGText(value, 
							$(this).css("left"),  $(this).css("top"),
							$(this).css("color"), $(this).css("opacity"),
							$(this).css("font-family"), $(this).css("font-size"),
							$(this).css("font-weight"), false));
					}
				});
				$("[id^="+cpref+"-label]").each(function(){
					$expDiv.children("svg").children("g").append(addSVGText($(this).html(), 
							$(this).css("left"),  $(this).css("top"), 
							$(this).css("color"), $(this).css("opacity"),
							$(this).css("font-family"), $(this).css("font-size"),
							$(this).css("font-weight"), true));
				});
				
				var svgContent = "data:image/svg+xml;base64," + window.btoa($expDiv.html()),
					encodedUri = encodeURI(svgContent);
				$('#'+cpref+'-data-export').attr("href", encodedUri);
				$('#'+cpref+'-data-export').attr("download", "jVenn_chart.svg");
				$('#'+cpref+'-data-export')[0].click();
				
			});
			$("#"+cpref+"-format-csv").click(function(event) {
				$("#"+cpref+"-canvasExport").css('background', ceColorOri);
				$("#"+cpref+"-canvasExport").css('border-color', "white");
				$("#"+cpref+"-canvasExport").hover(function() {
					$(this).css('background', 'linear-gradient(to bottom, white, #AECEFF) repeat scroll 0 0 transparent');
					$(this).css('border', '1px solid #6688AA');
					$(this).css('border-radius', '3px');
				}, function() {
					$(this).css('background', ceColorOri);
					$(this).css('border-color', "white");
				});
				select_form.hide();
				var	rawData = new Array(),
					comma = false;
				$("*[id^="+cpref+"-resultC]").each(function(){
					if (!this.empty) {
						var	currentRow = new Array(),
							tmpline = this.listnames.join("|");
						
						if (tmpline.indexOf(',') >= 0) {
							comma = true;
							tmpline = tmpline.replace(/,/g, '_');
						}
						currentRow.push(tmpline);
						
						for (var i=0 ; i < this.list.length ; i++) {
							tmpline = this.list[i];
							if (this.list[i].indexOf(",") >= 0) {
								comma = true;
								tmpline = this.list[i].replace(/,/g, '_');
							}
							currentRow.push(tmpline);
						}
						rawData.push(currentRow);
					}
	            });
				var csvContent = "data:text/csv;charset=utf-8,";
				if(comma) {
					csvContent += "##\n## Warning: comma(s) have been replaced by underscore(s)\n##\n";
				}
				transpose(rawData).forEach(function(infoArray, index){
					csvContent += infoArray.join(",") + "\n";
				});
				var encodedUri = encodeURI(csvContent);
				$('#'+cpref+'-data-export').attr("href", encodedUri);
				$('#'+cpref+'-data-export').attr("download", "jVenn.csv");
				$('#'+cpref+'-data-export')[0].click();
			});
			$("#"+cpref+"-canvasExport").click(function (event) {
				$(this).css('background', 'linear-gradient(to bottom, #AECEFF, white) repeat scroll 0 0 transparent');
				$(this).css('border', '1px solid #6688AA');
				$('#'+cpref+'-canvasExport').unbind('mouseenter mouseleave');
				select_form.show();
				if (select_form.is(":visible")){
					$(document).mouseup(function (event){
						if (select_form.has(event.target).length === 0){
							$("#"+cpref+"-canvasExport").css('background', ceColorOri);
							$("#"+cpref+"-canvasExport").css('border-color', "white");
							$("#"+cpref+"-canvasExport").hover(function() {
								$(this).css('background', 'linear-gradient(to bottom, white, #AECEFF) repeat scroll 0 0 transparent');
								$(this).css('border', '1px solid #6688AA');
								$(this).css('border-radius', '3px');
							}, function() {
								$(this).css('background', ceColorOri);
								$(this).css('border-color', "white");
							});
							select_form.hide();
						}
					});
				}
				var colorOrig=$("div[id^="+cpref+"-format-]").css('background');
				$("div[id^="+cpref+"-format-]").hover(function() {
					$(this).css('background', 'rgba(69,114,165,0.75)');
					$(this).css('color', 'white');
				}, function() {
					$(this).css('background', colorOrig);
					$(this).css('color', '');
				});
			});
			$("#"+cpref+"-canvasExport").hover(function() {
				$(this).css('background', 'linear-gradient(to bottom, white, #AECEFF) repeat scroll 0 0 transparent');
				$(this).css('border', '1px solid #6688AA');
				$(this).css('border-radius', '3px');
			}, function() {
				$(this).css('background', ceColorOri);
				$(this).css('border-color', "white");
			});
		}
		
		function unselect( unselected_name, init ) {
			var init = (init == null ? true : init);
			var unselected_idx = null ;
			// Find group index
			$("*[id^="+cpref+"-label]").each(function(){
				if( $(this).html() == unselected_name ) {
					$(this).css('opacity', 0.1);
					$(this).removeClass("is-selected");
					var id = $(this).attr('id');
					unselected_idx = id.charAt(id.length-1);
				}
			});
			// Change legend button
			legend_button = $("#"+cpref+"-item-" + unselected_idx) ;
			legend_button.children("span").css('margin-left', '0px');
			legend_button.css('opacity', 0.5);
			legend_button.children("span").text('off');
			// Change visibility
			var visible_id = cpref+"-resultC000000";
			$(".is-selected").each(function(){
				var id = $(this).attr('id');
				var replace_pos = cpref.length+7 + parseInt(id.charAt(id.length-1));
				visible_id = visible_id.substr(0, replace_pos) + 1 + visible_id.substr(replace_pos+1);
			});
			if (visible_id == cpref+"-resultC000000" && init) { // 0 selected element and init
				$("*[id^="+cpref+"-resultC]").each(function(){
					$(this).css('opacity', 1);
				});
				$("*[id^="+cpref+"-label]").each(function(){
					$(this).css('opacity', 0.6);
				});
			} else { // At least one selected element or non-init
				$("*[id^="+cpref+"-resultC]").each(function(){
					if( $(this).attr('id') == visible_id ) {
						$(this).css('opacity', 1);
					} else {
						$(this).css('opacity', 0.1);
					}
	        	});
			}
			return( visible_id );
		}
		
		function select( selected_name ) {
			var selected_idx = null ;
			// Find group index and show label
			$("*[id^="+cpref+"-label]").each(function(){
				if( $(this).html() == selected_name ) {
					$(this).css('opacity', 0.6);
					$(this).addClass("is-selected");
					var id = $(this).attr('id');
					selected_idx = id.charAt(id.length-1);
				} else if (!$(this).hasClass("is-selected")) {
					$(this).css('opacity', 0.1);
				}
			});
			// Change legend button
			legend_button = $("#"+cpref+"-item-" + selected_idx);
			legend_button.children("span").css('margin-left', '13px');
			legend_button.css('opacity', 0.75);
			legend_button.css('color', 'black');
			legend_button.children("span").text('on');
			// Change count visibility
			var visible_idx = cpref+"-resultC000000";
			$(".is-selected").each(function(){
				var id = $(this).attr('id');
				var replace_pos = cpref.length+7 + parseInt(id.charAt(id.length-1));
				visible_idx = visible_idx.substr(0, replace_pos) + 1 + visible_idx.substr(replace_pos+1);
			});
			$("*[id^="+cpref+"-resultC]").each(function(){
				if( $(this).attr('id') == visible_idx ) {
					$(this).css('opacity', 1);
				} else {
					$(this).css('opacity', 0.1);
				}
        	});
			return( visible_idx );
		}

		function search( val, min_size ) {
			var min_size = (min_size == null ? 1 : min_size)
			var groups_status = new Array();
			var visible_id = cpref+"-resultC000000";
			var nb_find = 0;
			$("*[id^="+cpref+"-label]").each( function() {
				if ($(this).html() != "") {
					groups_status[$(this).html()] = 'unselected';
				}
			});
			
			if (val == "" || min_size > val.length) { // Search is empty
				for (var group_name in  groups_status) {
					visible_id = unselect( group_name );
				}
			} else { // Search element
				// Find selected
				var listnames = null ;
				var perfect_match = false;
				$("*[id^="+cpref+"-resultC]:not(.number-empty)").each( function() {
					for (var idx = 0 ; idx < this.list.length && !perfect_match ; idx++) {
						if (this.list[idx].indexOf(val) != -1) {
							if (this.list[idx] == val) {
								perfect_match = true;
							}
							listnames = this.listnames;
							nb_find++ ;
						}
					}
				});
				if (nb_find == 1  || perfect_match) {
					for ( var idx_2 = 0 ; idx_2 < listnames.length ; idx_2++ ) {
						groups_status[listnames[idx_2]] = 'selected' ;
					}
					nb_find = 1;
				}
				// Update if search is not ambiguous
				if (nb_find <= 1) { // Find or unfind
					for (var group_name in  groups_status) {
						if( groups_status[group_name] == 'selected' ) {
							visible_id = select( group_name );
						} else {
							visible_id = unselect( group_name, false );
						}
					}
				} else { // Ambiguous
					for (var group_name in  groups_status) {
						visible_id = unselect( group_name );
					}
				}
			}
			
			// Update status
			if( opts.searchStatus != null ) {
				if (val == "" || min_size > val.length) {
					opts.searchStatus.text( "" );
				} else if (nb_find == 0) {
					opts.searchStatus.text( "not found" );
				} else if (nb_find == 1) {
					opts.searchStatus.text( "found" );
				} else if (nb_find > 1) {
					opts.searchStatus.text( "ambiguous" );
				} 
			}
			
			// Display
			clearCanvas();
        	if (opts.displayMode == 'edwards') {
               	placeEdwardsVenn( getVennType()[1] );
            } else {
               	placeClassicVenn( getVennType()[1] );
            }
        	if (opts.displayStat) {
        		placeStat( getVennType()[1] );
            }
        	if (visible_id != cpref+"-resultC000000") { 
        		$("#" + visible_id).show();
        	}
		}
		
		// Container prefix need for multiple jvenn on the same page
		var cpref = $(this).attr('id');
        
		this.each(function() {
            var $t = $(this);
            var extraheight = 0;
            if (opts.displayStat) {
            	extraheight = 280;
            	__canvasHeight += extraheight;
            }
            $t.css({"width": "500px", "height": ""+(450+extraheight)+"px"});
            $('<style>.number-black{font-weight:bold;color:#000000;text-decoration:none;font-size:12px;}.number-over{font-weight:bold;cursor:pointer;color:#0000FF;text-decoration:underline;}.number-empty{font-weight:normal;font-size:12px;}</style>').appendTo('body');
           
            var	legleft = 434,
            	legtop  = 130;
            if (opts.displayMode == 'edwards') {
            	legleft = 10;
            	legtop  = 425;
            }
            $('<style>.module-legend{border:1px solid lightgrey;border-radius:5px;position:relative;left:'+legleft+'px;top:-'+(legtop+extraheight)+'px;width:35px;padding-bottom:3px}</style>').appendTo('body');
            $('<style>.leg-items{padding-top:1px;margin:3px 3px 0px 3px;cursor:pointer;border: 1px solid grey;border-radius:2px;width:27px;height:11px;font-size:0.65em;line-height:10px;opacity:0.75}</style>').appendTo('body');
            
            var div_content = '<div id="'+cpref+'-frame" style="position: relative; left: 0pt; top: 5pt; width: 500px; height: "'+(445+extraheight)+'px;">';
			div_content += '<div id="'+cpref+'-canvasEllipse" width="500px" height="'+(415+extraheight)+'px;"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC100000"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC010000"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC001000"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC000100"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC000010"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC110000"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC101000"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC100100"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC100010"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC011000"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC010100"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC010010"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC001100"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC001010"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC000110"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC111000"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC110100"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC110010"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC101100"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC101010"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC100110"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC011100"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC011010"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC010110"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC001110"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC111100"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC111010"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC110110"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC101110"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC011110"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC111110"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC000001"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC000011"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC000101"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC000111"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC001001"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC001011"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC001101"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC001111"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC010001"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC010011"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC010101"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC010111"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC011001"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC011011"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC011101"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC011111"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC100001"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC100011"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC100101"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC100111"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC101001"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC101011"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC101101"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC101111"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC110001"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC110011"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC110101"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC110111"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC111001"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC111011"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC111101"></div>';
			div_content += '<div class="number-black" style="position: absolute; left: -1000px; top: -2200px;" id="'+cpref+'-resultC111111"></div>';
			div_content += '<div style="position: absolute; left: -1000px; top: -1000px; opacity:0.5;" id="'+cpref+'-label1"></div>';
			div_content += '<div style="position: absolute; left: -1000px; top: -1000px; opacity:0.5;" id="'+cpref+'-label2"></div>';
			div_content += '<div style="position: absolute; left: -1000px; top: -1000px; opacity:0.5;" id="'+cpref+'-label3"></div>';
			div_content += '<div style="position: absolute; left: -1000px; top: -1000px; opacity:0.5;" id="'+cpref+'-label4"></div>';
			div_content += '<div style="position: absolute; left: -1000px; top: -1000px; opacity:0.5;" id="'+cpref+'-label5"></div>';
			div_content += '<div style="position: absolute; left: -1000px; top: -1000px; opacity:0.5;" id="'+cpref+'-label6"></div>';
			div_content += '</div>';
            $t.html(div_content);
           
            // define internal variable
            __context = new C2S(__canvasWidth,__canvasHeight);
			
            // init with 0
            $("*[id^="+cpref+"-resultC]").each(function(){
				$(this).html(0);
				$(this).addClass("number-empty");
				$(this).css("font-size", opts.fontSize);
				$(this).css("font-family", opts.fontFamily);
			});
            $("*[id^="+cpref+"-label]").each(function(){
				$(this).css("font-family", opts.fontFamily);
            });
            var type = getVennType(); 
            if (type[0] == "list") {
            	fillListVenn();
            } else if (type[0] == "count") {
            	fillCountVenn();
            }
            
			if( opts.searchInput != null ) {
				opts.searchInput.keyup( function() {
					search( opts.searchInput.val(), opts.searchMinSize );
				});
			}
			
        	if (opts.displayMode == 'edwards') {
        		placeEdwardsVenn(type[1]);
        	} else {
        		placeClassicVenn(type[1]);
        	}
            if (opts.displayStat) {
            	placeStat(type[1]);
            }
            
            // if the exporting modul is requested
            if (opts.exporting === true){ addExportModule($t, extraheight, type[0]); }
            // if min 4 classes diagram is requested add legend
            if (opts.displaySwitch === true && type[1] >= 4) { addLegend($t, type[1]); }
            // number hover action
            $("*[id^="+cpref+"-resultC]").hover(
            	function(){
            		var activeleg = false;
            		$("*[id^="+cpref+"-item]").each(function(){
            			if($(this).children("span").text() === 'on') { activeleg = true; }
            		});
            		var activesearch = false;
            		if(opts.searchInput != null && opts.searchInput.val() != '') { activesearch= true; }
            		if(!activeleg && !activesearch &&  $(this).text() !== "") {
		            	var labels  = this.listnames;
		            	var current = this;
	            		$("*[id^="+cpref+"-label]").each(function(){
	            			if (labels.indexOf($(this).text()) < 0) {
	            				$(this).css('opacity', 0.1);
	            			} else {
	            				$(this).css('opacity', 0.6);
	            			}
	            		});
	            		$("*[id^="+cpref+"-resultC]").each(function(){
	                		if(this != current) {
	                			$(this).css('opacity', 0.1);
	                		}
	                	});
		            	clearCanvas();
		            	if (opts.displayMode == 'edwards') {
		                   	placeEdwardsVenn(type[1]);
		                } else {
		                   	placeClassicVenn(type[1]);
		                }
		            	if (opts.displayStat) {
		            		placeStat(type[1]);
		                }
            		}
            	},
            	function(){
            		var activeleg = false;
            		$("*[id^="+cpref+"-item]").each(function(){
            			if($(this).children("span").text() === 'on') { activeleg = true; }
            		});
            		var activesearch = false;
            		if(opts.searchInput != null && opts.searchInput.val() != '') { activesearch= true; }
            		if(!activeleg && !activesearch && $(this).text() !== "") {
	            		var labels = this.listnames;
	                	$("*[id^="+cpref+"-label]").each(function(){
	            			$(this).css('opacity', 0.5);
	                	});
	                	$("*[id^="+cpref+"-resultC]").each(function(){
	                		$(this).css('opacity', 1);
	                	});
	                	clearCanvas();
	                	if (opts.displayMode == 'edwards') {
	                       	placeEdwardsVenn(type[1]);
	                    } else {
	                       	placeClassicVenn(type[1]);
	                    }
	                	if (opts.displayStat) {
	                		placeStat(type[1]);
	                    }
            		}
            	}
            );
            $("#"+cpref+"-canvasEllipse").html(__context.getSvg());
        });
        return this;
	};
})(jQuery);
