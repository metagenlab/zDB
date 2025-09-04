;(function($) {

  $.fn.hottie = function(options) {

    var plugin = $(this);
    
    var init = function() {
      var max = settings.preMax;
      var min = settings.preMin;
      
      plugin.each(function() {
        var val = parseFloat(settings.readValue($(this)));
        if (val != null && !isNaN(val)) {
          min = Math.min( min, val );
          max = Math.max( max, val );
        }
        $(this).data('val', val);
      });

      var range = max - min;
      plugin.each(function() {
        var val = $(this).data('val');
        var c = settings.nullColor;
        if (val != null && !isNaN(val) && range != 0) {
          var adj = val - min;
          var pct = 1.0 * adj / range;
          c = getColor(pct);
        }
        $(this).css('background-color', c);
      });
    };

    var settings = {
      readValue : function(e) {
        return parseFloat(e.html());
      },
      preMin : 0,
      preMax : 0,
      colorArray : [ 
        "#2E3436",
        "#384471",
        "#4a538a",
        "#0094fe",
        "#00c7c9",
        "#6ad963",
        "#fff44e",
        "#fda659",
        "#ff4a53"
      ],
      nullColor : "#EFFBFB"
    };

    if ( options ) { 
      $.extend( settings, options );
    }

    function hex2num(hex) {
      if(hex.charAt(0) == "#") hex = hex.slice(1); //Remove the '#' char - if there is one.
      var triplet = new Array();
      var int1,int2;
      for(var i=0;i<6;i+=2) {
          var part = hex.charAt(i) + hex.charAt(i+1); 
          triplet.push(parseInt(part,16));
      }
      return triplet;
    }

    function num2hex(triplet) {
      var result = "#";
      for(var i=0;i<3;i++) {
          var hex = Math.round(triplet[i]).toString(16);
          while (hex.length < 2) {
            hex = "0" + hex;
          }
          result += hex; 
      }
      return result;
    }

    function rgbToHsv(triplet){
        var r, g, b;
        r = triplet[0];
        g = triplet[1];
        b = triplet[2];

        r = r/255, g = g/255, b = b/255;
        var max = Math.max(r, g, b), min = Math.min(r, g, b);
        var h, s, v = max;

        var d = max - min;
        s = max == 0 ? 0 : d / max;

        if(max == min){
            h = 0; // achromatic
        }else{
            switch(max){
                case r: h = (g - b) / d + (g < b ? 6 : 0); break;
                case g: h = (b - r) / d + 2; break;
                case b: h = (r - g) / d + 4; break;
            }
            h /= 6;
        }
        h = Math.max(0, Math.min(h, 1));
        s = Math.max(0, Math.min(s, 1));
        v = Math.max(0, Math.min(v, 1));
        return [h, s, v];
    }

    function hsvToRgb(triplet){
      var h = triplet[0];
      var s = triplet[1];
      var v = triplet[2];

      var r, g, b;

      var i = Math.floor(h * 6);
      var f = h * 6 - i;
      var p = v * (1 - s);
      var q = v * (1 - f * s);
      var t = v * (1 - (1 - f) * s);

      switch(i % 6){
          case 0: r = v, g = t, b = p; break;
          case 1: r = q, g = v, b = p; break;
          case 2: r = p, g = v, b = t; break;
          case 3: r = p, g = q, b = v; break;
          case 4: r = t, g = p, b = v; break;
          case 5: r = v, g = p, b = q; break;
      }
      r = Math.min(r, 1);
      g = Math.min(g, 1);
      b = Math.min(b, 1);
      return [r * 255, g * 255, b * 255];
      }

      var getColor = function(percent){
        
      if (percent == null)
        return nullColor;

      var colors = settings.colorArray.length;

      var colorPosition = (percent * (colors - 1));
      var sIndex = Math.floor(colorPosition);
      sIndex = Math.min(colors - 2, sIndex);

      var s = settings.colorArray[sIndex];
      var e = settings.colorArray[sIndex+1];

      var sHSL = rgbToHsv(hex2num(s));
      var eHSL = rgbToHsv(hex2num(e));

      var interiorPercent = (percent * (colors - 1)) - sIndex;

      var hsvResult = transition3(interiorPercent, 1, sHSL, eHSL);

      var dispRGB = hsvToRgb(hsvResult);
      return num2hex(dispRGB);
    };

    function transition(value, maximum, start_point, end_point){
      var r = start_point + (end_point - start_point)*value/maximum
      return r;
    }

    function transition3(value, maximum, s, e) {

      //handle situation where greyscale colors have red hue
      if (s[1] == 0) s[0] = e[0];
      if (e[1] == 0) e[0] = s[0];


      // handle black saturation issue
      if (s[2] == 0) s[1] = e[1];
      if (e[2] == 0) e[1] = s[1];
      
      
      var p = value / maximum;
      // hue has to wrap correctly around zero

      var distCCW = (s[0] >= e[0]) ? s[0] - e[0] : 1 + s[0] - e[0];
      var distCW  = (s[0] >= e[0]) ? 1 + e[0] - s[0] : e[0] - s[0];

      var hue = (distCW <= distCCW) ? s[0] + (distCW * p) : s[0] - (distCCW * p);
      if (hue < 0) hue += 1;

      var saturation = transition(value, maximum, s[1], e[1]); // s
      var value= transition(value, maximum, s[2], e[2]); // v
      return [hue, saturation, value];
    }

      init();

    }

})(jQuery);




