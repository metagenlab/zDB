

import matplotlib as mpl
import matplotlib.cm as cm


def to_rgb_str(rgba):
    """
    Encodes a rgba tuple in the form of a #rrggbb string in hexadecimal
    notation. This is mostly used to convert the rgba tuple from matplotlib
    colormap to ete tree drawing engine. As ete3 requires the precise #rrggbb
    format, add an extra 0 when number are lower than 16.
    """

    string_vals = []
    for val in rgba[0:3]:
        if val==0:
            string_vals.append("00")
        elif val<=15:
            string_vals.append("0"+hex(val).lstrip("0x"))
        else:
            string_vals.append(hex(val).lstrip("0x"))
    return "#"+"".join(string_vals)


def get_luminance(rgba):
    return (0.299*rgba[0]+0.587*rgba[1]+0.114*rgba[2])/255


def get_spaced_colors(n):
    
    '''
    Return n spaces colors color
    '''

    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]

    return ['#%02x%02x%02x' % (int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors]


def get_categorical_color_scale(value_list):

    '''
    Given input variable list,
    return a dictionnary of value2color
    '''
    
    value2col = dict(zip(value_list, get_spaced_colors(len(value_list))))
    
    return value2col


def get_continuous_scale(value_list, 
                         cm_scale="OrRd"):
        '''
        Generate color scale best of a list of value using maplotlib
        Return the scale and the max value.
        Example of scales:
        - OrRd
        - YlGnBu
        - PuBu
        '''
        
        cm_scale = getattr(cm, cm_scale)

        if min(value_list) == max(value_list):
            min_val = 0
            max_val = 1.5*max(value_list)
        else:
            min_val = min(value_list)
            max_val = max(value_list)

        norm = mpl.colors.Normalize(vmin=min_val*0.9,vmax=max_val) # *1.1
        m = cm.ScalarMappable(norm=norm, cmap=cm_scale)
        return [m, float(max_val)]
