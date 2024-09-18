#!/usr/bin/env python
import ete3
from ete3 import NodeStyle, StackedBarFace, TextFace, Tree, TreeStyle
from matplotlib.colors import rgb2hex

from lib import colors


class EteTree:

    DEFAULT_COLORS = ['#fc8d59', '#91bfdb', '#99d594',
                      '#c51b7d', '#f1a340', '#999999']
    RED = "#ff0000"
    BLUE = "#58ACFA"
    GREEN = "#99d594"
    ORANGE = "#fc8d59"

    def __init__(self, nwck_tree, **drawing_params):
        self.tree = nwck_tree
        self.params = drawing_params
        self.columns = []
        self.leaves_name = None

    def default_tree(nwck, **drawing_params):
        t = Tree(nwck)
        mid_point = t.get_midpoint_outgroup()
        if mid_point is not None:
            t.set_outgroup(mid_point)
        t.ladderize()
        return EteTree(t, **drawing_params)

    def pruned_tree(nwck, nodes, **drawing_params):
        t = Tree(nwck)
        t.prune([str(i) for i in nodes])
        mid_point = t.get_midpoint_outgroup()
        t.set_outgroup(mid_point)
        t.ladderize()
        return EteTree(t, **drawing_params)

    def add_column(self, face):
        self.columns.append(face)

    # Default style, may be overriden in a children class to change
    # the behaviour
    def get_style(self):
        tss = TreeStyle()
        tss.draw_guiding_lines = True
        tss.guiding_lines_color = "gray"
        tss.show_branch_support = False
        tss.show_leaf_name = False
        tss.margin_top = 5
        tss.margin_left = 5
        tss.margin_left = 5
        tss.margin_right = 5
        return tss

    # May be a good idea to be able to give custom parameters to the node names
    def get_leaf_name(self, index):
        if self.leaf_name_type is int:
            idx = int(index)
        elif self.leaf_name_type is str:
            idx = index
        else:
            raise Exception("Unsupported indexing type ", self.leaf_name_type)
        fgcolor = "black"
        if self.highlight_leaves is not None:
            if idx in self.highlight_leaves:
                fgcolor = "red"
        label = self.leaves_name.get(idx, self.default_val)
        t = TextFace(label, fgcolor=fgcolor, fsize=7)
        t.margin_right = 13
        return t

    def rename_leaves(self, hsh_names, default_val="-", leaf_name_type=int,
                      highlight_leaves=None):
        self.default_val = default_val
        if not isinstance(hsh_names, dict):
            raise Exception("Expects dict type for hsh_names")
        self.highlight_leaves = highlight_leaves
        self.leaves_name = hsh_names
        self.leaf_name_type = leaf_name_type

    def render(self, destination, **kwargs):
        for leaf in self.tree.iter_leaves():
            if self.leaves_name is not None:
                leaf.add_face(self.get_leaf_name(leaf.name), 0, "branch-right")

            for col_no, column in enumerate(self.columns):
                # Note: this assumes that only bioentries (integer)
                # are used as leaf names
                leaf.add_face(column.get_face(leaf.name), col_no, "aligned")

        tree_style = self.get_style()
        for col_no, column in enumerate(self.columns):
            header = column.get_header()
            if header is not None:
                tree_style.aligned_header.add_face(header, col_no)

        self.tree.render(destination, tree_style=tree_style, **kwargs)


class Column:
    def __init__(self, header=None, face_params=None, header_params=None):
        self.header = header
        self.header_params = header_params
        self.face_params = face_params

    def get_header(self):
        if self.header is None:
            return None
        face = TextFace(self.header)
        face.rotation = 270
        face.hz_align = 1
        face.vt_align = 1
        face.fsize = 7
        face.margin_right = 10
        # Put after the default values to erase any default value
        # in favor of a new one
        if self.header_params is not None:
            for param, value in self.header_params.items():
                setattr(face, param, value)
        return face

    def set_custom_header_params(self, header_params):
        self.header_params = header_params

    def set_default_params(self, text_face):
        text_face.margin_top = 2
        text_face.margin_right = 2
        text_face.margin_left = 2
        text_face.margin_bottom = 2
        text_face.hz_align = 1
        text_face.vt_align = 1
        text_face.border.width = 3
        text_face.border.color = "#ffffff"
        text_face.fsize = 7
        if self.face_params is not None:
            for name, value in self.face_params.items():
                setattr(text_face, name, value)


class SimpleColorColumn(Column):
    # should really be refactored.
    # Separation of concern sometimes broken...

    def __init__(self, values, header=None, use_col=True,
                 face_params=None, header_params=None, col_func=None,
                 default_val=0, default_val_is_num=False, color_gradient=False,
                 gradient_value_range=None, is_str_index=False):
        super().__init__(header, face_params, header_params)
        self.values = values
        self.header = header
        self.use_col = use_col
        self.default_val = default_val
        self.default_val_is_num = default_val_is_num
        self.is_str_index = is_str_index
        if face_params is None or "color" not in face_params:
            self.col = EteTree.BLUE
        else:
            self.col = face_params["color"]
        self.col_func = col_func
        self.color_gradient = color_gradient
        if color_gradient:
            if gradient_value_range is not None:
                self.cm, _ = colors.get_continuous_scale(gradient_value_range)
            else:
                my_values = list(values.values())
                if default_val_is_num:
                    my_values.append(default_val)
                self.cm, _ = colors.get_continuous_scale(my_values)

    def fromSeries(series, header=None, cls=None, **args):
        values = series.to_dict()
        if cls is None:
            return SimpleColorColumn(values, header=header, **args)
        else:
            return cls(values, header=header, **args)

    def get_face(self, index):
        if not self.is_str_index:
            index = int(index)
        val = self.values.get(index, self.default_val)

        italic = "normal"
        if self.face_params is not None:
            if self.face_params.get("italic", False):
                italic = "ITalic"

        text_face = TextFace(str(val), fstyle=italic)
        if (val == self.default_val and self.default_val_is_num) or \
                (val != self.default_val and self.color_gradient):
            rgba = self.cm.to_rgba(val, bytes=True)
            text_face.inner_background.color = colors.to_rgb_str(rgba)
            luminance = colors.get_luminance(rgba)
            if luminance >= .5:
                text_face.fgcolor = "#000000"
            else:
                text_face.fgcolor = "#ffffff"
        elif self.use_col and val != 0 and index in self.values:
            if self.col_func is None:
                text_face.inner_background.color = self.col
            else:
                text_face.inner_background.color = self.col_func(index)

        self.set_default_params(text_face)
        return text_face


class ValueColoredColumn(Column):
    """Gets the color by applying col_func to the value.
    """
    def __init__(self, values, col_func, header=None,
                 face_params=None, header_params=None):
        self.values = values
        self.col_func = col_func
        self.header = header
        self.face_params = face_params
        self.header_params = header_params

    def get_color(self, index, val):
        return self.col_func(val)

    def get_face(self, index):
        index = int(index)
        val = self.values.get(index)
        text_face = TextFace(str(val))
        text_face.inner_background.color = self.get_color(index, val)
        self.set_default_params(text_face)
        return text_face


class MatchingColorColumn(ValueColoredColumn):
    """Will color a face only if its value matches the value in
    to_match. In that case it gets the color by applying col_func
    to the value.
    """
    def __init__(self, values, to_match, col_func, header=None,
                 face_params=None, header_params=None):
        self.to_match = to_match
        super(MatchingColorColumn, self).__init__(
            values, col_func, header=header, face_params=face_params,
            header_params=header_params)

    def get_color(self, index, val):
        if val == self.to_match.get(index):
            return self.col_func(val)


class ModuleCompletenessColumn(Column):
    """
    Straightforward class:
    * the values contains the number of missing KO for a complete module
    * the text_face are colored according to the number of missing KOs:
    *  - none missing: green
    *  - missing KOs: orange
    """
    def __init__(self, values, header=None, add_missing=True):
        super().__init__(header)
        self.values = values
        self.add_missing = add_missing

    def get_face(self, index):
        index = int(index)
        val = self.values.get(index, 0)

        if self.add_missing:
            text_face = TextFace(val)
        elif val == 0:
            text_face = TextFace("C")
        elif val >= 1:
            text_face = TextFace("I")

        if val == 0:
            text_face.inner_background.color = EteTree.GREEN
        else:
            text_face.inner_background.color = EteTree.ORANGE
        self.set_default_params(text_face)
        return text_face


class KOAndCompleteness(Column):
    """
    This class should be used when displaying the number of KO
    in a module. The text faces are then colored according to the
    number of missing KOs in a module.
    If the module is complete (0 missing): text face is green
    If the module lacks some KOs: text face is orange
    """

    def __init__(self, n_kos, n_missing_kos, module):
        """
        The values should have two columns:
        * one for the number of ko in the module
        * one for the number of missing kos
        """
        super().__init__(header=module)
        self.values = n_kos
        self.n_missing = n_missing_kos

    def get_face(self, index):
        index = int(index)
        val = self.values.get(index, "-")
        if val == 0:
            val = "-"
        n_missing = self.n_missing.get(index, None)
        text_face = TextFace(val)

        if n_missing is not None and val != "-":
            color = EteTree.GREEN if n_missing == 0 else EteTree.ORANGE
            text_face.inner_background.color = color

        super().set_default_params(text_face)
        return text_face


class ReferenceColumn(Column):
    def __init__(self, values):
        pass

    def get_face(self, index):
        pass


class EteTool():

    '''
    Plot ete3 phylogenetic profiles.

    - self.add_simple_barplot: add a barplot face from taxon2value dictionnary
    - self.add_text_face: add text face
    - self.add_heatmap: add column with cells with value + colored background
    - self.rename_leaves: rename tree leaves from a dictionnary (old_name2new_name)
    '''

    def __init__(self,
                 tree_file):

        self.column_count = 0

        self.default_colors = ['#fc8d59', '#91bfdb', '#99d594', '#c51b7d',
                               '#f1a340', '#999999']

        self.color_index = 0

        self.rotate = False

        # if not tree instance, considfer it as a path or a newick string
        print("TREE TYOE:", type(tree_file))
        if isinstance(tree_file, Tree):
            self.tree = tree_file
        elif isinstance(tree_file, ete3.phylo.phylotree.PhyloNode):
            self.tree = tree_file
        else:
            self.tree = Tree(tree_file)
        # Calculate the midpoint node
        R = self.tree.get_midpoint_outgroup()
        # and set it as tree outgroup
        try:
            self.tree.set_outgroup(R)
        except Exception:
            pass

        self.tree.ladderize()

        self.tss = TreeStyle()
        self.tss.draw_guiding_lines = True
        self.tss.guiding_lines_color = "gray"
        self.tss.show_leaf_name = False

    def add_stacked_barplot(self, taxon2value_list, header_name,
                            color_list=False):

        pass

    def rename_leaves(self, taxon2new_taxon, keep_original=False,
                      add_face=True):
        for i, lf in enumerate(self.tree.iter_leaves()):
            if not keep_original:
                if lf.name in taxon2new_taxon:
                    label = taxon2new_taxon[lf.name]
                else:
                    label = 'n/a'
            else:
                if lf.name in taxon2new_taxon:
                    label = '%s (%s)' % (taxon2new_taxon[lf.name], lf.name)
                else:
                    label = 'n/a'
            if add_face:
                n = TextFace(label, fgcolor="black", fsize=12, fstyle='italic')
                lf.add_face(n, 0)
            lf.name = label

    def add_heatmap(self, taxon2value, header_name, continuous_scale=False,
                    show_text=False):

        from lib.colors import get_continuous_scale

        self._add_header(header_name)

        if continuous_scale:
            color_scale = get_continuous_scale(taxon2value.values())

        for i, lf in enumerate(self.tree.iter_leaves()):

            if lf.name not in taxon2value:
                n = TextFace('')
            else:
                value = taxon2value[lf.name]

                if show_text:
                    n = TextFace('%s' % value)
                else:
                    n = TextFace('    ')

                n.margin_top = 2
                n.margin_right = 3
                n.margin_left = 3
                n.margin_bottom = 2
                n.hz_align = 1
                n.vt_align = 1
                n.border.width = 3
                n.border.color = "#ffffff"
                if continuous_scale:
                    n.background.color = rgb2hex(color_scale[0].to_rgba(float(value)))
                n.opacity = 1.
                i += 1

            if self.rotate:
                n.rotation = 270
            lf.add_face(n, self.column_count, position="aligned")

        self.column_count += 1

    def _add_header(self, header_name, column_add=0):

        n = TextFace(f'{header_name}')
        n.margin_top = 1
        n.margin_right = 1
        n.margin_left = 20
        n.margin_bottom = 1
        n.hz_align = 2
        n.vt_align = 2
        n.rotation = 270
        n.inner_background.color = "white"
        n.opacity = 1.
        # add header
        self.tss.aligned_header.add_face(n, self.column_count-1+column_add)

    def _get_default_barplot_color(self,):

        col = self.default_colors[self.color_index]

        if self.color_index == 5:
            self.color_index = 0
        else:
            self.color_index += 1

        return col

    def add_simple_barplot(self, taxon2value, header_name, color=False,
                           show_values=False, substract_min=False,
                           highlight_cutoff=False, highlight_reverse=False,
                           max_value=False):

        if not show_values:
            self._add_header(header_name, column_add=0)
        else:
            self._add_header(header_name, column_add=1)

        values_lists = [float(i) for i in taxon2value.values()]

        min_value = min(values_lists)

        if substract_min:
            values_lists = [i-min_value for i in values_lists]
            for taxon in list(taxon2value.keys()):
                taxon2value[taxon] = taxon2value[taxon]-min_value

        if not color:
            color = self._get_default_barplot_color()

        for i, lf in enumerate(self.tree.iter_leaves()):

            try:
                value = taxon2value[lf.name]
            except KeyError:
                value = 0

            if show_values:
                barplot_column = 1
                if substract_min:
                    real_value = value + min_value
                else:
                    real_value = value
                if isinstance(real_value, float):
                    a = TextFace(" %s " % str(round(real_value,2)))
                else:
                    a = TextFace(" %s " % str(real_value))
                a.margin_top = 1
                a.margin_right = 2
                a.margin_left = 5
                a.margin_bottom = 1
                if self.rotate:
                    a.rotation = 270
                lf.add_face(a, self.column_count, position="aligned")
            else:
                barplot_column = 0
            if not max_value:
                fraction_biggest = (float(value)/max(values_lists))*100
            else:
                fraction_biggest = (float(value)/max_value)*100
            fraction_rest = 100-fraction_biggest

            if highlight_cutoff:
                if substract_min:
                    real_value = value + min_value
                else:
                    real_value = value
                if highlight_reverse:
                    if real_value > highlight_cutoff:
                        lcolor = "grey"
                    else:
                        lcolor = color
                else:
                    if real_value < highlight_cutoff:
                        lcolor = "grey"
                    else:
                        lcolor = color
            else:
                lcolor = color

            b = StackedBarFace([fraction_biggest, fraction_rest], width=100,
                               height=15, colors=[lcolor, 'white'])
            b.rotation = 0
            b.inner_border.color = "grey"
            b.inner_border.width = 0
            b.margin_right = 15
            b.margin_left = 0
            if self.rotate:
                b.rotation = 270
            lf.add_face(b, self.column_count + barplot_column, position="aligned")

        self.column_count += (1 + barplot_column)

    def add_barplot_counts(self,):
         # todo
        pass

    def remove_dots(self,):

        nstyle = NodeStyle()
        nstyle["shape"] = "sphere"
        nstyle["size"] = 0
        nstyle["fgcolor"] = "darkred"

        # Applies the same static style to all nodes in the tree. Note that,
        # if "nstyle" is modified, changes will affect to all nodes
        for n in self.tree.traverse():
            n.set_style(nstyle)

    def add_text_face(self, taxon2text, header_name, color_scale=False):

        from lib.colors import get_categorical_color_scale

        if color_scale:
            value2color = get_categorical_color_scale(taxon2text.values())

        self._add_header(header_name)

        # add column
        for i, lf in enumerate(self.tree.iter_leaves()):
            if lf.name in taxon2text:
                n = TextFace('%s' % taxon2text[lf.name])
                if color_scale:
                    n.background.color = value2color[taxon2text[lf.name]]
            else:
                print(lf.name, "not in", taxon2text)
                n = TextFace('-')
            n.margin_top = 1
            n.margin_right = 10
            n.margin_left = 10
            n.margin_bottom = 1
            n.opacity = 1.
            if self.rotate:
                n.rotation = 270
            lf.add_face(n, self.column_count, position="aligned")

        self.column_count += 1


class EteToolCompact():

    '''
    Plot ete3 phylogenetic profiles.

    - self.add_simple_barplot: add a barplot face from taxon2value dictionnary
    - self.add_heatmap: add column with cells with value + colored background
    - self.rename_leaves: rename tree leaves from a dictionnary (old_name2new_name)
    - self.add_categorical_colorscale_legend: add legend
    - self.add_continuous_colorscale_legend: add legend
    '''

    def __init__(self, tree_file):
        self.column_count = 0

        self.rotate = False

        self.tree = Tree(tree_file)

        self.tree_length = len([i for i in self.tree.iter_leaves()])

        self.text_scale = (self.tree_length)*0.01  # math.log2

        self.default_colors = ['#fc8d59', '#91bfdb', '#99d594', '#c51b7d',
                               '#f1a340', '#999999']

        self.color_index = 0

        # Calculate the midpoint node
        R = self.tree.get_midpoint_outgroup()
        # and set it as tree outgroup
        self.tree.set_outgroup(R)

        self.tss = TreeStyle()
        self.tss.draw_guiding_lines = True
        self.tss.guiding_lines_color = "gray"
        self.tss.show_leaf_name = False
        self.tss.branch_vertical_margin = 0

    def _get_default_barplot_color(self,):
        col = self.default_colors[self.color_index]

        if self.color_index == 5:
            self.color_index = 0
        else:
            self.color_index += 1

        return col

    def _add_header(self, header_name, column_add=0):
        n = TextFace(f'{header_name}')
        n.margin_top = 1
        n.margin_right = 1
        n.margin_left = 20
        n.margin_bottom = 1
        n.hz_align = 2
        n.vt_align = 2
        n.rotation = 270
        n.inner_background.color = "white"
        n.opacity = 1.
        # add header
        self.tss.aligned_header.add_face(n, self.column_count-1+column_add)

    def rename_leaves(self, taxon2new_taxon):
        for i, lf in enumerate(self.tree.iter_leaves()):
            n = TextFace(taxon2new_taxon[lf.name], fgcolor="black", fsize=12,
                         fstyle='italic')
            lf.add_face(n, 0)

    def add_continuous_colorscale_legend(self, title, min_val, max_val, scale):
        self.tss.legend.add_face(
            TextFace(f"{title}", fsize=4 * self.text_scale), column=0)

        if min_val != max_val:
            n = TextFace(" " * int(self.text_scale), fsize=4 * self.text_scale)
            n.margin_top = 1
            n.margin_right = 1
            n.margin_left = 10
            n.margin_bottom = 1
            n.inner_background.color = rgb2hex(scale[0].to_rgba(float(max_val)))

            n2 = TextFace(" " * int(self.text_scale), fsize=4 * self.text_scale)
            n2.margin_top = 1
            n2.margin_right = 1
            n2.margin_left = 10
            n2.margin_bottom = 1
            n2.inner_background.color = rgb2hex(scale[0].to_rgba(float(min_val)))

            self.tss.legend.add_face(n, column=1)
            self.tss.legend.add_face(TextFace(
                f"{max_val} % (max)", fsize=4 * self.text_scale), column=2)
            self.tss.legend.add_face(n2, column=1)
            self.tss.legend.add_face(TextFace(
                f"{min_val} % (min)", fsize=4 * self.text_scale), column=2)
        else:
            n2 = TextFace(" " * int(self.text_scale),
                          fsize=4 * self.text_scale)
            n2.margin_top = 1
            n2.margin_right = 1
            n2.margin_left = 10
            n2.margin_bottom = 1
            n2.inner_background.color = rgb2hex(scale[0].to_rgba(float(min_val)))

            self.tss.legend.add_face(n2, column=0)
            self.tss.legend.add_face(
                TextFace(f"{max_val} % Id", fsize=4 * self.text_scale),
                column=1)

    def add_categorical_colorscale_legend(self, title, scale):

        self.tss.legend.add_face(
            TextFace(f"{title}", fsize=4 * self.text_scale), column=0)

        col = 1
        for n, value in enumerate(scale):

            n2 = TextFace(" " * int(self.text_scale), fsize=4 * self.text_scale)
            n2.margin_top = 1
            n2.margin_right = 1
            n2.margin_left = 10
            n2.margin_bottom = 1
            n2.inner_background.color = scale[value]

            self.tss.legend.add_face(n2, column=col)
            self.tss.legend.add_face(
                TextFace(f"{value}", fsize=4 * self.text_scale), column=col+1)

            col += 2
            if col > 16:
                self.tss.legend.add_face(
                    TextFace("    ", fsize=4 * self.text_scale),
                    column=0)
                col = 1

    def add_simple_barplot(self, taxon2value, header_name, color=False,
                           show_values=False, substract_min=False,
                           max_value=False):
        if not show_values:
            self._add_header(header_name, column_add=0)
        else:
            self._add_header(header_name, column_add=1)

        values_lists = [float(i) for i in taxon2value.values()]

        min_value = min(values_lists)

        if substract_min:
            values_lists = [i-min_value for i in values_lists]
            for taxon in list(taxon2value.keys()):
                taxon2value[taxon] = taxon2value[taxon]-min_value

        if not color:
            color = self._get_default_barplot_color()

        for i, lf in enumerate(self.tree.iter_leaves()):

            try:
                value = taxon2value[lf.name]
            except Exception:
                value = 0

            if show_values:
                barplot_column = 1
                if isinstance(value, float):
                    a = TextFace(" %s " % str(round(value, 2)))
                else:
                    a = TextFace(" %s " % str(value))
                a.margin_top = 1
                a.margin_right = 2
                a.margin_left = 5
                a.margin_bottom = 1
                if self.rotate:
                    a.rotation = 270
                lf.add_face(a, self.column_count, position="aligned")
            else:
                barplot_column = 0
            if not max_value:
                fraction_biggest = (float(value)/max(values_lists))*100
            else:
                fraction_biggest = (float(value)/max_value)*100
            fraction_rest = 100-fraction_biggest

            b = StackedBarFace([fraction_biggest, fraction_rest],
                               width=100 * (self.text_scale/3),
                               height=18,
                               colors=[color, 'white'])
            b.rotation = 0
            b.margin_right = 10
            b.margin_left = 10
            b.hz_align = 2
            b.vt_align = 2
            b.rotable = False
            if self.rotate:
                b.rotation = 270
            lf.add_face(b, self.column_count + barplot_column, position="aligned")

        self.column_count += (1 + barplot_column)

    def add_heatmap(self, taxon2value, header_name, scale_type="continuous",
                    palette=False):
        from lib.colors import (get_categorical_color_scale,
                                get_continuous_scale)

        if scale_type == "continuous":
            scale = get_continuous_scale(taxon2value.values())
            self.add_continuous_colorscale_legend("Closest hit identity",
                                                  min(taxon2value.values()),
                                                  max(taxon2value.values()),
                                                  scale)
        elif scale_type == "categorical":
            scale = get_categorical_color_scale(taxon2value.values())
            self.add_categorical_colorscale_legend("MLST",
                                                   scale)
        else:
            raise IOError("unknown type")

        for i, lf in enumerate(self.tree.iter_leaves()):
            n = TextFace("   " * int(self.text_scale))
            if lf.name in taxon2value:
                value = taxon2value[lf.name]
                n = TextFace("   " * int(self.text_scale))
                if scale_type == "categorical":
                    n.inner_background.color = scale[value]
                if scale_type == "continuous":
                    n.inner_background.color = rgb2hex(scale[0].to_rgba(float(value)))

            n.margin_top = 0
            n.margin_right = 0
            n.margin_left = 10
            n.margin_bottom = 0
            n.opacity = 1.
            if self.rotate:
                n.rotation = 270
            lf.add_face(n, self.column_count, position="aligned")

        self.column_count += 1

    def remove_labels(self,):
        for i, lf in enumerate(self.tree.iter_leaves()):
            n = TextFace("")
            lf.add_face(n, 0)


def get_newick(node, newick, parentdist, leaf_names):
    '''
    convert hierarchical clustering to newick format
    '''
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), newick, node.dist, leaf_names)
        newick = get_newick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick
