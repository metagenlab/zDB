from ete4 import Tree
from ete4.smartview.faces import TextFace
from ete4.smartview.layout import Layout
from ete4.treeview import NodeStyle
from ete4.treeview import TreeStyle
from ete4.treeview.faces import StackedBarFace as treeview_StackedBarFace
from ete4.treeview.faces import TextFace as treeview_TextFace
from lib import colors
from matplotlib.colors import rgb2hex


class EteTree:
    DEFAULT_COLORS = ["#fc8d59", "#91bfdb", "#99d594", "#c51b7d", "#f1a340", "#999999"]
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
        if mid_point is not None and mid_point is not t.root:
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
        t = treeview_TextFace(label, fgcolor=fgcolor, fsize=7)
        t.margin_right = 13
        return t

    def draw_leaf_name(self, leaf):
        index = leaf.name
        if self.leaf_name_type is int:
            idx = int(index)
        elif self.leaf_name_type is str:
            idx = index
        else:
            raise Exception("Unsupported indexing type ", self.leaf_name_type)
        fgcolor = "black"
        if self.highlight_leaves is not None and idx in self.highlight_leaves:
            fgcolor = "red"
        label = self.leaves_name.get(idx, self.default_val)
        style = {"fill": fgcolor}
        t = TextFace(label, style=style)
        return [t]

    def rename_leaves(
        self, hsh_names, default_val="-", leaf_name_type=int, highlight_leaves=None
    ):
        self.default_val = default_val
        if not isinstance(hsh_names, dict):
            raise Exception("Expects dict type for hsh_names")
        self.highlight_leaves = highlight_leaves
        self.leaves_name = hsh_names
        self.leaf_name_type = leaf_name_type

    def render(self, destination, **kwargs):
        for leaf in self.tree.leaves():
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

    def get_layouts(self):
        layouts = []
        if self.draw_leaf_names:
            layouts.append(Layout("labels", draw_node=self.draw_leaf_name))
        return layouts


class Column:
    def __init__(self, header=None, face_params=None, header_params=None):
        self.header = header
        self.header_params = header_params
        self.face_params = face_params

    def get_header(self):
        if self.header is None:
            return None
        face = treeview_TextFace(self.header)
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

    def __init__(
        self,
        values,
        header=None,
        use_col=True,
        face_params=None,
        header_params=None,
        col_func=None,
        default_val=0,
        default_val_is_num=False,
        color_gradient=False,
        gradient_value_range=None,
        is_str_index=False,
    ):
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

        text_face = treeview_TextFace(str(val), fstyle=italic)
        if (val == self.default_val and self.default_val_is_num) or (
            val != self.default_val and self.color_gradient
        ):
            rgba = self.cm.to_rgba(val, bytes=True)
            text_face.inner_background.color = colors.to_rgb_str(rgba)
            luminance = colors.get_luminance(rgba)
            if luminance >= 0.5:
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
    """Gets the color by applying col_func to the value."""

    def __init__(
        self, values, col_func, header=None, face_params=None, header_params=None
    ):
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
        text_face = treeview_TextFace(str(val))
        text_face.inner_background.color = self.get_color(index, val)
        self.set_default_params(text_face)
        return text_face


class MatchingColorColumn(ValueColoredColumn):
    """Will color a face only if its value matches the value in
    to_match. In that case it gets the color by applying col_func
    to the value.
    """

    def __init__(
        self,
        values,
        to_match,
        col_func,
        header=None,
        face_params=None,
        header_params=None,
    ):
        self.to_match = to_match
        super().__init__(
            values,
            col_func,
            header=header,
            face_params=face_params,
            header_params=header_params,
        )

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
            text_face = treeview_TextFace(val)
        elif val == 0:
            text_face = treeview_TextFace("C")
        elif val >= 1:
            text_face = treeview_TextFace("I")

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
        text_face = treeview_TextFace(val)

        if n_missing is not None and val != "-":
            color = EteTree.GREEN if n_missing == 0 else EteTree.ORANGE
            text_face.inner_background.color = color

        super().set_default_params(text_face)
        return text_face
