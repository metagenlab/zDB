from views.utils import (format_amr, format_cog, format_ko, format_orthogroup,
                         format_pfam, optional2status)


class BaseObjectMetadata():

    @property
    @staticmethod
    def is_enabled(self):
        return optional2status.get(self.object_type, False)

    @property
    def object_name_plural(self):
        return f"{self.object_name}s"

    @property
    def object_name_singular_or_plural(self):
        return f"{self.object_name}(s)"


class AmrMetadata(BaseObjectMetadata):

    object_type = "amr"
    object_name = "AMR gene"

    @staticmethod
    def format_entry(gene, to_url=False):
        return format_amr(gene, to_url=to_url)


class CogMetadata(BaseObjectMetadata):

    object_type = "cog"
    object_name = "COG entry"
    object_name_plural = "COG entries"
    object_name_singular_or_plural = "COG entry(ies)"

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_cog(entry, as_url=to_url)


class KoMetadata(BaseObjectMetadata):

    object_type = "ko"
    object_name = "Kegg Ortholog"

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_ko(entry, as_url=to_url)


class PfamMetadata(BaseObjectMetadata):

    object_type = "pfam"
    object_name = "Pfam domain"

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_pfam(entry, to_url=to_url)


class OrthogroupMetadata(BaseObjectMetadata):

    object_type = "orthogroup"
    object_name = "Orthologous group"

    is_enabled = True

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_orthogroup(entry, to_url=to_url)


class VfMetadata(BaseObjectMetadata):

    object_type = "vf"
    object_name = "Virulence factor"

    @staticmethod
    def format_entry(entry, to_url=False):
        """Will point to the details page as soon as it exists
        """
        if to_url:
            return f"<a href=\"/fam_vf/{entry}\">{entry}</a>"
        return entry
