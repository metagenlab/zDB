from urllib.parse import quote

from views.utils import format_amr
from views.utils import format_cog
from views.utils import format_genomic_island_cluster
from views.utils import format_ko
from views.utils import format_orthogroup
from views.utils import format_pfam
from views.utils import missing_mandatory
from views.utils import optional2status


class BaseObjectMetadata:
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

    @property
    def index_comp_url(self):
        return "/index_comp/{}".format(self.object_type)


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
        """Will point to the details page as soon as it exists"""
        if to_url:
            return f'<a href="/fam_vf/{entry}">{entry}</a>'
        return entry


class GiMetadata(BaseObjectMetadata):
    object_type = "gic"
    object_name = "Genomic island cluster"

    @property
    @staticmethod
    def is_enabled(self):
        return optional2status.get("gi", False)

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_genomic_island_cluster(entry, to_url=to_url)


class ModuleMetadata(BaseObjectMetadata):
    object_type = "module"
    object_name = "KEGG Module"

    @staticmethod
    def format_entry(entry, to_url=False):
        """Will point to the details page as soon as it exists"""
        if to_url:
            return f"<a href=/KEGG_module_map/{entry}>{entry}</a>"
        return entry

    @property
    def index_comp_url(self):
        return None


class PathwayMetadata(BaseObjectMetadata):
    object_type = "pathway"
    object_name = "KEGG Pathway"

    @staticmethod
    def format_entry(entry, to_url=False):
        """Will point to the details page as soon as it exists"""
        if to_url:
            return f"<a href=/KEGG_mapp_ko/{entry}>{entry}</a>"
        return entry

    @property
    def index_comp_url(self):
        return None


class GroupMetadata(BaseObjectMetadata):
    object_type = "group"
    object_name = "Genome Group"
    overview_description = "Overview of defined groups of genomes."

    @staticmethod
    def format_entry(entry, to_url=False):
        if to_url:
            return f"<a href=/groups/{quote(entry)}>{entry}</a>"
        return entry


class MetadataGetter:
    metadata_classes = [
        AmrMetadata,
        CogMetadata,
        KoMetadata,
        PfamMetadata,
        OrthogroupMetadata,
        VfMetadata,
        GiMetadata,
        ModuleMetadata,
        PathwayMetadata,
    ]

    _annotations = ["cog", "pfam", "ko", "amr", "vf", "gic"]
    _orthology = ["orthogroup"]

    def __init__(self):
        self.metadata_instances = [cls() for cls in self.metadata_classes]
        self.object_type_to_metadata = {
            obj.object_type: obj for obj in self.metadata_instances
        }

    def get_annotations_metadata(self):
        return (
            self.object_type_to_metadata[object_type]
            for object_type in self._annotations
            if self.object_type_to_metadata[object_type].is_enabled
        )

    def get_orthology_metadata(self):
        return (
            self.object_type_to_metadata[object_type] for object_type in self._orthology
        )

    @property
    def group_metadata(self):
        return GroupMetadata()


def my_locals(local_dico):
    local_dico["optional2status"] = optional2status
    local_dico["missing_mandatory"] = missing_mandatory
    local_dico["metadata"] = MetadataGetter()
    return local_dico
