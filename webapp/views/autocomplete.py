from dal.autocomplete import Select2ListView
from views.utils import AccessionFieldHandler


class AutocompleteTaxid(Select2ListView):
    def get_list(self):
        with_plasmids = self.forwarded.get("include_plasmids", False)
        exclude = self.forwarded.get("exclude", [])
        exclude_taxids_in_groups = self.forwarded.get("exclude_taxids_in_groups", [])
        return AccessionFieldHandler().get_choices(
            with_plasmids=with_plasmids,
            exclude=exclude,
            exclude_taxids_in_groups=exclude_taxids_in_groups,
        )


class AutocompleteNMissing(Select2ListView):
    def get_list(self):
        taxids, plasmids = AccessionFieldHandler().extract_choices(
            self.forwarded.get("included", []),
            self.forwarded.get("include_plasmids", False),
        )
        n_max = len(taxids)
        if plasmids:
            n_max += len(plasmids)
        choices = [(i, i) for i in range(n_max)]
        return choices
