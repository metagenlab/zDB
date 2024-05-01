from dal.autocomplete import Select2ListView

from views.utils import AccessionFieldHandler


class AutocompleteTaxid(Select2ListView):

    def get_list(self):
        with_plasmids = self.forwarded["include_plasmids"]
        exclude = self.forwarded["exclude"]
        exclude_taxids_in_groups = self.forwarded["exclude_taxids_in_groups"]
        return AccessionFieldHandler().get_choices(
            with_plasmids=with_plasmids,
            exclude=exclude,
            exclude_taxids_in_groups=exclude_taxids_in_groups)


class AutocompleteNMissing(Select2ListView):

    def get_list(self):
        taxids, plasmids = AccessionFieldHandler().extract_choices(
            self.forwarded["included"], self.forwarded["include_plasmids"])
        n_max = len(taxids) + len(plasmids)
        choices = [(i, i) for i in range(n_max)]
        return choices
