from dal.autocomplete import Select2ListView

from views.utils import AccessionFieldHandler


class AutocompleteTaxid(Select2ListView):

    def get_list(self):
        with_plasmids = self.forwarded["include_plasmids"]
        to_exclude = self.forwarded["to_exclude"]
        return AccessionFieldHandler().get_choices(with_plasmids=with_plasmids,
                                                   to_exclude=to_exclude)


class AutocompleteNMissing(Select2ListView):

    def get_list(self):
        n_max = len(self.forwarded["included"])
        choices = [(i, i) for i in range(n_max)]
        return choices
