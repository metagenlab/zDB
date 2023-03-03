#from models import Genome
#from models import Database
#from django import forms

"""
import autocomplete_light


class GenomeAutocomplete(autocomplete_light.AutocompleteModelBase):
    search_fields = ['^genome', 'database']
autocomplete_light.register(Genome, GenomeAutocomplete)


class DBForm(forms.ModelForm):
    person = forms.ModelChoiceField(Genome.objects.all(),
        widget=autocomplete_light.ChoiceWidget('GenomeAutocomplete'))

    class Meta:
        model = Database
"""