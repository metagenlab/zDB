from django import template

register = template.Library()


@register.filter
def keyvalue(dict, key):
    return dict.get(key)


@register.filter
def remove_spaces(value):
    return value.replace(' ', '+')


def debug_mode(request):
    from django.conf import settings
    return {'debug_mode': settings.DEBUG}
