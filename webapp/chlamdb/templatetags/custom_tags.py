from django import template

register = template.Library()


@register.filter
def keyvalue(dict, key):
    return dict.get(key)


@register.filter
def remove_spaces(value):
    return value.replace(' ', '+')


@register.filter
def return_item(lst, i):
    try:
        print(i)
        return lst[i]
    except Exception:
        return None


@register.filter
def slice_up_to(lst, i):
    return lst[:i]


def debug_mode(request):
    from django.conf import settings
    return {'debug_mode': settings.DEBUG}
