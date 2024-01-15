def safe_replace(string, search_string, replace_string):
    if string:
        return string.replace(search_string, replace_string)
    return string
