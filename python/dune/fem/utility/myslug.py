import re
# taken form 'myst_parser/mdit_to_docutils/base.py' but removing the 'lower case'
_SLUGIFY_CLEAN_REGEX = re.compile(r"[^\w\u4e00-\u9fff\- ]")
def heading_slug_func(title: str) -> str:
    """Default slugify function.
    This aims to mimic the GitHub Markdown format - but without the lower case
    """
    return _SLUGIFY_CLEAN_REGEX.sub("", title.replace(" ", "-"))
