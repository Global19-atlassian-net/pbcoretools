"""Very independent little utilities.

These are easy to test and have zero side-effects.
They should depend on nothing put Python standard library modules.
"""


def split_filtStr(filtStr):
    """
    ...doctest:
        >>> split_filtStr('a>b;c<=d; e == f; bc=[0,1]')
        ['a>b', 'c<=d', 'e == f', 'bc=[0,1]']
        >>> split_filtStr('a>b AND c<=d AND  e == f')
        ['a>b', 'c<=d', 'e == f']
        >>> split_filtStr('bc=[0,1]')
        ['bc=[0,1]']
        >>> split_filtStr('zm = [3,4,5] AND length >= 1000')
        ['zm = [3,4,5]', 'length >= 1000']
    """
    if ';' in filtStr:
        return [s.strip() for s in filtStr.split(';')]
    elif ' AND ' in filtStr:
        return [s.strip() for s in filtStr.split(' AND ')]

    if ',' in filtStr and '[' not in filtStr:
        msg = "You're doing it wrong! You have ',' in the filter-string '{}', but not '['. That means you are probably trying to use a comma to separate conditions, which we do not support. Please use ' AND ' or ';' (semicolon) to separate conditions.".format(
            filtStr)
        raise ValueError(msg)
    return [filtStr]


def get_base_parser(description, log_level="INFO"):
    from pbcommand.cli import get_default_argparser_with_base_opts
    from pbcoretools import __VERSION__
    return get_default_argparser_with_base_opts(
        version=__VERSION__,
        description=description,
        default_level=log_level)
