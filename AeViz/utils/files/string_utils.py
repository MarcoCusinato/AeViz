import re, numpy as np
from AeViz.units import u

def merge_strings(*args):
    """
    Marges the strings keeping in mind that they can have latex sintax
    in them
    """
    if len(args) == 0:
        return None
    if len(args) == 1:
        return args
    if any([ar is None for ar in args]):
        return None
    assert all([type(ar) == str for ar in args]), "Can only be concatenating strings"
    out_string = r''
    for ar in args:
        if out_string.endswith('$') and ar.startswith('$'):
            out_string = out_string[:-1] + ar[1:]
        else:
            out_string += ar
    return out_string

def apply_symbol(latex_str: str, symbol: str = "\\tilde"):
    if latex_str is None:
        return None
    # Check if the string starts and ends with $
    in_math_mode = latex_str.startswith('$') and latex_str.endswith('$')
    
    # Remove surrounding $ if present
    core_str = latex_str[1:-1] if in_math_mode else latex_str
    
    # Find the first letter or word before an underscore
    match = re.match(r"([^_]+)(.*)", core_str)
    
    if not match:
        return latex_str  # Return unchanged if no valid match
    
    first_part, rest = match.groups()
    
    # Apply the symbol
    modified = f"{symbol}{{{first_part}}}{rest}"
    
    # Restore math mode if needed
    return f"${modified}$" if in_math_mode else modified
        
import re

def split_number_and_unit(s):
    match = re.match(r'\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)[ ]*([a-zA-Zµ°/%]*)\s*$', s)
    if match:
        number = match.group(1)
        unit = match.group(2)
        return (float(number) * u.Unit(unit))
    else:
        raise TypeError(f'String {s} not supported')