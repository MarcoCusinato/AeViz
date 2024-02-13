import sys

"""
Check if the module is imported from the terminal or not.
"""
def is_imported_from_terminal():
    return hasattr(sys, 'ps1')  # Check if the Python terminal prompt is
                                # available

global TERMINAL 
if is_imported_from_terminal():
    TERMINAL = True
else:
    TERMINAL = False