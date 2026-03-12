import sys
import builtins
import atexit

imported_modules = set()

original_import = builtins.__import__

def tracking_import(name, globals=None, locals=None, fromlist=(), level=0):
    imported_modules.add(name.split('.')[0])
    return original_import(name, globals, locals, fromlist, level)

builtins.__import__ = tracking_import


@atexit.register
def save_imports():
    with open("used_modules.txt", "w") as f:
        for m in sorted(imported_modules):
            f.write(m + "\n")
