import pkgutil
import importlib
import paircars.utils as paircars

missing = []

for module in pkgutil.walk_packages(paircars.__path__, paircars.__name__ + "."):
    name = module.name
    try:
        importlib.import_module(name)
        print(f"OK: {name}")
    except Exception as e:
        print(f"FAIL: {name} -> {e}")
        missing.append((name, str(e)))

print("\nMissing or failed modules:")
for m in missing:
    print(m)
    
import paircars.pipeline as paircars

missing = []

for module in pkgutil.walk_packages(paircars.__path__, paircars.__name__ + "."):
    name = module.name
    try:
        importlib.import_module(name)
        print(f"OK: {name}")
    except Exception as e:
        print(f"FAIL: {name} -> {e}")
        missing.append((name, str(e)))

print("\nMissing or failed modules:")
for m in missing:
    print(m)
    
import paircars.clusterutils as paircars

missing = []

for module in pkgutil.walk_packages(paircars.__path__, paircars.__name__ + "."):
    name = module.name
    try:
        importlib.import_module(name)
        print(f"OK: {name}")
    except Exception as e:
        print(f"FAIL: {name} -> {e}")
        missing.append((name, str(e)))

print("\nMissing or failed modules:")
for m in missing:
    print(m)
