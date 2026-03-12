import sys
import pkgutil
import importlib
from importlib.metadata import packages_distributions, version


def load_all_paircars():

    import paircars.utils as u
    import paircars.pipeline as p
    import paircars.clusterutils as c

    for pkg in [u, p, c]:
        for module in pkgutil.walk_packages(pkg.__path__, pkg.__name__ + "."):
            try:
                importlib.import_module(module.name)
            except Exception:
                pass


def get_runtime_packages():

    loaded_modules = set(sys.modules.keys())

    mapping = packages_distributions()

    required = {}

    for mod in loaded_modules:

        root = mod.split(".")[0]

        if root in mapping:

            for pkg in mapping[root]:

                try:
                    required[pkg] = version(pkg)
                except Exception:
                    pass

    return required


def main():

    print("Loading paircars modules...\n")

    load_all_paircars()

    print("Detecting runtime dependencies...\n")

    deps = get_runtime_packages()

    sorted_deps = sorted(deps.items())

    for pkg, ver in sorted_deps:
        print(f"{pkg}=={ver}")

    with open("runtime_requirements.txt", "w") as f:
        for pkg, ver in sorted_deps:
            f.write(f"{pkg}=={ver}\n")

    print("\nSaved runtime_requirements.txt")


if __name__ == "__main__":
    main()
