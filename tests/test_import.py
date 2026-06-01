import importlib

def get_module_and_version(pkg):

    module = importlib.import_module(pkg)
    version = getattr(module, "__version__", None)
    print("Module {}, version {}".format(pkg, version))
    return module, version

def test_import_astropy():
    module, version = get_module_and_version("astropy")
    assert version is not None

def test_import_ctapipe():

    module, version = get_module_and_version("ctapipe")
    assert version is not None
