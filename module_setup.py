def modifyParams(params):
    params["packages"].append("dune.fem._tutorial.doc")
    params["package_dir"]["dune.fem._tutorial.doc"] = "_tutorial/doc"
    params["package_data"] = {
        "dune.fem._tutorial.doc": ["*_nb.ipynb", "*_script.py"],
    }
    params["include_package_data"] = False
    return params
