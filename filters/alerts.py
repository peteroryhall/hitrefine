from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


def build_pains_catalog():
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog(params)


def build_brenk_catalog():
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog(params)


def check_pains(mol, catalog):
    entry = catalog.GetFirstMatch(mol)
    if entry:
        return False, entry.GetDescription()
    return True, None


def check_brenk(mol, catalog):
    entry = catalog.GetFirstMatch(mol)
    if entry:
        return False, entry.GetDescription()
    return True, None
