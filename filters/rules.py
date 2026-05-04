def lipinski(props):
    violations = sum([
        props["MW"] > 500,
        props["logP"] > 5,
        props["HBA"] > 10,
        props["HBD"] > 5,
    ])
    return violations, violations <= 1


def veber(props):
    passes = props["RotBonds"] <= 10 and props["TPSA"] <= 140
    return passes


def egan(props):
    passes = props["logP"] <= 5.88 and props["TPSA"] <= 131.6
    return passes
