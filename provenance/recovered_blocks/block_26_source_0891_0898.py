from pathlib import Path
original = Path("PCG13_31_BI.nex")
texto = original.read_text()
anterior = "mcmcp ngen=1000000"
assert texto.count(anterior) == 1, "No se encontró una única línea mcmcp esperada"
texto = texto.replace(anterior, "mcmcp ngen=1000000 append=yes", 1)
Path("PCG13_31_BI_continuar.nex").write_text(texto)
print("Preparado: PCG13_31_BI_continuar.nex")
