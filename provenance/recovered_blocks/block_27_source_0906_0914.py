from pathlib import Path
filas = [x.split() for x in Path("PCG13_31_BI.mcmc").read_text().splitlines()
         if x.strip() and x.lstrip()[0].isdigit()]
print("Generación y AvgStdDev:", filas[-1][0], filas[-1][-1])
datos = [x.split() for x in Path("PCG13_31_BI.pstat").read_text().splitlines()[2:]
         if len(x.split()) >= 9]
print("ESS mínimo:", min(float(x[6]) for x in datos))
print("PSRF mínimo y máximo:", min(float(x[8]) for x in datos),
      max(float(x[8]) for x in datos))
