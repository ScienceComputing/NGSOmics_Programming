import matplotlib.pyplot as plt
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus

# Read into a spectrum from public data repositories by its Universal Spectrum Identifier
usi_id = "mzspec:PXD000966:CPTAC_CompRef_00_iTRAQ_05_2Feb12_Cougar_11-10-09.mzML:scan:12298"
spectrum = sus.MsmsSpectrum.from_usi(usi_id)

# Annotate the spectrum with its ProForma string
peptide = "EM[L-methionine sulfoxide]EVEES[O-phospho-L-serine]PEK"
fragment_tol_mass, fragment_tol_mode = 10, "ppm"
spectrum = spectrum.annotate_proforma(peptide, fragment_tol_mass, fragment_tol_mode)

# Plot the spectrum
fig, ax = plt.subplots(figsize=(12, 6))
sup.spectrum(spectrum, grid=False, ax=ax)
ax.set_title(peptide, fontdict={"fontsize": "xx-large"})
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
plt.show()
