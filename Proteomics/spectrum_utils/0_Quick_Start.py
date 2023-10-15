# pip install spectrum_utils

import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
import matplotlib.pyplot as plt

# Read into a spectrum from public data repositories by its Universal Spectrum Identifier.
# https://www.psidev.info/usi 
usi_id = "mzspec:PXD000966:CPTAC_CompRef_00_iTRAQ_05_2Feb12_Cougar_11-10-09.mzML:scan:12298:[iTRAQ4plex]-LHFFM[Oxidation]PGFAPLTSR/3"
peptide = "WNQLQAFWGTGK"
spectrum = sus.MsmsSpectrum.from_usi(usi_id)

# Perform the quality control
fragment_tol_mass, fragment_tol_mode = 10, "ppm"
spectrum = (
    # Focus the mass range to 100 – 1400 m/z to filter out irrelevant peaks
    spectrum.set_mz_range(min_mz=100, max_mz=1400)
    # Remove the precursor peak
    .remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
    # Remove noise peaks of low intensity
    # Keep peaks that are at at least 5% of the base peak intensity; restrict the total number of peaks to the 50 most intense peaks
    .filter_intensity(min_intensity=0.05, max_num_peaks=50)
    # Scale the peak intensities by their square root, to de-emphasize overly intense peaks
    .scale_intensity("root")
    # Annotate peaks corresponding to a, b, and y peptide fragments
    .annotate_proforma(
        peptide, fragment_tol_mass, fragment_tol_mode, ion_types="aby"
    )
)

# Visualize the spectrum with the annotated peaks 
fig, ax = plt.subplots(figsize=(12, 6))
sup.spectrum(spectrum, grid=False, ax=ax)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
plt.show()
# plt.savefig("trial.png", bbox_inches="tight", dpi=300, transparent=True)
# plt.close()
