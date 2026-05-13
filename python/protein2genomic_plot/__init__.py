"""Plotting backend for protein2genomic.

Invoked via the top-level shell wrapper as `protein2genomic plot ...` or
directly as `python -m protein2genomic_plot.plot ...`.

Reads `isoform_structure.tsv` (the plot-ready table emitted by the C++ binary)
and renders one figure per `input_id` to PDF/PNG/SVG (matplotlib) or to an
interactive HTML page (plotly, optional).
"""

__version__ = "2.2.0"
