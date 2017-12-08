from collections import OrderedDict
from io import BytesIO

from Bio import Restriction
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from dna_features_viewer import (GraphicRecord, GraphicFeature,
                                 BiopythonTranslator)

def create_cuts_map_graphic_record(record, enzymes):
    """Create a DnaFeaturesViewer graphic record showing the cuts.

    Parameters
    ----------
    record
      A Biopython record

    enzymes
      A list of enzyme names. e.g. ('EcoRI'. 'BamHI')

    """
    batch = Restriction.RestrictionBatch(enzymes)
    cuts_dict = batch.search(record.seq)
    all_cuts = sorted(
        set([0, len(record)] +
            [c for cc in cuts_dict.values() for c in cc])
    )
    bands = list(zip(all_cuts, all_cuts[1:]))
    if not record.__dict__.get('linear', True) and (len(bands) > 1):
        start, end = bands.pop()
        band0 = [-(end - start), bands[0][1]]
        if bands == []:
            bands = [band0]
        else:
            bands[0] = band0
    sorted_bands = sorted(
        bands, key=lambda start_end: start_end[0] - start_end[1])
    gr_cuts = GraphicRecord(len(record), [
        GraphicFeature(_start, _end,  strand=1, label=name,
                       color="#ede15c", thickness=5)
        for name, (_start, _end) in zip("abcdefghijklmn", sorted_bands)
    ])
    gr_cuts.split_overflowing_features_circularly()
    return gr_cuts, all_cuts

def plot_cuts_map(record, enzymes, ax):
    """Plot the record and its cuts one above the other.

    Parameters
    ----------
    record
      A Biopython record

    enzymes
      A list of enzyme names. e.g. ('EcoRI'. 'BamHI')

    ax
      Matplotlib ax on which to plot the figure.

    """

    gr_cuts, all_cuts = create_cuts_map_graphic_record(
        record, enzymes)
    gr_cuts.plot(ax, fontsize=7, with_ruler=False)
    for cut in all_cuts:
        ax.axvline(cut, ls=":", color="k", lw=0.5)

    def features_prop(f):
        return dict(label=f.qualifiers.get("source", [False])[0],
                    color="#93a8ea",
                    thickness=5)

    def features_filter(f):
        return features_prop(f)["label"]

    translator = BiopythonTranslator([features_filter], features_prop)
    gr_record = translator.translate_record(record)
    gr_record.plot(ax, fontsize=4, level_offset=7)

def plot_all_constructs_cuts_maps(record_digestion_pairs, target=None,
                                  figsize=(12, 4)):
    """Plot schemas of all constructs with cuts, in a multipage PDF.

    Parameters
    ----------
    target
      Either None (at which case the function returns raw data of the PDF)
      or 'base64' (the function returns base64-encoded data, e.g. for
      web transfer), or a file path to be written to
    figsize
      The size in inches of each figure (=page of the pdf).

    """
    pdf_io = BytesIO()
    with PdfPages(pdf_io) as pdf:
        for record, digestion in record_digestion_pairs:
            fig, ax = plt.subplots(1, figsize=figsize)
            plot_cuts_map(record, digestion, ax)
            title = record.id + '\n' + " + ".join(digestion)
            ax.set_title(title, fontdict=dict(weight='bold'))
            pdf.attach_note(record.id)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    pdf_data = pdf_io.getvalue()
    if target is None:
        return pdf_data
    elif target == 'base64':
        return 'data:application/pdf;base64,' + pdf_data.decode("utf-8")
    else:
        with open(target, 'wb') as f:
            f.write(pdf_data)
