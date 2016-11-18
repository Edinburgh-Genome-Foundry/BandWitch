import matplotlib.pyplot as plt

def plot_separating_digests(gel_simulator, digestions,
                            sequences_digestions_dict,
                            axes=None, band_thickness=2):
    sequences = list(sequences_digestions_dict.keys())

    if axes is None:
        fig, axes = plt.subplots(len(digestions), 1,
                                 figsize=(0.5 * len(sequences),
                                          2 * len(digestions)))
    if len(digestions) == 1:
        axes = [axes]
    for ax, digestion in zip(axes, digestions):
        patterns = [(name, dig_dict[digestion]["observed_bands"])
                    for name, dig_dict in sequences_digestions_dict.items()]
        ax.axhline(1, c="k", lw=0.5)
        gel_simulator.format_ax(ax)
        gel_simulator.plot_bands_patterns(
            patterns, ax=ax, color='black', with_size_labels=False,
            band_thickness=band_thickness, plot_ladder=False,
            plot_ladder_ticks=True, plot_labels=(ax == axes[0]),
            ylabel=" + ".join(digestion)
        )
    return axes
