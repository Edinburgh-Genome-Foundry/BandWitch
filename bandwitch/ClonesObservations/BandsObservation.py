from collections import OrderedDict
from bandwagon.ladders import ladder_from_aati_fa_calibration_table
from bandwagon import BandsPattern
from .band_patterns_discrepancy import band_patterns_discrepancy

try:
    from plateo.parsers import plate_from_aati_fragment_analyzer_zip

    PLATEO_AVAILABLE = True
except ImportError:
    PLATEO_AVAILABLE = False


class BandsObservation:
    """One observation of a bands pattern.

    Parameters
    ----------
    name : str
        Name of the observation (used for the top label in plots).
    bands : BandPattern
        A BandPattern object.
    ladder : BandPattern
        A BandPattern object representing a ladder.
    migration_image
        Optional RGB array (HxWx3) representing the gel "image" (which will be
        displayed on the side of the plot).
    """

    def __init__(self, name, bands, ladder, migration_image=None):
        """Initialize."""
        self.name = name
        self.bands = sorted(bands)
        self.ladder = ladder
        self.migration_image = migration_image

    @staticmethod
    def from_aati_fa_archive(
        archive_path,
        min_rfu_size_ratio=0.3,
        ignore_bands_under=None,
        direction="column",
    ):
        """Return a dictionary of all band observations in AATI output files.

        Parameters
        ----------
        archive_path : str
            A path to a ZIP file containing all files output by the AATI fragment
            analyzer.
        min_rfu_size_ratio : float
            Cut-off ratio to filter out bands whose intensity is below some
            threshold. The higher the value, the more bands will be filtered out.

        Returns
        -------
        dict
            A dictionary ``{'A1': BandsObservation(), 'A2': ...}`` containing the
            measured pattern information for a whole 96-well microplate.
        """
        if not PLATEO_AVAILABLE:
            raise ImportError("Plateo must be installed to parse AATI zips.")

        ignore_bands_under = ignore_bands_under or 0

        plate = plate_from_aati_fragment_analyzer_zip(archive_path)
        ladder_data = plate.data["ladder"]
        ladder = ladder_from_aati_fa_calibration_table(dataframe=ladder_data)
        # get the set of constructs sorted by order of columnwise appearance:

        def band_is_strong_enough(band):
            """Return True if the band's intensity is above the set level."""
            return 1.0 * band["RFU"] / band["Size (bp)"] > min_rfu_size_ratio

        ordered_list = []
        for well in plate.iter_wells(direction=direction):
            bands_list = []
            for band in well.data["bands"].values():
                if band_is_strong_enough(band) and (
                    band["Size (bp)"] > ignore_bands_under
                ):
                    bands_list += [band["Size (bp)"]]

            bands_observation = BandsObservation(
                name=well.name,
                ladder=ladder,
                bands=bands_list,
                migration_image=well.data["migration_image"],
            )

            entry = (well.name, bands_observation)
            ordered_list += [entry]

        ordered_dict = OrderedDict(ordered_list)

        return ordered_dict

    def patterns_discrepancy(
        self,
        other_bands,
        relative_tolerance=0.1,
        min_band_cutoff=None,
        max_band_cutoff=None,
    ):
        """Return the maximal discrepancy between two band patterns.

        The discrepancy is defined as the largest distance between a band in
        one pattern and the closest band in the other pattern.

        Parameters
        ----------
        other_bands : list of int
            A list of bands (integers) to be compared with the current bands.
        relative_tolerance : float
            Tolerance, as a ratio of the full ladder span. If =0.1, then the
            discrepancy will have a value of 1 when a band's nearest
            correspondent in the other pattern is more than 10% of the ladder
            span apart.
        min_band_cutoff : int, optional
            Discrepancies involving at least one band below this minimal band
            size will be ignored. By default, it will be set to the smallest
            band size in the ladder.
        max_band_cutoff : int, optional
            Discrepancies involving at least one band above this minimum band
            size will be ignored. By default, it will be set to the largest
            band size in the ladder.
        """
        ladder_min, ladder_max = self.ladder.dna_size_span
        if min_band_cutoff is None:
            min_band_cutoff = ladder_min
        if max_band_cutoff is None:
            max_band_cutoff = ladder_max
        return band_patterns_discrepancy(
            other_bands,
            self.bands,
            ladder=self.ladder,
            relative_tolerance=relative_tolerance,
            zone=[min_band_cutoff, max_band_cutoff],
            reference_and_gel=True,
        )

    def to_bandwagon_bandpattern(self, background_color=None, label="auto"):
        """Return a pattern version for the plotting library BandWagon.

        If label is left to 'auto', it will be the pattern's name.
        """
        if label == "auto":
            label = self.name
        return BandsPattern(
            self.bands,
            corner_note="Total: %d bp." % sum(self.bands),
            ladder=self.ladder,
            label=label,
            gel_image=self.migration_image,
            background_color=background_color,
        )
