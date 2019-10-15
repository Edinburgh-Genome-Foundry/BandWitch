from plateo.parsers import plate_from_platemap_spreadsheet
from bandwitch import BandsObservation, ClonesObservations, Clone, load_record
import flametree

# LOAD ALL THE DATA (constructs records, clones maps, fragment analysis)

data_dir = flametree.file_tree('.').example_data.band_validation_data

constructs_dict = {
    f._name_no_extension: load_record(
        f._path, topology='circular', id=f._name_no_extension)
    for f in data_dir.constructs._all_files
}

clones_plate = plate_from_platemap_spreadsheet(data_dir.clones_map_xlsx._path,
                                               data_field='construct')
clones_map = {
    well.name: well.data.construct
    for well in clones_plate.iter_wells()
    if 'construct' in well.data
    and str(well.data.construct) != 'nan'
}
bands_observations = BandsObservation.from_aati_fa_archive(
    data_dir.aati_files_zip._path)

clones = {
    well_name: Clone(name=well_name, construct_id=clones_map[well_name],
                     digestions={("NcoI",): bands_observations[well_name]})
    for well_name, construct_id in clones_map.items()
}


# VALIDATE ALL CLONES WITH BANDWITCH
clones_observations = ClonesObservations(clones, constructs_dict)
validations = clones_observations.validate_all_clones(relative_tolerance=0.03)
validations_summary = clones_observations.validations_summary(validations)


# CREATE A FOLDER WITH VALIDATION REPORTS
report_root = flametree.file_tree('.').examples_output._dir('band_validation')
report_root._file('validations.pdf').write(
    clones_observations.plot_all_validations_patterns(validations)
)
ax = clones_observations.plot_validations_plate_map(validations)
ax.figure.savefig(report_root._file('success_map.pdf').open('wb'),
                  format='pdf', bbox_inches='tight')
