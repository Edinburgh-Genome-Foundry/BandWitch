import os
import bandwitch.tools as tools

this_directory = os.path.dirname(os.path.realpath(__file__))


def test_set_record_topology():
    record_path = os.path.join(
        this_directory, "test_data", "test_sequences", "asm_00.gb"
    )
    record = tools.load_record(record_path, topology="circular")
    assert record.annotations["topology"] == "circular"
    assert not tools.record_is_linear(record, default=True)

    # Default only when no topology is set
    tools.set_record_topology(record, "default_to_circular")
    assert not tools.record_is_linear(record, default=True)
    record.annotations.pop("topology")
    assert not tools.record_is_linear(record, default=False)
    tools.set_record_topology(record, "default_to_linear")
    assert tools.record_is_linear(record, default=False)
