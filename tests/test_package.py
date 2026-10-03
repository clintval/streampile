from importlib.resources import files

import streampile


def test_package_is_typed() -> None:
    assert files(streampile).joinpath("py.typed").is_file()
