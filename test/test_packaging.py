"""Tests for package metadata and command-line entry points."""

import importlib

import prepmd
import pytest


def test_package_version_is_defined():
    """The package should expose a valid version string."""
    assert isinstance(prepmd.__version__, str)
    assert prepmd.__version__


@pytest.mark.parametrize(
    ("module_name", "function_name"),
    [
        ("prepmd.prep", "entry_point"),
        ("prepmd.run", "entry_point"),
        ("prepmd.add_modeller_license", "entry_point"),
    ],
)
def test_command_entry_points_are_importable(module_name, function_name):
    """Every command declared in pyproject.toml should be importable."""
    module = importlib.import_module(module_name)

    assert callable(getattr(module, function_name))