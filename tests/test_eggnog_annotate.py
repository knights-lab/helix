import pytest
import os

from helix.eggnog_annotate import get_taxamap, get_annotations


@pytest.fixture()
def archaea_annotation():
    return os.path.join("fixtures", "archaea.emapper.annotations")


@pytest.fixture()
def bacteria_annotation():
    return os.path.join("fixtures", "bacteria.emapper.1m.annotations")


@pytest.fixture()
def taxonomy_map() -> dict:
    # return get_taxamap(os.path.join("fixtures", "r95.gtdb.tax"))
    return get_taxamap(os.path.join("fixtures", "taxonomy.tax"))


def test_archaea_annotation(taxonomy_map, archaea_annotation, bacteria_annotation):
    bacteria_map = get_annotations(taxonomy_map, bacteria_annotation)
    # archaea_map = get_annotations(taxonomy_map, archaea_annotation)

    i = 0
    for key, value in taxonomy_map.items():
        if "Archaea" in value:
            i += 1


    print(archaea_map)
