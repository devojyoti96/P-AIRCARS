import pytest
import shutil
import os


@pytest.fixture
def dummy_caltables(tmp_path):
    """
    Provide a list of dummy CASA caltables for testing.
    They are shallow copies of a real reference caltable.
    """
    # Path to a known valid CASA caltable
    path = os.path.dirname(os.path.abspath(__file__))
    ref_caltable = os.path.dirname(path) + "/testdata/test_caltable.bcal"
    if not os.path.exists(ref_caltable):
        pytest.skip("Reference caltable is not found")
    # Create two dummy copies
    cal1 = tmp_path / "cal1.K"
    cal2 = tmp_path / "cal2.K"
    shutil.copytree(ref_caltable, cal1)
    shutil.copytree(ref_caltable, cal2)
    return [str(cal1), str(cal2)]


@pytest.fixture
def dummy_caltable(tmp_path):
    path = os.path.dirname(os.path.abspath(__file__))
    ref_caltable = os.path.dirname(path) + "/testdata/test_caltable.bcal"
    if not os.path.exists(ref_caltable):
        pytest.skip("Caltable is not found")
    return ref_caltable


@pytest.fixture
def dummy_quartical_table(tmp_path):
    path = os.path.dirname(os.path.abspath(__file__))
    ref_caltable = os.path.dirname(path) + "/testdata/test_caltable.qcal"
    if not os.path.exists(ref_caltable):
        pytest.skip("Caltable is not found")
    return ref_caltable


@pytest.fixture
def dummy_msname(tmp_path):
    path = os.path.dirname(os.path.abspath(__file__))
    ref_msname = os.path.dirname(path) + "/testdata/test_ms.ms"
    if not os.path.exists(ref_msname):
        pytest.skip("Reference ms is not found")
    return ref_msname


@pytest.fixture
def dummy_image():
    path = os.path.dirname(os.path.abspath(__file__))
    ref_imagename = os.path.dirname(path) + "/testdata/test_image.fits"
    if not os.path.exists(ref_imagename):
        pytest.skip("Reference image is not found.")
    return ref_imagename


@pytest.fixture
def dummy_metafits():
    path = os.path.dirname(os.path.abspath(__file__))
    ref_metafits = os.path.dirname(path) + "/testdata/test.metafits"
    if not os.path.exists(ref_metafits):
        pytest.skip("Reference metafits is not found.")
    return ref_metafits
