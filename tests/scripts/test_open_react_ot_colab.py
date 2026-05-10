from scripts import open_react_ot_colab as opener


def test_colab_url_default_branch():
    url = opener.colab_url()
    assert url == (
        "https://colab.research.google.com/github/PabloAMC/Maillard/blob/"
        "qm-barriers/notebooks/react_ot_colab_gpu.ipynb"
    )


def test_github_url_default_branch():
    url = opener.github_url()
    assert url == (
        "https://github.com/PabloAMC/Maillard/blob/qm-barriers/"
        "notebooks/react_ot_colab_gpu.ipynb"
    )


def test_colab_url_custom_branch():
    assert "blob/main/" in opener.colab_url(branch="main")
