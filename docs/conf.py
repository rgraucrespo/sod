project = "SOD"
author = "Ricardo Grau-Crespo, Said Hamad, and contributors"
release = "0.81"

extensions = []
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_static_path = ["_static"]

try:
    import sphinx_rtd_theme  # noqa: F401
except ImportError:
    html_theme = "alabaster"
else:
    html_theme = "sphinx_rtd_theme"
html_title = f"SOD documentation ({release})"
html_logo = "_static/sod_logo.png"
