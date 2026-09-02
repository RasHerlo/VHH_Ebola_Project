"""Python entry point reserved for a later rewrite of Figure 5 / S3.

MATLAB is the current plotting stack. Use the same config files:

- ``config/paths.json``
- ``config/figure_manifest.json``

Example::

    from ebola_vhh.paths import load_paths, load_manifest
"""

from .paths import load_manifest, load_paths, repo_root

__all__ = ["load_paths", "load_manifest", "repo_root"]
