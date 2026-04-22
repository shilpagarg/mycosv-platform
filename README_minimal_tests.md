Minimal additions inspired by minimap2 / WhatsHap, but kept very small.

Files:
- test_golden_smoke.py      : reproducibility + tiny golden/regression smoke
- test_sv_report_smoke.py   : visualization report smoke test
- tox.ini                   : minimal pytest runner

Drop these files into your project root, then run:

pytest -q test_golden_smoke.py test_sv_report_smoke.py

or

tox
