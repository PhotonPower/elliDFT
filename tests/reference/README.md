# Referenzdaten

Hier liegen die Skripte, die Referenzwerte erzeugen, und die erzeugten
JSON-Dateien (`data/`). Jede Referenzzahl der Tests stammt aus einem dieser
Skripte. Referenzcode ist PySCF; es wird nie in elliDFT gelinkt, und die Tests
laufen ohne installiertes PySCF.

## Neu erzeugen

    python3 -m venv .venv && .venv/bin/pip install -r tests/reference/requirements.txt
    tests/reference/regenerate.sh m0|m1|m2|m3|m4|m5|m6|all

| Skript | Ausgabe | Inhalt |
|---|---|---|
| `gen_m0.py` | `data/m0/` | Schalen aller Elemente und Basissätze, Geometrie und Basisgrößen aller Testmoleküle |
| `gen_boys.py` | `data/m1/boys.json` | Boys-Funktion mit 60 Dezimalstellen (ohne PySCF) |
| `gen_m1.py` | `data/m1/` | Einelektronenmatrizen, Zweielektronenintegrale, Fall mit l bis 6, Hartree-Fock |
| `gen_m2.py` | `data/m2/` | Kleine Molekülgitter, Basisfunktionen an Punkten, LDA- und PBE-Energien auf (99, 590) |
| `gen_m3.py` | `data/m3/` | Integrale mit zwei und drei Zentren, Energien mit und ohne Density Fitting |
| `gen_m4.py` | `data/m4/` | erf-gedämpfte Integrale, Energien mit PBE0, B3LYP und CAM-B3LYP |
| `gen_m5.py` | `data/m5/` | Singulett-Anregungen (volle Antwort und Tamm-Dancoff) und statische Polarisierbarkeit aus PySCFs Matrizen A und B |
| `gen_m6.py` | `data/m6/` | Gemessene Brechzahlen von Polystyrol und PMMA mit Quelle, Zielwerte (ohne PySCF) |

## Gleichgestellte Konventionen

Festgelegt in `common.py` und im Kopf von `gen_m1.py`:

- Geometrie: aus `data/geometries/*.xyz`, Umrechnung mit 0,529177210903 Å pro
  Bohr, Übergabe an PySCF in Bohr.
- Basis: dieselben Dateien aus `data/basis` wie elliDFT, gelesen mit
  `optimize=False` (Dateireihenfolge innerhalb eines Drehimpulses bleibt).
- Sphärisch oder kartesisch wird immer ausdrücklich gesetzt. 6-31G* ist
  kartesisch, alle anderen Basissätze sind sphärisch.
- Keine Punktgruppensymmetrie, kein Density Fitting.
- Jede Datei enthält unter `meta` die Versionen von PySCF, NumPy, SciPy und
  Python sowie alle Einstellungen.

Die Dateien enthalten kein Datum, damit ein erneuter Lauf ohne Änderung der
Zahlen keinen Unterschied im Repository erzeugt.
