# Hinweise zu Fremdmaterial

Der Code von elliDFT steht unter der MIT-Lizenz (`LICENSE`). Diese Datei nennt
das Fremdmaterial, das im Repository liegt oder aus dem Teile erzeugt wurden,
und die Bibliotheken, gegen die gebaut wird.

## Im Repository enthalten

### Basissätze (`data/basis/*.nw`)

Heruntergeladen von der Basis Set Exchange (https://www.basissetexchange.org),
Version 0.12, im NWChem-Format. Die Dateien sind unverändert; Quelle, Datum und
SHA-256 jeder Datei stehen in `data/basis/MANIFEST.json`.

Die Basis Set Exchange wird unter der folgenden Lizenz verbreitet:

```
Copyright (c) 2020 The Molecular Sciences Software Institute, Virginia Tech

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived
from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
```

Bitte um Zitat: B. P. Pritchard, D. Altarawy, B. Didier, T. D. Gibson,
T. L. Windus, „A New Basis Set Exchange: An Open, Up-to-date Resource for the
Molecular Sciences Community“, J. Chem. Inf. Model. 59, 4814–4820 (2019),
doi:10.1021/acs.jcim.9b00725. Die Originalarbeiten der einzelnen Basissätze
nennt die Basis Set Exchange je Basissatz.

### Lebedev-Winkelgitter (`src/grid/lebedev_tables.inc`)

Die Stützstellen und Gewichte sind mit `scripts/generate_lebedev_tables.py`
Ziffer für Ziffer aus `Lebedev-Laikov.F` übernommen. Diese Routine stammt von
D. N. Laikov (C-Original), wurde von C. van Wüllen nach Fortran übertragen und
über CCL verbreitet. Die Fortran-Datei selbst liegt nicht im Repository; Quelle
und Prüfsumme stehen in `third_party/lebedev/README.md`.

Bitte der Autoren um Zitat in Veröffentlichungen und Handbüchern:
V. I. Lebedev, D. N. Laikov, „A quadrature formula for the sphere of the 131st
algebraic order of accuracy“, Doklady Mathematics 59, 477–481 (1999).

## Bibliotheken (nicht im Repository)

Beim Bauen werden diese Bibliotheken vorausgesetzt. Keine davon liegt im
Repository.

| Bibliothek | Zweck | Lizenz |
|---|---|---|
| libxc | Austausch-Korrelations-Funktionale | MPL-2.0 |
| Eigen 3 | lineare Algebra | MPL-2.0 |
| nlohmann-json | JSON-Ein- und -Ausgabe | MIT |
| Catch2 3 | Tests | BSL-1.0 |
| OpenMP-Laufzeit (z. B. libomp) | Parallelisierung, optional | je nach Compiler |

Wer elliDFT als Binärpaket weitergibt, in das libxc oder Eigen eingebunden
ist, muss die Bedingungen der MPL-2.0 einhalten: Der Quellcode dieser
Bibliotheken muss zugänglich sein.

libxc bittet um Zitat: S. Lehtola, C. Steigemann, M. J. T. Oliveira,
M. A. L. Marques, „Recent developments in libxc — A comprehensive library of
functionals for density functional theory“, SoftwareX 7, 1–5 (2018).

## Nur als Werkzeug verwendet

Die folgenden Pakete dienen nur der Erzeugung von Referenzdaten und Geometrien
(`tests/reference`, `scripts/make_geometries.py`). Sie werden weder eingebunden
noch mitgeliefert. Die erzeugten Zahlen in `tests/reference/data` und
`data/geometries` sind Ausgaben dieser Programme.

| Paket | Lizenz |
|---|---|
| PySCF | Apache-2.0 |
| NumPy, SciPy | BSD-3-Clause |
| RDKit | BSD-3-Clause |
