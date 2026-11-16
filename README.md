# elliDFT

Eigenständiges Kohn-Sham-DFT-Programm mit Gauß-Basis für die Ellipsometrie:
Polarisierbarkeiten α(ω) und elektronische Anregungen von Polymer-Oligomeren,
daraus Brechzahl n(λ) und Extinktion k(λ).

- Hauptprogramm mit Steuerdatei und Protokoll: `docs/Steuerdatei.md`
- Entscheidungen des Auftraggebers: `docs/Entscheidungen.md`
- Auswertung von Ellipsometrie-Messungen: `docs/Ellipsometrie.md`
- Zeitbereichsverfahren (Echtzeit-TDDFT): `docs/Zeitbereich.md`
- Nichtlineare optische Größen: `docs/Nichtlineare-Optik.md`
- Stand und Abnahmekriterien je Meilenstein: `docs/progress/`
- Eingesehene Quellen der implementierten Formeln: `docs/sources/README.md`
- Lizenz: MIT (`LICENSE`); mitgelieferte Fremddaten und Abhängigkeiten: `THIRD_PARTY_NOTICES.md`

Verweise im Code der Form „specification, section N“ meinen das
Auftragsdokument des Projekts. Es wird nicht veröffentlicht; die Festlegungen,
auf die es ankommt, stehen jeweils an der Stelle im Code und in
`docs/Entscheidungen.md`.

## Bauen und testen

    scripts/build_and_test.sh          # kurze Tests, wenige Sekunden
    scripts/build_and_test.sh --slow   # zusätzlich die vollen Abnahmesätze (Kennung slow), für einen Server

Voraussetzungen: C++20-Compiler, CMake ≥ 3.20, Eigen 3, libxc, OpenMP (unter macOS `libomp`),
Catch2 3, nlohmann-json. Die Threadzahl steuert `OMP_NUM_THREADS`.
Die Tests laufen ohne PySCF; die Referenzdaten liegen unter `tests/reference/data`
und lassen sich mit `tests/reference/regenerate.sh` neu erzeugen.

## Programme

    build/ellidft-basisinfo <geometrie.xyz> <basis>   # Größe der Basis als JSON
    build/ellidft-scf <geometrie.xyz> <basis>         # Hartree-Fock, Ergebnis als JSON
    build/ellidft-scf <geometrie.xyz> <basis> --functional CAM-B3LYP --grid 99,590   # Kohn-Sham (LDA, PBE, PBE0, B3LYP, CAM-B3LYP)
    build/ellidft-scf <geometrie.xyz> <basis> --functional PBE --aux-basis def2-universal-JKFIT   # mit Density Fitting

    build/ellidft-response <geometrie.xyz> <basis> --wavelengths 400,589,1000 --states 10
                                                      # Polarisierbarkeit und Anregungen, JSON nach Abschnitt 8
    build/ellidft-optics --density 1.05 --molar-mass 104.15 einheit1.json einheit2.json ...
                                                      # Inkrement, n(λ) und k(λ) aus einer Oligomer-Reihe
    scripts/run_oligomer_series.sh                    # beides in einem Ablauf

Alle lesen XYZ-Dateien in Ångström (`--bohr` für Bohr) und die Basissätze aus
`data/basis` (oder aus dem Verzeichnis in `ELLIDFT_BASIS_DIR`).
