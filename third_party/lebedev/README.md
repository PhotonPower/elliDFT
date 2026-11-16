# Lebedev-Laikov-Routine

Die Winkelgitter in `src/grid/lebedev_tables.inc` stammen aus
`Lebedev-Laikov.F`, der Fortran-Fassung der Routine für Lebedev-Winkelgitter
(C-Original von D. N. Laikov, nach Fortran übertragen von C. van Wüllen), wie
sie über CCL verbreitet wird.

**Die Fortran-Datei selbst liegt nicht im Repository** (Entscheidung vom
16.11.2026, `docs/Entscheidungen.md`). Mitgeliefert werden nur die daraus
erzeugten Tabellen. Wer sie neu erzeugen will, lädt die Datei von der Quelle
unten in dieses Verzeichnis und prüft die Prüfsumme; `.gitignore` hält sie aus
der Versionsverwaltung heraus.

| | |
|---|---|
| Quelle | http://www.ccl.net/cca/software/SOURCES/FORTRAN/Lebedev-Laikov-Grids/Lebedev-Laikov.F |
| Geladen am | 14.11.2026 |
| SHA-256 | 8c4c594ea16160da117ba84d683f53341c36a99210c1aa3d91253aa416df0bb1 |
| Verwendung | `scripts/generate_lebedev_tables.py` liest die Datei und schreibt `src/grid/lebedev_tables.inc`. Die Fortran-Datei wird nicht übersetzt. |
| Entscheidung | Vom Auftraggeber am 14.11.2026 als Quelle bestimmt |

**Lizenz.** Die Datei enthält keinen Lizenztext. Laut dem README der
CCL-Verteilung hat D. N. Laikov erlaubt, die Routinen öffentlich zugänglich zu
machen. Im Kopf der Datei steht die Bitte, V. I. Lebedev, D. N. Laikov,
„A quadrature formula for the sphere of the 131st algebraic order of
accuracy“, Doklady Mathematics 59, 477–481 (1999) in Veröffentlichungen und in
Handbüchern zu zitieren, die den Code beschreiben. Siehe auch
`THIRD_PARTY_NOTICES.md`.
