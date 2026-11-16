# Hauptprogramm mit Steuerdatei und Protokoll

Der Auftrag verlangt eine Steuerdatei (JSON oder TOML) und ein Protokoll neben
der JSON-Ausgabe. Das Programm `ellidft` liefert beides; gewählt ist JSON, weil
die Bibliothek dafür schon freigegeben ist.

    ellidft <steuerdatei.json>

```json
{
  "name": "styrene_4mer",
  "geometry": "styrene_4mer.xyz", "bohr": false, "charge": 0,
  "basis": "def2-SVPD", "aux_basis": "def2-universal-JKFIT",
  "functional": "CAM-B3LYP", "grid": [75, 302],
  "scf": {"energy_tolerance": 1e-10, "error_tolerance": 1e-8, "max_iterations": 200, "guess": "atomic"},
  "response": {"wavelengths_nm": [589, 633], "gamma_hartree": 0.0, "states": 0, "tda": false},
  "threads": 0,
  "output_directory": "."
}
```

Pflicht sind nur `geometry` und `basis`. Kommentare mit `//` sind erlaubt,
unbekannte Schlüssel werden abgewiesen. Dateinamen gelten relativ zur
Steuerdatei.

| Ausgabe | Inhalt |
|---|---|
| `<name>.json` | Ergebnisse im Aufbau von Abschnitt 8 des Auftrags; das liest `ellidft-optics` |
| `<name>.log` | Protokoll: Programmversion mit Commit, Übersetzer, Bibliotheksversionen, Einstellungen, jede SCF-Iteration (Energie, Änderung, Kommutator), Energieanteile, Elektronenzahl auf dem Gitter, Antwortgleichungen je Frequenz, Anregungen, Zeiten |

Das Protokoll geht zugleich auf die Fehlerausgabe. Die Einzelprogramme
(`ellidft-scf`, `ellidft-response`) bleiben; `ellidft` rechnet dasselbe wie
`ellidft-response` (geprüft an H₂: gleiche Zahlen bis auf die letzte Stelle).

## Fortsetzen abgebrochener Läufe

`<name>.chk` enthält die Dichtematrix, neu geschrieben nach jeder
SCF-Iteration, und jede fertig gerechnete Polarisierbarkeit. Ein späterer Lauf
mit derselben Steuerdatei beginnt das SCF bei dieser Dichte und rechnet
fertige Frequenzen nicht noch einmal. Das geschieht von selbst, wenn die Datei
zu Molekül, Basis, Methode und Gitter passt; sonst wird sie mit einer Meldung
im Protokoll übergangen. `"restart": false` schaltet es ab.

Geprüft an Wasser (PBE0): Der zweite Lauf braucht 2 statt 8 SCF-Iterationen,
keine Anwendung des Antwortoperators, und liefert dieselben Zahlen (Energie
auf 10⁻¹³ Hartree, Polarisierbarkeiten bitgleich).

Die Datei wird zuerst unter anderem Namen geschrieben und dann umbenannt; ein
Abbruch während des Schreibens lässt die alte Datei unversehrt. Nicht
gespeichert werden Anregungsenergien und Zwischenstände der Antwortgleichungen
innerhalb einer Frequenz.

Nicht enthalten: die Zeitpropagation und die nichtlinearen Größen (eigene
Programme, ohne Checkpoint).
