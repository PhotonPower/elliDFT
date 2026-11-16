# Beispiel: Polymerfilm auf einem Substrat

Alle Zahlen hier sind berechnet oder erfunden, keine Messung.

| Datei | Inhalt |
|---|---|
| `simulate.json` | Vorwärtsrechnung: 132 nm Film (Sellmeier-Koeffizienten von Polystyrol) auf 2 nm Oxid |
| `example_data.txt` | damit berechnete Ψ und Δ bei 65°, 70° und 75° |
| `fit.json` | Fit von Dicke und Cauchy-Koeffizienten an diese Tabelle |
| `substrate_invented.txt` | **erfundene** Substrattabelle in der Größenordnung von Silizium |

```
ellidft-ellipsometry simulate simulate.json --table sim.txt > sim.json
ellidft-ellipsometry fit fit.json --table fit.txt > fit_result.json
```

Der Fit liefert 132,000 nm und eine Brechzahl, die über den Bereich höchstens
um 3·10⁻⁵ von der Sellmeier-Kurve abweicht (Cauchy- gegen Sellmeier-Form).
