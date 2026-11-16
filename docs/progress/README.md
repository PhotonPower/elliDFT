# Stand der Meilensteine

Stand 14.11.2026. Einzelheiten und Messwerte stehen in `M0.md` bis `M6.md`;
die Erweiterungen in `E1.md` (Auswertung der Messungen), `E2.md` und `E3.md`
(Zeitbereichsverfahren), `E4.md` und `E5.md` (nichtlineare Größen).
Die Zweige bauen aufeinander auf (`m6-optik`, `entscheidungen-1`,
`ellipsometrie-auswertung`, `zeitbereich`, `nichtlinear`, `nichtlineare-groessen`,
`photopolymere`, `steuerdatei`, `checkpoint`, `speicher-df`, `weiter`); der letzte enthält alles einschließlich der Umsetzung der
Entscheidungen vom 14.11.2026. Auf `main` ist
nichts zusammengeführt.

| Meilenstein | Zweig | Code | Abnahme |
|---|---|---|---|
| M0 Gerüst | `m0-geruest` | fertig | erfüllt |
| M1 Integrale, Hartree-Fock | `m1-integrale` | fertig | erfüllt (alle 60 HF-Fälle einmal vollständig gerechnet) |
| M2 Kohn-Sham | `m2-kohn-sham` | fertig | Lebedev, eigene LDA, Drehinvarianz erfüllt; Energien und Elektronenzahl (relatives Kriterium) für 8 von 32 Fällen gerechnet und erfüllt |
| M3 Skalierung | `m3-skalierung` | fertig | Density Fitting für 7 von 24 Fällen erfüllt; Styrol-Tetramer und volle Skalierung offen |
| M4 Hybride | `m4-hybride` | fertig | Integrale erfüllt; Energien für 10 von 36 Fällen erfüllt; Vergleich RI-K gegen direkt bei 66 Atomen offen |
| M5 Lineare Antwort | `m5-antwort` | fertig | alle referenzfreien Kriterien erfüllt; Referenzvergleich für 1 von 9 Fällen erfüllt |
| M6 Optik | `m6-optik` | fertig | nichts gerechnet (Oligomere, Server) |
| E1 Ellipsometrie-Auswertung | `ellipsometrie-auswertung` | erste Fassung | keine Kriterien gesetzt; Tests mit berechneten Daten bestanden, keine echte Messung |
| E2 Zeitbereichsverfahren | `zeitbereich` | erste Fassung | vorgeschlagene Kriterien erfüllt (H₂, H₂O gegen den Antwortcode ≤ 2·10⁻⁵); Energieerhaltung des Magnus-Schritts bei starken Pulsen schwach |
| E3 Zeitbereich: Energieerhaltung, nichtlineare Antwort | `nichtlinear` | erste Fassung | Energie nach starkem Puls auf 10⁻¹¹ erhalten; Antwort bis zur dritten Ordnung im statischen Grenzfall gegen das SCF im Feld geprüft; frequenzabhängige Größen fehlen |
| E4 Nichtlineare Größen | `nichtlineare-groessen` | erste Fassung | keine Kriterien gesetzt; zweite Ordnung gegen unabhängigen Weg auf 5·10⁻⁵, dritte auf 0,4 %; Betrag des Zwei-Photonen-Querschnitts ungeprüft |
| E5 Nichtlineare Größen für Photopolymere | `photopolymere` | erste Fassung | Orientierungsmittel, n₂ des Films, Programm; Proben gegen den Antwortcode auf 5·10⁻⁵ bis 10⁻³; keine Rechnung an einem echten Initiator |
| E6 Hauptprogramm, Server-Paket, Checkpoint | `weiter` | fertig | Proben an kleinen Molekülen bestanden; Messungen am Tetramer stehen aus |
| M7 GPU | — | nicht begonnen | auf diesem Rechner nicht möglich (keine NVIDIA-Karte) |

## Auf dem Server als Erstes

    python3 -m venv .venv && .venv/bin/pip install -r tests/reference/requirements.txt
    scripts/server_run.sh --series

Ein Aufruf baut, fährt die kurzen und die `slow`-Tests (alle offenen
Referenzvergleiche aus M1 bis M5 und die ausführlichen Proben der
Erweiterungen), misst am Styrol-Tetramer die Zeiten und Speicherspitzen für M3
und M4 (mit PySCF als Zeitreferenz, wenn die `.venv` vorhanden ist) und rechnet
mit `--series` die Oligomer-Reihen für M6. Alles landet in einem Verzeichnis
`results/server_<Datum>`; `report.md` darin ist die Zusammenfassung. Dieses
Verzeichnis bitte zurückgeben.

- `--with-direct` misst zusätzlich das Tetramer ohne Density Fitting (lang).
- `--quick` prüft nur den Ablauf, in Sekunden, mit Wasser statt des Tetramers.
- Für das Laufzeitkriterium von M3 müssen beide Programme gleich viele Threads
  bekommen (`OMP_NUM_THREADS` vor dem Aufruf setzen).
- Ein fehlgeschlagener Schritt bricht den Lauf nicht ab; der Bericht nennt ihn.

Auf dem Entwicklungsrechner ist nur `--quick` gelaufen. Dabei stimmten die
Energien von elliDFT und PySCF für Wasser (PBE, RI-J) auf 2·10⁻¹³ Hartree
überein. Die Zeiten und der Speicherbedarf am Tetramer sind nicht gemessen.

## Entscheidungen des Auftraggebers

Alle getroffenen Entscheidungen stehen in `docs/Entscheidungen.md`. Offen sind
noch: der Austauschaufbau für große Systeme (nach den Messungen auf dem
Server), der Zahlenwert der Dichte für Polystyrol und PMMA und die
gemischte Genauigkeit auf der GPU.
