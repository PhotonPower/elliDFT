# Entscheidungen des Auftraggebers

Der Auftrag (nicht veröffentlicht, siehe Eintrag vom 16.11.2026) behält sich mehrere Punkte vor. Hier steht, was
entschieden wurde, wann, und wo es umgesetzt ist. Wo eine Entscheidung den
Wortlaut des Auftrags ändert, gilt diese Liste.

| Datum | Punkt | Entscheidung | Umsetzung |
|---|---|---|---|
| 14.11.2026 | JSON in C++ | nlohmann-json | `CMakeLists.txt` |
| 14.11.2026 | libxc | per Homebrew installieren | Rechner; `find_package(Libxc)` |
| 14.11.2026 | Lebedev-Gitter | aus der Originalroutine (Lebedev-Laikov.F über CCL), Tabellen per Skript | `third_party/lebedev`, `scripts/generate_lebedev_tables.py` |
| 14.11.2026 | Tests auf dem Entwicklungsrechner | nur sehr kurze Tests, weiter implementieren; volle Abnahmeläufe später auf einem Server | Kennung `slow`, `scripts/build_and_test.sh --slow` |
| 14.11.2026 | Oligomer-Geometrien (Abschnitt 9) | gestreckte Ketten (all-trans-Rückgrat) entlang z, syndiotaktisch, RDKit mit MMFF94; isotaktisch für 1 bis 3 Einheiten als Kontrolle; gemeinsame Kettenachse | `scripts/make_geometries.py`, `data/geometries`, `ellidft-optics --frame lab` |
| 14.11.2026 | Elektronenzahl-Kriterium in M2 | relativ statt absolut: Fehler pro Elektron höchstens 10⁻⁷ auf (99, 590) und 10⁻⁶ auf (75, 302) | `tests/regression/test_m2_kohn_sham.cpp`, `tests/unit/test_kohn_sham.cpp` |
| 14.11.2026 | B3LYP-Variante | „B3LYP“ ist VWN-RPA (`HYB_GGA_XC_B3LYP`); VWN5 heißt `B3LYP5` | `src/xc/functional.cpp` |
| 14.11.2026 | BLAS/LAPACK | als Build-Option, erst nach Messung einschalten | `-DELLIDFT_USE_BLAS=ON`, Voreinstellung aus |
| 14.11.2026 | Argon | bleibt Testatom; Produktumfang bleibt H bis Cl | `src/core/elements.hpp` |
| 14.11.2026 | Dichte für M6 | der eigene Messwert der Proben oder das Datenblatt des Herstellers, kein Lehrbuchwert | Eingabe von `ellidft-optics`; **der Wert selbst steht noch aus** |
| 14.11.2026 | Lizenzen | MIT für den Code behalten; Basissatzdateien und Lebedev-Datei vor einer Veröffentlichung prüfen lassen | `LICENSE`; Prüfung steht aus |
| 14.11.2026 | Auswertung der Ellipsometrie-Messung (Abschnitt 9: bisher ausgeschlossen) | kommt in den Umfang: Vorwärtsmodell und Fit. Der Auftraggeber hat keine Auswertesoftware („auswerteroutinen könntest du auch hinzufügen“). | `src/ellipsometry`, `ellidft-ellipsometry`, `docs/Ellipsometrie.md` |
| 14.11.2026 | Zeitbereichsverfahren (im Auftrag nicht genannt) | als eigenes Modul bauen: Echtzeit-Propagation nach einem Feldpuls. Hauptzweck Gegenprüfung des Antwortcodes und Breitbandspektren; die Rechnung für n und k bleibt im Frequenzbereich. Kriterien: α(ω) gegen den Antwortcode auf 10⁻⁴ relativ (H₂, H₂O; HF, PBE, PBE0), Elektronenzahl auf 10⁻¹⁰, Energie nach dem Puls auf 10⁻⁷ Hartree. | `src/realtime`, `ellidft-realtime`, `docs/Zeitbereich.md` |
| 14.11.2026 | Zeitbereich: Propagator | den energieerhaltenden Propagator vierter Ordnung bauen („1 ja“) | `Propagator::Composition4`, `--propagator composition4` |
| 14.11.2026 | Zeitbereich: Zweck | nichtlineare Optik („2: nichtlineare optik“) | bisher: Felder mit Zeitverlauf und Antwort nach Ordnungen im Feld (`dipole_orders`); welche Größen gebraucht werden, ist offen |
| 14.11.2026 | Nichtlineare Größen | Zwei-Photonen-Absorption, Frequenzverdopplung, Kerr-Effekt. „Kerr“ als optischen Kerr-Effekt gelesen; der elektrooptische ist zusätzlich da. | `src/realtime/nonlinear`, `docs/Nichtlineare-Optik.md` |
| 15.11.2026 | Nichtlineare Optik: Anwendung | Photopolymere; gemeint ist der optische Kerr-Effekt | `docs/Nichtlineare-Optik.md` |
| 15.11.2026 | Zwei-Photonen-Querschnitt und Filmgrößen | Vorschlag angenommen: GM, je absorbiertem Photon, ein Strahl, linear polarisiert, Orientierungsmittel, Lorentz-Linie mit 0,1 eV halber Breite als Voreinstellung; n₂ des Films mit Lorentz-Lokalfeld; keine Filmgröße für die Frequenzverdopplung; Programm mit Steuerdatei; Tests mit PBE0 und CAM-B3LYP | `src/realtime/nonlinear`, `ellidft-nonlinear` |
| 15.11.2026 | Photopolymere: Anwendung | Zwei-Photonen-Lithographie (nicht Holographie) | `docs/Nichtlineare-Optik.md`, Abschnitt „Zwei-Photonen-Lithographie“ |
| 15.11.2026 | Push nach GitHub | Das Repository wird zuerst auf privat gestellt, dann werden alle Zweige gepusht. Grund: Es war öffentlich, und die Fremddaten (Basissätze, Lebedev-Datei) sind vor einer Veröffentlichung rechtlich zu prüfen. | Sichtbarkeit am 15.11.2026 auf privat gestellt |
| 16.11.2026 | Rechte und Lizenz | Das Projekt ist in der Freizeit des Auftraggebers entstanden; die Rechte liegen bei ihm. MIT bleibt. | `LICENSE` |
| 16.11.2026 | Fremddaten vor der Veröffentlichung | Hinweise auf Fremddaten und Abhängigkeiten in einer eigenen Datei; `Lebedev-Laikov.F` nicht mehr im Repository (nur die daraus erzeugten Tabellen); das Auftragsdokument nicht mehr im Repository. Beide liegen nur lokal und sind in `.gitignore` eingetragen. | `THIRD_PARTY_NOTICES.md`, `third_party/lebedev/README.md`, `.gitignore` |
| 16.11.2026 | Git-Historie | Die Veröffentlichung beginnt mit einer neuen Historie aus einem Commit; die Entwicklungsgeschichte bis dahin wird nicht veröffentlicht. Neues GitHub-Repository statt Überschreiben des alten. | Commit „elliDFT: erster öffentlicher Stand“ |
| 16.11.2026 | Veröffentlichung | Repository öffentlich; Änderungen gehen direkt auf `main` und werden ohne Rückfrage gepusht | Sichtbarkeit am 16.11.2026 auf öffentlich gestellt |

## Erläuterungen

**Elektronenzahl.** Der Auftrag nennt einen absoluten Fehler von 10⁻⁶ auf dem
Gitter (99, 590) und 10⁻⁵ auf (75, 302). Der absolute Fehler wächst mit der
Molekülgröße; schon für Benzol liefert der Referenzcode auf demselben Gitter
1,7·10⁻⁶. Nach den Referenzwerten ist der größte relative Fehler auf (99, 590)
5,2·10⁻⁸ (CH₄), Benzol liegt bei 4,0·10⁻⁸.

**Oligomer-Geometrien.** Jedes Oligomer ist ein Ausschnitt eines um eine
Einheit längeren Vorlage-Oligomers: Schweratome an dessen Positionen gehalten
(Abweichung höchstens 0,08 Å), Wasserstoffe relaxiert. Dadurch liegt Einheit k
in jedem Oligomer am selben Ort, und das Inkrement α(m) − α(m−1) beschreibt
das Hinzufügen einer Einheit, nicht eine Konformationsänderung. Das Rückgrat
liegt in der xz-Ebene mit 2,54 Å pro Einheit entlang z. MMFF94 ist ein
Kraftfeld; eine Optimierung auf dem Niveau der Rechnungen ist das nicht.
Die isotaktische all-trans-Kette ist sterisch gedrängt und nur als Kontrolle
gedacht.

**Ellipsometrie-Auswertung.** Die Entscheidung hebt nur diesen einen Ausschluss
aus Abschnitt 9 auf. Die übrigen (offene Schalen, schwerere Elemente,
periodische Rechnungen, Gradienten, Dispersion und Lösungsmittel,
Molekulardynamik) gelten weiter. Die Regeln des Auftrags gelten auch für den
neuen Teil; Abnahmekriterien hat der Auftraggeber dafür noch nicht gesetzt.

**Zeitbereichsverfahren.** Vorschlag und Kriterien stammen von mir; der
Auftraggeber hat auf die Frage „jetzt bauen oder erst nach Steuerdatei,
Ablaufskript und Server-Paket?“ mit „ja bitte“ geantwortet. Ich habe das als
„jetzt bauen“ gelesen.

## Noch offen

- Austauschaufbau für große Systeme (Abschnitt 9): nach den Messungen am
  Styrol-Tetramer auf dem Server.
- Gemischte Genauigkeit auf der GPU (M7).
- Zahlenwert der Dichte für Polystyrol und PMMA.
- Ellipsometrie-Auswertung: Abnahmekriterien; Quelle für Substratdaten
  (Silizium, Oxid, Glas); Gerät und Dateiformat der Messungen.
- Zwei-Photonen-Lithographie: welche Photoinitiatoren und Monomere, welche
  Schreibwellenlänge; Abgleich der Querschnittsfestlegung mit der Literatur.
