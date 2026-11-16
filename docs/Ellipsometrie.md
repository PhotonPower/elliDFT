# Auswertung von Ellipsometrie-Messungen

Stand: 15.11.2026, zweite Fassung (einachsige Schichten, Rauigkeit und
Mischschichten, Rückseite dicker Substrate, absorbierende Modelle). Das Programm `ellidft-ellipsometry` rechnet
für einen Schichtstapel Ψ und Δ aus (`simulate`) und passt Schichtdicken und
Materialparameter an gemessene Werte an (`fit`). Der Auftrag schloss das
zunächst aus; der Auftraggeber hat es am 14.11.2026 in den Umfang genommen
(`docs/Entscheidungen.md`).

**Was noch nie geprüft wurde:** Das Programm hat noch keine echte Messung
gesehen und ist mit keinem anderen Ellipsometrie-Programm verglichen. Alle
Tests rechnen mit Daten, die aus bekannten Modellen erzeugt sind. Siehe
„Grenzen“ unten.

## Aufruf

```
ellidft-ellipsometry simulate <modell.json> [--table <datei>]
ellidft-ellipsometry fit      <modell.json> [--table <datei>]
```

Das Ergebnis kommt als JSON auf die Standardausgabe, ein kurzes Protokoll auf
die Fehlerausgabe. `--table` schreibt gemessene und berechnete Werte als
Texttabelle. Dateinamen in der Modelldatei gelten relativ zu ihrem Verzeichnis.
Ein Beispiel liegt in `examples/ellipsometry`.

## Messdaten

Texttabelle, eine Zeile je Messpunkt, Trennzeichen Leerraum, Komma oder
Semikolon; Leerzeilen und Zeilen mit `#` werden übersprungen:

```
# wavelength_nm  angle_deg  psi_deg  delta_deg  [sigma_psi_deg  sigma_delta_deg]
500.0  70.0  12.51  170.25
```

Mehrere Einfallswinkel stehen einfach untereinander. Gerätedateien liest das
Programm nicht.

## Konventionen

| Größe | Festlegung |
|---|---|
| ρ | ρ = r_p / r_s = tan Ψ · e^{iΔ} |
| Vorzeichen von r_p | wie in den Fresnel-Formeln der eingesehenen Quelle: bei senkrechtem Einfall r_p = −r_s |
| Δ | in [0°, 360°); ein blankes durchsichtiges Substrat hat unterhalb des Brewster-Winkels Δ = 180°, darüber 0° |
| Brechzahl | N = n + i k, k ≥ 0 bedeutet Absorption |
| Wellenlänge | Vakuumwellenlänge in nm |

Liefert ein Gerät Δ mit der anderen Vorzeichenwahl für r_p, sind die Werte um
180° verschoben. Dafür gibt es `"delta_offset_deg": 180` beim Einlesen. Ob das
bei Ihrem Gerät nötig ist, zeigt eine Messung am blanken Substrat.

## Modelldatei

```json
{
  "ambient":   {"model": "constant", "n": 1.0},
  "layers": [
    {"name": "film",
     "thickness_nm": {"value": 100, "fit": true, "min": 0, "max": 1000, "scan": [50, 950, 10]},
     "material": {"model": "cauchy", "A": {"value": 1.5, "fit": true}, "B": {"value": 0.005, "fit": true}}},
    {"name": "oxide", "thickness_nm": 2.0, "material": {"model": "cauchy", "A": 1.45, "B": 0.0036}}
  ],
  "substrate": {"model": "table", "file": "substrat.txt"},
  "data":      {"file": "messung.txt", "delta_offset_deg": 0},
  "residuals": "ncs",
  "simulate":  {"wavelengths_nm": {"from": 400, "to": 900, "step": 10}, "angles_deg": [65, 70, 75]}
}
```

Die Schichten stehen von oben (an der Umgebung) nach unten. `ambient` darf
fehlen (dann n = 1). `data` und `residuals` braucht nur `fit`, `simulate` nur
`simulate`. Kommentare mit `//` sind erlaubt. Unbekannte Schlüssel werden
abgewiesen.

Jede Zahl eines Modells ist entweder fest (einfache Zahl) oder ein Objekt:

| Schlüssel | Bedeutung |
|---|---|
| `value` | Wert bzw. Startwert |
| `fit` | `true`: wird angepasst |
| `min`, `max` | Schranken (Voreinstellung: Dicke ≥ 0, Cauchy-A ≥ 1, Sellmeier-Koeffizienten ≥ 0, Dichte ≥ 0) |
| `scan` | `[von, bis, Anzahl]`: zusätzliche Startwerte; der Fit läuft von jedem und behält das beste Ergebnis |

`scan` ist für Schichtdicken gedacht: Bei dickeren Filmen hat der Fit viele
Nebenminima, und ein einzelner Startwert landet leicht im falschen.

### Materialmodelle

| `model` | Formel | Parameter |
|---|---|---|
| `constant` | n, k unabhängig von der Wellenlänge | `n`, `k` |
| `cauchy` | n = A + B/λ² + C/λ⁴, λ in µm, k = 0 | `A`, `B`, `C` |
| `sellmeier` | n² = 1 + Σ B_i λ²/(λ² − C_i), λ in µm, k = 0 | `B` und `C_um2` als Listen |
| `lorentz` | ε = ε∞ + Σ s_j/(E_j² − E² − iΓ_jE), Energien in eV | `eps_inf`, `oscillators`: Liste mit `strength_ev2`, `energy_ev`, `width_ev` |
| `tauc_lorentz` | ε₂ = A·E₀·C·(E − E_g)²/(E·[(E² − E₀²)² + C²E²]) oberhalb der Lücke, sonst 0; ε₁ über Kramers-Kronig | `eps_inf`, `A_ev`, `E0_ev`, `C_ev`, `Eg_ev` |
| `bruggeman` | Effektivmedium aus zwei Materialien | `a`, `b` (Materialien), `fraction_b` |
| `table` | Tabelle `wavelength_nm n k`, linear interpoliert, keine Extrapolation | `file`, wahlweise `columns`: [Spalte n, Spalte k] |
| `ellidft` | α′(λ) je Wiederholeinheit aus einer Ausgabe von `ellidft-optics`, Lorentz-Lorenz bei gegebener Dichte | `file`, `density_g_cm3` |

### Benannte Materialien

Unter `materials` lassen sich Materialien mit Namen anlegen; überall, wo ein
Material erwartet wird, darf dann der Name stehen. Ein mehrfach verwendetes
Material ist ein einziges Objekt mit einem Satz anpassbarer Zahlen. Das braucht
man, sobald dasselbe Material im Film und in einer Rauigkeitsschicht steckt.

```json
"materials": {"ps": {"model": "cauchy", "A": {"value": 1.57, "fit": true}, "B": 0.008}},
"layers": [
  {"name": "rau", "thickness_nm": {"value": 3, "fit": true},
   "material": {"model": "bruggeman", "a": "ps", "b": {"model": "constant", "n": 1.0}, "fraction_b": 0.5}},
  {"name": "film", "thickness_nm": {"value": 100, "fit": true}, "material": "ps"}
]
```

### Rauigkeit und Mischschichten

Eine raue Oberfläche wird als dünne Schicht aus Filmmaterial und Umgebung
beschrieben (üblich: je zur Hälfte), eine durchmischte Grenzfläche oder ein
poröser Film ebenso mit anderem Anteil. `fraction_b` ist anpassbar.

### Einachsige Schichten

Mit `material_extraordinary` ist eine Schicht einachsig, die optische Achse
steht senkrecht zur Oberfläche: `material` gilt dann für das Feld in der
Schichtebene (ordentlich), `material_extraordinary` für das Feld entlang der
Normalen. Das ist der Fall orientierter Polymerfilme und passt zu n_o und n_e
aus `ellidft-optics` (Tabelle mit `--table`, Spalten über `columns` wählen).
Andere Achsenrichtungen gibt es nicht.

### Dicke, durchsichtige Substrate

Mit `substrate_thickness_mm` ist das Substrat eine dicke Platte; das an der
Rückseite reflektierte Licht wird inkohärent hinzugezählt. `backside` ist das
Medium dahinter (ohne Angabe die Umgebung). Bei Glasträgern ändert das Ψ und Δ
erheblich; im Test liegt ein Modell ohne Rückseite um mehr als 0,5° daneben.

Bei teilweise depolarisiertem Licht sind Ψ und Δ hier festgelegt als
Ψ = ½ arccos N und Δ = atan2(S, C) mit den über alle Teilstrahlen gemittelten
Größen N, C, S. **Wie Ihr Gerät Ψ und Δ in diesem Fall ausgibt, ist am Handbuch
zu prüfen**; die Festlegung ist nicht aus einer Quelle übernommen.

Das Modell `ellidft` ist die Verbindung zur Rechnung: Die Dispersion kommt aus
der berechneten Polarisierbarkeit, angepasst werden nur Dichte und Dicke. Es
braucht eine Rechnung, die den ganzen gemessenen Wellenlängenbereich abdeckt
(`--wavelengths` in `ellidft-response`), weil nicht extrapoliert wird.

### Dispersionsmodell aus der Rechnung

`ellidft-optics` passt an das berechnete n(λ) ein Cauchy- und ein
Sellmeier-Modell (ein Term) an und gibt die Koeffizienten unter
`dispersion_models` aus, mit mittlerer und größter Abweichung. Sie lassen sich
als Material `cauchy` oder `sellmeier` in die Modelldatei übernehmen. Es zählen
nur Wellenlängen mit k < 10⁻⁴, mindestens drei; die Modelle gelten nur im
gerechneten Bereich. Probe: Aus Werten der Sellmeier-Kurven von Polystyrol und
PMMA kommen die Koeffizienten auf 10⁻⁷ bzw. 10⁻⁶ zurück.

## Was angepasst wird

| `residuals` | Verglichene Größen | Gewichte |
|---|---|---|
| `ncs` (Voreinstellung) | N = cos 2Ψ, C = sin 2Ψ cos Δ, S = sin 2Ψ sin Δ | keine |
| `psi_delta` | Ψ und Δ in Grad, Δ-Differenz in (−180°, 180°] | 1/σ aus den Spalten 5 und 6, sonst 1° |

`ncs` ist überall stetig, auch wo Δ um 360° springt oder für Ψ → 0 unbestimmt
wird. Minimiert wird die Summe der Fehlerquadrate mit dem
Levenberg-Marquardt-Verfahren; Ableitungen aus zentralen Differenzen.

Ausgegeben werden die Parameter mit Standardfehlern, die Korrelationsmatrix,
die mittleren und größten Abweichungen in Ψ und Δ und n, k jeder Schicht bei
den gemessenen Wellenlängen.

**Zu den Fehlerangaben.** Die Standardfehler stammen aus der Krümmung am
Minimum; die Streuung der Messpunkte wird aus dem Fit selbst geschätzt. Sie
enthalten keine Fehler des Modells und keine Fehler der festgehaltenen Größen
(Substrat, Oxiddicke, Einfallswinkel). Bei systematischen Abweichungen sind sie
zu klein. Sind Parameter nicht unabhängig bestimmbar, meldet das Programm das
und gibt keine Fehler an. Korrelationen nahe ±1 zeigen, dass die Messung zwei
Parameter kaum trennt.

## Nachweise

Ohne Referenzprogramm, alle in den kurzen Tests:

| Prüfung | Ergebnis |
|---|---|
| Fresnel-Koeffizienten gegen die Formeln mit Winkeln aus dem Brechungsgesetz | Abweichung ≤ 10⁻¹⁵ |
| Brewster-Winkel, senkrechter Einfall, Totalreflexion | erfüllt |
| Rückrechnung der Brechzahl eines absorbierenden Substrats aus Ψ und Δ | relativ ≤ 10⁻¹¹ |
| Einzelschicht gegen die ausgeschriebene Summe der Mehrfachreflexe | ≤ 10⁻¹³ |
| Glasplatte in Luft gegen die geschlossene Formel (Betrag) | ≤ 10⁻¹⁴ |
| Grenzfälle: Dicke null, Schicht aus Substrat- oder Umgebungsmaterial, Halbwellenschicht, geteilte Schicht | ≤ 10⁻¹⁴ |
| Vier Schichten: Rekursion gegen charakteristische Matrizen (zweiter, getrennt geschriebener Algorithmus) | ≤ 10⁻¹³ |
| Ausgleichsrechnung: Gerade mit geschlossener Lösung und Kovarianz, Rosenbrock, Exponentialabfall, Schranken, abhängige Parameter | erfüllt |
| Fit an exakte Daten: Dicke, A, B | 10⁻⁶ nm, 10⁻⁸ |
| Fit an verrauschte Daten (0,05°): Ergebnis innerhalb von 4 Standardfehlern, geschätzte Streuung trifft das Rauschen auf 25 % | erfüllt |
| Film von 1500 nm mit 25 Startwerten | gefunden |
| Dichte aus α′(λ) | 10⁻⁸ g/cm³ |

Zur zweiten Fassung:

| Prüfung | Ergebnis |
|---|---|
| Kramers-Kronig-Integral gegen den Realteil eines Lorentz-Oszillators | 1,5·10⁻⁹ |
| Tauc-Lorentz: numerisches ε₁ gegen die geschlossene Formel der Quelle | 2·10⁻¹² |
| Bruggeman: definierende Gleichung, Grenzfälle, vertauschte Bestandteile, verdünnter Grenzfall | 10⁻¹³ |
| Einachsige Schicht: gleiche Indizes, senkrechter Einfall, s-Polarisation unabhängig von n_e | 10⁻¹⁵ |
| Einachsige Schicht: Phasendicke π für p-Polarisation, Brewster-Winkel tan²θ = n_e²(n_o² − 1)/(n_e² − 1) | erfüllt |
| Glasplatte in Luft mit Rückseite: R = 2R₀/(1 + R₀) je Polarisation | 10⁻¹³ |
| Undurchsichtige Platte und angepasste Rückseite: nur der Vorderseitenreflex | 10⁻⁹ |
| Über das Programm: orientierter, rauer Film auf Glas, drei Parameter | 240 nm, 4,0000 nm, A = 1,58 zurückgewonnen |

Die geschlossene Tauc-Lorentz-Formel kam beim Abruf der Quelle ohne die
Wurzel im letzten Logarithmus an. Mit der Wurzel stimmt sie auf 10⁻¹² mit dem
Integral überein; verwendet wird im Programm nur das Integral.

Die Formeln der einachsigen Schicht sind aus der Dispersionsbeziehung der
Quelle und den Maxwell-Gleichungen hergeleitet; die Prüfungen sind Grenzfälle
und Folgerungen daraus, kein unabhängiger Vergleich.

Beide Algorithmen des Vorwärtsmodells stammen von derselben Hand. Ein
gemeinsamer Denkfehler, etwa in einer Konvention, fiele so nicht auf.

## Grenzen

- **Einachsig nur mit der Achse entlang der Normalen.** Keine verkippten
  Achsen, keine zweiachsigen Schichten.
- **Rauigkeit nur als Effektivmedium-Schicht.** Keine Streuung, keine
  Gradienten innerhalb einer Schicht.
- **Rückseite nur vollständig inkohärent**, isotropes Substrat. Wie ein Gerät
  mit Blenden einen Teil des Rückseitenreflexes ausblendet, ist nicht
  modelliert.
- **Tauc-Lorentz ist langsam** (ein Integral je Wellenlänge und Auswertung).
- **Ideales Gerät.** Keine Depolarisation, keine Bandbreite, keine
  Winkelstreuung, kein Fit des Einfallswinkels.
- **Keine Substratdaten im Repository.** Für Silizium, Oxid oder Glas braucht
  es Tabellen aus einer Quelle Ihrer Wahl. Die Tabelle im Beispiel ist erfunden.
- **Keine Gerätedateien.** Nur die Texttabelle oben.
