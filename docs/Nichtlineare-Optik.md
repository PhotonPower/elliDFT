# Nichtlineare Optik aus der Zeitpropagation

Stand: 15.11.2026, zweite Fassung (ungeordnete Materialien, Programm). Der Auftraggeber hat am 14.11.2026 drei
Größen genannt: Zwei-Photonen-Absorption, Frequenzverdopplung und Kerr-Effekt
(`docs/Entscheidungen.md`). Das Modul `src/realtime/nonlinear` liefert sie für
ein einzelnes Molekül. Für ungeordnete Materialien gibt es Orientierungsmittel, n₂ des Films und
das Programm `ellidft-nonlinear` (Abschnitt unten).

**Alle Zahlen unten sind Prüfwerte** an H₂ und Wasser in kleinen Basen ohne
diffuse Funktionen. Für Hyperpolarisierbarkeiten sind solche Basen ungeeignet
(γ von H₂ kommt hier sogar negativ heraus).

## Festlegungen

Für statische Felder gilt die Taylor-Reihe

    μ_i = μ⁰_i + α_ij F_j + ½ β_ijk F_j F_k + ⅙ γ_ijkl F_j F_k F_l .

Die frequenzabhängigen Größen sind so festgelegt, dass sie für ω → 0 in diese
Ableitungen übergehen. Das ist eine Festlegung dieses Programms; in der
Literatur sind mehrere üblich, die sich um Zahlenfaktoren unterscheiden. Für
das Feld F(t) = F₀ + F₁ cos ωt:

| Ordnung | Anteil | Größe | Effekt |
|---|---|---|---|
| 2 | cos 2ωt | ¼ β(−2ω; ω, ω) F₁² | Frequenzverdopplung |
| 2 | cos ωt | β(−ω; ω, 0) F₁F₀ | elektrooptischer (Pockels-)Effekt |
| 3 | cos ωt | ⅛ γ(−ω; ω, ω, −ω) F₁³ | optischer Kerr-Effekt |
| 3 | sin ωt | Imaginärteil derselben Größe | Zwei-Photonen-Absorption |
| 3 | cos ωt | ½ γ(−ω; ω, 0, 0) F₁F₀² | elektrooptischer Kerr-Effekt (statisches Feld) |
| 3 | cos 3ωt | 1/24 γ(−3ω; ω, ω, ω) F₁³ | Frequenzverdreifachung |

„Kerr“ habe ich als optischen Kerr-Effekt gelesen; den elektrooptischen gibt es
zusätzlich, er dient auch als Probe.

## Verfahren

**Harmonische Analyse** (`harmonic_response`). Das Feld λ·s(t)·(F₀ + F₁ cos ωt)
wird langsam eingeschaltet. Vier Propagationen mit λ = ±a, ±2a trennen die
Ordnungen (`dipole_orders`), eine Ausgleichsrechnung im Fenster nach dem
Einschalten trennt die Vielfachen von ω. Ausgegeben wird auch, was die
Anpassung übrig lässt (Einschwingreste, Resonanzen).

**Zwei-Photonen-Absorption** (`two_photon_absorption`). Ein Puls
F₀ sin²(πt/T) cos ω(t − T/2) hinterlässt die Energie W = W₂F₀² + W₄F₀⁴; W₄ ist
der Zwei-Photonen-Anteil. Dafür ist der energieerhaltende Propagator
`composition4` gedacht. Der Querschnitt ist hier über die Zahl absorbierter
Photonen festgelegt, dN/dt = σ₂ (I/ω)², mit der Intensität I = cF²/(8π):

    σ₂ = W₄ ω (8π/c)² / ∫ sin⁸(πt/T) dt ,   1 a.u. = 1,897 GM.

Das ist ein Mittel über das Zwei-Photonen-Spektrum des Pulses (Breite etwa
4π/T), für festgehaltene Moleküle und entlang der Feldrichtung polarisiertes
Licht.

## Nachweise

Unabhängige Rechenwege: der Antwortcode für das Molekül im statischen Feld
(dα(ω)/dF = β(−ω; ω, 0), d²α(ω)/dF² = γ(−ω; ω, 0, 0)); die Entwicklung bei
kleinen Frequenzen (für alle Felder entlang einer Richtung bleibt wegen der
Vertauschungssymmetrie ein einziger quadratischer Term, proportional zur Summe
der Frequenzquadrate; hier hergeleitet); und Anregungsenergien dunkler
Zustände aus dem Antwortcode.

| Prüfung | Zeitpropagation | unabhängig |
|---|---|---|
| Wasser (STO-3G, HF), ω = 0,08: α(ω) | 3,574116 | 3,574113 |
| Pockels-Effekt β(−ω; ω, 0) | −6,28082 | −6,28052 |
| Frequenzverdopplung β(−2ω; ω, ω); statisch −6,17134 | −6,51058 | Dispersion 3,11-mal die des Pockels-Effekts (Grenzwert 3) |
| H₂ (cc-pVDZ, HF), ω = 0,1: elektrooptischer Kerr-Effekt | −45,25 | −45,43 |
| optischer Kerr-Effekt γ(−ω; ω, ω, −ω); statisch −44,25 | −46,26 | Dispersion 1,70-mal die des elektrooptischen (Grenzwert 2) |
| Frequenzverdreifachung | −52,40 | Dispersion 6,9-mal (Grenzwert 6) |
| Zwei-Photonen-Absorption von H₂, dunkler Zustand bei 0,7835 Hartree | 22,5 GM bei ω = 0,3917 (halbe Energie); −7·10⁻⁵ GM bei ω = 0,3317 | Lage des Zustands aus dem Antwortcode |

Was das trägt und was nicht:

- **Zweite Ordnung: belastbar.** Der Pockels-Effekt stimmt auf 5·10⁻⁵, das
  Dispersionsverhältnis liegt beim Grenzwert.
- **Dritte Ordnung: auf etwa 0,4 % genau** (elektrooptischer Kerr-Effekt gegen
  den unabhängigen Weg). Die Dispersion ist hier nur 3 bis 5 % der Größe; das
  Verhältnis für den optischen Kerr-Effekt (1,70 statt 2) bestätigt deshalb nur
  grob. Die dritte Ordnung ist eine Differenz geteilt durch 6a³ und empfindlich
  gegen die fünfte Ordnung (Anteil 5a²) und gegen Einschwingreste.
- **Zwei-Photonen-Absorption: die Lage ist geprüft, der Betrag nicht.** Die
  Absorption erscheint bei der halben Energie des dunklen Zustands und
  verschwindet 0,06 Hartree daneben. Für den Zahlenwert des Querschnitts gibt
  es keinen unabhängigen Vergleich, und die Festlegung (je Photon, feste
  Orientierung, Pulsmittel) ist meine. Für die Lichtgeschwindigkeit steht
  137,036 im Code; diese Ziffern stammen nicht aus einer eingesehenen Quelle.

## Ungeordnete Materialien (Photopolymere)

Entscheidungen des Auftraggebers vom 15.11.2026: Anwendung sind Photopolymere,
gemeint ist der optische Kerr-Effekt, und der Zwei-Photonen-Querschnitt wird
wie vorgeschlagen festgelegt.

| Punkt | Festlegung |
|---|---|
| Orientierung | Mittel über alle Orientierungen, linear polarisiertes Licht, alle Felder parallel |
| Zwei-Photonen-Querschnitt | je absorbiertem Photon, ein Strahl, in GM |
| Linienform | Lorentz-Linie, halbe Breite als Eingabe (Voreinstellung 0,1 eV) |
| Film | Lorentz-Lokalfeld (n² + 2)/3 für jedes Feld; nur der elektronische Anteil |

**Orientierungsmittel.** u·γ·uuu und |u·S·u|² sind Polynome vierten Grades in
der Feldrichtung u. Sieben Richtungen (drei Achsen mit Gewicht 2/15, vier
Raumdiagonalen mit 3/20) mitteln sie exakt.

**Zwei-Photonen-Stärke je Zustand** (`two_photon_strengths`). Für jeden
Zustand mit Energie ω_f aus dem Antwortcode läuft ein Puls bei ω_f/2. Aus dem
Anteil vierter Ordnung der hinterlassenen Energie folgt die Stärke |e·S·e|²
(Störungsrechnung zweiter Ordnung, hier hergeleitet); überlappen Zustände
innerhalb der Pulsbreite, wird das lineare Gleichungssystem gelöst. Das
Verfahren verlangt Photonen unterhalb jeder Ein-Photonen-Linie: Zustände, deren
halbe Energie näher als 0,05 Hartree an der ersten erlaubten Anregung liegt,
lässt das Programm aus.

**Querschnitt.** σ₂(ω) = 4π³ α² ω² Σ_f Stärke_f · g(2ω − ω_f) mit normierter
Lorentz-Linie g. Hergeleitet aus der Übergangsrate und dem Photonenfluss; die
Form stimmt mit N·π³·α·a₀⁵·ω²/c und N = 4 überein, die eine Websuche für
Beerepoot et al., PCCP 17, 19306 (2015) nennt. Die Arbeit selbst konnte nicht
gelesen werden; ein Abgleich steht aus.

**n₂ des Films** (`nonlinear_index`). n = n₀ + n₂I mit
n₂ = 2π² N L⁴ γ/(n₀² c) in atomaren Einheiten, L = (n₀² + 2)/3; hergeleitet
aus Lorentz-Lorenz mit der Polarisierbarkeit α + γF²/8 bei der Lichtfrequenz.
Für ein Polymer ist γ je Wiederholeinheit einzusetzen, also die Differenz
aufeinanderfolgender Oligomere wie bei α. `ellidft-nonlinear` nimmt das
gerechnete Molekül als Einheit; für Polymere bildet `ellidft-kerr-film` die
Differenzen:

    ellidft-kerr-film --density 1.18 --molar-mass 100.12 --optics pmma.optics.json \
        mma_1.json mma_2.json mma_3.json mma_4.json

Eingaben sind Ausgaben von `ellidft-nonlinear` (Aufgabe `kerr`) für 1, 2, 3, …
Einheiten bei denselben Wellenlängen. Ausgegeben werden die Reihe der
Inkremente, ihre Änderung zwischen den beiden letzten Gliedern und n₂ aus dem
letzten Inkrement. Die lineare Brechzahl kommt aus der Ausgabe von
`ellidft-optics` oder als feste Zahl (`--refractive-index`). Geprüft mit
konstruierten Eingaben gegen die Handrechnung; an echten Oligomeren nicht.

**Frequenzverdopplung.** In einem ungeordneten Material verschwindet sie im
Volumen. Das Programm gibt nur die molekulare Größe entlang einer Richtung aus.

### Zwei-Photonen-Lithographie

Anwendung laut Auftraggeber (15.11.2026). Welche Größe wofür gebraucht wird:

| Frage | Größe | Woher |
|---|---|---|
| Wie stark absorbiert der Initiator bei der Schreibwellenlänge? | σ₂(λ) in GM | `two_photon` am Initiator |
| Bleibt das Harz bei dieser Wellenlänge durchsichtig? | Ein-Photonen-Anregungen, k(λ) | `ellidft-response` am Initiator und am Monomer |
| Wie bricht das Harz, vor und nach der Belichtung? | n(λ) von Monomer und Polymer | Oligomer-Reihe und `ellidft-optics` |
| Ändert die hohe Intensität im Fokus die Brechzahl? | n₂ | `kerr` mit `film` am Monomer oder an der Wiederholeinheit |

Die Zwei-Photonen-Auswertung setzt voraus, dass die Schreibwellenlänge
unterhalb jeder Ein-Photonen-Absorption liegt. Für Initiatoren, die im UV
absorbieren, und Schreiben im Nahinfraroten ist das erfüllt. Eine Vorlage für
die Steuerdatei liegt in `examples/nonlinear/zwei_photonen_lithographie.json`;
sie ist nicht gerechnet.

Nicht enthalten: was nach der Absorption geschieht (Radikalbildung,
Polymerisationsschwelle, Diffusion), thermische Effekte und die Form des Fokus.

### Programm

    ellidft-nonlinear <steuerdatei.json>

```json
{
  "geometry": "molekuel.xyz", "basis": "def2-SVPD", "functional": "CAM-B3LYP",
  "grid": [75, 302], "aux_basis": "def2-universal-JKFIT",
  "time_step": 0.2, "fock_tolerance": 1e-11,
  "kerr": {"wavelengths_nm": [780], "amplitude": 0.01, "ramp_time": 300, "periods": 2,
           "film": {"density_g_cm3": 1.18, "molar_mass_g_mol": 100.12, "refractive_index": 1.49}},
  "two_photon": {"states": 6, "amplitude": 0.01, "duration": 300, "half_width_ev": 0.1,
                 "wavelengths_nm": {"from": 600, "to": 1000, "step": 10}},
  "second_harmonic": {"wavelengths_nm": [1064], "direction": [0, 0, 1]}
}
```

Fertige Teile (jede Wellenlänge von `kerr` und `second_harmonic`, die Stärken
von `two_photon`) schreibt das Programm sofort nach
`<steuerdatei>.partial.json`. Ein späterer Lauf mit demselben Molekül,
derselben Methode und denselben Einstellungen übernimmt sie von dort;
`"restart": false` schaltet das ab. Innerhalb einer Wellenlänge (sieben
Richtungen) wird nicht zwischengespeichert.

Jede der drei Aufgaben ist freigestellt. Ergebnis als JSON, Protokoll auf der
Fehlerausgabe, mit Warnungen bei schlechter Anpassung und bei negativen
Stärken. Ein kleinstes Beispiel liegt in `examples/nonlinear`.

### Weitere Nachweise

| Prüfung | Zeitpropagation | unabhängig |
|---|---|---|
| Ein-Photonen-Absorption eines Pulses (H₂, STO-3G) gegen Oszillatorstärke mal Pulsspektrum | 616,12 | 615,93 |
| Absorption eines dunklen Zustands im statischen Feld 0,02 (H₂, cc-pVDZ) gegen den Antwortcode im Feld | 8,33980 | 8,33941 |
| Kerr-Mittel eines gedrehten H₂ (6-31G*, nur s-Funktionen) gegen ein Fünftel des Werts entlang der Bindung | −24,819 | −24,802 |
| Mittlere Zwei-Photonen-Stärke, H₂ gedreht gegen ausgerichtet (cc-pVDZ) | 99,68 | 99,55 |
| Pockels-Effekt an Wasser mit PBE0 und CAM-B3LYP gegen den Antwortcode im Feld | −5,65072; −5,86315 | −5,65044; −5,86283 |
| n₂-Formel gegen Lorentz-Lorenz mit feldabhängiger Polarisierbarkeit | relativ 10⁻⁵ | |
| Sieben Richtungen gegen die geschlossenen Mittelwerte | 10⁻¹⁴ | |

Die Absorption im statischen Feld ist der Grenzfall der
Zwei-Photonen-Absorption, in dem eines der Photonen die Frequenz null hat. Sie
prüft den Betrag der Absorption zweiter Ordnung gegen einen unabhängigen Weg.
Für zwei gleiche Photonen gibt es weiter keinen unabhängigen Vergleich des
Betrags; dort stützt sich der Vorfaktor auf die Herleitung.

Der Spitzenwert für den dunklen Zustand von H₂ mit 0,1 eV halber Breite ist
16,6 GM. Das ist ein Prüfwert.

## Grenzen

- Nur molekulare Größen, nur Hartree-Fock geprüft (Dichtefunktionale laufen
  über denselben Code, sind hier aber nicht gerechnet).
- Unterhalb der ersten Anregung. Nahe Resonanzen gibt es keinen eingeschwungenen
  Zustand; die harmonische Analyse meldet das nur über den Rest der Anpassung.
- Kosten: vier Propagationen je Frequenz und Feldrichtung. Die Tests brauchten
  16 000 (Wasser) bis 65 000 (H₂, dritte Ordnung) Fock-Matrizen.
- Kein Orientierungsmittel, keine Tensorkomponenten außer den gerechneten
  Richtungen, keine Schwingungsbeiträge.
- Zwei-Photonen-Absorption nur unterhalb der ersten erlaubten Anregung. Liegt
  die halbe Energie eines Zustands darüber, liefert die Auswertung Unsinn
  (negative Stärken); das Programm lässt solche Zustände aus und warnt.
- Nach einem Suchergebnis unterschätzt CAM-B3LYP Zwei-Photonen-Stärken um das
  1,5- bis 3-Fache gegenüber genaueren Verfahren. Nicht selbst geprüft.
- Kosten der Mittelwerte: sieben Richtungen. Der Test der Orientierungsmittel
  an H₂ brauchte 510 000 Fock-Matrizen.
