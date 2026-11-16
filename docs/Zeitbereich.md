# Zeitbereichsverfahren (Echtzeit-TDDFT)

Stand: 14.11.2026, zweite Fassung (Kompositionspropagator, Felder mit
Zeitverlauf, Antwort nach Ordnungen im Feld). Das Modul `src/realtime` propagiert die
besetzten Orbitale eines Moleküls mit abgeschlossenen Schalen nach einem kurzen
Feldpuls in der Zeit. Aus dem Dipolmoment folgt die Polarisierbarkeit
α(ω + iγ). Vom Auftraggeber am 14.11.2026 in den Umfang genommen
(`docs/Entscheidungen.md`).

**Wofür es da ist.** Erstens als Gegenprüfung des Moduls `response`: Für
schwache Pulse muss dieselbe Polarisierbarkeit herauskommen, auf einem völlig
anderen Rechenweg (Zeitschritte mit dem Potential statt linearer Gleichungen
mit dem Kern). Zweitens für breite Spektren aus einer einzigen Propagation.
Drittens, nach der Entscheidung des Auftraggebers vom 14.11.2026, als Grundlage
für die nichtlineare Optik: Die Propagation enthält alle Ordnungen im Feld.
Für n(λ) bei einzelnen Wellenlängen bleibt `ellidft-response` der billigere Weg.

## Aufruf

```
ellidft-realtime <geometrie.xyz> <basis> --gamma <g>
    [--functional <name>] [--grid <radial>,<winkel>] [--aux-basis <name>]
    [--time-step <dt>] [--time <T>] [--kick <kappa>] [--one-sided]
    [--directions <xyz>] [--propagator magnus4|composition4|midpoint]
    [--frequencies <w,...>] [--wavelengths <nm,...>] [--spectrum <von>,<bis>,<n>]
    [--compare] [--signal <datei>] [--bohr] [--charge <n>]
```

Alle Größen in atomaren Einheiten (Zeit: 1 a.u. = 24,19 as; Frequenz und
Dämpfung in Hartree). Die Ausgabe ist JSON; `polarizability_dynamic` hat
denselben Aufbau wie bei `ellidft-response`. `--compare` löst zusätzlich die
Antwortgleichungen bei denselben Frequenzen und gibt die Abweichung aus.
`--signal` schreibt dμ_i/dκ_j(t) als Texttabelle. Felder mit Zeitverlauf und
die Antwort nach Ordnungen gibt es bisher nur in der Bibliothek (unten), noch
nicht im Programm.

| Option | Voreinstellung | Bedeutung |
|---|---|---|
| `--gamma` | keine, Pflicht | Dämpfung γ; dieselbe Größe wie bei `ellidft-response` |
| `--time-step` | 0,1 | Zeitschritt |
| `--time` | 12/γ | Länge der Propagation; der abgeschnittene Rest des Signals hat die Größe e^(−γT) |
| `--kick` | 10⁻⁴ | Stärke κ des Pulses |
| `--one-sided` | aus | eine statt zwei Propagationen je Richtung |
| `--propagator` | `magnus4` | siehe unten |

## Verfahren

| Schritt | Formel | Quelle |
|---|---|---|
| Puls | H(t) = H₀ − μ̂·u F(t), F(t) = κ δ(t); daraus ψ(0⁺) = e^(−iκ u·r) ψ(0⁻) | Hauge et al. 2023, Gl. 3, 4, 21; der zweite Teil ist hier hergeleitet |
| Fock-Matrix | F = h + J[Re P] + V_xc[Re P] − ½ (c_x K[P] + c_lr K_LR[P]), P = 2CC⁺ komplex | wie im SCF; der Imaginärteil von P geht nur in den Austausch ein (hier hergeleitet) |
| Schritt, `midpoint` | C(t+Δt) = e^(−iΔt F_M) C(t), F_M = ½ (F(t) + F(t+Δt)) | Herbert 2022, Gl. 78–80 |
| Schritt, `magnus4` | C(t+Δt) = e^(−iH) C(t), H = Δt/2 (F₁+F₂) + i √3 Δt²/12 [F₁, F₂], F₁, F₂ an den Gauß-Legendre-Punkten des Schritts | Blanes et al. 2008, Gl. 252–254 |
| Schritt, `composition4` | drei Mittelpunktsschritte der Längen w₁Δt, w₀Δt, w₁Δt mit w₁ = 1/(2 − 2^(1/3)), w₀ = 1 − 2w₁ < 0 | Wikipedia, „Symplectic integrator“ (Dreisprung) und „Leapfrog integration“ (Koeffizienten nach Yoshida) |
| Feld mit Zeitverlauf | H(t) = H₀ − μ̂·F(t); im Mittelpunktsschritt am Mittelpunkt, im Magnus-Schritt an den Gauß-Punkten ausgewertet | Hauge et al. 2023, Gl. 3, 4 |
| Polarisierbarkeit | α_ij(ω+iγ) = ∫₀ᵀ dμ_i/dκ_j(t) e^(i(ω+iγ)t) dt, Simpson-Regel | Hauge et al. 2023, Gl. 8; Herbert 2022, Gl. 82 |

Die Fock-Matrix innerhalb eines Schritts hängt von den Orbitalen an seinem
Ende ab. Das ist hier gelöst und nicht aus einer Quelle übernommen. Der
Mittelpunktsschritt wird wiederholt, bis die Fock-Matrix am Ende zu den
Orbitalen passt. Im Magnus-Schritt wird F(t)
durch das Polynom durch die Fock-Matrizen am Schrittende (Versuchswert), am
Schrittanfang und an zwei früheren Zeiten dargestellt, und der Schritt wird
wiederholt, bis sich die Fock-Matrix am Ende nicht mehr ändert. Direkt nach dem
Puls gibt es keine früheren Zeiten; die ersten drei Schritte werden deshalb
gemeinsam nachverfeinert. Ohne diesen Anlauf konvergierte das Verfahren nur
mit dritter Ordnung (gemessen).

Gerechnet wird in der Basis der Grundzustandsorbitale; die Exponentialfunktion
entsteht durch Diagonalisierung und ist für jeden Zeitschritt unitär.

Mit zwei Propagationen je Richtung (+κ und −κ, Voreinstellung) ist die Antwort
die halbe Differenz geteilt durch κ. Das entfernt alle geraden Ordnungen im
Feld und die Restbewegung eines nicht ganz konvergierten Grundzustands.

## Nachweise

Vergleich mit dem Antwortcode: größte Abweichung eines Tensorelements bei
ω = 0, 0,08, 0,2 und 0,45 Hartree, bezogen auf das größte Element.
Vorgeschlagenes Kriterium: höchstens 10⁻⁴.

| Fall | Einstellungen | Abweichung | Lauf |
|---|---|---|---|
| H₂, cc-pVDZ, HF | Δt 0,2, T 90, γ 0,15, zwei Richtungen | 1,3·10⁻⁵ | `slow` |
| H₂, cc-pVDZ, PBE | wie oben | 4,3·10⁻⁶ | `slow` |
| H₂, cc-pVDZ, PBE0 | wie oben | 5,8·10⁻⁶ | `slow` |
| H₂, cc-pVDZ, CAM-B3LYP | wie oben | 8,3·10⁻⁶ | `slow` |
| H₂O gedreht, STO-3G, HF, voller Tensor | Δt 0,1, T 45, γ 0,3 | 1,9·10⁻⁵ | kurz |
| H₂O gedreht, STO-3G, PBE | Δt 0,1, T 90, γ 0,15 | 1,7·10⁻⁵ | `slow` |
| H₂O gedreht, STO-3G, PBE0 | wie oben | 1,9·10⁻⁵ | `slow` |
| H₂, verkleinerte Fälle (HF, PBE, PBE0) | Δt 0,2, T 45, γ 0,3 | ≤ 1,9·10⁻⁵ | kurz |

Die Dichtefunktionale laufen in den Tests auf kleinen Gittern (30, 110) bzw.
(20, 50). Beide Rechenwege nutzen dasselbe Gitter, die Übereinstimmung hängt
also nicht von dessen Güte ab.

Weitere Prüfungen:

| Prüfung | Ergebnis |
|---|---|
| Fehler bei halbem Zeitschritt (H₂, HF, Δt 0,4 gegen 0,2): `midpoint` | Faktor 4,0 (zweite Ordnung) |
| dasselbe mit `magnus4` | Faktor 15,6 (vierte Ordnung wäre 16) |
| dasselbe mit `composition4` | Faktor 15,2 |
| `composition4` gegen den Antwortcode: H₂ (Δt 0,2), H₂O (Δt 0,1) | 1,1·10⁻⁵ und 4,5·10⁻⁵ |
| Fock-Matrix und Energie der Grundzustandsdichte gegen das SCF (HF, PBE, PBE0, CAM-B3LYP) | 10⁻⁹ und 10⁻¹⁰ Hartree |
| Ohne Puls: Dipol und Energie zeitlich konstant | 10⁻⁹ und 10⁻¹⁰ |
| Fock-Matrix einer komplexen Dichte hermitesch | 10⁻¹² |
| Orthonormalität der Orbitale, also die Elektronenzahl | ≤ 2·10⁻¹² nach 900 Schritten |
| Energie nach schwachem Puls | Schwankung ≤ 1,1·10⁻⁹ Hartree |
| Transformation eines bekannten Signals (zwei Sinuslinien) | relativ ≤ 10⁻⁵ |

## Energieerhaltung und Wahl des Propagators

Nach dem Puls ist die Energie eine Erhaltungsgröße. Gemessen an Wasser
(STO-3G, HF) nach einem starken Puls (κ bis 0,05, zugeführte Energie
2,2·10⁻³ Hartree):

| Propagator | Δt | größte Abweichung der Energie |
|---|---|---|
| `midpoint` | 0,1 | 4·10⁻¹³ |
| `composition4` | 0,2 bis 0,05 | ≤ 4·10⁻¹¹ (HF), ≤ 8·10⁻¹⁰ (PBE0) |
| `magnus4` | 0,1 | 1,6·10⁻⁴ (über T = 20) |
| `magnus4` | 0,05 | 7·10⁻⁵ (über T = 20) |
| `magnus4` | 0,025 | 1,7·10⁻⁶ (über T = 10) |
| `magnus4` | 0,0125 | 6·10⁻⁸ (über T = 10) |

Der Mittelpunktsschritt erhält die Hartree-Fock-Energie exakt: Sie ist
quadratisch in der Dichtematrix, damit gilt E(t+Δt) − E(t) = Sp[F_M (ΔP)], und
das verschwindet, weil der Schritt mit F_M vertauscht. `composition4` erbt das
und ist zugleich von vierter Ordnung. Mit Dichtefunktional gilt es nicht mehr
exakt, die Abweichung bleibt aber klein.

Der Magnus-Schritt dieses Moduls erhält die Energie nicht, und die Abweichung
wächst mit der Zeit. Er benutzt Fock-Matrizen früherer Schritte; die schnellen
Anteile aus dem 1s-Orbital des Sauerstoffs (Frequenzen um 20 Hartree) gibt das
Polynom erst bei kleinen Schritten richtig wieder. Bei schwachen Pulsen bleibt
die Schwankung absolut winzig (≤ 10⁻⁹), relativ zur zugeführten Energie liegt
sie aber auch dort im Prozentbereich.

| Zweck | Empfehlung | Fock-Matrizen je Schritt |
|---|---|---|
| lineare Antwort, schwacher Puls | `magnus4` (Voreinstellung) | etwa 3 |
| starke Pulse, Felder mit Zeitverlauf, nichtlineare Antwort | `composition4` | etwa 12, bei enger Toleranz oder starkem Puls bis 22 |
| Kontrolle | `midpoint` | etwa 3 |

Die Zahlen sind gemessen (H₂ und Wasser, schwache Pulse, Voreinstellung der
Toleranz). Bei gleichem Zeitschritt ist `composition4` in der linearen Antwort
zwei- bis dreimal ungenauer als `magnus4` und etwa viermal teurer.

## Felder mit Zeitverlauf und nichtlineare Antwort

In der Bibliothek (`RealTimeOptions::field`) kann während der Propagation ein
homogenes Feld F(t) mit beliebigem Zeitverlauf wirken, mit derselben
Vorzeichenwahl wie das statische Feld im SCF. `dipole_orders` zerlegt das
Dipolmoment nach Ordnungen in der Feldstärke: Für das Feld λ·f(t) ist

    μ(t; λ) = μ₀ + λ μ₁(t) + λ² μ₂(t) + λ³ μ₃(t) + …

μ₁ ist die lineare Antwort auf das Profil f, μ₂ und μ₃ sind die Antworten
zweiter und dritter Ordnung (im statischen Grenzfall α·u, β:uu/2 und
γ:uuu/6). Sie entstehen aus vier Propagationen mit λ = ±a und ±2a durch
Differenzen; der Restfehler ist von der Ordnung a⁴ bzw. a² mal der fünften
Ordnung.

Nachweise, jeweils gegen einen unabhängigen Rechenweg:

| Prüfung | Ergebnis |
|---|---|
| Aufgenommene Energie eines Pulses (H₂, F₀ = 0,02, ω = 0,3, Dauer 40) gegen die Arbeit −∫F′(t) μ(t) dt | 2,2234·10⁻³ Hartree; Differenz 1,3·10⁻⁹ (`magnus4`, Δt 0,1), 2,1·10⁻⁸ (`composition4`, Δt 0,2), 2,4·10⁻⁷ (`midpoint`, Δt 0,1) |
| Langsam eingeschaltetes Feld 0,1 an H₂ gegen das SCF im statischen Feld | μ = 0,634834 (linear wären 0,642018, also 1,1 % nichtlinearer Anteil); größte Abweichung 7·10⁻⁶ |
| Wasser, Feld entlang einer allgemeinen Richtung, a = 0,02: μ₁ gegen α·u aus dem Antwortcode | relativ 1·10⁻⁶ |
| μ₂ gegen dieselben Differenzen aus vier SCF-Rechnungen im Feld | relativ 4·10⁻⁵ |
| μ₃ ebenso | relativ 2·10⁻³ |

Die dritte Ordnung ist eine Differenz von Dipolmomenten geteilt durch 6a³;
deshalb ist sie am empfindlichsten gegen die Iterationstoleranz und gegen das
nicht ganz adiabatische Einschalten.

**Noch nicht gebaut:** frequenzabhängige nichtlineare Größen (Frequenz-
verdopplung β(−2ω; ω, ω), Frequenzverdreifachung, intensitätsabhängige
Brechzahl) aus μ₂(t) und μ₃(t). Dafür fehlen die Wahl der Pulsform, die
Auswertung im Frequenzbereich und ein unabhängiger Vergleich; bisher sind die
höheren Ordnungen nur im statischen Grenzfall geprüft.

## Wahl von Zeitschritt und Dauer

- **Zeitschritt.** Er muss die höchste Frequenz auflösen, die interessiert.
  Für Valenzspektren genügten in den Tests 0,2 (H₂) und 0,1 (Wasser, mit
  Rumpforbital). Bei Wasser mit Δt = 0,2 lag die Abweichung bei 2,7·10⁻⁴.
- **Dauer.** e^(−γT) ist der abgeschnittene Rest. Kleine Dämpfung (scharfe
  Linien) verlangt lange Läufe: γ = 0,1 eV braucht T ≈ 3300 a.u.
- **Kosten.** Etwa drei Fock-Matrizen je Schritt mit `magnus4` (zwölf mit
  `composition4`), mal Zahl der Richtungen, mal zwei bei zweiseitigem Puls.

## Grenzen

- **Nur abgeschlossene Schalen**, dieselben Funktionale wie im SCF,
  adiabatische Näherung (das Potential hängt nur von der momentanen Dichte ab).
- **Felder mit Zeitverlauf nur in der Bibliothek**, nicht im Programm.
- **Kein Verlassen der Basis.** Ionisation wird nicht beschrieben, es gibt
  keine absorbierenden Ränder; für starke Felder ist das Modul nicht geprüft.
- **Höhere Ordnungen nur im Zeitbereich und nur statisch geprüft.** Die
  Umrechnung in frequenzabhängige Hyperpolarisierbarkeiten fehlt.
- **Keine Signalverlängerung** (Padé oder Anpassung), die kürzere Läufe
  erlaubte.
- **Diagonalisierung in jedem Schritt.** Für große Basen wäre das zu ersetzen.
- **Kein Referenzprogramm.** Der Nachweis ist die Übereinstimmung mit dem
  eigenen Antwortcode, der seinerseits gegen PySCF geprüft ist. Rechnungen an
  größeren Molekülen gibt es noch nicht.
