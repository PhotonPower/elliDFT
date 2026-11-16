# Eingesehene Quellen

Regel 9 des Auftrags: Jede implementierte Formel trägt im Code einen Kommentar
mit Quelle, und zitiert wird nur, was tatsächlich eingesehen wurde. Diese Liste
hält fest, welche Quellen eingesehen wurden und wofür sie dienen. Die
Lehrbücher und Originalarbeiten aus Abschnitt 10 des Auftrags waren nicht
zugänglich und werden deshalb im Code nicht zitiert.

| Kürzel im Code | Quelle | Eingesehen | Verwendet für |
|---|---|---|---|
| Hernández Vera 2025 | M. Hernández Vera, „FSIM: A Pedagogical and Extensible HPC Framework for the Hartree–Fock Method“, arXiv:2510.24968v1 (2025) | 14.11.2026, PDF S. 3–17 | Hermite-Gauß-Funktionen und Entwicklungskoeffizienten (Gl. 11–20), Überlappung (Gl. 35, 36), kinetische Energie (Gl. 38–43), Kernanziehung und Hermite-Coulomb-Integrale (Gl. 46–51), Zweielektronenintegrale (Gl. 53–64), RHF-Gleichungen (Gl. 21–24), SCF-Ablauf mit symmetrischer Orthogonalisierung und Energieausdruck (Abschnitt 3.2), DIIS (Gl. 65–69) |
| Vikhamar-Sandberg & Repisky 2025 | R. Vikhamar-Sandberg, M. Repisky, „Efficient Boys function evaluation using minimax approximation“, arXiv:2512.10059v2 (2025) | 14.11.2026, PDF S. 1–5 | Definition der Boys-Funktion (Gl. 3), asymptotische Formel und ihr Fehler (Gl. 11, 12), Aufwärts- und Abwärtsrekursion (Gl. 14, 15) |
| Ribaldone & Desmarais 2025 | C. Ribaldone, J. K. Desmarais, „Spherical to Cartesian Coordinates Transformation for Solid Harmonics Revisited: Construction of the Hartree Potential“, arXiv:2412.16733v2 (2025) | 14.11.2026, PDF S. 1–6 | Reelle Kugelflächenfunktionen als Polynome (Gl. 1, 2, 5, 14), Koeffiziententabellen bis l = 6 (Tabellen I und II) als Testwerte |
| Lebedev-Laikov.F | Fortran-Routine der Lebedev-Gitter (Laikov, van Wüllen), über CCL verbreitet; nicht im Repository, Herkunft und Prüfsumme in `third_party/lebedev/README.md` | 14.11.2026, Kopfkommentare und Aufbau der Routinen | Stützstellen und Gewichte der 32 Winkelgitter, Bahnen der Oktaedergruppe (Kopf von GEN_OH) |
| PySCF-Quelltext | PySCF 2.15.0, `dft/radi.py`, `dft/gen_grid.py`, `dft/LebedevGrid.py`, `data/radii.py` | 14.11.2026 | Als Beschreibung des Gitters der Referenz gelesen, nicht als Code übernommen: Treutler-Ahlrichs-Radialgitter (M4) mit Parametern ξ, Becke-Funktion, Größenkorrektur nach Treutler/Ahlrichs, Bragg-Radien, Nenngrade der Lebedev-Gitter. Die dort genannten Originalarbeiten (Becke 1988, Treutler und Ahlrichs 1995) wurden nicht eingesehen. |
| Wikipedia, Legendre | „Associated Legendre polynomials“, en.wikipedia.org | 14.11.2026 | Rekursionen und Orthogonalität der zugeordneten Legendre-Funktionen, nur im Test der Lebedev-Gitter |
| libxc vwn.mpl | libxc 7.0.0, `maple/vwn.mpl` (gitlab.com/libxc/libxc) | 14.11.2026 | Parameter und Formel der VWN5-Korrelation für die eigene LDA |
| Wikipedia, LDA | „Local-density approximation“, en.wikipedia.org | 14.11.2026 | Slater-Austausch −¾(3/π)^{1/3}∫ρ^{4/3}, Definition von r_s |
| refractiveindex.info | M. N. Polyanskiy, Sci. Data 11, 94 (2024); Seiten „PS“ und „PMMA“ mit Sellmeier-Anpassungen an N. Sultanova, S. Kasarova, I. Nikolov, Acta Phys. Pol. A 116, 585 (2009) | 14.11.2026, Datenbankseiten; die Originalarbeit nicht | Gemessene Brechzahlen von Polystyrol und PMMA als Richtwerte für M6 (`tests/reference/gen_m6.py`) |
| Wikipedia, Fresnel | „Fresnel equations“, en.wikipedia.org | 14.11.2026 | Reflexions- und Transmissionskoeffizienten r_s, r_p, t_s, t_p für nichtmagnetische Medien und ihre Vorzeichenwahl, Brechungsgesetz, Brewster-Winkel, Reflexionsgrad bei senkrechtem Einfall (`src/ellipsometry/stack.cpp` und Tests) |
| Wikipedia, Ellipsometrie | „Ellipsometry“, en.wikipedia.org | 14.11.2026 | Definition ρ = r_p/r_s = tan Ψ e^{iΔ} |
| Wikipedia, Transfermatrix | „Transfer-matrix method (optics)“, en.wikipedia.org | 14.11.2026 | Reflexion einer Glasplatte in Luft bei senkrechtem Einfall, nur als Testwert; N = n + ik mit abklingender Welle für k > 0 |
| Wikipedia, Levenberg-Marquardt | „Levenberg–Marquardt algorithm“, en.wikipedia.org | 14.11.2026 | Gedämpfte Normalgleichungen mit der Diagonale von JᵀJ, Regel für den Dämpfungsparameter (Faktor 2 hoch, Faktor 3 herunter) (`src/ellipsometry/least_squares.cpp`) |
| Herbert 2022 | J. M. Herbert, „Density Functional Theory for Electronic Excited States“, arXiv:2204.10135v1 (2022) | 14.11.2026, Abschnitt 4.1, PDF S. 34–36 | Echtzeit-TDDFT: Propagator U = exp(−iΔt F) am Mittelpunkt des Schritts in der orthonormalen Basis (Gl. 78–80), Polarisierbarkeit als ∂μ(ω)/∂E(ω) aus dem Dipolsignal (Gl. 81–83), übliche Zeitschritte (`src/realtime`) |
| Hauge et al. 2023 | E. Hauge, H. E. Kristiansen, L. Konecny, M. Kadek, M. Repisky, T. B. Pedersen, „Cost-Efficient High-Resolution Linear Absorption Spectra Through Extrapolating the Dipole Moment from Real-Time Time-Dependent Electronic-Structure Theory“, arXiv:2307.01511v2 (2023) | 14.11.2026, PDF S. 6–8, 12–13, 16–17 | Feldkopplung V(t) = −μ̂·u F(t) (Gl. 3, 4), Delta-Puls F(t) = κδ(t) (Gl. 21), gedämpfte Fourier-Transformation als Laplace-Transformation mit ω + iγ (Gl. 8), lineare Antwortfunktion (Gl. 22) als Vorzeichenprobe, Bedingung Δt < π/ω_max, dort verwendete Pulsstärken und Zeitschritte |
| Blanes et al. 2008 | S. Blanes, F. Casas, J. A. Oteo, J. Ros, „The Magnus expansion and some of its applications“, arXiv:0810.5488v1 (2008) | 14.11.2026, Abschnitt um Gl. 243–254 | Magnus-Integrator vierter Ordnung mit zwei Gauß-Legendre-Punkten (Gl. 252–254) (`src/realtime/propagator.cpp`) |
| Wikipedia, Dreisprung | „Symplectic integrator“ und „Leapfrog integration“, en.wikipedia.org | 14.11.2026 | Aus einem symmetrischen Verfahren zweiter Ordnung entsteht durch drei Schritte eines vierter Ordnung („triple jump“); Koeffizienten w₁ = 1/(2 − 2^(1/3)), w₀ = −2^(1/3)/(2 − 2^(1/3)) (Yoshida). Die Originalarbeit von Yoshida (1990) wurde nicht eingesehen. (`src/realtime/propagator.cpp`) |
| Wikipedia, nichtlineare Optik | „Two-photon absorption“, „Hyperpolarizability“, „Hartree atomic units“, en.wikipedia.org | 14.11.2026 | Einheit Göppert-Mayer (10⁻⁵⁰ cm⁴ s/Photon), Auswahlregel (ein- und zweiphotonenerlaubte Übergänge schließen sich bei Inversionssymmetrie aus), atomare Zeiteinheit 2,4188843265864·10⁻¹⁷ s, c = 1/α ≈ 137 a.u. Die Seite zur Hyperpolarisierbarkeit nennt keine Vorfaktoren; die Festlegungen in `src/realtime/nonlinear.hpp` sind eigene. |
| Wikipedia, Kerr-Effekt | „Kerr effect“, en.wikipedia.org | 15.11.2026 | Optischer Kerr-Effekt: P = ε₀(χ⁽¹⁾ + ¾χ⁽³⁾\|E\|²)E cos ωt und n = n₀ + 3χ⁽³⁾\|E\|²/(8n₀), als Gegenprobe zur hier hergeleiteten Formel für n₂ (`nonlinear_index`) |
| nicht eingesehen | M. T. P. Beerepoot, D. H. Friese, N. H. List, J. Kongsted, K. Ruud, Phys. Chem. Chem. Phys. 17, 19306 (2015) | 15.11.2026: Abruf abgelehnt, nur der Text eines Suchergebnisses | Nicht als Quelle verwendet. Das Suchergebnis nennt die Form N·π³·α·a₀⁵·ω²/c des Vorfaktors für den Zwei-Photonen-Querschnitt, einen strittigen Faktor 4 in der Literatur und dass CAM-B3LYP die Stärken um das 1,5- bis 3-Fache unterschätzt. Die Formel im Code ist selbst hergeleitet; der Abgleich mit der Arbeit steht aus. |
| Wikipedia, Ellipsometrie-Modelle | „Lorentz oscillator model“, „Tauc–Lorentz model“, „Effective medium approximations“, „Birefringence“, en.wikipedia.org | 15.11.2026 | Lorentz-Oszillator ε∞ + Σ s/(ω₀² − ω² − iΓω) mit exp(−iωt); Tauc-Lorentz: ε₂ und Kramers-Kronig-Integral, die geschlossene Form für ε₁ nur als Gegenprobe im Test (beim Abruf unvollständig wiedergegeben); Bruggeman-Formel mit dem Hinweis zur Wahl der Wurzel; Wellenvektorfläche des einachsigen Mediums (`src/ellipsometry`) |

Hinweis zu Hernández Vera 2025, Gl. 17: Der mittlere Term ist dort als
E^t_ij gesetzt; gemeint ist E_t^ij, wie Gl. 18 und die Herleitung aus Gl. 14
zeigen.

Hinweis zur Ellipsometrie: Die Rekursion für Schichtstapel (Airy-Summe je
Schicht) und die charakteristischen Matrizen im Test sind aus den
Fresnel-Formeln selbst hergeleitet und im Code so bezeichnet; die Herleitung
steht in den Kommentaren. Fachbücher der Ellipsometrie wurden nicht eingesehen.
Die Größen N, C, S sind im Code definiert und keiner Quelle zugeschrieben.

Hinweis zum Zeitbereichsverfahren: Aus den Quellen stammen der Propagator, der
Puls und die Transformation. Selbst hergeleitet und im Code so bezeichnet sind
die Wirkung des Pulses auf die Orbitale, die Fock-Matrix für komplexe Dichten
und die Behandlung der Selbstkonsistenz im Zeitschritt (Polynom durch frühere
Fock-Matrizen, gemeinsame Verfeinerung der ersten drei Schritte). Abschnitt 5.9
von Blanes et al. über nichtlineare Gleichungen wurde nicht gelesen.
