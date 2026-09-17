# Ableitungen per AD auf einem Ausdrucksgraphen (cppDE, dMod2)

## Kontext

Der Codegen bildet jede Ableitung Eintrag für Eintrag mit `sp.diff` auf dem ganzen Ausdruck. Danach
folgen `sp.cse` über alle Einträge und die C++-Ausgabe jedes Eintrags. Betroffen sind Jacobimatrix,
df/dp, die Reverse-Kontraktionen, die Rosenbrock-Terme zweiter Ordnung, Events und CVODE.

Der Aufwand wächst mit Σ_i (Symbole in f_i) × Größe(f_i), bei gekoppelten Modellen also quadratisch.
Beispiel Fe-Spinwellenmodell: 481 Zustände, rund 155 000 Einträge zu je 13 ms, über eine Stunde
Codegen. Der `ThreadPoolExecutor` bringt wegen des GIL nichts.

**Ziel:** Codegen- und Compile-Zeit wachsen nur linear mit der Modellgröße, bei regelmäßiger Struktur
gar nicht mit N. Die Ableitungen bleiben exakt (bis auf Rundung), und die linearen Löser bleiben
direkt (dichte LU, KLU).

## Entscheidungen

**Python und Branch**
- Der Codegen bleibt in Python (reticulate), mit einem eigenen schlanken Ausdrucksgraphen.
- SymPy wird nur noch gebraucht für den Parser-Fallback, für gecachte Ableitungsvorlagen seltener
  Funktionen und als Referenz im Dev-Harness.
- Gearbeitet wird weiter auf `devel-reverseAD`, in beiden Repos. Simon committet selbst.

**Umfang**
- Gebaut werden: A (Graph, Parser, Drucker), B (JVP und VJP), C (Jacobimatrix aus Einträgen oder
  gefärbten JVPs), D (Klassengraph, Schleifen), E (lineare Kombinationen als Daten).
- F (Aufteilen in mehrere Übersetzungseinheiten, Codegen-Cache) kommt vorerst nicht.
- CVODE wird als letzte Phase migriert.

**Was wegfällt**
- `derivSymb()`
- der R-Fallback in cppFUN
- `derivMode = "symbolic"` in cppFUN, samt `_jacobian`, `_hessian`, `_chain_jac`, `_chain_hess` und
  `cppde_chain_blas.hpp`
- `attr(model, "jacobian")` in cppODE und in CVODE

**Folgen für cppFUN und dMod2**
- Ein unkompiliertes cppFUN-Objekt bricht beim Auswerten mit einer klaren Meldung ab.
- dMod2 passt sich an: `plotFluxes()` rechnet die Werte direkt in R, `Pimpl()` läuft auf `forward`,
  und die Tests kompilieren vorher.
- Alles Neue, was `cppDE::compile` braucht, muss auch `dMod2::compile` können. Geplant ist nichts
  davon: Es bleibt eine `.cpp`-Datei mit unveränderten Flags, der neue Header liegt im selben
  Include-Verzeichnis, und der Objekt-Cache hängt an `.headerStamp()`.
- In `dMod2/R/compile.R` sind nur die zwei veralteten BLAS-/chain-Kommentare (551-556, 693-697)
  anzupassen.

**Nicht geändert, aber dokumentiert**
- `dfdp_t_vec_axpy` ignoriert weiterhin die Ableitung nach `s_0`, genau wie heute.
- Verhaltensänderungen kommen in `NEWS.md`:
  - `abs` kompiliert im Forward-Modus;
  - die Ableitung von `floor`/`ceiling` ist 0 statt eines Fehlers;
  - die zweite Ableitung von `Heaviside` bei 0 ist 0.

**Beim Erkunden gefundene Fehler**
- **Vorzeichenfehler, bestätigt.** Das in R zusammengesetzte `compute_ydd` (`R/cppODE.R` 1056-1081)
  addiert `J_init·f`, obwohl der Funktor −J schreibt. Das betrifft nur den Startschritt von rb4.
- **Dedup-Loop-Pfad** (`codegen_cppODE.py` 863-882): Ein Eintrag, der nur von `t` abhängt, wird als
  konstant eingefroren.

## Verifikation gegen den Code (2026-09-16)

Geprüft auf `devel-reverseAD` bei e55c599. Der Plan trägt; Korrekturen und Ergänzungen:

**Umgebung**
- `dev/python/` gibt es noch nicht. `CPPDE_PY_DIR` ist neu: Die vier Loader in `R/zzz.R` lesen fest
  `system.file("python")`. reticulate 1.46 setzt `sys.path` schon nur während des Imports.
- `py_require("sympy")` trägt keine Python-Version. `DESCRIPTION` sagt Python ≥ 3.8 und wird 3.9.

**Zeilen und Namen**
- Der Symbolik-Test in `test-piecewise.R` steht in 63-88. `test-Pexpl.R` endet bei 150.
- `data_code` und `codegen_stats` sind neue Schlüssel. `data_code` kommt in `R/cppODE.R` direkt hinter
  `"namespace {"`. CVODE schreibt die Datei in Python selbst, dort gehört es in `_render_source`.
- Attribut `jacobian`: In `cppODE.R` betrifft es 1594-1604 und 1614, in `cvode.R` 247-253 und 264.
  `cvode.R` liest `jac_nnz_rows` für eine verbose-Meldung, also bleiben `jac_nnz_rows`/`jac_nnz_cols`.
- `dfdt_dot` wird für **jedes** Reverse-Modell erzeugt, weil die Sprünge es brauchen
  (`cppde_adjoint_step.hpp`). Nur die vier jvp/dfdt-Paare sind rb4-exklusiv.
- Färbung und Farbseeds gibt es noch nicht; „Farbseeds als `double`“ ist eine neue Festlegung.
- `_init_consts` steht heute nur im dünnen Funktor und wird nur im Dedup-Pfad benutzt. Das bleibt so:
  nur der dünne Pfad, denn dort lebt `m_W_sparse` einen Solve lang.
- `inst/examples/example_cppFUN.R` und `example_fun.R` lesen `jacobian.symb`/`hessian.symb`.
- `cppDE/CLAUDE.md` existiert nicht; gemeint ist die Workspace-`CLAUDE.md`.
- In dMod2 tragen auch `R/utils.R`, `R/parameters.R` (175, 269, 280-285), `R/prediction.R`
  (985-994, 1129, 1196, 1284), `man/Pexpl.Rd` und `man/Y.Rd` die Symbolik.

**Fehler**
- `compute_ydd`: Der Code rechnet `f_t − J f` statt `f_t + J f`, dicht wie dünn (`csc_matvec_add`).
  Die Wirkung ist nur der Startschritt von rb4 bei `hini == 0` und `f_t ≠ 0`, denn
  `weighted_sup_norm` nimmt Beträge.
- Dedup: Mit einem 80-Zustands-Modell reproduziert. Dort zeigt außerdem `x0_0` auf `params[n+k]`,
  im Standardpfad auf `params[i]`. `LEAF(INIT)` zeigt künftig einheitlich auf `params[i]`.
- Neu: `_generate_root_gradient_lambdas` ersetzt `DiracDelta` nicht. Ein `Heaviside` in der Wurzel
  eines nicht-terminalen Events scheitert deshalb im Drucker. Graph-AD behebt das.

**Fallback**
- Der R-Wert-Fallback greift in jedem Modus. Ohne `compile()` werten aus:
  - `test-cppFUN.R` 7-22, 24-46, 87-115, 119-144 und 148-161;
  - in dMod2 `test-steadystates.R` 28, `test-Pequil-Pimpl.R` (rund acht Pimpl), `test-Pexpl.R`
    111-149 und `test-Y.R` 141-157;
  - `inst/workshops/Pimpl_Pexpl.Rmd`.
  
  Alle diese Stellen kompilieren künftig vorher. Entschieden: Der Fallback fällt strikt weg, auch für
  Werte.
- In `Y()` prüft der AD-Zweig nicht, ob der Eintrag geladen ist.
- `_write_vjp_impl` differenziert ohne `jacobian` selbst. `_eval_one<T>` gibt es schon.

**Funktionen**
- Für duale Typen fehlen `fabs`, `floor`, `ceil`, `sign`, `heaviside`, `erf`, `atan2`, `cbrt` und
  `tgamma`. Deshalb:
  - Der Emitter senkt `sign`/`Heaviside` auf `cppde::select` ab (H(0) = 1/2) und `abs` auf `abs`.
  - `floor`/`ceiling` bekommen einen Helfer mit Tangente 0.
  - Andere Funktionen ohne duale Überladung sind bei dualem T ein Codegen-Fehler.
- CVODE schreibt +J, ohne CSE und in `double`, und sein dünnes Muster jedes Mal neu. Phase 8 behält
  das.

**Benchmarks**
- `benchmarks/R/harness.R` 277 ist seit 5154227 nicht parsebar.
- `benchmarks/cache/` fehlt und wird per `fetch-models.R` geklont.
- `max_states` ist 200/400.
- Das 481er-Fe-Modell liegt in keinem Repo. Die Phasen 6 und 7 nutzen das synthetische LLG-Dipolmodell.

**Anschluss (Phase 10):** Die großen PEtab-Modelle (≥ 100 Zustände), Brusselator/FHN mit großem N
und das 481er-LLG laufen einmal durch Codegen, `cppDE::compile` und `dMod2::compile`. Gemessen wird
forward und reverse, ohne Vergleich mit dem alten Codegen.

## Stand

**Reihenfolge ab Phase 8 (Entscheidung 2026-09-17):**

| Schritt | Inhalt | Stand |
|:--|:--|:--|
| 0–5 | Streichungen, Graph/Parser/AD, Reverse-Kontraktionen, Jacobimatrix, Events, cppFUN | fertig |
| 6 | E: lange lineare Summen als Tabellen | fertig |
| 7 | D: Schleifen über Klassen gleicher Struktur, Summen als innere Schleifen | fertig |
| 8 | CVODE auf dem Graphen | fertig |
| F | C++: festes Event mit Zeitparameter und uhrlesender rechter Seite, zweite Ordnung | fertig (1.8e-1 → 7e-14) |
| 10a | Neuinstallation, komplette cppDE-Benchmark-Suite lokal | fertig (42 Probleme, 879 Zellen) |
| 10b | große Modelle: Codegen- und Compile-Zeit mit `cppDE::compile` und `dMod2::compile` | fertig |
| 10c | dMod2: PEtab-Roundtrip einmal verifizieren | fertig: 32 von 32 laufen durch, 2 wie dokumentiert ausgelassen |
| 9 | erst danach: alten Codegen entfernen, Doku, Neuinstallation, beide testthat-Suiten nacheinander | fertig: cppDE 180 Tests, dMod2 382 Tests, 0 Fehler |

Ergebnisse 10a (Teil 1 `benchmarks/results/20260917-023711_full_c1_nosens-sens1`, Teil 2
`…062619…`; neben anderen Läufen, daher absolute Zeiten überhöht):
- Gegen den eingecheckten Lauf vom 2026-09-01 (Linux, alter Generator), gemeinsame Zellen:
  - Schritte neu/alt 1,00, Fehler 0,97 bis 1,03.
  - Dieselben vier Zellen scheitern (Oliveira dreimal, Weber, sens1 bei 1e-10).
- Neu ist Lang_PLOSComputBiol2024: CVODES scheitert mit Sensitivitäten bei 1e-4 und 1e-7 (s. u.).
- A/B neuer gegen alten Generator im selben Prozess:
  - `cppODE()`: Laufzeit 0,99 ohne und 1,01 mit Sensitivitäten.
  - `cvode()`: 0,97 und 0,87 bei gleicher Schrittzahl; die Spanne reicht von 0,28 bis 1,36.
    Die Sensitivitäts-RHS als ein JVP rechnet den Primalteil je Parameter neu.
- Die Head-to-head-Zahl mit Sensitivitäten in `Methods.Rmd` (1,84) stammt vom alten Generator
  und dürfte jetzt kleiner sein. Vor dem Zitieren auf einer ruhigen Maschine neu messen.

Smoke-Phasen 6 und 7:
- Alle Fälle grün; mit und ohne E bzw. Schleifen gleich bis 1e-14 / 1.6e-10 / 6e-17.
- LLG 481 läuft jetzt auf bdf mit FD-Schritt 1e-6. rb4 zerlegt je Schritt eine duale Matrix,
  und h = 1e-4 lag im Abbruchfehler.
- Zweite Ordnung (S5s) mit allen 27 Richtungen: Die Arena wächst im Solve binnen Sekunden um GB
  (neu 8,6 GB, alt 16,9 GB). Der Fall läuft jetzt über die Parameter.
  - Ergebnis: Gradient 4.1e-12, Hessematrix 4.9e-12.
  - Mit und ohne E gleich bis 3.6e-15; Spitzenspeicher 1,5 GB.
- Roundtrip-Nachläufe: Raimundez, Isensee, Bachmann und Schwen ändern ihren Wert beim Reimport.
  Mit dem alten Generator ist das identisch (Isensee 33649,2 → 1990350).

Tests und Compiles laufen einzeln und unter einer Speichergrenze.

Ergebnisse 10b (`benchmarks/results/compile-large-20260917-035729.csv`, neben drei anderen Läufen):

| Modell | Zustände | Codegen fwd/rev | `cppDE::compile` fwd/rev | `dMod2::compile` fwd/rev | Prüfung fwd/rev |
|:--|--:|--:|--:|--:|:--|
| Brusselator 1D | 1000 | 5.4 / 0.8 s | 13 / 13 s | 27 / 29 s | 3.5e-9 / 3.3e-8 |
| FitzHugh-Nagumo | 1000 | 0.6 / 0.7 s | 14 / 12 s | 19 / 28 s | 1.2e-6 / 5.2e-7 |
| LLG-Dipol | 481 | 105 / 95 s* | 16 / 13 s | 22 / 33 s | 2.8e-7 / 4.1e-7 |
| Chen_MSB2009 | 504 | 2.1 / 1.8 s | 222 / 22 s | 232 / 228 s | Solve scheitert** |
| Froehlich_CellSystems2018 | 1228 | 4.6 / 4.8 s | 179 / 27 s | 186 / 210 s | 9.2e-3 (FD) / 6.1e-8 |
| Lang_PLOSComputBiol2024 | 124 | 1.0 / 1.0 s | 97 / 16 s | 103 / 118 s | 5.2e-6 / 1.0e-10 |

- \* Davon 65 s in `checkSymbolNames()`: 35 TRE-Durchläufe über 20 MB Text. Jetzt ein
  PCRE-Durchlauf, 0,05 s. Python allein braucht 5,7 s.
- \*\* Mit dem alten Generator identisch (gleiche Zahl angenommener und verworfener Schritte bei
  1e-6, 1e-8 und 1e-10, mit und ohne Sensitivitäten). Der Benchmark-Aufbau von Chen scheitert in
  cppDE unabhängig vom Codegen.
- dMod2 baut in beiden Modi das Wert- und das Vorwärtsmodell, seine Zeit folgt dem Vorwärts-Compile.
- Ein Vorwärtsmodell (dual) verbringt unter `-O2` 59 % der Compile-Zeit in GCCs RTL-GCSE
  (Chen: 211 s; ohne GCSE 93 s; ohne Schleifen 394 s). Ab 40 kB dualem Modellcode schaltet ein
  GCC-Pragma im Quelltext den Pass ab. Lang: 90 → 55 s bei gleicher Laufzeit (6,56 / 6,32 s).
- Zum Vergleich, nebenbei gemessen: Der alte Generator braucht für Chen 47 s ohne und 414 s mit
  Sensitivitäten, der neue 19 s und 223 s (vor dem Pragma).
- Nach der Neuinstallation (Pragma, Schlüsselwortprüfung), mit weniger Last, jeweils forward:

  | Modell | Codegen | `cppDE::compile` | `dMod2::compile` |
  |:--|--:|--:|--:|
  | LLG | 15 s | 12 s | 23 s |
  | Chen | 0,7 s | 53 s | 59 s |
  | Lang | 0,2 s | 34 s | 39 s |

  Die Prüfungen sind unverändert.

Ergebnisse 10c:
- `test-petab`: 40 Tests, 0 Fehler, 0 übersprungen.
- Roundtrip über die Benchmark-Models-Sammlung: 25 von 34 in einem langen Arbeitsverzeichnis,
  2 ausgelassen wie dokumentiert.
  - Bachmann, Smith, Fiedler, Isensee und Raimundez scheiterten an Windows-Grenzen
    (Pfadlänge, Kommandozeile).
  - Alkan scheiterte an einer Parser-Regression (s. u.).
  - Im kurzen Verzeichnis laufen Alkan, Bachmann und Fiedler durch.
  - Isensee und Raimundez brauchten die dMod2-Korrektur der Kommandozeilengrenze (s. u.).
- Bachmann (−838,26 → −533,60) und Schwen (1901,98 → 2306,99) ändern ihren Wert beim Reimport,
  mit dem alten Generator identisch. Das liegt am Export von dMod2, nicht am Codegen.
- Lang: CVODES scheitert mit Sensitivitäten bei rtol 1e-4 und 1e-7 in 10a, ebenso die Referenz.
  Keine Regression:
  - Die neuen Rückrufe stimmen an Punkten der Trajektorie mit SymPy überein.
  - Über neun Toleranzen scheitert der neue Generator sechsmal, der alte siebenmal.
  - Ein Bit Unterschied in atol entscheidet über Erfolg.

Module: `cppde_graph.py` (Graph, Parser, AD), `cppde_emit.py` (Scheduler, Drucker), `cppde_model.py`
(Funktionen des Modells, Events, CVODE-Rümpfe), `cppde_struct.py` (E und D; ersetzt die geplanten
`cppde_linear.py` und `cppde_vector.py`), `cppsympy.py` (nur SymPy). Der Harness liegt in
`dev/python/` (`harness.py`, `run_checks.py`, `pyrt.py`, `smoke.R`, `corpus.json`,
`legacy/` = e55c599).

Ergebnisse:
- Parität zum alten Generator in `run_checks.py`: 1e-12. Die Summationsreihenfolge langer Summen
  verschiebt die letzten Stellen, deshalb nicht 1e-13.
- Vergleichspunkte, an denen schon der Referenzwert nicht endlich ist oder ein Zwischenwert
  1e12 übersteigt, werden gezählt, aber nicht verglichen. Werte, die bis auf Rundung auslöschen,
  gelten bis 1e-12 des größten Ausgabewerts derselben Funktion als gleich.
- Im gefärbten Pfad macht ein singulärer Eintrag über `0·∞` die ganze Zeile nicht endlich.
- `smoke.R` vergleicht jeden Fall zusätzlich mit dem alten Generator. Die Abweichung liegt bei
  höchstens 1e-10, meist unter 1e-15.
- E: jede Funktion mit und ohne Tabellen gleich (1,8 Mio. Vergleiche). Das LLG-Modell mit 481
  Zuständen und 20 MB Modelltext wird in etwa 5 s erzeugt (3,7 s Parsen, 1–2 s Code).
- D: mit und ohne Schleifen bitgleich; Summen als Schleife bis 1e-13. Die Zahl der Anweisungen
  ist bei N und 2N gleich; Brusselator 2D mit 4608 Zuständen: 6 s Codegen.
- **Abweichung bei D (bewusst):** Die Ableitungen entstehen weiter skalar auf dem Graphen, die
  Klassen und Schleifen werden danach über den fertigen Anweisungen gebildet
  (`cppde_struct.vector_block`). Codegröße und Compile-Zeit hängen damit nicht von N ab, die
  Codegen-Zeit wächst linear (rund 1 ms pro Zustand für rb4 zweiter Ordnung; 10^5 Zustände
  brauchen Minuten und GB in Python). Dafür teilen alle Ableitungen einen geprüften AD-Pfad, und
  die Schleifen sind bitgleich zum skalaren Code. AD auf dem Klassengraphen lässt sich später vor
  `vector_block` setzen, ohne die Ausgabe zu ändern.
- CVODE: Rümpfe gegen SymPy grün; C++-Syntax von zehn Varianten grün.

Gefunden:
- `normalise_logic` las das Fakultäts-`!` als Not.
- `cppde::max`/`min` haben keine Überladung für Expression-Templates. Der alte Generator
  kompilierte `max(x, 2*y)` im AD-Modus deshalb nicht; die Argumente werden jetzt Temporaries.
- Konstruktionszeit-`fixed` wirkte in `cppFUN(derivMode = "forward")` nicht; behoben.
- Summen ab einigen tausend Termen brachen `ast.parse` mit RecursionError ab; sie werden vorher in
  einen flachen Aufruf umgeschrieben.
- Hash-Consing legte `0.5` (float) und `Fraction(1, 2)` im selben Knoten ab; der Typ hing von der
  Reihenfolge ab. Floats mit kleinem Zweiernenner werden jetzt als Fraction gespeichert.
- Ohne E kompiliert ein LLG-Modell mit 19 Zuständen in zweiter Ordnung etwa 9 Minuten bei
  1,7 GB Speicher; sechs solche Compiles parallel sind zu viel.
- Der alte CVODE-Generator las `<state>_0` als Anfangswert, differenzierte aber nach der
  gleichnamigen Parameterzeile. Der neue bleibt aus Kompatibilität dabei.
- **C++, unabhängig vom Codegen:** Ein festes Event mit parameterabhängiger Zeit und eine rechte
  Seite, die die Uhr liest, ergeben in forward-reverse einen falschen Eintrag d²/dθ² für den
  Zeitparameter (Beispiel: `t_dose*2` und `-k1*A + k2*B*time`, Abweichung 0.9 gegen ff und FD).
  Ursache: `apply_fixed_jump_adjoint` sammelt das Heun-Sandwich geprunt und lässt die Beiträge
  `df/dt` von `f(x_e, t_e)` und `f(x_a, t_e)` zur Kotangente des Shifts weg, das Gegenstück zu
  955ecbe für Wurzel-Events. Behoben.
- Eine Jacobimatrix erster Ordnung hielt ihre Temporaries bis zum Ende des Solves in der Arena
  (LLG 481: 18 GB). Sie bindet jetzt ihre Ausgaben und öffnet einen eigenen Scope. Verschachtelte
  Duale haben weiter keinen Scope pro Aufruf; ein langer Solve zweiter Ordnung wächst mit der
  Schrittzahl.
- SBML-`piecewise(1, x > 0, 0, x <= 0)` hat keinen otherwise-Zweig, aber erschöpfende Bedingungen.
  SymPy macht daraus `True`, der Graph-Parser lehnte ab (Alkan_SciSignal2018). Ohne
  otherwise-Zweig fällt der Parser jetzt auf SymPy zurück; Prüfung in `run_checks.py fallback`.
- dMod2 rechnete unter Windows mit 24000 Zeichen Kommandozeile, `R CMD` läuft aber durch cmd.exe
  (8191). Jetzt 8000, und ein `ar`-Block zählt den Aufruf davor mit.

## Compile-Budget (gilt für alle Phasen)

- **Während der Entwicklung wird nur in Python geprüft:**
  - Ein Python-Backend des Emitters erzeugt dieselben Anweisungen, Temporaries, Schleifen und
    Tabellen wie das C++-Backend, nur mit anderer Blattsyntax. Der Code läuft mit `exec`,
    mit `float`, einer `Dual`-Klasse oder einer verschachtelten `Dual`-Klasse.
  - Referenz sind die eingefrorenen SymPy-Pfade und finite Differenzen.
  - Aufruf in Sekunden (Windows; unter Linux liegt uv unter `~/.cache/R/reticulate/uv/bin/uv`):
    `$LOCALAPPDATA/R/cache/R/reticulate/uv/bin/uv.exe run --python 3.12 --with sympy==1.14 python dev/python/run_checks.py`
- **Am Phasenende** kompiliert `dev/python/smoke.R` höchstens eine Handvoll Modelle, parallel.
- **Die volle testthat-Suite** läuft in beiden Paketen genau einmal ganz am Ende, im Hintergrund.
- **Die Benchmark-Suite** läuft in Schritt 10a lokal, vor dem Entfernen des alten Codegens.
- **Neuinstallation** gibt es nur an den markierten Stellen. Alternativ liest `R/zzz.R` optional
  `CPPDE_PY_DIR`, dann nutzen die Smoke-Tests den Quellbaum ohne Reinstall.

## Architektur (`cppDE/inst/python`)

| Modul | Inhalt |
|:--|:--|
| `cppde_ir.py` | Graph mit Hash-Consing: ganzzahlige IDs, `(op, args, attr)` als Schlüssel, deterministische Reihenfolge (keine Iteration über Mengen, wegen dMod2s md5-Wiederverwendung). Knoten: `NUM` (Fraction/float), `NAMED`, `BOOL`, `LEAF(STATE, PARAM, INIT, TIME, FORCING, FRATE, VEC, LINROW, TABLE, LOOPVAR)`, n-äres `ADD` und `MUL` mit Koeffizienten, `POW`, `CALL`, `SELECT`, `CMP`/`AND`/`OR`/`NOT`. Abhängigkeitsflags und Bitsets (`sbits`, `dbits`, `pbits`), Funktionstabelle. |
| `cppde_parse.py` | `normalise_logic`, `^`→`**`, `ast.parse`, SymPy-kompatible Normalform (Zusammenfassen, Zahlen falten wie SymPys Auto-Auswertung). Fallback bei Syntax oder exotischen Funktionen: `cppsympy.safe_sympify` und dann `from_sympy`. Unbekannter nackter Name ergibt einen Fehler. |
| `cppde_ad.py` | `jvp` (symbolische Seeds), `vjp` (eine n-äre Summe pro Adjoint), zweite Ordnung als `vjp(jvp(...))`, `time_derivative` mit Forcing-Raten, Lie-Ableitung, Muster aus Bitsets, symbolische Einträge per dünner Vorwärtsakkumulation, Distanz-2-Färbung mit Scatter-Tabellen, Regel „Einträge oder Färbung“ (siehe unten). |
| `cppde_emit.py` | Scheduler pro Funktion: materialisiert bei ≥ 2 Nutzungen oder großer Teilausdrücke; ein einmal genutztes `select` bleibt inline (`test-piecewise.R`). Erst Parameter-Knoten, dann Zeit-Knoten, dann der Rest topologisch, Invarianten aus Schleifen heraus. Anweisungen `Decl`, `Store`, `For`, `If`, `Tag`, `LinApply`, `Table`. C++-Backend für festes T und Template-Modus; Python-Backend. |
| `cppde_linear.py` | E: Summen aus mindestens 12 Termen `c·x_j` (c numerisch, gruppiert nach zustandsfreiem Faktor) werden zu `s·LINROW(map, r)`. `LinMap` dicht oder CSR, gleiche Zeilen zusammengelegt, optional mit Wertepool. |
| `cppde_vector.py` | D: Klassenschlüssel aus Skelett, Schnittpunkten und Slot-Tabellen; Aufspaltung nach Gleichheitsmustern; unter 16 Instanzen wieder skalar. Reduktionsklassen für Summen mit variabler Länge (Kuramoto). AD läuft über dieselbe Knotenschnittstelle auf dem Klassengraphen, Akkumulation als Scatter-Add. |
| `inst/include/cppde/cppde_linmap.hpp` | `linmap_csr`/`linmap_dense`, `apply`, `apply_t_add`, `row_dot`, `axpy_row_dense`, `axpy_row_idx`, generisch im Elementtyp. Wird von `cppde.hpp` eingebunden. |
| `cppsympy.py` | Die einzige Stelle für SymPy-Hilfen: `safe_sympify`, `sbml_piecewise`, Parse-Dictionary (ersetzt die Kopien in `codegen_cppODE.py`, `codegen_cppFUN.py` und `derivSymb.py`), `from_sympy`. |

**Regeln**
- Alle Importe zwischen Modulen stehen auf Modulebene. `import_from_path` setzt `sys.path` nur
  während des Imports.
- Die öffentliche API bleibt: `generate_ode_cpp`, `generate_event_code`, `generate_rootfunc_code`,
  `generate_forcing_init_code`, `fixed_event_time_exprs`, `decide_sparse`, `analyze_klu_settings`,
  `generate_fun_cpp`, `generate_cvode_cpp`.
- Neue Rückgabeschlüssel: `data_code` (Tabellen auf Namespace-Ebene, von R direkt nach
  `namespace {` eingesetzt) und `codegen_stats`.
- Wegfallende Rückgabeschlüssel: `jac_nnz_exprs` und `time_derivs`.

**Invarianten, die erhalten bleiben**
- **Jacobi-Funktor:**
  - schreibt −J;
  - dicht mit Dirty-Reset über `_dr`/`_dc`;
  - dünn besetzt in CSC-Ordnung mit Auffüllen fehlender Diagonalen und `build_pattern`;
  - `_init_consts` nur für Einträge ohne Abhängigkeit von Zustand, Zeit oder Forcing (behebt den
    Dedup-Fehler);
  - alle Anweisungen im Skalartyp T, damit Dual- und Dual2nd-Tangenten für die IFT-LU stimmen;
  - Farbseeds als `double`.
- **Arena-Scope** nur bei `ad_level == 1`, nur in `ode_system` und `adjoint_terms`, nie in
  `jacobian`, `dfdt_dot` oder Events.
- **Kontraktionen** mit unveränderten Signaturen und Ausgabeformen:
  - `jac_t_vec`, `dfdp_t_vec_axpy` (Slots `[0, n_states)` unberührt);
  - für rb4: `jvp_x_t_vec`, `jvp_p_t_vec_axpy`, `dfdt_x_t_vec`, `dfdt_p_t_vec_axpy`, `dfdt_dot`.
- **Event-Lambdas** im Kontext `full_params`; G_tt als `double` mit `.val()`-Peeling; `nullptr` bei
  Forcings.

**Einträge oder Färbung** (Startwerte, per `CPPDE_JAC` erzwingbar, später per Benchmark justiert)
1. Einträge, wenn `W_ent ≤ max(4·E_s, 50 000)`.
2. Sonst Färbung, wenn χ ≤ n/4.
3. Sonst Einträge.

Mit E zählen `LINROW`-Zeilen als Pseudovariablen, und es gilt J = A + G·C, wobei C per
`axpy_row_*` als Daten eingeht.

## Phasen

### Phase 0: Streichungen und dMod2-Anpassung

**cppDE**
- **`derivSymb` entfernen:**
  - `R/symbolics.R`: `derivSymb` löschen (`getSymbols` und `checkSymbolNames` bleiben);
  - `R/zzz.R`: `get_derivSymb_py` löschen;
  - löschen: `inst/python/derivSymb.py`, `man/derivSymb.Rd`, `cppde_chain_blas.hpp`;
  - `NAMESPACE` per roxygen neu erzeugen.
- **`R/cppFUN.R`:**
  - `derivMode` nur noch `c("forward", "reverse")`;
  - löschen: `safeParse`, `parsed_*`, `sym_*`, `.raw_*_sym`, `.chain_*_sym`, `.buildSeedMatrix`,
    `.buildSeedTensor2`, die symbolischen Zweige in `.jac_impl`, `.hess_impl` und `.evaluate_impl`,
    die Attribute `jacobian.symb`/`hessian.symb`;
  - `fwd <- emit_deriv && use_ad`;
  - nicht geladene Einträge laufen auf das neue `.notCompiled(st)` („… is not compiled; call
    compile()“);
  - die `.C`-Pfade für Altobjekte bleiben.
- **`R/tools.R`, `matchDerivMode`:** eigene Fehlermeldung für `"symbolic"`.
- **`codegen_cppFUN.py`:** `jacobian`/`hessian`-Argumente, symbolische Schreiber, Chain-Wrapper und
  `n_sym_*` entfernen. `_write_vjp_impl` bleibt bis Phase 5.
- **`R/cppODE.R` (1594-1614, Doku Zeile 76) und `R/cvode.R` (245-264):** Attribut `jacobian` entfernen.
- **Tests:**
  - `test-cppFUN.R`: Symbolik-Vergleiche ersetzen durch geschlossene Formen oder finite Differenzen;
    neuer Test, dass ein unkompiliertes Objekt einen Fehler wirft;
  - `test-dual2nd-primitives.R` und `test-piecewise.R` 63-88: Referenz über `stats::D()` in R;
  - `test-no-parameters.R` 137-154: auf `forward` umstellen.
- **`NEWS.md`.**

**dMod2**
- `R/utils.R`: `.matchDerivMode` wie in cppDE.
- `R/parameters.R`:
  - `Pexpl`: Auswahl ohne `"symbolic"`; der Rückfallzweig in `.Pexpl_p2p` liefert nur Werte bei
    `deriv = FALSE`, sonst „compile first“;
  - `Pimpl()`: `derivMode = "forward"`, dokumentiert als erst nach `compile()` auswertbar.
- `R/prediction.R`: `Y()`-Auswahl und Doku; der symbolische Zweig in `X2Y` gibt nur Werte oder
  „compile first“ zurück.
- `R/plots.R`, `plotFluxes()`: `parse()`/`eval()` pro Bedingung mit kleiner Umgebung (`Heaviside`,
  `exp10`, `piecewise`).
- `R/compile.R`: die zwei Kommentare.
- Tests: `test-Pequil-Pimpl.R`, `test-Pexpl.R` 78-149, `test-Y.R` 57-63 und 143-160,
  `test-constructor-deriv-flag.R`, `test-deriv2.R` („symbolic“ wird „forward“).

**Prüfung:** Reinstall beider Pakete, dann nur die umgeschriebenen Testdateien.

**Ende der Phase:** Diese Dateien laufen grün, und `grep` nach `symbolic`/`derivSymb` findet in beiden
`R/`-Bäumen nur noch Fremdtreffer.

### Phase 1: Graph, Parser, Drucker, Python-Backend, Harness

**Dateien:**
- neu: `cppde_ir.py`, `cppde_parse.py`, `cppde_emit.py`;
- `cppsympy.py`: gemeinsame Hilfen;
- beide Codegen-Module importieren diese Hilfen;
- `DESCRIPTION`: Python ≥ 3.9 (wegen `ast.unparse`);
- optional `CPPDE_PY_DIR` in `R/zzz.R`.

**Harness `cppDE/dev/python/`** (per `^dev$` vom Build ausgeschlossen):

| Datei | Inhalt |
|:--|:--|
| `pyrt.py` | `select`, `heaviside`, `sign`, `linmap_*`, Fake-`PchipForcing`, `Dual` |
| `oracle.py` | eingefrorene Kopien von `_safe_sympify`, `_compute_ode_jacobian_serial`, `_compute_ode_dfdp`, `_replace_dirac_delta`, Jv/dfdt, Event-Mathematik, CVODE-Ausdrücke |
| `models.py` | Modellsammlung |
| `corpus.json` | Korpus, einmalig erzeugt von `collect_corpus.R` aus Tests und Beispielen |
| `harness.py`, `run_checks.py` | Prüflauf |
| `smoke.R` | Compile-Stichproben am Phasenende |

**Prüfungen (nur Python):**
- Paritätstest über den Korpus: Parser, Python-Backend und Auswertung gegen `lambdify`, auf beiden
  Seiten jeder Verzweigung, relativer Fehler ≤ 1e-13.
- Der Fallback wird mit `x!`, `re` und `beta` geprüft.
- Die Ausgabe ist unter zwei `PYTHONHASHSEED` identisch.

**Compile:** kein Modell. Nur ein `g++ -fsyntax-only -I inst/include` über eine erzeugte Datei mit dem
Korpus in `double`, `dual` und `dual2nd`.

**Ende der Phase:** 100 % Parität, und der Syntaxcheck läuft sauber.

### Phase 2: Reverse-Kontraktionen aus dem VJP (einschließlich rb4 zweiter Ordnung)

**Dateien:**
- `cppde_ad.py`: `jvp`, `vjp`, `time_derivative`.
- `codegen_cppODE.py`: neuer `adjoint_terms`-Builder. Die Jacobimatrix, die bei `skip_jacobian` nur für
  die Kontraktionen berechnet wurde, entfällt.
- `R/cppODE.R`: `data_code` einsetzen.

**Prüfungen (nur Python):**
- Alle sieben Funktionen gegen die SymPy-Kontraktionen des Referenzpfads, bei zufälligen `x, v, lam,
  t, sc`.
- Abgleich mit finiten Differenzen.
- `T = Dual` (Forward über Reverse) gegen finite Differenzen.
- Modelle mit Forcings, Piecewise, `abs`/`sign`/`min`/`max`, `s_0`, ohne Parameter, mit expliziter
  Zeitabhängigkeit.

**Compile:**
- S1: 3 Zustände, Forcing, zeitabhängiges Piecewise; gebaut als reverse bdf und forward-reverse rb4.
- S2: 40er-Kette, dünn besetzt, reverse.
- Abgleich: reverse gegen `sens1`, `adjoint2` gegen `sens2`.

**Ende der Phase:** Übereinstimmung auf 1e-8.

### Phase 3: Jacobimatrix (Einträge und Färbung), `ode_system`, `dfdt`, alle Modi

**Dateien:**
- `cppde_ad.py`: Einträge, Muster, Färbung, Auswahlregel.
- `codegen_cppODE.py`:
  - `ode_system` und `jacobian` in dicht, dünn und No-op;
  - Muster für `decide_sparse`/`analyze_klu_settings`;
  - `_try_template_dedup` entfällt; große Modelle laufen bis Phase 7 über Einträge oder Färbung;
  - `jac_nnz_exprs` und `time_derivs` fallen weg.
- `R/cppODE.R`: Vorzeichenkorrektur in `compute_ydd` als eigene Änderung.

**Prüfungen (nur Python):**
- Beide Strategien werden auf jedem Modell erzwungen.
- Emulierte `dense_matrix` (mit Reset) und `csc_matrix` (mit `build_pattern`) gegen den
  Referenzpfad (−J und `dfdt`).
- Die Tangenten von `Dual` und verschachteltem `Dual` gegen SymPy d/dp.
- Muster, `klu_settings` und `decide_sparse` gleich wie bisher (Modelle aus
  `test-sparse-autodetect.R`).
- Ein nur zeitabhängiger Eintrag im dünnen Pfad bei zweitem `t`.

**Compile:**
- S1 forward bdf (dicht, IFT).
- S2 forward-forward rb4 (dünn) und mit `CPPDE_JAC=colour` (dual).
- S1 tsit5.
- `test-sparse-autodetect.R`.

**Ende der Phase:**
- `sens1`/`sens2` stimmen mit finiten Differenzen auf etwa 1e-6 überein.
- Färbung und Einträge stimmen auf 1e-10 überein.

### Phase 4: Events und Wurzeln

**Dateien:** `codegen_cppODE.py`:
- `generate_event_code` in beiden Modi;
- `generate_rootfunc_code`;
- `fixed_event_time_exprs`;
- Lie-Ableitungen;
- G_tt als `double` mit Peeling.

Danach ruft der ODE-Codegen SymPy nur noch über Fallback und Vorlagen auf.

**Prüfungen (nur Python):** Lambda-Rümpfe und Switch-Fälle gegen die Event-Mathematik des
Referenzpfads (dg/dx, dg/dt mit Forcing-Termen, G_tt, `_event_grad_case`, `_root_gdot_sym`) und gegen
finite Differenzen.

**Compile, S3:**
- Modell: nicht-terminales, zeitabhängiges, nichtlineares Wurzel-Event plus festes Event mit
  parametrisierter Zeit.
- Gebaut als forward-forward und forward-reverse.
- Abgleich: Hessematrizen gegeneinander, Gradient gegen finite Differenzen.

**Ende der Phase:** Übereinstimmung auf 1e-8.

### Phase 5: cppFUN

**Dateien:** `codegen_cppFUN.py`:
- `_eval_one<T>` im Template-Modus;
- `_vjp_impl` aus dem Graph-VJP, mit Seed-Schleife und Hoisting;
- restlicher SymPy-Code raus.

**Prüfungen (nur Python):**
- `eval` gegen SymPy.
- `vjp` gegen die Kontraktion der Jacobimatrix.
- `vjp` mit `Dual` (entspricht `vjp2`) gegen finite Differenzen.
- Randfälle: keine Variablen, keine Parameter.

**Compile, S4:**
- Modell: Observable mit Piecewise, `pow`, `abs` und `fixed`; forward und reverse mit `deriv2`.
- Abgleich: `jac`/`hess` gegen `D()`, `vjp` gegen `jac`, `vjp2` gegen die Hessematrix.

**Ende der Phase:** Übereinstimmung auf 1e-10.

### Phase 6: E (lineare Kombinationen als Daten)

**Dateien:**
- `cppde_linear.py`;
- `LINROW`/`LinApply` und die G·C-Assemblierung in AD und Emitter;
- `cppde_linmap.hpp` plus Include;
- `dev/cxx/test_linmap.cpp` in `run.sh`.

**Prüfungen (nur Python):**
- Jede Funktion wird mit E und ohne E erzeugt, relativer Abstand ≤ 1e-13.
- Kleine n zusätzlich gegen den Referenzpfad.
- Einmal die Codegen-Zeit des synthetischen Dipol-Modells mit 481 Zuständen ausgeben.

**Compile, S5:** Dipol-artiges Modell mit n ≈ 60, forward dicht und reverse dünn, jeweils mit und ohne E.

**Ende der Phase:**
- Übereinstimmung auf 1e-10.
- Das 481er-Modell wird in Sekunden erzeugt.

### Phase 7: D (Klassengraph, Schleifen)

**Dateien:**
- `cppde_vector.py`;
- Klassengraph-Schnittstelle in `cppde_ad.py`;
- Schleifen, Reduktionen und Scatter-Add in `cppde_emit.py`.

**Prüfungen (nur Python):**
- Modelle: 1D-Reaktion-Diffusion (periodisch und nicht), 2D-5-Punkt-Stencil, Kuramoto, LLG mit E.
- Jede Funktion mit und ohne Vektorisierung (`CPPDE_VECTORISE`) identisch.
- Die Zahl der Anweisungen ohne Tabellen ist bei N und 2N gleich.

**Compile, S6:**
- 2D-MOL mit N ≈ 10⁴, forward dünn und reverse rb4.
- Das 481er-LLG-Dipolmodell, forward.
- Codegen- und Compile-Zeit einmal notieren.

**Ende der Phase:**
- Übereinstimmung auf 1e-10.
- Die Codegröße hängt nicht von N ab.

### Phase 8: CVODE

**Dateien:** `codegen_cvode.py`:
- Jacobi-Einträge ohne Diagonal-Auffüllung;
- `sens_rhs1_fn` als ein JVP mit Zustands- und Parameter-Seeds;
- `adj_rhs_fn` und `adj_quad_fn` als VJP mit negativem Vorzeichen;
- Event- und Wurzelableitungen;
- `linmap`-Include, wenn E aktiv ist.

**Prüfungen (nur Python):** gegen die eingefrorenen SymPy-Ausdrücke des alten Generators.

**Compile:** je ein CVODE-Modell mit Sensitivitäten und Events, reverse, und dünn besetzt, jeweils
gegen das native Backend.

**Ende der Phase:** Übereinstimmung im Rahmen der Solver-Toleranz.

### Schritt F: fester Event-Zeitpunkt und uhrlesende rechte Seite (C++)

**Datei:** `inst/include/cppde/cppde_adjoint_step.hpp`, `apply_fixed_jump_adjoint`: beide
Auswertungen der rechten Seite am Eventzeitpunkt geben dem Shift ihr `df/dt` (`dfdt_dot`), wie
der Wurzelpfad seit 955ecbe.

**Prüfung:** `dev/cxx/test_reverse_events2.cpp` mit parameterabhängiger Eventzeit und einem Reset
auf der uhrlesenden Komponente; ohne Fix 1.8e-1, mit Fix 7e-14. `smoke.R` Phase 4 wieder mit
`-k1*A + k2*B*time`.

### Schritt 10: Benchmarks und PEtab, vor dem Aufräumen

- **10a:** cppDE neu installieren, dann `Rscript benchmarks/run-benchmarks.R --tier full` lokal
  und seriell.
- **10b:** `Rscript benchmarks/compile-large.R`: PEtab-Modelle ab 100 Zuständen, Brusselator und
  FHN mit 1000 Zuständen, LLG mit 481; Codegen, `cppDE::compile`, `dMod2::odemodel` +
  `dMod2::compile`, Stichprobe gegen FD bzw. vorwärts.
- **10c:** dMod2 neu installieren; `test-petab.R` mit `DMOD_PETABTESTS`; Import, Auswertung mit
  Ableitungen, Export nach v2 und Reimport über die Benchmark-Models-Sammlung (33 von 35 laut
  NEWS).

### Phase 9: Aufräumen, Doku, Gesamtlauf

**Aufräumen:**
- verbleibende SymPy-Pfade in `codegen_cppODE.py` entfernen (`_to_cpp`, `_safe_sympify`, alte
  Emitter und Helfer);
- `CppdePrinter` entfernen, falls ungenutzt.

**Doku:**
- `NEWS.md`;
- `dev/methods/Methods.Rmd` (57-59, 527, 1942, 2011, 2077-2226, 2258) und neu rendern;
- Description-Feld in `DESCRIPTION`;
- `cppDE/CLAUDE.md`: veraltete Namen (`derivSymb.py`, `funCpp`), Anzahl der Testdateien (12), die
  Aussage, dass `compile()` Objekte wiederverwendet;
- Workspace-`CLAUDE.md`: `funCpp()`;
- roxygen.

**Gesamtlauf:**
- Reinstall cppDE, dann dMod2.
- Beide testthat-Suiten einmal im Hintergrund (`TESTTHAT_CPUS=6` bzw. `10`).
- Optional `Rscript dev/check-like-ci.R`.

**Ende der Phase:** Beide Suiten laufen grün.

## Verifikation (Zusammenfassung)

1. `dev/python/run_checks.py` nach jeder Änderung. Das ist schnell und kompiliert nichts.
2. `dev/python/smoke.R` am Ende jeder Phase, mit wenigen Modellen und parallel.
3. `dev/cxx/run.sh` für die neuen Laufzeit-Helfer. Die Dateien sind klein, das dauert Sekunden.
4. Die volle testthat-Suite in cppDE und dMod2 einmal in Phase 9. Die Benchmark-Suite folgt
   danach als eigenes Vorhaben.

## Risiken

- **Abweichungen gegenüber SymPys automatischer Vereinfachung**, zum Beispiel bei Kürzungen oder
  0·inf an Singularitäten. Abgefangen durch die Korpusparität; bewusste Abweichungen werden
  dokumentiert.
- **Arena-Wachstum im gefärbten Pfad unter Dual-Typen.** Wird in Phase 3 gemessen; wenn nötig, fällt
  die Wahl auf Einträge zurück.
- **Große Zahlentabellen**, zum Beispiel eine dichte Map. Gegenmittel: gleiche Zeilen zusammenlegen
  und einen Wertepool nutzen. F bleibt zurückgestellt.
- **Summationsreihenfolge ändert sich durch E und D** (Unterschiede um 1e-16). Tests, die auf exakte
  Gleichheit prüfen, müssen geprüft werden.
- **Pimpl auf forward** macht einen Dual-`.Call` pro nleqslv-Jacobimatrix. Unkritisch für kleine
  Systeme.
- **Kuramoto-artige nichtlineare dichte Kopplung** wird erst mit den Reduktionsklassen aus D linear.
