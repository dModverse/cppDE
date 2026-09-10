# Der geschlossene Adjungierte: das Tape aus dem Solver nehmen

**Nachfolgeplan zu `dev/adjoint-plan.md`, angelegt 2026-09-10.** Jener Plan hat den
diskreten Adjungierten gebaut und ist abgearbeitet. Dieser nimmt ihm das Tape ab.

## Kontext

Der diskrete Adjungierte steht und liefert richtige Gradienten, ist aber langsamer als
CVODES ASA. Auf Bachmann, ganze Kette bis `normL2`, 36 Bedingungen, 113 Parameter,
`atol = 1e-11`, `rtol = 1e-9`, ein Kern:

| Route | ms | in Wertläufen |
|---|---|---|
| Wertlauf | 36,5 | 1,00 |
| Vorwärts | 1026,4 | 28,1 |
| diskreter Adjungierter | 335,7 | 9,19 |
| CVODES ASA | 216,6 | 5,93 |

Der Profiler sagt, wo die Zeit liegt. Bedingung `long`, 721 Schritte, ein Solve:

| | ms | je Schritt |
|---|---|---|
| Vorwärtsintegration, ganzer Stepper | 0,506 | 0,70 µs |
| `rev_replay`, Schritt auf dem Tape-Typ nachspielen | 4,356 | 6,04 µs |
| `rev_sweep`, Tape rückwärts | 1,978 | 2,74 µs |
| `rev_prepare` + `rev_solve`, transponierte Algebra | 0,665 | 0,92 µs |

Der Rückwärtslauf kostet das 12,5-fache der Vorwärtsintegration, und 90 Prozent davon
sind das Tape. Die transponierten Solves sind zehn Prozent. Eine Auswertung der rechten
Seite dauert 26 ns; Replay und Sweep zusammen sind 8,8 µs je Schritt, also rund
dreihundert Auswertungen der rechten Seite für einen Schritt, dessen Vorwärtsvariante
siebenundzwanzig kostet.

Das Tape differenziert die Schrittimplementierung operationsweise, obwohl die Ableitung
feststeht. Zwei Dinge, die es bei seiner Konstruktion nicht gab, machen die geschlossene
Form heute möglich:

- **Das Gitter ist eingefroren.** Schrittweite, Ordnung und Fehlernorm werden seit dem
  2026-09-09 nicht mitdifferenziert, und das bleibt so. Zu adjungieren ist eine feste
  Schrittfolge.
- **Die Kontraktionen existieren bereits.** `codegen_cvode.py:144-155` leitet `df/dp`
  symbolisch ab, `:457-470` emittiert `λ' = −J'λ` und `q' = −(df/dp)'λ` als flache
  Akkumulationen. Beides für die ASA gebaut, beides direkt übertragbar. Der native
  Emitter hat kein `df/dp`, aber dieselbe Maschinerie.

**Das Ziel:** der Adjungierte wird erzeugt und geschrieben statt aufgezeichnet. Danach
trägt `cppde::codual` nichts mehr und der Typ verschwindet, im Solver wie in `cppFUN`.

## Was das bringt, und was nicht

Ein geschlossener Schritt-Adjungierter kostet den transponierten Solve, ein bis zwei
Kontraktionen und Historien-Bookkeeping, geschätzt gut eine Mikrosekunde statt 8,8. Der
Rückwärtslauf läge dann beim rund 1,7-fachen der Vorwärtsintegration statt beim
12,5-fachen. Auf Kettenebene fiele die Zielfunktion von 336 ms auf grob 165 ms, gegen
ASAs 217 ms: aus 1,55-mal langsamer wird ungefähr 1,3-mal schneller.

**Eine Größenordnung über ASA ist nicht erreichbar, und das ist Arithmetik.** Ein
Wertlauf der Kette kostet 36,5 ms, ASA 217 ms. Ein Zehntel davon wäre weniger als ein
einziger Wertlauf, und ein Gradient braucht mindestens die Trajektorie plus einen
Rückwärtslauf. Der Boden liegt bei zwei bis drei Wertläufen. Eine Größenordnung ist auf
dem Tape zu holen, nicht auf dem Gradienten.

Danach ist der nächste Posten die Kette darüber: bei `long` sind von 10,76 ms rund
3,8 ms weder Vorwärtslauf noch Sweep, sondern die vjps von `g` und `p` und der
R-seitige Aufwand. Heute 36 Prozent, danach die Mehrheit. Das ist die Nachfolge dieses
Plans, nicht sein Inhalt.

## Eine Festlegung fällt

`dev/adjoint-plan.md` legt fest: *"Der Weg ist ein Reverse-Skalartyp, nicht
handgeschriebene Adjoint-Gleichungen pro Stepper."* Diese Festlegung wird umgedreht, aus
zwei Gründen, die es bei ihrer Niederschrift nicht gab: die Steuerungskette ist raus,
und die Kosten sind gemessen. Was von ihr richtig bleibt, ist ein gemeinsames Gerüst
über alle Stepper statt eines Sonderfalls je Verfahren.

Das Abnahmekriterium aus `dev/adjoint-plan.md` gilt unverändert: **jede Kombination, die
vorwärts eine Ableitung liefert, liefert auch rückwärts eine.**

---

## Die Architektur

### Was der Solver wirklich braucht

**bdf und adams** (`cppde_multistepper.hpp`, ein Template, zwei Koeffizientensätze).
Getapet ist genau **ein Residuum**:

    res = (y − zn0) + rl1·zn1 − γ·f(y, t_new)          multistepper.hpp:1348-1350

davor der Nordsieck-Rescale und der Pascal-Shift (`ndfRescale`, `ndfPredict`,
`:1318-1327`), danach der Tail `zn[j] += l[j]·acor` (`:1493-1617`). **Alle Koeffizienten
sind bereits plain double** (`rl1`, `γ`, `l[]`, `tq[]`, `m_h`, `:1310-1341`), weil
`time_type` skalarisiert ist. Damit ist der Adjungierte:

- `∂res/∂y` ist `W = I − γJ`, genau die Matrix, die `equation_solver` bereits
  faktorisiert und transponiert löst (`cppde_reverse_step.hpp:320-355`).
- `∂res/∂zn_j` sind feste ganzzahlige beziehungsweise diagonale Abbildungen.
- `∂res/∂θ = −γ·(df/dp)`. **Die einzige neue erzeugte Kontraktion, die bdf und adams
  brauchen.** `J'λ` wird nicht separat gebraucht, es steckt im Solve.
- Der Tail ist ein Rang-1-Update mit bekannten Koeffizienten.
- Der Dense-Output-Adjungierte ist **schon heute geschlossen**: Horner über `zn_dense`
  mit skalarer Zeit (`multistepper.hpp:2011-2013`), eine feste `[n × (q+1)]`-Gewichtung.

**rb4** (`cppde_rosenbrock4.hpp`). Sechs Stufenlösungen gegen dieselbe Faktorisierung,
`t` und `dt` bereits skalarisiert (`:415-417`). Der Adjungierte ist die transponierte
Stufenrekursion: sechs `W'`-Solves, dazu `J'λ` und `(df/dp)'λ` je Stufe. Der Interpolant
ist kubisch in `s` mit `cont3`/`cont4` aus `g1..g5`.

**tsit5** (`cppde_tsit5.hpp`). Sieben explizite Stufen, Lehrbuch-Rückwärtsrekursion
`λ_i = b_i λ_out + h Σ_{j>i} a_ji J(x_j)' λ_j`, also `J'λ` und `(df/dp)'λ` je Stufe.
Hermite-Interpolant explizit. Hier ist `h` heute als einziges symbolisch
(`ad_traits::step_coef`, `cppde_ad_traits.hpp:110-139`); mit dem eingefrorenen Gitter
fällt das weg und `step_coef` samt Fallunterscheidung verschwindet.

**Fazit: zwei erzeugte Kontraktionen tragen alles.**

    jac_t_vec(x, t, params, λ)  ->  J'λ          rb4 und tsit5
    dfdp_t_vec(x, t, params, λ) ->  (df/dp)'λ    alle vier

Beide existieren in `codegen_cvode.py` als `adj_rhs_fn` (`:1929-1949`) und `adj_quad_fn`
(`:1951-1969`). Sie werden nach `codegen_cppODE.py` gehoben, damit beide Backends
dieselbe Ableitung benutzen.

### Wo das im Baum sitzt

Die Schrittalgebra ist modellunabhängig und gehört in die Header-Bibliothek, das Modell
liefert nur die Kontraktionen. Der Aufhängepunkt existiert: `step_checkpoint<Stepper, T>`
(`cppde_reverse_step.hpp:81`, spezialisiert bei `:142` tsit5, `:152` rb4, `:180`
multistepper) ist die per-Verfahren-Policy. Daneben tritt `adjoint_step<Stepper>`, das
dieselbe Schnittstelle bedient, die `step_recorder` heute nach außen zeigt:

    begin / independent / load / attempt* / interpolate / seed / seed_carry /
    sweep / accumulate / wx / whistory / wt / wdt / wout / error_scale / xerr / stepper

Damit bleibt `cppde_reverse_trajectory.hpp` zu etwa 85 Prozent unverändert: der
Checkpoint-Speicher, die Beobachtungsliste, die äußere Schleife, die
Parameter-Akkumulation, das Lambda- und eta-Tracing, `clamp_to_step`. Verfahrensspezifisch
sind nur `replay_one` (22 Zeilen, `:533-555`), die `snapshot_family` des Collectors
(`:206-263`) und der Nordsieck-Restart nach einem Ereignis (`:615-632`).

### Ereignisse

Heute wird die Sprungarithmetik in `cppde_saltation.hpp` unter `codual` nachgerechnet;
**eine Saltationsmatrix wird nie gebildet**, ihre Transponierte entsteht implizit. Die
Ereigniszeit kommt aus einem IFT-Quotienten `compute_dt_star` (`:36-91`) mit dem Trick
eines numerisch verschwindenden, aber tape-lebendigen Shifts (`:59-61`).

Geschlossen heißt: die Saltationsmatrix wird gebildet und transponiert angewandt, aus
Größen, die `compute_dt_star` bereits ausrechnet. Der Zero-Shift-Trick entfällt mit dem
Tape.

---

## Stufen

Geliefert wird in einem Zug: alles bauen, am Ende einmal umschalten. Die Stufen sind
Entwicklungsreihenfolge, nicht Auslieferungsschnitte.

### Stufe 0. Zwei Messungen, bevor eine Zeile entsteht

- *Ausdrucksgröße.* `df/dp` und die beiden Kontraktionen für Bachmann (25 Zustände,
  113 Parameter) und Lang (124 Zustände, 294 Parameter) erzeugen, ohne sie zu benutzen,
  und Quellgröße und Übersetzungszeit gegen den heutigen Stand messen. Das einzige
  Risiko, das den Entwurf kippen kann.
- *Kostenmodell prüfen.* Am `long`-Profil nachrechnen, welcher Anteil der 8,8 µs auf das
  Replay der Nordsieck-Arithmetik entfällt und welcher auf die eine RHS-Auswertung. Die
  Schätzung "gut eine Mikrosekunde" steht und fällt damit.

### Stufe 1. Die Verifikationslöcher schließen, vor dem Umbau

Der Umbau ist genau dort am riskantesten, wo heute nicht geprüft wird. Drei Löcher, alle
in `dev/cxx/`:

- **Ereignisse decken nur bdf und tsit5 ab** (`test_reverse_events.cpp:436-437`). rb4 und
  adams haben in C++ **keine** Ereignisabdeckung. Genau die beiden bekommen in Stufe 4
  eine neu geschriebene Saltation.
- **Forcings haben null C++-Abdeckung.** Nur ein R-Test, auf bdf. `J'λ` ist betroffen,
  sobald eine Forcing multiplikativ auftritt.
- **Der dünne Pfad läuft nie durch den Reverse-Stepper.** Der transponierte Kern ist
  geprüft (`test_sparse_transpose.cpp`), ein KLU-Modell rückwärts nicht.

Dazu die Kleinigkeiten: `bench_codual.cpp:12` und `bench_revmem.cpp:14` nennen
`run.sh`-Flags, die es nicht gibt (`run.sh:21-71`); sie werden verdrahtet, weil der Bench
die Wirkung des Umbaus belegen muss. Und `run.sh --record` läuft einmal über alle
Prüfstände, um den Ist-Stand als Referenzausgabe festzuhalten.

### Stufe 2. Die erzeugten Kontraktionen

`df/dp` und die beiden Kontraktionen aus `codegen_cvode.py` nach `codegen_cppODE.py`
heben, mit CSE (`_cse_temps`, `:982-1006`), über `ScalarType` templatisierbar. Neuer
Schlüssel im Rückgabe-Dict (`:893-906`), Splice in `R/cppODE.R:1550-1557`. Der
CVODE-Emitter importiert danach aus derselben Quelle statt eine eigene Ableitung zu
führen. Dabei fällt der fehlende `DiracDelta`-Guard auf (siehe Risiken).

### Stufe 3. `adjoint_step<Stepper>`, Verfahren für Verfahren

- 3a bdf und adams: Residuum, Nordsieck-Rescale, Pascal-Shift, Tail, Dense-Output.
- 3b rb4: sechs Stufenlösungen gegen dieselbe Faktorisierung.
- 3c tsit5: explizite Rückwärtsrekursion und Hermite-Interpolant.

Jede gegen den Vorwärtsmodus geprüft, bevor die nächste beginnt.

### Stufe 4. Ereignisse und Wurzeln

Saltationsmatrix bilden und transponiert anwenden, Ereigniszeit-Ableitung aus
`compute_dt_star` ohne Zero-Shift-Trick, feste Ereignisse und Wurzelereignisse getrennt.

### Stufe 5. Der Nordsieck-Restart

`replay_boundary`s Multistep-Zweig (`cppde_reverse_trajectory.hpp:615-632`) läuft heute
über `initialize()` auf dem codual-Stepper. `initialize` baut jeden Slot aus einem
Zustand, der Adjungierte kollabiert also auf diesen Zustand, und das ist hinschreibbar.

### Stufe 6. `cppFUN` bekommt einen symbolischen vjp

Heute instanziiert `_write_vjp_impl` (`codegen_cppFUN.py:946-1018`) den Ausdruckskörper
ein zweites Mal über `codual`. Ersetzt durch eine erzeugte Kontraktion `w' ∂f/∂x` und
`w' ∂f/∂p` aus `derivSymb.jac_hess_symb` (`derivSymb.py:325-404`), mit CSE (`_cse_exprs`,
`:438-449`). Die R-Schnittstelle `$vjp(x, p, W)` bleibt zeichengleich, dMod2 merkt nichts.

### Stufe 7. Umschalten und löschen

`cppde_codual.hpp`, `cppde_codual_math.hpp`, `cppde_codual_tape.hpp` entfallen. Mit ihnen
entfällt die **zweite Erzeugung des Modellrumpfs**: `R/cppODE.R:296-312` und `:346-359`
rufen den Generator heute ein zweites Mal mit `num_type = "cppde::codual<double>"`; der
geschlossene Adjungierte braucht `f` und `J` nur in `double`, also fällt `namespace rev_`
(`:1545-1548`) weg. Nebenbei weniger Übersetzungszeit je Reverse-Modell.

### Stufe 8. Zweite Ordnung

Der geschlossene Adjungierte wird über den Skalartyp instanziiert, `dual<double,N>` statt
`double`, und liefert Hesse-Vektor-Produkte über die vorhandenen Dual-Zahlen. Damit
entfällt `codual<dual<...>>` als Weg dorthin. Offener Punkt bleibt der KLU-Solve mit
dual-wertiger Iterationsmatrix: entweder ein dual-fähiger dünner Solve oder ein dichter
Rückfall für zweite Ordnung.

---

## Verifikation

Das Orakel ist der Vorwärtsmodus. Der Wechsel kostet fast nichts, denn **die
C++-Prüfstände vergleichen schon heute gegen den Vorwärtsmodus**, nicht gegen das Tape:
sie bauen mit `cppde::dual` die volle Schrittmatrix `S` und prüfen `w'S` gegen `S'w`.

| Prüfstand | prüft | Toleranz | Verfahren |
|---|---|---|---|
| `test_reverse_step.cpp` | ein tsit5-Schritt | 1e-14 | tsit5 |
| `test_reverse_step_rb4.cpp` | ein rb4-Schritt | 1e-11 | rb4 |
| `test_reverse_step_multistep.cpp` | ein Schritt auf dem Nordsieck-Carry, jeder Slot einzeln geseedet | 1e-9 | bdf, adams |
| `test_reverse_trajectory_methods.cpp` | ganze Trajektorie, dazu finite Differenzen als dritte Quelle | 1e-6 bis 1e-10 | alle vier |
| `test_reverse_events.cpp` | Trajektorie mit Sprüngen | 1e-6 / 1e-9 | bdf, tsit5 |

`test_reverse_step_multistep.cpp` ist das schärfste Orakel im Baum und die einzige Stelle,
an der bdf und adams bei 1e-9 statt bei 1e-6 geprüft werden. Stufe 3a wird daran
gemessen.

`run.sh --record F` / `--against F` diffed zwei Stände byteweise. Der geschlossene
Adjungierte rechnet in anderer Reihenfolge, wird also nicht bitgleich sein; der Diff
bleibt trotzdem das Werkzeug, das jede *unbeabsichtigte* Änderung zeigt.

Der ganze Reverse-Satz kostet zwei bis drei Minuten Wandzeit, fast ausschließlich
Übersetzung, dazu der ASan-Durchgang.

Auf R-Ebene bleiben `cppDE/tests/testthat/test-reverse.R` und
`dMod2/tests/testthat/test-reverse.R` unverändert bestehen, ebenso
`test-cppFUN.R:244-320` für Stufe 6. Auf Kettenebene läuft
`dMod2/inst/benchmarks/bench_adjointPerCondition.R` auf Bachmann: der Winkel `1 − cos`
gegen den Vorwärtsgradienten darf sich nicht verschlechtern, heute 3,55e-07.

**Abnahme, in dieser Reihenfolge:**

1. Der Rückwärtslauf kostet höchstens das Doppelte der Vorwärtsintegration, je Schritt,
   auf allen vier Verfahren, gemessen mit `cppODE(profile = TRUE)`. Heute 12,5x. Der
   Bericht geht nur nach stderr, das reicht.
2. Die Bachmann-Kette über alle 36 Bedingungen ist schneller als CVODES ASA, heute 336
   gegen 217 ms. Fällt sie nicht, liegt es an der Kette darüber, und das ist ein eigener
   Befund, kein Fehlschlag dieses Plans.

## Risiken

- **Ausdrucksschwellung.** Die Kontraktion für Lang ist `124 × 294` Partialableitungen.
  Stufe 0 misst das zuerst. Rückfall: die Kontraktion nicht ausmultiplizieren, sondern
  das dünne `df/dp` als Datentabelle emittieren und zur Laufzeit kontrahieren, wie es die
  dünne Jacobi schon tut.
- **Nicht glatte Konstrukte sind besser gestellt als gedacht.** `piecewise`, `min`, `max`,
  `abs` und die logischen Operatoren werden bereits symbolisch abgeleitet
  (`cppsympy.py:37-46`, `:81-93`), nicht über AD. "Symbolisch oder Fehler" ist damit
  größtenteils erfüllt. Eine Lücke: `Heaviside` bricht im ODE-Emitter ab, weil
  `DiracDelta` keinen Drucker hat; `codegen_cppFUN.py:277-291` hat den Guard,
  `codegen_cppODE.py` nicht. Gehört in Stufe 2.
- **Forcings sind der wunde Punkt.** Sie werden nie symbolisch abgeleitet; die Kettenregel
  liegt an zwei Stellen von Hand (`codegen_cppODE.py:1129-1137`, `:1533-1548`).
  `(df/dp)` ist nicht betroffen, `J'λ` schon, sobald eine Forcing multiplikativ auftritt.
  Deshalb steht die Forcing-Abdeckung in Stufe 1 und nicht später.
- **`step_coef` und `CPPDE_SYMBOLIC_STEPSIZE` entfallen.** Damit fällt auch der
  Untersuchungspunkt "ein symbolisches `t` in den Saltationskorrekturen" aus
  `dev/adjoint-plan.md` weg: ohne Tape gibt es keinen symbolischen Zeittyp mehr. Und
  `test_reverse_step.cpp` verliert seine `wt`/`wdt`-Prüfung, die einzige im Baum.
- **Zweite Ordnung mit KLU.** Siehe Stufe 8.

## Nebenbefunde, die nicht in diesen Plan gehören

Zwei echte Fehler, bei der Erkundung gefunden, beide unabhängig vom Umbau:

1. **Eine multiplikativ auftretende Forcing erzeugt nicht übersetzbaren C++.**
   `_generate_jac_code_plain` reicht `forcings_list = []` an `_to_cpp` für die
   Jacobi-Einträge (`codegen_cppODE.py:1118`, `:1123`), aber die echte Liste für `dfdt`
   (`:1128`, `:1134`). Bei `A' = u*A − k*A` mit `forcings = "u"` entsteht
   `J(0,0) = -(-params[1]+u);` mit einem nirgends deklarierten `u`. Latent, weil jedes
   Beispiel und jeder Test Forcings additiv verwendet.
2. **`Heaviside` bricht den ODE-Generator ab**, `PrintMethodNotImplementedError` für
   `DiracDelta`.

## Dokument

Der Plan lebt als `cppDE/dev/closed-adjoint-plan.md`, mit einem Verweis aus
`dev/adjoint-plan.md`. Das bestehende Dokument bleibt als Aufzeichnung dessen stehen, was
gebaut wurde, und bekommt an der Festlegung "Reverse-Skalartyp statt handgeschriebener
Gleichungen" eine Zeile, dass sie gefallen ist und warum.

`dev/methods/Methods.Rmd` Abschnitt 2.9 (`:1090-1290`) und die Zeile über den
`cppFUN`-vjp (`:1438`) beschreiben das Tape und werden mit Stufe 7 neu geschrieben.
