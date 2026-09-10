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

**Erledigt am 2026-09-10. Beide fallen günstiger aus als der Entwurf angenommen hat.**

*Ausdrucksgröße.* `df/dp` und die beiden Kontraktionen symbolisch abgeleitet und
emittiert, ohne sie zu benutzen:

| | Zustände | Parameter | `df/dx` | `df/dp` | Ableitung |
|---|---|---|---|---|---|
| Bachmann | 25 | 29 | 89 | 67 | 0,1 s |
| Lang | 124 | 218 | 1175 | 1394 | 2,4 s |

Quellgröße, gegen das, was das Modell heute ohnehin emittiert:

| | rhs + Jacobi heute | `J'λ` + `(df/dp)'λ` neu | Verhältnis |
|---|---|---|---|
| Bachmann | 7,5 KB | 11,8 KB | 1,6 |
| Lang | 132,4 KB | 202,4 KB | 1,5 |

**Die Übersetzungszeit fällt, statt zu steigen.** Auf Lang, eine
Übersetzungseinheit, Minimum aus drei Läufen:

| | Quelle | `g++ -O2` |
|---|---|---|
| die zweite Erzeugung über `codual`, wie heute | 171,3 KB | 5,1 s |
| die beiden Kontraktionen, in `double` | 202,8 KB | 1,9 s |

Neunzehn Prozent mehr Quelle, 2,7-mal schneller übersetzt. Die
Expression-Templates des Tape-Typs sind teurer als der Text, den sie ersetzen.
Bei der Messung ist aufgefallen, dass eine nicht instanziierte codual-Struktur
968 Byte Objektcode erzeugt statt 210 KB: die Elementfunktionen sind implizit
inline, ein Vergleich ohne Treiber misst nichts.

*Kostenmodell.* Die beiden Kontraktionen auf Bachmann, je Aufruf, `-O2`:

| | ns |
|---|---|
| `J'λ`, 89 Einträge | 10,1 |
| `(df/dp)'λ`, 67 Einträge | 7,9 |
| beide zusammen | 18,0 |

Achtzehn Nanosekunden gegen 180 ns für den transponierten Solve und 740 ns für
`rev_prepare`. **Der geschlossene Schritt-Adjungierte ist damit von der linearen
Algebra dominiert, nicht von den Kontraktionen**, und landet bei etwa einer
Mikrosekunde gegen 8,8 heute. Das ist Faktor neun auf dem Rückwärtslauf und
setzt ihn auf rund das 1,3-fache der Vorwärtsintegration, besser als die
angesetzten 1,7.

**Der nächste Hebel danach steht damit auch fest:** `rev_prepare` kostet je
Schritt das Vierfache des Solves, obwohl sich die transponierte Struktur
zwischen zwei Faktorisierungen nicht ändert.

### Stufe 1. Die Verifikationslöcher schließen, vor dem Umbau

**Erledigt am 2026-09-10.** Drei Löcher, alle in `dev/cxx/`, alle zu. Der heutige
Reverse-Pfad besteht die neue Abdeckung überall dort, wo er sie überhaupt
annehmen kann, und an der einen Stelle, wo er es nicht kann, war es ein echtes
Loch im Vollständigkeitskriterium.

**Ereignisse decken jetzt alle vier Verfahren ab** statt bdf und tsit5.
`test_reverse_events.cpp` bekommt eine `pipeline`-Spezialisierung für rb4 und
zwei Läufe mehr. adams, rb4 und tsit5 bestehen auf Anhieb. Tape-Breite je
Schritt, nebenbei gemessen: bdf 249 Knoten, adams 306, rb4 576, tsit5 463.

**Forcings haben eine eigene Abdeckung**, `test_reverse_forcing.cpp`, alle vier
Verfahren, mit einer Forcing multiplikativ in zwei Gleichungen und additiv in
einer dritten. Das ist die Form, in der sie in `df/dx` eingeht statt nur in
`dfdt`, und der Reverse-Pfad trägt sie. Die Carry-Handoff-Prüfung steht dort
auf 1e-7 statt 1e-9, und der Grund ist das Orakel: unter `dual` maximiert die
Abbruchregel des Korrektors über die Sensitivitätsspalten und verlässt die
Iteration an einer leicht anderen Stelle als der Doppelt-Lauf, die Forcing
verstärkt das auf etwa 3e-9. Keine einzige Adjungierten-Prüfung ist davon
betroffen.

**Der dünne Pfad läuft jetzt rückwärts**, `test_reverse_sparse.cpp`, dasselbe
Modell und dieselben Toleranzen wie der dichte Prüfstand, nur mit
`sparse_lu_tag`. bdf und adams bestehen.

**Ein Fund dabei: rb4 hat auf einer dünnen Jacobi keinen Rückwärtspfad.**
`rosenbrock4::replay_step` baut sich eine eigene dichte Jacobi und bildet das
Residuum über jeden Eintrag (`cppde_rosenbrock4.hpp:430-446`), also übersetzt
ein dünnes rb4-Reverse-Modell nicht. An der R-Oberfläche bestätigt: bdf dünn
rückwärts übersetzt, rb4 dünn rückwärts bricht im Compiler ab. `cppODE()`
verweigert die Kombination jetzt vorher, mit Test, und Stufe 3b hebt die
Verweigerung wieder auf, weil sie das Replay ohnehin ersetzt.

Dazu die Kleinigkeiten: `--bench-codual` und `--bench-revmem` sind verdrahtet,
sie waren in beiden Benches dokumentiert und existierten nicht. Und
`run.sh --record` hat den Ist-Stand aller elf Prüfstände festgehalten.

### Stufe 2. Die erzeugten Kontraktionen

**Erledigt am 2026-09-10.** `generate_ode_cpp(emit_contractions = TRUE)` emittiert
`struct adjoint_terms` mit `jac_t_vec` und `dfdp_t_vec`, beide mit CSE, beide
über `ScalarType` templatisierbar, beide nach Ausgabeslot gruppiert statt
akkumuliert. `R/cppODE.R` reicht das für jedes Reverse-Modell durch, in `double`,
neben dem Wertkörper. Ein explizites Verfahren emittiert keine Jacobi und braucht
für die Integration auch keine, leitet sie für die Kontraktion aber trotzdem ab.

Nachgerechnet, nicht behauptet: `J'λ` gegen die emittierte Jacobi auf 1e-13,
`(df/dp)'λ` gegen zentrale Differenzen auf 1e-6, der Anfangswertblock exakt null.

**Zwei Generatorfehler dabei behoben, beide älter als dieser Plan.**

*Eine multiplikative Forcing erzeugte nicht übersetzbaren C++.*
`_generate_jac_code_plain` reichte `forcings_list = []` an den Drucker für die
Jacobi-Einträge und die CSE-Temps, aber die echte Liste für `dfdt`. Eine Forcing,
die die Differentiation überlebt, kam damit als undeklarierter Bezeichner heraus.
Latent, weil jedes Beispiel und jeder Test sie additiv verwendet, wo sie aus
`df/dx` verschwindet. Regressionstest in `test-ode-methods.R`.

*`Heaviside` brach den Generator ab.* `DiracDelta` hat keinen Drucker. Jede
Ableitung des ODE-Emitters geht jetzt durch denselben Guard, den
`codegen_cppFUN.py` seit je hat. Regressionstest in `test-piecewise.R`.

**Offen geblieben:** der CVODE-Emitter führt seine `df/dp`-Ableitung weiterhin
selbst. Das Zusammenlegen gehört in Stufe 7, wenn ohnehin an beiden Emittern
gearbeitet wird.

### Stufe 3. `adjoint_step<Stepper>`, Verfahren für Verfahren

**3a bdf und adams: erledigt am 2026-09-10.** `cppde_adjoint_step.hpp`.

Der Entwurf hat sich beim Bauen an einer Stelle gedreht, und zum Besseren. Ein
Schritt zerfällt in

    zn_pred = A zn_in                    Rescale und Pascal-Shift
    res(y, zn_pred, theta) = 0           der Korrektor, über die IFT
    zn_out  = B zn_pred + c acor         der Schwanz, acor = y - zn_pred[0]

und `A`, `B`, `c` wirken **auf den Slot-Index allein**, mit Koeffizienten aus dem
Carry, für jede Zustandskomponente gleich. Sie werden deshalb nicht von Hand
transponiert, sondern **aus dem Stepper selbst gelesen**: ein Probe-Stepper,
dessen Zustände die Slots sind, bekommt die Einheitsmatrix ins
`[Slot x Komponente]`-Feld und läuft durch dieselben Routinen wie der
Vorwärtsschritt. Damit kann hier nichts von der Implementierung abweichen, weil
hier nichts wiederholt wird: geschrieben sind nur die Transponierten und die
IFT.

Dafür ist `replay_predict` aus `replay_residual` herausgezogen worden, die
Vorstufe ohne die Auswertung der rechten Seite. Eine Definition, zwei Aufrufer.

**Geprüft** im vorhandenen `test_reverse_step_multistep.cpp`, dem schärfsten
Orakel im Baum: derselbe Schritt, dieselbe Dual-Referenz, jeder Nordsieck-Slot
einzeln geseedet, Aufwärmlängen 3/8/20/45/80, Ordnung hoch und runter, Rescale,
bdf und adams. Der geschriebene Adjungierte trifft die Referenz überall auf
1e-9.

**Gemessen**, `dev/cxx/bench_adjoint_step.cpp`, derselbe Schritt:

| | Vorwärtsschritt | Tape | geschrieben |
|---|---|---|---|
| bdf | 0,141 µs | 1,603 µs (11,4x) | 0,320 µs (2,27x) |
| adams | 0,157 µs | 2,185 µs (13,9x) | 0,319 µs (2,03x) |

Drei Dinge haben den Weg dorthin gebaut, jedes gemessen statt geraten. Alle
Spalten in einem Durchlauf statt einer je Spalte, weil die Operatoren
komponentenweise gleich wirken. Der Probe-Stepper einmal angelegt statt je
Schritt, denn sein Konstruktor allokiert jeden Slot, den er je brauchen könnte,
und das war zuerst 78 Prozent des ganzen Pfades. Und ein Zwischenspeicher auf
dem Carry: über lange Strecken hält ein Lauf Ordnung und Schrittweite, und dann
sind es dieselben Matrizen.

Ohne Zwischenspeicher liegt der geschriebene Pfad bei 0,64 µs, also immer noch
2,5-mal schneller als das Tape. Die Trefferquote auf einem echten Lauf ist noch
zu messen.

**3c tsit5: erledigt am 2026-09-10.** Die Rückwärtsrekursion eines expliziten
Verfahrens steht in seiner Tableau und sonst nirgends, also gibt `tsit5` seine
Koeffizienten heraus statt dass der Adjungierte eine zweite Abschrift trägt.
Die Stufenzustände sind nicht gecheckpointet, also läuft der Schritt einmal in
`double` vorwärts, um sie zu holen, und die Rekursion kontrahiert danach `J'`
und `(df/dp)'` je Stufe. Stufe 7 ist FSAL und trägt keine Lösung, die Rekursion
läuft also über sechs.

Geprüft in `test_reverse_step.cpp` über einen und über vier Schritte, gegen die
eingefrorene Dual-Referenz. Der Prüfstand beißt: eine relative Störung von 1e-6
an einem einzigen Koeffizienten bricht vierzehn Zusicherungen.

**3d. Die Verdrahtung in die Trajektorie. Fehlte im Plan.**

Die Stufenliste sprang von den Schritt-Adjungierten zu den Ereignissen, aber
dazwischen liegt das Stück, ohne das kein Solve den geschriebenen Adjungierten
je benutzt und die zweite Abnahmezahl gar nicht messbar ist.

`closed_trajectory<Stepper>` läuft denselben `trajectory_store` rückwärts wie
`trajectory_recorder::sweep`, aber ohne Tape. Je Schritt: die Beobachtungen im
Schritt über den Dense-Output einseeden, den Schritt-Adjungierten anwenden, den
Carry an den vorigen Schritt reichen, `wp` aufsummieren.

Der Dense-Output-Adjungierte fällt dabei mit ab, statt hingeschrieben zu werden:
der Tail-Probe hat nach seinem Durchlauf einen gültigen Interpolanten über `B`,
also gibt `eval_dense_into(t_obs, .)` auf ihm unmittelbar die Zeile
`d x_interp / d (zn_pred, acor)`. Eine Auswertung je Beobachtung, exakt, und
wieder nichts wiederholt.

Zuerst ohne Ereignisse, was für die Abnahmemessung reicht, denn Bachmann hat
keine. Die Ereignisgrenzen kommen mit Stufe 4.

- 3b rb4: sechs Stufenlösungen gegen dieselbe Faktorisierung. Hebt zugleich die
  Verweigerung von dünn plus rb4 plus reverse aus Stufe 1 wieder auf.

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
