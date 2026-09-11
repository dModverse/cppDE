# Der geschlossene Adjungierte: das Tape aus dem Solver nehmen

**Nachfolgeplan zu `dev/adjoint-plan.md`, angelegt 2026-09-10.** Jener Plan hat den
diskreten Adjungierten gebaut und ist abgearbeitet. Dieser nimmt ihm das Tape ab.

**Stand 2026-09-11: die Stufen 0 bis 7 sind abgearbeitet.** Ein Reverse-Modell
zeichnet nichts mehr auf, auf keinem Verfahren und mit keiner Art von
Intervention; das Wort `codual` kommt in cppDE nicht mehr vor. Offen ist Stufe 8,
die zweite Ordnung.

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

## Wo es nach 3d steht

**Gemessen am 2026-09-11**, gleiche Maschine, gleicher Kern, gleiche Toleranzen.
Bachmann, ganze Kette bis `normL2`, 36 Bedingungen, 113 Parameter:

| Route | ms | in Wertläufen | vorher |
|---|---|---|---|
| Wertlauf | 35,4 | 1,00 | 36,5 |
| Vorwärts | 979,1 | 27,6 | 1026,4 |
| Adjungierter, geschrieben | 119,6 | 3,37 | 335,7 mit Tape |
| CVODES ASA | 205,5 | 5,80 | 216,6 |

Der Gradient kostet jetzt 3,37 Wertläufe statt 9,19 und liegt **um das
1,72-fache unter ASA**, wo er vorher um das 1,55-fache darüber lag. Der Winkel
gegen den Vorwärtsgradienten steht unverändert bei 3,55e-07.

Je Schritt, Bedingung `long`, 721 Schritte, `cppODE(profile = TRUE)`:

| | ms | je Schritt | Anteil |
|---|---|---|---|
| Vorwärtsintegration, ganzer Stepper | 0,530 | 0,735 µs | |
| `rev_prepare`, Jacobi und Faktorisierung | 0,411 | 0,570 µs | 39% |
| `rev_adjoint`, die Schrittalgebra | 0,415 | 0,576 µs | 39% |
| `rev_solve`, der transponierte Solve | 0,126 | 0,175 µs | 12% |
| `rev_operators`, die Operatoren lesen | 0,101 | 0,140 µs | 10% |
| **Rückwärtslauf** | **1,056** | **1,465 µs** | **1,99x** |

Damit ist Abnahmezahl 1 erreicht, gegen 12,5x zu Beginn, und Abnahmezahl 2 auch.
Die lineare Algebra ist jetzt die Hälfte des Rückwärtslaufs, die Schrittalgebra
die andere. Das Tape ist aus dem Profil verschwunden.

**Die Trefferquote des Operator-Zwischenspeichers lag bei null**, und der Grund
war eine Zeile: `qwait` stand im Schlüssel. Es zählt herunter bis zur nächsten
Ordnungsentscheidung und wechselt darum fast jeden Schritt, erreicht aber nur
die Fehlerkonstanten für hoch und runter. Adams nimmt es entgegen und liest es
nie, NDF benutzt es für `tq[1]` und `tq[3]`. Beides gehört dem Regler, der nicht
differenziert wird. Ohne die Zeile: 142 Neuaufbauten auf 721 Schritte, und
`rev_operators` fällt von 0,317 auf 0,101 ms.

**Zwei Fehler haben die Messung selbst gefunden, beide durch die Abdeckung
gerutscht.**

*Die erzeugte Kontraktion dimensionierte ihre Ausgabe nicht, die handgeschriebene
schon.* Zwei Verträge für dieselbe Funktion, und die Prüfstände sahen nur den
einen. Die Startgrenze der Trajektorie reichte `jac_t_vec` den Puffer der
Dense-Output-Zeile, der einen Eintrag je Nordsieck-Slot lang ist. Bei zwei
Zuständen passt das, bei fünfundzwanzig nicht. Der Generator dimensioniert jetzt
selbst, und ein Prüfstand mit mehr Zuständen als Slots steht in `test-reverse.R`.

*Der geschriebene Pfad hatte kein Profil.* Der Berichtaufruf stand nur im
getapeten Zweig, also war Abnahmezahl 1 auf ihm gar nicht messbar. Zwei
Kategorien dazu, `rev_operators` und `rev_adjoint`, und der Zeitnehmer sitzt am
Neuaufbau statt am Aufruf, damit die Aufrufzahl die Fehlzugriffe zählt.

**Der nächste Hebel ist damit benannt und beziffert:** `rev_prepare`, 39 Prozent
des Rückwärtslaufs.

## Wo es nach der Leistungsrunde steht

**Gemessen am 2026-09-11**, nach 3d, dem Einschritt-Teil und einer Runde auf den
Kosten. Dieselbe Kette, dieselbe Maschine:

| Route | ms | in Wertläufen | mit Tape |
|---|---|---|---|
| Wertlauf | 34,8 | 1,00 | 36,5 |
| Vorwärts | 974,6 | 28,0 | 1026,4 |
| Adjungierter, geschrieben | 106,6 | 3,06 | 335,7 |
| CVODES ASA | 201,1 | 5,77 | 216,6 |

**Der Gradient kostet drei Wertläufe und liegt um das 1,89-fache unter ASA.** Der
Winkel gegen den Vorwärtsgradienten steht unverändert bei 3,55e-07: keine der
Änderungen rechnet anders, sie rechnen dasselbe weniger oft.

Je Schritt, Bedingung `long`, 721 Schritte, Minimum über zwanzig Läufe im selben
Prozess, weil die Maschine zwischen zwei Läufen um ein Drittel driftet:

| | je Schritt | Anteil |
|---|---|---|
| Vorwärtsintegration, ganzer Stepper | 0,639 µs | |
| `lu_factor`, KLUs Refaktorisierung | 0,444 µs | 42% |
| `rev_adjoint`, die Schrittalgebra | 0,280 µs | 26% |
| `rev_solve`, der transponierte Solve | 0,155 µs | 15% |
| `rev_operators`, die Operatoren lesen | 0,094 µs | 9% |
| `jac_eval` und `w_build` | 0,086 µs | 8% |
| **Rückwärtslauf** | **1,065 µs** | **1,67x** |

Was die Runde gebracht hat, und woran sie gescheitert wäre:

**Die Kontraktion wird akkumulierend emittiert.** Jeder Aufrufer addierte ihr
Ergebnis in etwas, das er schon hatte, also hat die zuweisende Form einen Vektor
in der Breite des ganzen Parametersatzes genullt, für jeden Parameter einen Slot
geschrieben, ob er einen Term trägt oder nicht, und dem Aufrufer einen zweiten
Durchgang zum Addieren gelassen. Ein Drittel der Schrittalgebra.

**Der Zusammenbau von W ging durch einen BLAS-Aufruf**, für ein paar hundert
Doubles, dessen Vermittlung mehr kostet als die Kopie. Jetzt `memcpy`, und was
von `lu_factor` bleibt, ist zu 92 Prozent KLU selbst.

**Gegen KLUs eigenen Aufwand ist mit Iteration nichts zu holen, und das ist
Arithmetik.** Ein gehaltener Faktor plus Nachkorrektur konvergiert linear mit der
relativen Änderung von W, und die ist über einen Schritt zu groß, um in weniger
als drei Durchgängen auf 1e-13 zu kommen. Drei Durchgänge kosten mehr als die
Refaktorisierung. Bei 25 Zuständen sind 0,44 µs für 150 Flops fast nur Aufwand,
kein Rechnen; das hebt nur eine erzeugte LU mit fester Pivotfolge, und die
tauscht Genauigkeit gegen Zeit. Bleibt draußen.

**Die Kette über dem Solver ist inzwischen die Hälfte.** `.Call` ist 50 Prozent
einer Reverse-Zielfunktion, der Rest ist R, und davon war ein Viertel
Mengenoperationen auf Namen. `intersect` und `setdiff` aus base nullen,
vereinheitlichen und suchen; ein `match` reicht, weil Parameternamen eindeutig
sind. `.pickCotangent` fällt damit von 10 auf 3 Prozent, die ganze Zielfunktion
um zwölf.

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

**3d. Die Verdrahtung in die Trajektorie. Fehlte im Plan. Mehrschritt erledigt
am 2026-09-10, tsit5 am 2026-09-11, rb4 wartet auf 3b.**

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

*Einschritt, erledigt am 2026-09-11.* `closed_onestep_trajectory` läuft denselben
Speicher, mit einfacherem Carry: ein Einschrittverfahren reicht nur den Zustand
über eine Schrittgrenze, es gibt also keine Historie zu kollabieren und keine
Startgrenze. Der Dense-Output wird hingeschrieben statt probiert: tsit5
interpoliert Hermite-kubisch über `x_alt`, `x_neu`, `k1` und `k7` und gibt seine
vier Gewichte heraus, wie es schon seine Tableau herausgibt, also steht die Basis
einmal da und wird zweimal gelesen. Zwei der vier Gewichte sitzen auf `f`, eine
Beobachtung kostet damit je eine `J'`- und eine `(df/dp)'`-Kontraktion an beiden
Enden ihres Schritts. Der Probe-Stepper bleibt der Sonderfall des
Mehrschrittverfahrens, dessen Nordsieck-Operatoren man sonst abschreiben müsste.

rb4 hat noch keinen Schritt-Adjungierten und bleibt bis 3b auf dem Tape.

**3b rb4: erledigt am 2026-09-11. Es brauchte zweite Ableitungen, und der Plan
hat das nicht gesehen.**

Die sechs Stufen lösen gegen `W = I/(γh) − J(x, t)` mit `J` am Schrittanfang.
Jede Stufenlösung hängt damit über `W` von `x` und `θ` ab, nicht nur über ihre
rechte Seite: aus `g = W⁻¹ r` folgt `dg = W⁻¹(dr + (dJ) g)`, und der Adjungierte
dieses Terms ist die Kontraktion von λ mit der Ableitung von `J g`. Das ist eine
zweite Ableitung, die das Tape stillschweigend geliefert hat.

**Die Messung aus Stufe 0, auf der zweiten Ordnung: sie fällt kleiner aus als
die erste.**

| | Zustände | Parameter | nnz `J` | nnz `∂(Jv)/∂x` | Quelle 1. Ord | 2. Ord | `g++` |
|---|---|---|---|---|---|---|---|
| Bachmann | 25 | 54 | 89 | 55 | 13,4 KB | 19,5 KB | 0,16 / 0,16 s |
| Lang | 124 | 342 | 1175 | 1108 | 260,8 KB | 298,2 KB | 0,26 / 0,28 s |

Eine rechte Seite zweimal abzuleiten tötet Terme, statt welche zu erzeugen. Das
Risiko der Ausdrucksschwellung ist damit erledigt, und der Rückfall auf eine
Datentabelle wird nicht gebraucht.

**Der Generator emittiert vier weitere Kontraktionen, aber nur für rb4**, denn
nur dort werden sie gebraucht und ableiten ist nicht umsonst: `jvp_x_t_vec` und
`jvp_p_t_vec_axpy` über `J(x,p,t)·v`, und `dfdt_x_t_vec` und
`dfdt_p_t_vec_axpy` über `df/dt`. Letztere sind nötig, weil vier der Stufen ein
Vielfaches von `df/dt` addieren. Die Emission läuft durch dieselbe Routine wie
die erste Ordnung, mit einem anderen Körper: `_emit_contraction_pair`.

*Nebenbei:* `df/dt` trägt für jede Forcing einen Kettenterm, den der
Jacobi-Emitter als Text anhängt. Damit dessen Ableitung nach Zustand und
Parameter stimmt, wird die Zeitableitung der Forcing als erfundenes Symbol
geführt und differenziert wie die Konstante, die sie an dieser Stelle ist.

**Die transponierten Solves gehen durch die Faktorisierung des Steppers
selbst**, nicht durch eine zweite. Ein Korrektorverfahren darf das nicht, seine
Iterationsmatrix ist absichtlich veraltet; eine Rosenbrock-Stufe ist ein direkter
Solve, und die Matrix, mit der sie gelöst hat, ist die Matrix, die ihre Ableitung
braucht. Nachgesehen, dass das trägt: der Regler gibt `reuse_jacobian` nur nach
einem verworfenen Versuch heraus, und der begann am selben Zustand.

**Geprüft** in `test_reverse_step_rb4.cpp` über einen, zwei und vier Schritte
gegen die Dual-Referenz auf 1e-11, und in `test_reverse_trajectory_methods.cpp`
über die ganze Trajektorie samt Beobachtungen im Schritt. Der Prüfstand beißt:
die sechste Stufe wertet an `X_5 + g_5` aus, und `X_5` trägt die a-Koeffizienten
der fünften; mit einer eigenen Zeile für die sechste lag der Gradient um 1e-3
daneben.

**Die Verweigerung von dünn plus rb4 plus reverse fällt** auf die eine
Kombination zusammen, die noch nachspielt: dünn plus rb4 plus Ereignis. Dafür
wird der getapete Sweep nur noch dort emittiert, wo der geschriebene nicht
greift. Wo er greift, ist er nicht die Vorgabe, sondern der einzige Pfad.

### Stufe 4. Ereignisse und Wurzeln, mit dem Restart darin. Erledigt am 2026-09-11

**Der Sprung ist im Wert ein Sprung.** Die Maschine trägt eine Unstetigkeit über
eine Heun-Klammer: eine Schicht vorwärts auf die Ereignisfläche, die
Rücksetzungen, eine Schicht zurück auf die Gitterzeit. Geschoben wird um den
Rest der Ereigniszeit, und von dem wird der Wert abgezogen, bevor er benutzt
wird. Er ist also **numerisch exakt null**, und nur seine Ableitung überlebt.

Damit ist die Klammer im Wert die Identität und kollabiert auf die klassische
Saltation:

    dx = R'(dx_before + f_before s) − f_after s

mit `R` den Rücksetzungen und `s` dem Differential des Zeitrests. Für eine
Wurzel ist das `−(∇g·dx + ∂g/∂p·dp)/ġ`, für ein festes Ereignis das Differential
seiner Zeit. Die Zweitordnungskorrektur multipliziert `dt*` mit sich selbst und
fällt mit ihm weg.

**Ein Sprung-Adjungierter braucht darum keine Jacobi.** Er braucht die rechte
Seite an beiden Enden und die Ableitungen der Ereignisausdrücke; jeder Term, den
die Schichten beigetragen hätten, trägt einen Faktor, der null ist. Der Plan
hatte hier `J'` je Schicht angesetzt, also vier Kontraktionen, die nicht
gebraucht werden.

**Die Modellseite** steht in `event_adjoint_terms`, erzeugt neben den
Kontraktionen von `f`: je Ereignis die Ableitung seines Werts nach Zustand und
Parametern, für eine Wurzel die ihrer Bedingung nach den Parametern, für ein
festes Ereignis die seiner Zeit. Indiziert nach der Stelle, an der der
Vorwärtslauf das Ereignis einreiht, damit ein Index auf beiden Seiten dasselbe
meint.

**Der Restart ist ein Aufrufer geworden**, wie vorgesehen: `collapse_restart`
steht einmal da und wird von der Startgrenze und von jeder Ereignisgrenze
gerufen.

**Der Kotangens landet nicht immer auf einem Carry.** Der Zustand, der in den
Sprung geht, wurde vom Dense-Output des Schritts darunter gelesen, also geht
sein Kotangens dorthin und der Carry wird gar nicht gelesen: der Restart hat ihn
weggeworfen. Beide Trajektorien führen das als anstehende Interpolation, genau
wie der getapete Sweep.

**Damit fällt die Verweigerung von dünn plus rb4 plus reverse ganz**, denn der
getapete Sweep wird für ein Reverse-Modell überhaupt nicht mehr emittiert.

Geprüft in `test-reverse.R`: Ereignisse auf allen vier Verfahren, mit einem
festen Ereignis über eine Parameterzeit, einem Wurzelereignis multiplikativ und
einer Forcing daneben, dazu ein dünnes Modell mit Sprung auf bdf und rb4.

### Eine Beschränkung, die bleibt

**Die Mehrschrittverfahren laufen nur mit Dense-Interpolation**, vorwärts wie
rückwärts. Der Nordsieck-Interpolant ist dort nicht eine Ausgabeform neben
anderen, sondern das Verfahren selbst: die Historie *ist* der Interpolant.
`useDenseOutput = FALSE` wird für bdf und adams deshalb verworfen, mit einer
Warnung, und war es schon vor diesem Plan. Es gibt keinen Grund, dem
Rückwärtspfad einen Weg zu bauen, den der Vorwärtspfad nicht hat.

### Stufe 5. Der Nordsieck-Restart. Geht in Stufe 4 auf.

`initialize()` baut jeden Slot aus einem Zustand, der Adjungierte kollabiert also
auf diesen Zustand. Das steht seit 3d geschrieben, als Startgrenze der
Trajektorie (`cppde_adjoint_step.hpp`, das Ende von `sweep`). Nach einem Ereignis
gilt dieselbe Abbildung, Stufe 5 ist damit ein Aufrufer und kein Baustein.

### Stufe 5b. `rev_prepare`. Abgeschlossen am 2026-09-11, ohne den Vorschlag

Der Vorschlag war eine gehaltene Faktorisierung plus iterative Nachkorrektur.
**Er trägt nicht**, und das steht in der Leistungsrunde oben: die Nachkorrektur
konvergiert linear mit der relativen Änderung von `W` über einen Schritt, braucht
damit drei Durchgänge für die Genauigkeit eines Gradienten, und drei Durchgänge
kosten mehr als neu zu faktorisieren.

Was geblieben ist: die Montage von `W` ging durch einen BLAS-Aufruf für ein paar
hundert Doubles. Jetzt `memcpy`, und `lu_factor` ist zu 92 Prozent KLU selbst.
Bei dieser Modellgröße sind 0,44 µs für 150 Flops fast nur Aufwand; das hebt nur
eine erzeugte LU mit fester Pivotfolge, und die tauscht Genauigkeit gegen Zeit.
Bleibt draußen, solange das Kriterium lautet: keine Genauigkeit einbüßen.

### Stufe 6. `cppFUN` bekommt einen symbolischen vjp. Erledigt am 2026-09-11

`_write_vjp_impl` instanziierte den Ausdruckskörper ein zweites Mal über `codual`
und legte dafür ein Tape an. Jetzt kontrahiert es die symbolische Jacobi: dieselben
Einträge, die `_jacobian` emittiert, mit CSE, einmal je Beobachtung berechnet und
je Seed kontrahiert. Ein Aufrufer im Dual-Modus übergibt keine symbolische Jacobi,
also wird sie dort abgeleitet; die Rückwärtsrichtung ist so oder so ihre
Kontraktion.

Die R-Schnittstelle `$vjp(x, p, W)` bleibt zeichengleich, dMod2 merkt nichts.
`codegen_cppFUN.py` nennt `codual` nicht mehr.

### Stufe 7. Umschalten und löschen

**Vorbedingung, aus 3d gelernt: erledigt am 2026-09-11.** Die Lambda-Spur hing am
Tape, und `adjointGrid` wählte darum die Implementierung. dMod2 schaltet die Spur
ein, sobald `optionsReverse$gradtol` gesetzt ist, also schlossen sich der
adjoint-gewichtete Regler und der geschriebene Pfad aus. Beide Trajektorien
halten `lambda` und `eta` jetzt selbst; `wt` und `wdt` sind weg, sie waren die
Ableitungen nach Schrittzeit und Schrittweite, die das eingefrorene Gitter aus
der Kettenregel nimmt. Auf Bachmann kostet die gewichtete Zielfunktion 131 ms
gegen 122 ungewichtet, statt auf das Tape zurückzufallen.

**Die zweite Vorbedingung ist seit Stufe 4 erfüllt:** ein Reverse-Modell
emittiert den getapeten Sweep überhaupt nicht mehr. Was noch steht, steht
ungenutzt.

**Die Bestandsaufnahme**, damit das Löschen mechanisch bleibt:

| weg | bleibt, in derselben Datei |
|---|---|
| `cppde_codual.hpp`, `_math`, `_tape` | |
| `step_recorder` und seine Replay-Weichen | `equation_solver`, `step_checkpoint` |
| `trajectory_recorder`, `replay_one`, `replay_boundary`, `apply_jump` | `trajectory_store`, `event_record`, `observation` |
| `rosenbrock4::replay_step`, `replay_inv_gamma_dt` | die Tableau-Zugänge aus 3b |
| `multistepper::replay_residual`, `replay_outputs` | `replay_predict`, das der Probe-Stepper braucht |
| `step_checkpoint::load` und `finish` | `capture`, `apply_tail`, `start_state` |
| die codual-Zweige in `cppde_ad_traits.hpp` und `cppde_profiler.hpp` | |
| `dev/cxx/test_codual.cpp`, `bench_codual.cpp`, `bench_revmem.cpp` | |
| der getapete Vergleichspfad in den Reverse-Prüfständen | ihr Orakel, der Vorwärtsmodus |

Mit dem Tape entfällt die **zweite Erzeugung des Modellrumpfs**: `R/cppODE.R`
ruft den Generator ein zweites Mal mit `num_type = "cppde::codual<double>"` und
emittiert `namespace rev_` samt einem zweiten Satz Forcings. Der geschlossene
Adjungierte braucht `f` und `J` nur in `double`. Nebenbei weniger
Übersetzungszeit je Reverse-Modell, und nach der Messung aus Stufe 0 ist das der
größere Posten: die Expression-Templates des Tape-Typs kosteten auf Lang 5,1
Sekunden gegen 1,9 für mehr Text in `double`.

**Das Wort `codual` kommt in cppDE nicht mehr vor**, in keinem Header, keinem
Generator, keinem Prüfstand und keiner Zeile Vignette.

Zwei Dinge waren beim Löschen nicht offensichtlich. `replay_outputs` sieht nach
Tape aus, ist aber die Art, wie der Probe-Stepper `acor` bildet: es bleibt, mit
einem Kommentar, der sagt warum. Und `has_step_snapshot` stand im gelöschten
Teil, wird aber vom Checkpoint-Sammler gebraucht.

**Die Prüfstände sind mitgezogen worden, nicht weggeworfen.** Ihr Orakel war nie
das Tape, sondern der Vorwärtsmodus, also verlieren sie nur den zweiten Weg
dorthin. Drei von ihnen, Forcings, Ereignisse und der dünne Pfad, hatten bis
dahin gar keinen geschriebenen Vergleich; sie haben jetzt einen, mit von Hand
geschriebenen Kontraktionen, und `test_reverse_events.cpp` prüft damit den
Sprung-Adjungierten auf allen vier Verfahren gegen die Dual-Referenz.

Was ersatzlos entfällt: der Schritt-Adjungierte einer Beobachtung mitten im
Schritt, denn der Dense-Output-Adjungierte gehört der Trajektorie und wird dort
geprüft; die Kotangenten nach Schrittzeit und Schrittweite; und der Vergleich
der Lambda-Spur gegen eine zweite Quelle, die es nicht mehr gibt.

*Offen geblieben:* der CVODE-Emitter führt seine `df/dp`-Ableitung weiterhin
selbst. Das Zusammenlegen ist kein Teil des Löschens und wartet.

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
   auf allen vier Verfahren, gemessen mit `cppODE(profile = TRUE)`. Der Bericht geht nur
   nach stderr, das reicht. **Auf bdf erreicht am 2026-09-11, 1,99x gegen 12,5x zu
   Beginn.** Die anderen drei stehen aus, bis 3b und der Einschritt-Teil von 3d sie
   überhaupt auf den geschriebenen Pfad bringen.
2. Die Bachmann-Kette über alle 36 Bedingungen ist schneller als CVODES ASA.
   **Erreicht am 2026-09-11, 120 gegen 206 ms.**

`rev_prepare` wird dabei getrennt ausgewiesen: es ist Jacobi und Faktorisierung, also
dieselbe Arbeit, die der Vorwärtslauf seltener tut, und es skaliert mit `n³` statt mit
dem Adjungierten. Stufe 5b misst sich daran.

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
