# End-to-end diskreter Adjoint für cppDE und dMod2

## Kontext

**Der Grund ist die Skalierung der Vorwärts-Duals mit der Parameterzahl.** dMod2 bezieht Gradient
und Gauss-Newton-Hesse heute aus einem einzigen Objekt, den Vorwärtssensitivitäten. Ein
`dual<double, N>` trägt `1 + N` Doubles je Skalar, also skaliert jede Zustandsgröße, jede
RHS-Auswertung und das ganze Nordsieck-Array mit `n_theta`.

Am 2026-09-07 gemessen, und die Zahlen sind moderater, als die Struktur vermuten lässt:

- **Der Schrittweitenzuschlag ist real und begrenzt.** `wrms_max_ewt` (`cppde_newton.hpp:158-218`)
  nimmt in `finish()` das Maximum über die Zustandsnorm und jede einzelne Sensitivitätsspalte, und
  ein Maximum über eine wachsende Menge kann nur wachsen. Auf einer linearen Kette mit 40 Zuständen
  steigt die akzeptierte Schrittzahl von 105 ohne Sensitivitäten über 138 bei 2 Richtungen und 217
  bei 10 auf 305 bei 40, und bleibt dort für 60 und 80 stehen. Der Zuschlag ist also Faktor 3 und
  sättigt, sobald die ungünstigste Richtung im Maximum steht.
- **Die Arbeit je Richtung ist klein und fällt.** Eine Dual-Auswertung liefert Wert und alle
  Tangenten in einem Durchlauf mit geteilten Teilausdrücken.
- **Netto ist der Preis stark unterlinear.** `factor_grad` aus
  `dMod2/inst/benchmarks/bench_gradientCost.R`, der Preis eines Objective-Gradienten in reinen
  Wertlösungen, auf Bachmann: 3.5 bei 10 Parametern, 6.7 bei 25, 7.2 bei 50, 8.0 bei 75, 8.5 bei
  113. Elfmal so viele Parameter kosten das 2.4-fache. Auf Boehm mit 9 Parametern liegt er bei 2,
  dort aber an der Auflösungsgrenze der Uhr.

Der diskrete Adjoint hat keinen dieser Terme. Sein Vorwärtslauf ist eine reine Wertlösung in
`double`, also mit der Schrittfolge und der Arbeitsmenge einer Wertlösung, und sein Rückwärtslauf
ist eine Neurechnung je Schritt unter `codual` mit einem Zustand von `n_states`. Die Kosten sind
unabhängig von `n_theta`.

**Erster Ordnung ist der Gewinn damit ein Faktor und keine Größenordnung.** Gegen `factor_grad = 8.5`
auf Bachmann steht ein Adjoint bei geschätzt 2 bis 4 Wertlösungen, also Faktor 2 bis 4. Das trägt das
Tor aus dem Abschnitt unten nicht allein: `sr1` aus dem Identitäts-Seed braucht 5.8-mal so viele
Auswertungen wie `gn`, und 5.8 geteilt durch 3 ist noch immer größer als 1. Der erwartete Wert liegt
stattdessen bei größeren Modellen, wo die Kurve noch nicht gemessen ist, und beim Speicher.

**Zweiter Ordnung ändert sich die Ordnung, und dort liegt der Fall.** `deriv2` propagiert
`n_theta^2` Richtungen: auf Bachmann etwa 2.5 MB, auf Lang (124 Zustände, 294 Parameter) etwa 86 MB
im integrierten Zustand, also unbenutzbar. Forward-over-reverse ersetzt das durch `O(n_theta)`. Der
Adjoint erster Ordnung ist das Fundament dafür und wird deshalb gebaut, auch wenn er für sich
genommen nur einen Faktor bringt. Dazu kommt die Festlegung im Abschnitt "Vollständigkeit", die
ohnehin nicht aus einer Kostenrechnung folgt.

**Eine Verwechslung, die der Plan an jeder Stelle ausschließen muss.** Reverse gewinnt nur beim
Skalar. Die volle Ableitungsmatrix rückwärts zu füllen kostet einen Durchlauf je Ausgabe, auf Lang
9600 gegen 294 vorwärts, also Faktor 32 in die falsche Richtung (Lang, gemessen 2026-09-04). Was
den Adjoint billig macht, ist die Akkumulation auf einen Vektor an jeder Beobachtungszeit. Ein
`seed` ist für vorwärts optional und für rückwärts zwingend, und er kommt immer aus einer skalaren
Reduktion.

Ein diskreter Adjoint liefert den Gradienten unabhängig von `n_theta` und macht über
forward-over-reverse die exakte Hesse in `O(n_theta)` statt `O(n_theta^2)` erreichbar. Der
Konsument steht seit `e1a88bc` und `45261ba` auf master: `trust()` nimmt seinen Modell-Hesse aus
einer austauschbaren Quelle (`hessianMethod`, dazu `hessianFallback` und `fallbackLimit`), und der
quasi-Newton-Zweig braucht nur Gradienten. `hessian` ist dabei ein Aufrufargument, das von
`trust_driver.h:181` über `objClass.R:197` und `.evalProd` (`classes.R:317`) bis in die Kerne
durchgereicht wird. Der gradientenreine Pfad existiert also bereits durchgehend und wird benutzt.

Dieser Plan baut den Produzenten: reverse-mode Ableitungen durch die vollständige Kette, vom
Solver einschließlich seiner Schrittweiten- und Ordnungssteuerung über die algebraischen
Schichten bis zum Skalar von `normL2`.

## Festlegungen

- Der Adjoint geht durch die **vollständige Schrittweiten- und Ordnungssteuerung**.
- **Am Ende trägt alles reverse.** Kein Verfahren, kein Schalter und kein Modellmerkmal bleibt
  draußen. Was das heißt, steht im Abschnitt "Vollständigkeit" und ist das Abnahmekriterium des
  Plans.
- Ein gemeinsames Gerüst über **alle Stepper** von Anfang an, auch wenn zuerst nur einer verdrahtet
  wird.
- Der erste Meilenstein ist ein **`obj()`-Aufruf im Reverse-Modus auf einem kleinen Modell**, gegen
  den Vorwärtsaufruf geprüft. Siehe den eigenen Abschnitt weiter unten.
- **Kein eingefrorener Modus** im Auslieferungsstand. Ein Schalter dafür existiert nur im
  Reverse-Pfad als Testartefakt, siehe unten.
- Der Weg ist ein **Reverse-Skalartyp**, nicht handgeschriebene Adjoint-Gleichungen pro Stepper.
  Was dafür zusätzlich templatisiert werden muss, wird templatisiert.
- cppDE bekommt einen Branch `devel-reverseAD`; dMod2 hat ihn bereits. Beide stehen am 2026-09-07
  auf ihrem master, die Optimierer-Arbeit ist dort gelandet. `sync-devel.yaml` merged in dMod2 jeden
  master-Push nach `devel-reverseAD` und legt bei Konflikt ein Issue an, statt den Branch zu
  überschreiben; die dMod2-Schicht aus Stufe 7 bleibt deshalb schmal und merge-freundlich.

## Wer den Adjoint aufruft, und was er heute kostet

Die Messungen aus `dMod2/dev/optimiser-plan.md`, Stand 2026-09-07, beantworten die Frage, welcher
Optimierer den Adjoint überhaupt aufruft. Die Antwort ist nicht die, von der dieser Plan
ausgegangen ist.

**Der Fallback-Zweig ist kein Aufrufer.** Auf Boehm sind in `gn -> sr1` 14 von 142 Auswertungen
quasi-Newton, der Arm gewinnt keinen Start, den `gn` nicht auch erreicht, und ist auf 21 von 100
Starts bitgleich mit `gn`. Auf Bachmann kostet er 308 Auswertungen je Plateau-Start gegen 276 für
`gn`. Ein Adjoint, der nur diesen Zweig bedient, lohnt sich nicht.

**Der Aufrufer ist `sr1` aus einem Identitäts-Seed über den ganzen Abstieg.** Auf Bachmann mit 113
Parametern ist das der einzige Arm, der das Plateau verlässt: bestes Ergebnis -393.96 gegen -389.79
für beide Gauss-Newton-Arme, zwei Starts unter -390, kein Parameter an einer Schranke. Er ist zu
100 Prozent gradientenrein und zahlt heute 1590 Auswertungen je Plateau-Start gegen 276.

**Damit hat Stufe 7 ein Tor, das vor der ersten C++-Zeile ablesbar ist.** Der Arm braucht auf
Bachmann den Faktor 5.8 an Auswertungen. `inst/benchmarks/bench_gradientCost.R` misst `factor_grad`,
den Preis eines Objective-Gradienten in reinen Wertlösungen, über `n_theta` bis 113. Liegt
`factor_grad` dort über 5.8, macht der Adjoint den Arm, der das bessere Optimum findet, in Wandzeit
billiger als `gn`, nicht nur konkurrenzfähig. Diese Zahl gehört gemessen, bevor Stufe 1 beginnt; das
Skript steht bereits und braucht nur eine ruhige Maschine.

**Ein zweiter Aufrufer ist bereits benannt.** Sampled quasi-Newton (Berahas, Jahani, Takac) ist im
Optimierer-Plan ausdrücklich hierher verschoben: zusätzliche Gradienten entlang zufälliger
Richtungen sind unter Vorwärtssensitivitäten von `gn` dominiert und werden erst mit einem
konstanten Gradientenpreis sinnvoll. Er gehört in die Nachfolge dieses Plans, nicht in ihn.

## Tragende Befunde aus der Erkundung

**Die Stepper sind bereits über den Skalartyp templatisiert.**
`multistepper<Method, Value, JacobianPattern, Resizer>` (`cppde_multistepper.hpp:525-531`),
`rosenbrock4<Value, ...>` (`cppde_rosenbrock4.hpp:124-130`), `tsit5<Value, Resizer>`
(`cppde_tsit5.hpp:42-46`), beide Controller ebenso. Ein zusätzlicher Skalartyp ist damit die
vorgesehene Erweiterungsachse, nicht ein Fremdkörper.

**Die Zeit ist die Ausnahme.** `time_type = ad_lu::scalar_type_t<Value>`
(`cppde_multistepper.hpp:545`, `cppde_rosenbrock4.hpp:137`, `cppde_tsit5.hpp:53`) ist immer
`double`, und `do_step` skalarisiert `t` und `dt` in der ersten Zeile
(`cppde_multistepper.hpp:772-773`, `:1037-1038`; `cppde_rosenbrock4.hpp:244-245`;
`cppde_tsit5.hpp:187`). Damit sind `m_h`, `m_hscale`, `m_eta`, `m_tau`, `m_l`, `m_gamma`,
`m_tq[]` tangentenfrei.

**Die Fehlerkontrolle sieht die Tangenten, das Ergebnis ist trotzdem reell.** `m_ewt` trägt einen
Eintrag je Tangente (`cppde_multistepper.hpp:812-829`), `wrms_max_ewt`
(`cppde_newton.hpp:158-218`) nimmt das Maximum über Zustand und jede Sensitivitätsscheibe. `dsm`
ist ein `double`. Die Schrittfolge hängt also von den Sensitivitäten ab, `h` trägt aber keine
Ableitung, und der Term `dy/dh_k · dh_k/dtheta` fehlt heute überall.

**Der Newton-Korrektor ist bereits per IFT differenziert, nicht per Replay.** `cppde_ad_lu.hpp`
macht "AD-aware implicit-function-theorem peeling" im `solve`. Das Muster, das der Adjoint
braucht, existiert also schon in Vorwärtsrichtung und muss nur gespiegelt werden.

**Ausgabe erfolgt interpoliert, nie an Schrittenden.** `process_dense`
(`cppde_event_engine.hpp:431-597`) ruft `calc_state` an den angeforderten Zeiten; der
Interpolant ist das Nordsieck-Horner-Polynom (`cppde_multistepper.hpp:1791-1874`). Der
Adjoint-Seed sitzt damit an interpolierten Punkten, und die Interpolation gehört mit adjungiert.
Sie ist linear in `zn`, also unproblematisch.

**Forcings tragen keine Parametersensitivität** (`cppde_pchip_forcing.hpp`, Knoten und
Koeffizienten sind `double`). Ein Block weniger.

**Die Naht kennt den Seed-Begriff bereits.** `sens1ini` ist Phi'(theta), `[n_phi_rows, n_theta]`,
und ihre Spaltenzahl setzt die Sensitivitätsbreite zur Laufzeit (`cppDE/R/cppODE.R:443-445`,
gespiegelt in `cppde_r_batch.hpp:592`), gesetzt in `dMod2/R/prediction.R:245-271`. Ein Argument
an Position 15, vor `dimnames`, ist mit den vorhandenen positionellen Zugriffen verträglich, und
der Batch-Pfad ist über `cond_elt` (`cppde_r_batch.hpp:253`) in beide Richtungen
abwärtskompatibel.

**Es gibt heute nichts Reverses.** Weder `adjoint` noch `vjp` noch `costate` kommt in cppDE vor.
`cppde::chain_jac` ist `J · S`, also Rechtsmultiplikation. `derivSymb.py` liefert nur
algebraische Partialableitungen.

**Namenskonflikt.** `derivMode` ist in `cppDE/R/funCpp.R:52-55` bereits vergeben und bedeutet
`"dual"` gegen `"symbolic"`; dMod2 reicht es durch. Der neue Modus braucht einen anderen Namen,
Vorschlag `sweep = c("forward", "reverse")`.

**Der Modus ist keine Laufzeitfahne, sondern ein viertes kompiliertes Objekt.** `solveODE()`
(`cppDE/R/solveODE.R:583-590`) hat überhaupt kein Ableitungsargument; es liest
`attr(model, "deriv")` und `attr(model, "deriv2")` (`:20-21`) und verweigert `sens1ini` an einem
Modell, das nicht dafür übersetzt wurde (`:34-37`). `funCpp`s `derivMode` friert ebenso beim Bauen
ein (`funCpp.R:60`, `:113-119`, gestempelt `:139`). Das etablierte Muster steht bereits dreifach in
`odemodel()`: `func`, `extended`, `extended2` (`odeClass.R:206`, `:214`, `:220`), ausgewählt je
Aufruf über `pickModel(deriv, deriv2)` (`prediction.R:298`). Reverse ist damit ein viertes Objekt
plus ein Zweig in `pickModel`, kein Durchreichen einer Fahne durch die ganze Kette. Das verkleinert
die Stufen 6 und 7 erheblich.

**dMod2 wählt den Stepper nie.** `odemodel()` (`odeClass.R:102`) hat kein `method`-Argument, die
drei `cppODE()`-Aufrufe setzen keines, und `Xs.cppDE`s `optionsOde` warnt bei unbekannten Optionen
(`prediction.R:206-211`), kann es also nicht nachreichen. Der Produktionspfad ist damit ausnahmslos
`method = "bdf"`, der Nordsieck-Multistepper. Ein anderer Stepper ist nur über die Punkte von
`odemodel()` erreichbar, die `pick(cppDE::cppODE, dots)` (`odeClass.R:204-205`) durchreicht, und
`method` ist ein Formal von `cppODE`. In `dMod2/R` und `dMod2/tests` kommt `tsit5` oder `rb4` kein
einziges Mal vor.

**`codual` ist kein weiterer `is_ad`-Typ.** `cppde_ad_traits.hpp` trägt neben `is_ad`,
`inner_type` und `scalar_value` auch `extract_derivs`, `max_deriv_size`, `any_deriv` und
`bulk_inject_results`, die alle Tangentenplätze voraussetzen, und `is_ad` zieht einen Typ in die
Sensitivitätsschleife von `wrms_max_ewt` und in das Peeling von `cppde_ad_lu.hpp`. Ein `codual`
trägt keinen Tangentenplatz, sondern einen Tape-Index. Also ein eigenes Merkmal `is_reverse`,
`is_ad` bleibt falsch, und die wenigen Stellen, die verzweigen müssen, verzweigen darauf.

## Architektur

### Der Skalartyp `cppde::codual<T>`

Neben `cppde::dual<T,N>` tritt `cppde::codual<T>`. Ein `codual` hält einen Wert und einen Index
in ein Tape; jede Operation schreibt einen Eintrag mit den lokalen Partialableitungen und den
Indizes ihrer Operanden. Der Rückwärtslauf über das Tape akkumuliert Kotangenten.

Dateien nach vorhandener Konvention: `cppde_codual.hpp` (Typ), `cppde_codual_math.hpp`
(Operatoren und Mathematik), `cppde_codual_tape.hpp` (Band und Rückwärtslauf),
`cppde_codual_traits.hpp` oder Erweiterung von `cppde_ad_traits.hpp`.

**Expression Templates gelten für `codual` nicht.** `cppde_dual_expr.hpp` existiert, um
Temporaries zu kollabieren; ein Tape braucht genau diese Zwischenwerte. Das ist kein Konflikt,
solange `codual` keine ET-Schicht bekommt: die vorhandenen ETs sind auf `dual` spezialisiert,
`codual` geht den eager Pfad. Der Verlust ist hinnehmbar, weil der Rückwärtslauf ohnehin von der
Bandbreite des Tapes dominiert wird.

### Tape pro Schritt, Checkpoint pro Schritt

Ein Tape über die ganze Trajektorie ist nicht tragbar und auch nicht nötig. Stattdessen:

- **Vorwärts** läuft der Solver unverändert in `double` oder `dual`. Nach jedem akzeptierten
  Schritt wird der vollständige Solver-Zustand als Checkpoint abgelegt.
- **Rückwärts** wird jeder Schritt aus seinem Checkpoint mit `codual` neu gerechnet, wobei das
  Tape entsteht, sofort rückwärts abgearbeitet und wieder verworfen wird.

Kosten: etwa eine zusätzliche Vorwärtsrechnung plus Tape-Overhead je Schritt. Speicher: ein
Schritt-Tape, das in einer Arena leben kann, die pro Schritt zurückgesetzt wird, plus die
Checkpoints.

**Checkpointgröße.** Für den reinen Gradienten trägt der Zustand keine Forward-Tangenten mehr,
also ist das Nordsieck-Array `(max_order+2) * n_x` Doubles, 7 Slots bei BDF, 14 bei Adams, dazu
die Skalare (`m_h`, `m_hscale`, `m_eta`, `m_q`, `m_qwait`, `m_tau`, `m_l`, `m_tq`, `m_crate`,
`m_gamma`, `m_gammap`, `m_nst`, `m_nstlp`, `m_nstlj`, `m_ewt`). Auf Lang sind das etwa 7 KB je
Schritt, bei tausend Schritten etwa 7 MB. Faktorisierungen fallen wegen `MSBP = 20` und
`MSBJ = 51` nur selten an und werden separat zwischengespeichert.

Das ist der Grund, warum kein Revolve nötig ist: alles speichern ist bezahlbar.

### Der Newton-Korrektor wird nicht getapet

`ndf_newton_solve` (`cppde_newton.hpp:231-384`) konvergiert auf `F(y) = 0`. Sein Adjoint ist eine
Lösung mit der transponierten, bereits faktorisierten Matrix, kein Replay der Iterationen. Der
`codual`-Pfad muss dort dieselbe IFT-Abkürzung nehmen, die `cppde_ad_lu.hpp` in Vorwärtsrichtung
schon nimmt. Das ist die eine Stelle, an der von Hand adjungiert wird, und sie ist klein.

### Semantik an den nicht differenzierbaren Punkten

Die Steuerung entscheidet an sehr vielen Stellen über Gleitkommavergleiche, darunter der
Akzeptanzschalter `dsm <= 1.0` (`cppde_multistepper_controller.hpp:265`), die Ordnungsauswahl
über exakte Gleichheit `etam == etaq` und `etam == etaqm1` (`:425`, `:429`), die
`callSetup`-Disjunktion (`cppde_multistepper.hpp:858-862`) und der Newton-Abbruch `dcon <= 1.0`
(`cppde_newton.hpp:356`).

Festlegung: **das Tape differenziert den tatsächlich genommenen Zweig.** Kontrollentscheidungen,
also Akzeptanz, Ordnung, Zahl der Schritte, Zahl der Newton-Iterationen und
Jacobi-Neuberechnungen, sind stückweise konstant und werden nicht differenziert. Alles, was
innerhalb eines festen Kontrollpfads glatt ist, wird differenziert, einschließlich der Werte von
`h`, `eta`, `gamma` und der Nordsieck-Koeffizienten.

Das ist die Standardsemantik von Operator-Overloading-AD, sie ist fast überall korrekt, und sie
ist genau die Ableitung derjenigen stückweise glatten Funktion, die der Optimierer tatsächlich
sieht. Sie gehört in die Roxygen und in `Methods.Rmd`.

### Was "exakt" hier heißt

Zwei Genauigkeitsbegriffe sind zu trennen, und der Unterschied ist der eigentliche Ertrag des
diskreten Wegs.

**Zur berechneten Funktion.** Der Solver liefert `y_h(theta)`, die diskretisierte Lösung. Der
diskrete Adjoint gibt deren Ableitung auf Maschinengenauigkeit, nicht auf Solvertoleranz. Wert und
Gradient gehören damit exakt zusammen, was Liniensuchen und Profile trägt.

**Zur wahren Lösung.** Davon weicht er um `O(tol)` ab, aber das ist der Fehler von `y_h` selbst und
nicht der der Ableitung. Er ist unvermeidbar, egal wie differenziert wird.

| | zur berechneten Funktion | zur wahren Lösung |
|---|---|---|
| voller diskreter Adjoint | Maschinengenauigkeit | `O(tol)` |
| Forward-Sensitivitäten heute | `O(tol)`, es fehlt der Schrittweitenterm | `O(tol)` |
| CVODES ASA | `O(tol)` | `O(tol)` |

Daraus folgt, warum das Orakel scharf ist: der eingefrorene Adjoint und die Forward-Sensitivitäten
berechnen dieselbe Größe, einmal vorwärts und einmal rückwärts. Zwischen ihnen liegt kein
Toleranzfehler, nur Rundung. Finite Differenzen hätten diese Eigenschaft nicht.

Zwei Einschränkungen, damit "Maschinengenauigkeit" nicht mehr verspricht als es hält: Rundung
akkumuliert wie in jeder Rechnung, und der Newton-Korrektor wird per Implicit-Function-Theorem
adjungiert, liefert also die Ableitung der exakt gelösten Korrektorgleichung und nicht die der bei
`dcon <= 1` abgebrochenen Iteration. Der Vorwärtspfad nimmt in `cppde_ad_lu.hpp` dieselbe
Abkürzung, weshalb beide übereinstimmen und das Orakel gültig bleibt.

### Der Verifikationsschalter sitzt im Reverse-Pfad, nicht im Forward

Der Vorwärtsmodus bleibt unangetastet. `h` dual zu machen hieße, Fehlernorm, Kontrollgesetz und
Stage-Alphas dual zu rechnen, denn `onestep_controller::error()` und
`multistepper::error_norm()` liefern beide `double` und `dt *= factor` kann `dt` daher keine
Tangente geben, unabhängig von `time_type`. Das würde den Vorwärtspfad dauerhaft verteuern, für
einen Nutzen, den nur die Verifikation hat.

Stattdessen entscheidet der **Reverse-Pfad**, ob die Schrittweiten- und Ordnungskette mit ins Tape
geht:

- **ohne** ist der Adjoint exakt gegen die heutigen Forward-Sensitivitäten prüfbar, bis auf
  Rundung;
- **mit** unterscheidet er sich davon um genau den Schrittweitenterm.

Das ist ein Testartefakt, kein zweiter ausgelieferter Modus, und es kostet im Normalbetrieb
nichts. Der Auslieferungsstand ist weiterhin der volle Adjoint.

Damit verteilt sich die Verifikationslast: die Primitiven gegen `dual` (Stufe 1), der eingefrorene
Schritt und die eingefrorene Trajektorie gegen die Forward-Sensitivitäten (Stufen 3 und 4), und
der volle Term als dieselbe geprüfte Maschinerie auf mehr Operationen, plausibilisiert über den
`rtol`-Sweep aus Stufe 0.

Sollte sich ein exaktes Orakel für den vollen Term später doch als nötig erweisen, ist `h` dual
hinter einem Compile-Makro nach dem Muster von `CPPDE_STEP_TRACE` die Form, die nichts kostet.

## Vollständigkeit

Das Abnahmekriterium des Plans: **jede Kombination, die vorwärts eine Ableitung liefert, liefert am
Ende auch rückwärts eine.** Der Reverse-Modus ist dann kein Sonderfall mit Kleingedrucktem, sondern
eine zweite Richtung durch dieselbe Fläche. Was das aufzählt:

**Verfahren.** `method = c("bdf", "adams", "rb4", "tsit5")` (`cppODE.R:77`), vier Werte über drei
Steppertypen, weil `bdf` und `adams` dasselbe `multistepper`-Template mit anderem
`multistep_method`-Enum sind. Dazu `useNDF` in beiden Multistep-Varianten.

**Lineare Algebra.** `dense_lu_tag` und `sparse_lu_tag` (`cppODE.R:784-799`), also auch der
KLU-Pfad, wo er zur Verfügung steht.

**Modellmerkmale.** `events`, `rootfunc`, `forcings`, `fixed`, `includeTimeZero`, `nStack`. Für
`useDenseOutput = FALSE` mit Wurzelereignissen bleibt es bei einem Fehler statt einer stillen
Falschrechnung, siehe Stufe 5; das ist die eine benannte Ausnahme und sie steht in der Roxygen.

**Einstiegspunkte.** `solveODE()` und `solveODEBatch()`, also auch der Bedingungs-Batch über OpenMP,
und die algebraische Schicht `funCpp()` in beiden `derivMode`-Varianten.

**Backend.** `cvode()` bekommt seinen Reverse-Modus über die ASA von CVODES, Stufe 8. Damit trägt
auch der zweite Backend reverse, und zwar mit der Mathematik, die dort hingehört.

**Ordnung.** Erste Ordnung überall, danach zweite über forward-over-reverse. Das ist die Nachfolge
dieses Plans und steht unter "Was nicht", nicht weil es entfällt, sondern weil es den verifizierten
Adjoint erster Ordnung als Fundament braucht.

Was am Ende nicht existiert, existiert auch vorwärts nicht.

## Der Meilenstein

Die Langstrecke steht, und der Endzustand ist der Abschnitt darüber. Die erste Frage ist trotzdem
eine kleinere: **trägt reverse überhaupt, und können wir damit `obj()` aufrufen?** Das ist der
Meilenstein, und verifiziert wird er mit Vorwärtsaufrufen auf kleinen Modellen.

**Er sitzt auf `tsit5`.** Das ist eine Entscheidung über die Reihenfolge, nicht über den Umfang: der
Produktionspfad ist ausnahmslos `bdf` (siehe Befunde), und der Multistepper kommt unmittelbar danach
als Stufe 3b. Die Begründung ist Fehlerlokalisierung. Der Meilenstein verbindet zum ersten Mal sechs
neue Teile: `codual`, das Tape, die Checkpoint-Struktur, den Rückwärtslauf über die Trajektorie, die
Naht und die dMod2-Kette. Steht er auf `tsit5`, sind das 415 Zeilen expliziter Runge-Kutta ohne
Newton, ohne Ordnungswahl und ohne Nordsieck-Historie, und ein falsches Ergebnis liegt sicher in
einem dieser sechs Teile. Steht er auf dem Multistepper, kommen 2381 Zeilen mit Newton-Korrektor,
variabler Ordnung, Historie über Schrittgrenzen und dem `qwait`-Fenster dazu, und der Fehler kann
überall liegen. Der Adjoint eines Mehrschrittverfahrens mit variabler Ordnung ist deutlich
schwerer korrekt herzuleiten als der eines Einschrittverfahrens, weil die Abhängigkeit eines
Schritts ein Fenster über mehrere frühere Schritte plus die Ordnungs- und Schrittweitenwechsel
umfasst.

Die Stufen 4, 6 und 7, also Trajektorie, Naht und dMod2-Kette, sind stepperunabhängig. Der Umweg
kostet deshalb nur die Checkpoint-Struktur von `tsit5`, und die wird ohnehin gebraucht.

**Was dazugehört.** Stufen 1, 2, 3 (nur `tsit5`), 4, 6 und 7 auf den kleinen Fixtures.

**Stufe 2 ist dabei, und zwar wegen `Pexpl`, nicht wegen `Y`.** `Y` ist vermeidbar,
`normL2(data, Xs(m) * Pexpl(...))` mit zustandsbenannten Daten braucht keine Beobachtungsfunktion
(`test-Xs.R:20`, `inst/examples/normL2.R`). `Pexpl` liegt auf jedem Pfad und ist per Vorgabe ein
kompilierter AD-Aufruf.

**Was nicht.** Stufe 5, Events und Wurzeln: kleine Testmodelle haben keine. Stufe 8 ohnehin nicht.

**Das Fixture.** `odemodel(reactions, method = "tsit5", ...)` erreicht `tsit5` über die Punkte, die
`pick(cppDE::cppODE, dots)` (`odeClass.R:204-205`) durchreicht; `Xs.cppDE`s `optionsOde` kann es
nicht, es warnt bei unbekannten Optionen. Ansonsten `fx_decay_compiled()`
(`helper-fixtures.R:54`) unverändert.

**Warum die Verifikation scharf ist.** Der eingefrorene Reverse-Pfad und die heutigen
Vorwärtssensitivitäten berechnen dieselbe Größe aus entgegengesetzter Richtung. Zwischen ihnen liegt
kein Toleranzfehler, nur Rundung, siehe "Was exakt hier heißt". Der Meilenstein ist damit auf
Rundungsniveau geprüft und nicht auf Solvertoleranz, und das gilt für jede seiner Stufen einzeln.
Das ist eine schärfere Aussage als jeder heutige dMod2-Gradiententest: `test-normL2.R:449` und `:481`
prüfen gegen `numDeriv::grad` bei Toleranz 1e-3.

## Stufen

Jede Stufe endet auf einem Test, der ohne die folgende Stufe läuft.

### Stufe 0. Branch und die Größe des Schrittweitenterms messen

```sh
cd ~/Documents/Projects/dModverse/cppDE && git checkout -b devel-reverseAD   # erledigt
```

Kein C++. Die Forward-Sensitivitäten enthalten den Schrittweitenterm nicht, eine zentrale
Differenz mit laufender adaptiver Steuerung enthält ihn. Die Differenz der beiden, aufgetragen
über `rtol`, trennt ihn vom FD-Rauschen: der Term skaliert mit `rtol`, das Rauschen folgt `eps`
und der Kondition.

Skript `dev/stepsize-term.R`, linear zum Durchklicken. Modelle mit arbeitender
Schrittweitensteuerung, also Robertson und ein Modell aus `tests/benchmarks/bench_stiff_suite.R`,
über alle vier Methoden, `rtol` von 1e-4 bis 1e-11 bei festem `eps`.

Diese Zahl entscheidet, wie viel Aufwand die Schrittweitenkette in den Stufen 3 und 4 rechtfertigt,
und sie ist die Plausibilisierung, die dort an die Stelle eines exakten Orakels tritt.

**Ergebnis, 2026-09-06.** Das Skript läuft, die Methode ist nicht die geplante finite Differenz,
sondern der Diskretisierungsfehler `e_h(theta)` gegen einen Lauf bei `rtol = 1e-13`; FD kann den Term
grundsätzlich nicht auflösen. Die Kennzahl sieht nach Sprung-Artefakt aus, nicht nach glattem
Beitrag: auf Robertson bei `rtol = 1e-4` ist sie das Zweihundertfache des Gradienten, und auf Decay
bei `rtol = 1e-6`, wo die Schrittzahl über den ganzen Sweep konstant bei 99 bleibt, fällt sie auf
`6e-6`. Der Kopf von `dev/stepsize-term.R` beschreibt die Größe und die Lesart.

**Gegenprobe, 2026-09-07: bestätigt.** Dieselbe Kennzahl bei `N_SWEEP` 200, 400 und 800. Sie
verdoppelt sich mit der Sweep-Dichte, statt stehen zu bleiben:

| Modell | rtol | N=200 | N=400 | N=800 | x400/200 | x800/400 |
|---|---|---|---|---|---|---|
| Decay | 1e-4 | 1.38 | 2.77 | 5.54 | 2.01 | 2.00 |
| Decay | 1e-6 | 4.0e-6 | 6.0e-6 | 1.0e-5 | 1.50 | 1.76 |
| Decay | 1e-8 | 7.4e-4 | 1.5e-3 | 3.0e-3 | 2.01 | 2.00 |
| Decay | 1e-10 | 7.9e-4 | 1.5e-3 | 3.2e-3 | 1.90 | 2.11 |
| Robertson | 1e-4 | 292 | 791 | 1431 | 2.71 | 1.81 |
| Robertson | 1e-6 | 36.7 | 73.1 | 148 | 1.99 | 2.02 |
| Robertson | 1e-8 | 0.31 | 1.24 | 3.13 | 3.94 | 2.53 |
| Robertson | 1e-10 | 0.048 | 0.087 | 0.386 | 1.83 | 4.41 |

Kein Verhältnis liegt bei 1, acht von zwölf liegen zwischen 1.8 und 2.1, und die Ausreißer liegen
darüber, weil ein dichterer Sweep Sprünge trifft, die ein gröberer übersprungen hat. Die Kennzahl ist
Sprunghöhe geteilt durch `dtheta` und damit **keine Ableitung**. Die Decay-Zeile bei `rtol = 1e-6`
bestätigt es von der anderen Seite: dort ist die Schrittzahl über den ganzen Sweep konstant bei 99,
es wird nirgends umgeschaltet, und die Kennzahl liegt bei 1e-5, ist also gar nichts.

**Was daraus folgt.** `e_h(theta)` ist stückweise glatt mit Sprüngen an den Umschaltpunkten der
Schrittfolge, und innerhalb eines Stücks ist seine theta-Ableitung vernachlässigbar. Das ist genau,
was Fehlerkontrolle leisten soll: bei fester Ordnung hängt `y_h` nur über den lokalen Fehler von `h`
ab, also ist `dy/dh * dh/dtheta` von dessen Größe. Über einen Sprung hinweg ist `y_h` in theta
unstetig, und das fängt keine Ableitung ein; der Optimierer sieht dort Rauschen der Größe `O(tol)`.

Drei Konsequenzen, alle an einer Zahl statt an einem Argument:

1. `h` dual bleibt draußen, dauerhaft. Nicht weil es zu teuer wäre, sondern weil es nichts misst.
2. Der Schwerpunkt liegt vollständig auf dem Adjoint selbst.
3. **Das Orakel wird schärfer, als der Plan angenommen hat.** Eingefrorener und voller Reverse-Pfad
   müssen nicht um den Schrittweitenterm auseinanderliegen, sondern auf Höhe des lokalen Fehlers
   übereinstimmen. Der Test in Stufe 4 ist damit eine Gleichheit und keine Größenordnungsprüfung.

### Stufe 1. `codual` und das Tape

`cppde_codual.hpp`, `cppde_codual_math.hpp`, `cppde_codual_tape.hpp`. Vollständiger
Operatorensatz spiegelbildlich zu `cppde_dual_math.hpp`, Arena-gestütztes Band, Rückwärtslauf,
Erweiterung von `cppde_ad_traits.hpp` um `is_reverse`, `inner_type` und `scalar_value` für den neuen
Typ. `is_ad` bleibt für `codual` falsch, siehe die Befunde oben.

**Stand 2026-09-07: steht.** Der Test ist `dev/cxx/test_codual.cpp`, gefahren über
`dev/cxx/run.sh --codual`, nicht die geplante testthat-Datei: `codual` ist vor Stufe 2 von R aus
nicht erreichbar, weil kein Emitter ihn ausgibt. Das Orakel ist trotzdem dasselbe, `dual` gegen
`codual` auf derselben Formel, und der C++-Weg ist schärfer, weil er auf 1e-14 vergleicht statt auf
die 1e-10 der R-Suite.

Abgedeckt: die vier Operatoren mit ihren gemischten Skalar-Überladungen, fünfzehn unäre Funktionen,
`abs` beidseitig der Null, `pow` in allen drei Formen, `min`, `max`, `clamp`, dazu geteilte
Teilausdrücke und Verkettungen. Dazu drei strukturelle Prüfungen: ein Ausdruck aus lauter Konstanten
schreibt null Knoten, wiederholtes Seeden akkumuliert, und `scope` stellt den Knotenstand wieder her.

Die testthat-Fassung nach dem Muster von `tests/testthat/test-dual2nd-primitives.R` kommt mit
Stufe 2 dazu, wenn `funCpp` den Typ emittiert; `run_modes` und `expect_modes_agree` (`:21`, `:38`)
sind dann die Vorlage.

**Offen auf dieser Maschine:** `libasan` und `libubsan` sind nicht installiert, deshalb bricht
`dev/cxx/run.sh` im Sanitizer-Schritt ab, für den bestehenden Test genauso wie für den neuen. Der
`-O2`-Build läuft mit `-Wall -Wextra` warnungsfrei.

### Stufe 2. `funCpp` bekommt ein vjp

Die algebraische Schicht zuerst, weil dort kein Zeitintegrationsproblem sitzt und weil sie die
Konventionen für alles Weitere festlegt. Der Name ist `funCpp`, nicht `cppFUN`; die von Simon
erwogene Umbenennung ist offen und gehört nicht in diesen Plan.

`codegen_funCpp.py` emittiert `<model>_vjp` und `<model>_vjp_c`, die `w' J` bilden. Der bestehende
`_eval_ad`-Pfad bleibt unverändert. `R/funCpp.R` bekommt den passenden Einstieg als weiteres
kompiliertes Objekt, nach dem Muster oben.

**Der Kunde ist `Pexpl`, nicht `Y`.** Beide bauen ein `funCpp` (`parameters.R:255`,
`prediction.R:805`), aber `Y` ist vermeidbar: `normL2(data, Xs(m) * Pexpl(...))` mit
zustandsbenannten Daten läuft ohne Beobachtungsfunktion, wie `test-Xs.R:20` und
`inst/examples/normL2.R` zeigen. `Pexpl` liegt dagegen auf jedem Pfad, und sein `derivMode` ist
per Vorgabe `"dual"`, also ein kompilierter AD-Aufruf. Ausweichen ginge nur über eine
Identitätstrafo oder `derivMode = "symbolic"` (`parameters.R:137-143`), und beides ist als
Dauerlösung nichts wert.

**Stand 2026-09-07: steht.** `codegen_funCpp.py` emittiert `<model>_vjp` und `<model>_vjp_c`,
`R/funCpp.R` gibt den Einstieg als `f$vjp(vars, params, w)` heraus, und drei Tests in
`tests/testthat/test-funCpp.R` prüfen ihn gegen den symbolischen Jacobi-Pfad bei 1e-12.

Signatur und Formen. `w` ist `[n_obs, n_out]` oder `[n_obs, n_out, n_seed]`, zurück kommen `y`,
`wx` `[n_obs, n_vars, n_seed]` und `wp` `[n_params, n_seed]`. `wp` summiert über die
Beobachtungen, weil die Parameter über sie geteilt sind, `wx` nicht.

Zwei Eigenschaften, die die Kette später braucht und die geprüft sind:

- **Eine Aufzeichnung trägt beliebig viele Seeds.** Die Partialableitungen hängen nicht vom Seed ab,
  also wird je Beobachtung einmal aufgezeichnet und danach nur noch gesweept. Ein Identitäts-Seed
  über die Ausgänge liefert damit die volle Jacobi aus einem Durchlauf.
- **Die Form ohne Variablen funktioniert.** Das ist `Pexpl`s Fall, `variables = NULL` und eine
  Beobachtung.

**Fallstrick für Stufe 7.** `funCpp(compile = FALSE)` ist die Vorgabe, und der Reverse-Pfad hat
keinen interpretierten Rückfall, anders als `func` und der symbolische `jac`. Ohne Kompilat gibt es
einen Fehler, der genau das sagt. `derivMode = "symbolic"` liefert kein `vjp`.

**Die eigentliche Lehre aus dieser Stufe, und sie gilt für Stufe 3 genauso.** Das vjp instanziiert
den Modellrumpf mit `codual`, also muss die codual-Oberfläche **alles** abdecken, was ein Generator
emittieren kann, nicht nur die Rechenoperationen. `test-piecewise.R` ist genau daran gescheitert:
`cppde::select`, die Form, in die der piecewise-Codegen den Bedingungsoperator übersetzt, gab es für
`codual` nicht, und dazu `cppde::value_of`. Beides ist nachgetragen, `select` gibt den gewählten
Operanden ganz zurück und schreibt nichts aufs Band, weil der Zweig eine Kontrollentscheidung ist.
Vor Stufe 3 gehört dieselbe Prüfung über die Stepper-Header: jeder Aufruf, den ein Schritt auf dem
Zustandstyp macht, braucht seine codual-Form.

Die testthat-Fassung nach `test-Pexpl.R:78-94` und `test-Y.R:69-78` bleibt die Vorlage, wenn dMod2
in Stufe 7 seine eigene Zweiwege-Parität bekommt.

**Warum nicht der analytische Weg.** `derivSymb.py` könnte `(df/dp)^T lambda` direkt emittieren, und
das ist das stärkste Argument dafür, den Adjoint diskret statt über CVODES ASA zu bauen: kein Tape,
keine Reverse-Infrastruktur, und der Zustandsanteil braucht nur `J^T lambda`, wobei `J` für Newton
ohnehin assembliert wird. Für die rechte Seite stimmt das auch. Es trägt nur nicht weiter: der Rückwärtslauf muss auch die Arithmetik des Steppers
selbst differenzieren, die Nordsieck-Kombination, die Fehlernorm, das Kontrollgesetz und die
Interpolation, und die emittiert kein Codegen. Deshalb `codual`.

### Stufe 3. Ein Schritt rückwärts

Checkpoint-Struktur je Steppertyp, Neurechnung eines Schritts unter `codual`, Rückwärtslauf,
Kotangens des Schrittanfangs aus dem Kotangens des Schrittendes. Der Newton-Korrektor nimmt die
IFT-Abkürzung. Die Stepper-Templates werden dort erweitert, wo `codual` noch nicht durchgeht.

**Zerfällt in zwei Teile, und die Reihenfolge ist die Entscheidung aus dem Meilenstein-Abschnitt.**

- **3a, `tsit5`.** Expliziter Runge-Kutta, sieben Stufen mit FSAL. Kein Newton, keine Ordnungswahl,
  keine Historie über Schrittgrenzen. Der Checkpoint ist `x`, `t`, `dt`; `k7` gehört nicht hinein,
  siehe Befund 2 unten. Die Abstraktion, die hier entsteht, muss die beiden anderen tragen, also
  wird sie danach gebaut und nicht nach `tsit5` allein.
- **3b, `multistepper` und `rosenbrock4`.** Nach dem Meilenstein, also nach Stufe 7. Newton-IFT,
  variable Ordnung, Nordsieck-Historie, `qwait` mit `saved_tq5`, dazu beide LU-Tags. Das ist der
  schwerste Block des Plans, und er trifft dann auf eine stehende, gegen den Vorwärtsmodus geprüfte
  Maschinerie statt auf sechs gleichzeitig neue Teile. `bdf` und `adams` sind dasselbe Template,
  also ein Stück Arbeit für zwei der vier Methoden.

Hier entsteht auch der Verifikationsschalter: die Neurechnung nimmt die Fehlernorm und das
Kontrollgesetz wahlweise mit ins Tape oder behandelt sie als Konstanten. Beides läuft über
denselben Code, der Unterschied ist, ob die Steuergrößen als `codual` oder als `double` geführt
werden, also ein `if constexpr` und kein zweiter Pfad.

Test, mit ausgeschalteter Steuerungskette: für einen einzelnen Schritt ist
`w' * (dy_{k+1}/dy_k)` aus dem Rückwärtslauf gleich der Kontraktion der Vorwärtssensitivität
desselben Schritts, auf Rundungsniveau. Das ist prüfbar, ohne dass die Trajektorie steht.

**Stand 3a, 2026-09-08: steht.** `inst/include/cppde/cppde_reverse_step.hpp` trägt
`step_checkpoint<Stepper, T>` und `step_recorder<Stepper, T>`; der Test ist
`dev/cxx/test_reverse_step.cpp`, gefahren über `dev/cxx/run.sh --reverse-step`. Geprüft ist
`w' S` gegen `S' w` bei 1e-14 relativ, für die drei Einheitsseeds und für einen gemischten,
über einen, zwei und vier verkettete Schritte, dazu der neu gerechnete Schrittwert bitgleich
gegen den Vorwärtslauf.

Sechs Befunde:

1. **`codual` fehlten die zusammengesetzten Zuweisungen.** `vec_axpy` schreibt
   `y[i] += alpha * x[i]`, und `+=` gab es nicht; ohne sie instanziiert kein Stepper. Nachgetragen
   nach dem Muster von `dual`, also frei definiert im Mathe-Header. Die skalaren `+=` und `-=`
   behalten den Slot statt einen Knoten zu schreiben, weil `d(x + c)/dx` gleich eins ist. Das ist
   dieselbe Lehre wie bei `select` in Stufe 2: die codual-Oberfläche muss alles tragen, was ein
   Aufrufer auf dem Zustandstyp macht, nicht nur die Rechenoperationen. Sonst ging `tsit5`
   unverändert durch.
2. **Der Checkpoint von `tsit5` enthält `k7` nicht, anders als oben angenommen.** Das
   FSAL-Recycling ist eine Ersparnis, keine Abhängigkeit: das übernommene `k1` ist `f(x, t)` am
   Checkpoint, die Neurechnung erzeugt es bitgleich. Ein gespeichertes `k7` als Konstante
   einzuspielen würde dem Band gerade die Abhängigkeit `dk1/dx` nehmen, also den Gradienten
   verfälschen. Der Checkpoint ist damit `x`, `t`, `dt`. Der Vier-Schritt-Fall im Test deckt genau
   das ab.
3. **Die Abstraktion heißt `carry`.** Ein Schritt ist `(carry, theta) -> x_out`, wobei `carry`
   alles ist, was er aus einem früheren Schritt liest. Für ein Einschrittverfahren ist das der
   Zustand, für den Multistepper das Nordsieck-Array; `load()` legt beides als Tape-Unabhängige an,
   der zweite Teil im Ausgabeparameter `history`. Damit steht die Naht für 3b, ohne dass sie heute
   geraten werden musste. Dazu `rebind_value` an allen drei Steppern, ein Alias je Klasse, weil nur
   der Stepper seine übrigen Template-Argumente kennt.
4. **`h` ist im Reverse-Pfad eine Tape-Unabhängige, keine Konstante.** Das ist Festlegung 1 und
   der Auslieferungsstand; der Verifikationsschalter schaltet die Kette *ab*, nicht *an*.
   `ad_traits::step_coef` entscheidet den Typ der Stage-Koeffizienten: `double` für einen
   arithmetischen oder vorwärts-AD-Zeittyp, symbolisch unter `is_reverse`. Der Vorwärtspfad ist
   damit unverändert, und zwar nicht aus Sparsamkeit, sondern weil `dt` dort nie eine Tangente
   bekommt: `onestep_controller::error()` gibt ein `double` zurück und das Kontrollgesetz ist
   `dt *= factor` mit `double`-Faktor. Ein duales `h` schleppte dort eine Nulltangente durch jede
   Stage-AXPY. Wichtig dabei: `time_type` des Controllers ist `stepper_type::value_type`
   (`cppde_onestep_controller.hpp:87`) und **nicht** `double`, ein Umschalten auf "symbolisch für
   jeden AD-Zeittyp" träfe also sehr wohl den Produktionspfad.

   Der Schalter ist `step_recorder::tape_stepsize(false)`. Er nimmt `t` und `h` aus den
   Unabhängigen; der Rest des Bandes ist unberührt, weil eine nicht abhängige `codual` ohnehin
   keinen Knoten schreibt. Der Test prüft beide Stellungen gegeneinander: eingefroren sind beide
   Zeitkotangenten exakt null und der x/theta-Block ist identisch mit dem vollen Lauf.

5. **Das Orakel für den Schrittweitenterm ist das Compile-Makro, das der Plan vorgesehen hat.**
   `CPPDE_SYMBOLIC_STEPSIZE` lässt `step_coef` auch für einen Vorwärts-AD-Zeittyp symbolisch
   werden. Damit trägt die Vorwärtsreferenz `h` als eigene Richtung und der Vergleich ist wieder
   eine Gleichheit auf Rundungsniveau statt einer Plausibilisierung. Nur die Test-TU definiert es.
   Damit sind es im Test `n_x + n_theta + 2` Richtungen: Zustand, Parameter, `t0` und `h`.
6. **Wer `do_step` direkt fährt, muss `prepare_sensitivities()` selbst rufen.** Sonst nimmt der
   FSAL-Zweig unter `dual` den Wert von `k7` nach `k1` und lässt dessen Tangenten stehen, weil
   `m_K.slot_stride()` bei ungeprimtem Stagemakel null ist. Das hat den ersten Testlauf über zwei
   Schritte um 1e-3 danebenliegen lassen; der Fehler saß im Vorwärtsorakel, nicht im Adjoint.
   Betrifft nur Testtreiber, im Produktionspfad ruft der Controller es.

Nebenher: `dev/cxx/run.sh` kennt `--reverse-step`, findet die R-Bibliothek auch unter Windows, wo
sie nicht in `R.home("lib")` liegt, und überspringt den Sanitizer-Schritt statt an ihm abzubrechen,
wenn die Runtimes fehlen. Der `-O2`-Build ist mit `-Wall -Wextra` warnungsfrei; dafür hat der
No-op-`scoped_timer` in `cppde_profiler.hpp` einen eigenen Destruktor bekommen.

Berührt für das symbolische `h`: `ad_traits::step_coef` und `step_coef_of` als der eine Ort, an
dem die Entscheidung fällt; `vec_axpy_stage` in `cppde_dual_slab.hpp`, das bei einem
`double`-Koeffizienten den Slab-Pfad nimmt und sonst elementweise geht, weil der Reverse-Replay
keinen Slab hat; in `tsit5::do_step` die fünfzehn Stage-AXPYs, die Lösung und die Fehlerschätzung.
Der Zeitpunkt `t` ist aus demselben Grund eine Unabhängige: `t_k = t_0 + sum h_j` trägt eine
Ableitung, sobald die `h_j` eine tragen, und ein nicht-autonomer Modellrumpf liest sie.

**Was 3b von hier aus noch braucht:** eine `step_checkpoint`-Spezialisierung je Stepper, die
`history` füllt, plus die IFT-Abkürzung im Newton-Korrektor. `step_recorder` selbst ist
stepperfrei und sollte unverändert tragen. Für `rosenbrock4` und den Multistepper müssen die
Stage-Koeffizienten dieselbe `step_coef`-Behandlung bekommen wie `tsit5`; im Multistepper hängen
zusätzlich die Nordsieck-Koeffizienten `m_l` und `m_tq` an der Schrittweitenhistorie und gehören
damit ebenfalls aufs Band.

### Stufe 4. Die Trajektorie rückwärts

Checkpointspeicher über die Schrittfolge, Rückwärtsschleife, Adjungierung der
Dense-Output-Interpolation, Einspeisung der Seeds an den Beobachtungszeiten. Quadratur für
`dL/dtheta` als Akkumulator ohne Zustandsdimension.

Test, mit ausgeschalteter Steuerungskette: `w' * S` aus dem Vorwärtslauf gegen den
Adjoint-Gradienten, auf Rundungsniveau. Für den Meilenstein auf dem Zerfallsmodell unter `tsit5`;
mit Stufe 3b kommen die steifen Modelle aus `tests/benchmarks/` und die übrigen drei Methoden dazu.
Die Vorlage für die Testform ist `test-ode-methods.R:52`, die vorhandene Sensitivitätsprüfung über
`methods_all <- c("bdf", "adams", "rb4", "tsit5")` (`:15`), nur mit dem Vorwärtsmodus als Orakel
statt finiter Differenzen.

Mit eingeschalteter Steuerungskette, also im Auslieferungsstand, ist die Differenz zum
ausgeschalteten Lauf der Schrittweitenterm. Nach der Gegenprobe aus Stufe 0 muss sie **auf Höhe des
lokalen Fehlers verschwinden** und nicht mit `rtol` skalieren, denn innerhalb eines Kontrollpfads
trägt die Schrittweite nichts bei. Eine große Differenz ist damit ein Fehler im Rückwärtslauf und
kein wiedergewonnener Term. Das macht Stufe 4 zu einer Gleichheitsprüfung und ersetzt die
Plausibilisierung, die hier ursprünglich vorgesehen war.

### Stufe 5. Events und Wurzeln

Transponierte Saltation. `cppde_saltation.hpp` hat die Vorwärtsvariante inklusive
`compute_dt_star` mit IFT und Korrektur zweiter Ordnung; der Rückwärtsfall ist deren
Transponierte. Der Neustart nach einem Ereignis verwirft die Nordsieck-Historie
(`restart_from_order1`, `cppde_multistepper.hpp:1504-1518`), was den Rückwärtslauf an dieser
Stelle abschneidet und die Buchführung vereinfacht.

Nur der `process_dense`-Pfad wird unterstützt. `localize_root_controlled`
(`cppde_event_engine.hpp:600-627`) re-integriert innerhalb der Lokalisierung und mutiert dabei
den Stepper; das ist rückwärts nicht sinnvoll zu rekonstruieren. Für Modelle mit Wurzelereignissen
und `useDenseOutput = FALSE` erhebt der Reverse-Modus einen Fehler statt still falsch zu rechnen.

Test: gegen die Vorwärtssensitivitäten auf den Ereignismodellen aus `test-ode-methods.R`.

### Stufe 6. Die Naht

Zwei Teile, und der erste ist kleiner als gedacht.

**Das Modellobjekt.** `cppODE(..., sweep = "reverse")` übersetzt ein viertes Objekt und stempelt
`attr(modelname, "sweep")` neben die vorhandenen `"deriv"`, `"deriv2"` und `"method"`
(`cppODE.R:1232-1236`). `solveODE()` liest das wie heute schon `attr(model, "deriv")`
(`solveODE.R:20-21`) und verweigert einen `seed` an einem Vorwärtsmodell so, wie es heute
`sens1ini` an einem `deriv = FALSE`-Modell verweigert (`:34-37`). Keine Modusfahne wandert durch die
Kette.

**Das Argument.** `seed` als 15. positionelles Argument vor `dimnames`. Betroffen sind
`cppDE/R/solveODE.R:348-352` und die Validierung ab `:3`, die Aufrufstellen `:583-590`,
`:661-670`, `:697-728`, `:831-837`, `:872-874`, dazu `read_solve_args`
(`cppde_r_batch.hpp:205`), `read_cond_args` (`:256`, braucht `cond_elt(cond, 14)`) und beide
Emitter (`R/cppODE.R:1127-1134`, `codegen_cvode.py:2008-2019`). Der Batch-Einstieg behält seine
drei Argumente.

`read_solve_args` nimmt heute vierzehn SEXPs, `read_cond_args` reicht `cond_elt(cond, 0)` bis
`cond_elt(cond, 13)` durch, und `cond_elt` (`cppde_r_batch.hpp:253`) fällt für fehlende Indizes auf
`R_NilValue` zurück. Der Seed ist damit Index 14 und der Batch-Pfad bleibt in beide Richtungen
kompatibel; nur der Einzelpfad bricht gegen alte kompilierte Modelle.

Offen und bewusst nicht in dieser Stufe: `sens1ini` und `sens2ini` in ein gemeinsames `seed`-Argument
zu falten. Beide sind nach der Implementierung benannt, nicht nach ihrer Rolle, und die Rolle ist in
allen Fällen dieselbe, die Ableitung der anschließenden Abbildung am Rand. Das wäre ein Umbau der Positionen 2 und 3
und damit ein Bruch, während das Anhängen an Position 14 additiv ist. Erst reverse zum Laufen
bringen, dann über die Umbenennung entscheiden.

`cvode()` bleibt in dieser Stufe unberührt; sein Reverse-Modus ist Stufe 8.

### Stufe 7. Die Kette in dMod2

Auf `devel-reverseAD` in dMod2. Jede Schicht bekommt neben ihrer Vorwärtsform eine
Rückwärtsform:

- `normL2` erzeugt den Seed. Er ist nicht der Residuenvektor: das Fehlermodell hängt selbst von
  theta ab, `src/residual_kernel.cpp:144-158` koppelt über
  `dwr[k] = inv_s * dx - wr * inv_s * ds`. Der Kotangens wirkt auf Prediction und Fehlermodell.
  Die BLOQ-Zweige haben ihre eigenen Formen (`:325-372`).
- `Y` (`R/prediction.R:731`): statt `dG/dX %bmm% dX + dG/dP %bmm% dP` zwei vjps.
- `Xs.cppDE` (`R/prediction.R:172-389`): statt `prep1` (`:245`) den Seed nach hinten geben.
- `Pexpl` (`R/parameters.R:232`): `w' Jac` statt `Jac %*% dP`.
- `Pimpl` (`R/parameters.R:903`): transponierter IFT-Solve. Kompiliert sich ein eigenes
  `cppODE`-Modell mit `rootfunc = "equilibrate"` und löst darin eine geschachtelte ODE
  (`:1329-1338`, `:1392`), fällt also unter Stufe 6 mit.
- `Pequil` (`R/parameters.R:1590`): transponierte Endpunktsensitivität, ebenso geschachtelt
  (`:1632-1641`, `:1731`).
- `R/classes.R`: `.evalProd` (`:317`) wertet strikt vorwärts aus. Ein Rückwärtsdurchlauf
  durch denselben Baum ist dort neu und ist die einzige strukturelle Änderung an der
  Kompositionsalgebra.

Zwei Vorlagen, beide vorhanden, keine davon zu erfinden:

- **Die Auswahl je Aufruf** ist `pickModel(deriv, deriv2)` (`prediction.R:298`), das heute zwischen
  `func`, `extended` und `extended2` wählt. Reverse ist ein vierter Zweig, und `odemodel()` baut
  entsprechend ein viertes Objekt neben `odeClass.R:206`, `:214`, `:220`.
- **Das Durchreichen der Absicht** ist `hessian`: seit `e1a88bc` läuft es von `trust_driver.h:181`
  über die Objective-Closure und `.evalProd` bis in die Kerne, mit `build_hessian` als abgeleitetem
  Flag (`objClass.R:208`). `sweep` tritt als dritte Achse daneben und nimmt denselben Weg.

Bedingungen bleiben unabhängig, die vorhandene Parallelisierung ist nicht betroffen.

Test: `normL2(deriv = TRUE)$gradient` vorwärts gegen rückwärts, auf den Fixtures
`fx_decay_compiled()` und `fx_decay_multicond_compiled()`
(`tests/testthat/helper-fixtures.R:54`, `:108`), einschließlich Log-Trafo und
bedingungsspezifischer Theta-Teilmengen. Toleranz ist Rundungsniveau, nicht die 1e-3 der
FD-Prüfungen in `test-normL2.R:449` und `:481`. Der Multicond-Fixture nimmt zusätzlich den
Batch-Pfad `P2Xbatch` (`prediction.R:330`), der erst ab zwei Bedingungen greift
(`classes.R:247`).

### Stufe 8. CVODES ASA als Vergleich

`cvode()` bekommt einen Reverse-Modus über die ASA von CVODES. SUNDIALS wird von `./configure`
bereits geprüft und in `inst/cvodeConfig.dcf` gemeldet, `codegen_cvode.py` emittiert schon die
analytischen Vorwärtssensitivitäten. Zu bauen sind die Adjoint-Rechte-Seite, der Term
`(df/dp)^T lambda` und die Checkpointkonfiguration.

**Dies ist kein Orakel.** ASA ist der kontinuierliche Adjoint, also erst differenzieren, dann
diskretisieren, mit eigener Schrittweitensteuerung im Rückwärtslauf. Der diskrete Adjoint aus den
Stufen 3 bis 5 differenziert die Diskretisierung. Beide stimmen nur bis auf `O(tol)` überein, und
das ist genau die Größenordnung des Schrittweitenterms aus Stufe 0. Ein auf `O(tol)` genauer
Vergleich kann eine Abweichung von `O(tol)` nicht beurteilen. Die Verifikation bleibt beim
Vorwärtsmodus.

Was der Vergleich liefert:

1. Plausibilität durch eine unabhängige Implementierung fremder Mathematik.
2. Wandzeit und Speicher gegen den diskreten Adjoint, auf Bachmann und Lang.
3. Die Größe von ASAs Inkonsistenz zwischen Wert und Gradient. Das ist das zweite Argument der
   Begründung für den diskreten Weg und bisher nur behauptet. Messbar als Differenz zwischen dem
   ASA-Gradienten und dem Vorwärtsgradienten desselben Laufs, aufgetragen über `rtol`.

Die Messung aus Punkt 3 gehört mit der aus Stufe 0 in dieselbe Auswertung: beide betreffen die
Frage, wie viel Gradientenkonsistenz auf dieser Modellklasse überhaupt wert ist.

## Verifikation

Das Orakel ist durchgehend der Vorwärtsmodus, gegen den eingefrorenen Reverse-Pfad gestellt. Die
beiden berechnen dieselbe Größe aus entgegengesetzter Richtung, also gilt Rundungstoleranz und
nicht Solvertoleranz. Kein Test ruht auf finiten Differenzen, außer den bereits vorhandenen, und
keiner ruht auf ASA.

Vorlagen, alle vorhanden, alle nach demselben Muster "zwei Wege, dieselbe Zahl":

| Ebene | vorhandener Test | Form |
|---|---|---|
| Primitive | `cppDE test-dual2nd-primitives.R:21`, `:38` | `run_modes()` über `dual` und `symbolic`, `expect_modes_agree()` bei 1e-12 auf `y`, `dy`, `d2y` |
| Stepper | `cppDE test-ode-methods.R:15`, `:52` | `methods_all` über alle vier Methoden, Sensitivitäten gegen FD bei `eps = 1e-5` |
| algebraisch | `dMod2 test-Pexpl.R:78-94`, `test-Y.R:69-78` | `derivMode` dual gegen symbolic, Wert und Jacobi, 1e-10 |
| Objective | `dMod2 test-deriv2.R:439`, `:470`, `:504` | `deriv2`-Gradient gegen `deriv1`-Gradient, 1e-7 |

Der Reverse-Modus tritt in jede dieser Zeilen als weiterer Weg ein. Die Stepper-Zeile ist die
einzige, die heute auf finiten Differenzen ruht, und genau dort ersetzt der Vorwärtsmodus sie als
schärferes Orakel.

- cppDE: `TESTTHAT_CPUS=6 NOT_CRAN=true Rscript -e 'testthat::test_dir("tests/testthat",
  package = "cppDE", load_package = "installed", reporter = "summary")'`
- dMod2 nach `R CMD INSTALL` von cppDE zuerst, dann dMod2, seriell wegen des kaputten
  Parallel-Runners unter `load_all()`.
- `inst/benchmarks/bench_gradientCost.R` misst bereits `factor_grad` über `n_theta` und ist damit
  der Bench dieses Plans, nicht ein neu zu schreibender: sein eigener Kopf nennt `factor` "the number
  an adjoint has to beat". Er bekommt eine Adjoint-Spalte, sobald Stufe 7 steht. Braucht eine ruhige
  Maschine und `OMP_NUM_THREADS=1`, `min` über Wiederholungen.
- Danach `inst/benchmarks/bench_hessianSource.R` mit `hessianMethod = "sr1"`, `hessianInit =
  "identity"` unter Adjoint-Gradienten gegen `gn`, in `evalPerHit` und in Wandzeit. Das ist das Tor
  aus dem Abschnitt oben, diesmal gemessen statt hochgerechnet.

## Risiken

- **Der Schrittweitenterm ist vernachlässigbar, und das ist gemessen.** Stufe 0 hat es am
  2026-09-07 entschieden: die Kennzahl skaliert mit der Sweep-Dichte, ist also Sprunghöhe und keine
  Ableitung. Das Risiko ist erledigt und in einen Vorteil umgeschlagen, weil das Orakel für Stufe 4
  dadurch schärfer wird. Was bleibt, ist die Unstetigkeit von `y_h` in theta selbst, und die trifft
  jedes Ableitungsverfahren gleichermaßen.
- **Der variable-Ordnungs-Multistepper ist das schwerste Stück.** Die Ordnungsauswahl über exakte
  Gleitkommagleichheit und das `qwait`-Fenster mit `saved_tq5`, das eine Schrittweite Verzögerung
  trägt (`cppde_multistepper.hpp:1419-1423`), sind Zustand über Schrittgrenzen hinweg und gehören
  vollständig in den Checkpoint.
- **Falsche Gradienten sind unsichtbar.** Die Vorwärtstrajektorie bleibt korrekt, egal wie falsch
  der Adjoint ist. Deshalb steht in jeder Stufe ein Vergleich gegen den Vorwärtsmodus, und
  deshalb ist Stufe 3 auf Einzelschrittebene und nicht erst auf Trajektorienebene abgesichert.
- **`inst/COPYRIGHTS` muss mitwachsen**, falls Teile aus CVODES ASA abgeleitet werden. Der
  Nordsieck-Stepper ist bereits als SUNDIALS-Portierung ausgewiesen.
- **Die ABI-Änderung in Stufe 6 bricht Einzelaufrufe** gegen alte kompilierte Modelle. Der
  Batch-Pfad ist über `cond_elt`s Rückfall auf `R_NilValue` abwärtskompatibel, der Einzelpfad nicht.
- **Der Meilenstein berührt den Produktionspfad nicht.** Er steht auf `tsit5`, dMod2 kompiliert
  ausnahmslos `bdf`. Das ist der Preis der Fehlerlokalisierung und er ist bewusst bezahlt, aber
  daraus folgt: Stufe 3b ist keine Kür, sondern die Fortsetzung, und der Meilenstein darf nicht als
  "reverse läuft" verbucht werden, solange sie aussteht.
- **Vollständigkeit ist die teuerste Festlegung dieses Plans.** Vier Methoden, zwei LU-Tags, Events,
  Wurzeln, Forcings, Batch und der CVODES-Backend ergeben eine Fläche, die vorwärts über Jahre
  gewachsen ist. Jede Stufe trägt deshalb ihren Anteil der Fläche mit, statt sie am Ende
  nachzureichen; die Kombinationsmatrix gehört von Stufe 3b an in die Tests und nicht in eine
  Schlussstufe.

## Noch zu untersuchen

- **Ein symbolisches `t` in den Saltationskorrekturen, und was es vorwärts kosten würde.**
  `cppde_saltation.hpp` leitet die Ereigniszeit heute von Hand her: `compute_dt_star` löst die
  Wurzelbedingung per IFT und trägt eine Korrektur zweiter Ordnung nach. Seit Stufe 3a kann der
  Zeittyp selbst eine Ableitung tragen, `dual` wie `codual`, und `ad_traits::step_coef` ist der
  eine Ort, an dem das umgestellt würde. Drei Fragen, in dieser Reihenfolge:

  1. Bekommt `t*` seine Ableitung direkt aus der Differentiation der Wurzelbedingung, und wie viel
     der handgeschriebenen Herleitung entfällt dadurch? Betrifft beide Richtungen, nicht nur den
     Adjoint.
  2. Falls ja: ist das auf Ereignismodellen auch schneller, oder nur kürzer? Ein duales `t*` ersetzt
     handgeschriebene Formeln durch AD-Arithmetik, und das ist nicht automatisch der billigere Weg.
  3. Falls es sich lohnt: sind `t` und `dt` dual im **Normalfall** verkraftbar? Heute stehen sie
     bewusst draußen, weil `dt` vorwärts ohnehin keine Tangente bekommt, die Nulltangente also durch
     jede Stage-AXPY liefe. Ob das messbar ist, ist offen — die Stage-AXPYs sind BLAS-gestützt und
     der Zuschlag ist eine Spalte auf `n_theta`, es kann also gut sein, dass es im Rauschen liegt.
     Zu messen auf einem Modell ohne Ereignisse, sonst zahlt man den Preis für einen Nutzen, den
     nur die Ereignismodelle haben. Wenn es nichts kostet, entfällt die Fallunterscheidung in
     `step_coef` und der Vorwärtsmodus wird einfacher statt komplizierter.

  Gehört ans Ende: erst wenn Stufe 5 steht und die vorhandene Saltation gegen den Adjoint geprüft
  ist, gibt es einen Vergleichsmaßstab für Punkt 1 und ein Modell für Punkt 2.

## Was nicht

- ASA nicht als Verifikationsorakel, nur als Vergleich und Benchmark. Der Grund steht in Stufe 8.
- Keine zweite Ordnung. Forward-over-reverse ist der nächste Schritt und braucht den verifizierten
  Adjoint erster Ordnung als Fundament; das heutige `deriv2` ist dann das Orakel auf kleinen
  Modellen. Wenn es soweit ist, heißt zweite Ordnung hier die **volle Matrix aus `n_theta` Seeds**
  und kein matrixfreies Krylov-Teilproblem: der Optimierer-Plan schließt Steihaug-Toint und GLTR
  strukturell aus, weil ein schlaffes Spektrum keinen effektiven Rang hat und die Krylov-Zahl damit
  auf der Parameterzahl landet. Beides zusammen ist widerspruchsfrei, denn `n_theta`
  Hesse-Vektor-Produkte sind genau die volle Matrix, und die kostet rückwärts `O(n_theta)` statt
  `O(n_theta^2)`.
- Kein Revolve-Checkpointing. Alles speichern ist bei diesen Größen bezahlbar.
- Keine Expression-Template-Schicht für `codual`.
- Kein `localize_root_controlled` im Reverse-Modus.
