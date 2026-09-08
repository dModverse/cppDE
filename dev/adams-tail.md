# Adams: der replizierte Nachlauf trifft Nordsieck-Slot 1 nicht

Eingrenzung zu dem einen offenen Punkt aus `dev/adjoint-plan.md`, Stufe 4. Kurz: der Rückwärtslauf
für `adams` stimmt auf Schrittebene, auf Trajektorienebene nicht, und die Ursache sitzt im
replizierten Nachlauf eines Schritts, nicht in der Tape-Maschinerie.

## Der Befund

`dev/cxx/test_reverse_trajectory_methods.cpp` fährt die aufgezeichnete Schrittfolge in reinen
Doubles nach: Carry laden, `do_step`, `complete_step`, `set_tn_current`, `prepare_dense_output`,
`save_acor_to_zn_qmax`, `set_order_for_next_step`, `rescale`. Was dieser Nachlauf am Ende eines
Schritts im Nordsieck-Array stehen hat, muss der Carry sein, den der Lauf für den nächsten Schritt
aufgezeichnet hat.

Für `bdf`, `rb4` und `tsit5` stimmt das exakt. Für `adams` gilt:

- **Slot 0 stimmt**, über alle Schritte.
- **Slot 1 weicht ab**, ab Schritt 2, also ab dem ersten Ordnungswechsel, relativ um 4e-4 und
  wachsend. Beispiel: aufgezeichnet `-1.9428344005227227e-04`, nachgerechnet
  `-6.0568351513154341e-04`, Verhältnis 3.117.
- Die Abweichung ist deterministisch und reproduzierbar.

Auf der Ableitungsseite schlägt das als 1e-3 bis 2e-3 zwischen Vorwärts- und Rückwärtsmodus durch,
weil die Referenz vorwärts verkettet und ein einmal abweichender Carry alles Folgende mitnimmt.

## Was ausgeschlossen ist

- **Der Adjoint selbst.** `dev/cxx/test_reverse_step_multistep.cpp` fährt beide Multistep-Verfahren
  bei 1e-9: jeder Nordsieck-Slot einzeln geseedet, Ordnung 2 bis 5, Ordnung gehalten, erhöht und
  gesenkt, Rescale angewandt, und der Seed wahlweise auf dem Carry oder durch die Dense-Output-
  Interpolation. Die Gleichung, der transponierte Solve und der Nachlauf stimmen dort für `adams`.
- **Die Tape-Maschinerie.** Slots sind monoton über `rewind()`, ein Wert aus einem früheren Schritt
  liest sich als Konstante. Der Befund liegt ohnehin in reinen Doubles, ganz ohne Band.
- **Der Korrektor.** Die Abweichung ist Slot 1, nicht Slot 0; der Korrektor bestimmt Slot 0.
- **Verworfene Versuche.** Die Prüfung überspringt Schritte, deren Nachfolger einen Versuch verworfen
  hat, weil dessen Wiederholung `reload_zn1_from_f` fährt und damit `zn[1]` neu setzt. Die
  Abweichung bleibt.
- **`set_tn_current`, `save_acor_to_zn_qmax`, `prepare_dense_output`.** Alle drei fehlten und sind
  nachgetragen; sie haben `bdf`, `rb4` und `tsit5` in Ordnung gebracht und `adams` nicht.

## Wo weiterzusuchen wäre

Slot 1 wird auf dem Weg von drei Stellen berührt, und nur diese kommen in Frage:

1. **`complete_step`**, `zn[1] += l[1] * acor`. `l` kommt aus `adams_set_coefficients`, das an `q`,
   `qwait`, `h` und `tau` hängt. Alle vier stehen im Carry, `qwait` wird aber vom Controller in
   `prepare_next_step` gesetzt (`set_qwait(2)` oder `set_qwait(L)`) und im Nachlauf nicht. Für den
   Rückwärtslauf ist das folgenlos, weil jeder Schritt seinen eigenen Carry lädt; für den
   verketteten Nachlauf ist es die erste Stelle, an der die beiden auseinanderlaufen können.
2. **`rescale(eta)`**, `zn[j] *= eta^j`. `eta` wird im Sammler als `hscale_nachher / dt` abgeleitet.
   Das ist richtig, solange der Controller genau einmal skaliert und nichts anderes `m_hscale`
   anfasst. Zu prüfen, ob das für `adams` gilt.
3. **`adamsIncreaseOrder`**, `zn[m_L] = 0`. Rührt Slot 1 nicht an, verschiebt aber, welche Slots
   danach skaliert werden, weil `q` wächst.

Das Verhältnis 3.117 gehört gegen `eta` und gegen `l[1]` des betreffenden Schritts gehalten; es ist
weder `eta` noch eine Potenz davon, was gegen (2) allein spricht.

## Wie es zu reproduzieren ist

```sh
dev/cxx/run.sh --reverse-trajectory-methods
```

Die Prüfung heißt `adams carry handoff step <k> slot <j>.<i>` und steht in `replay_values()`. Sie
läuft nur im ungestörten Lauf, nicht in den beiden gestörten Läufen der Finite-Differenzen-Probe,
wo sie zu Recht abweichen würde.

Solange der Punkt offen ist, steht die Trajektorien-Prüfung für `adams` auf 1e-2 und trägt nichts;
`bdf` steht auf 1e-6, die beiden Einschrittverfahren auf 1e-10.
