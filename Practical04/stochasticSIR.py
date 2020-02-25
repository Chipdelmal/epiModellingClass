
import model
import functions as fun

(bet, gamm, iota, N, dt) = [0.1, 0.05, 0.01, 1000.0, 0.1]
parms = (bet, gamm, iota, N, dt)
sir_out = model.simulateSIR(parms, 500, 2001, 1)
(t, S, I, R) = [list(sir_out[i]) for i in ['t', 'S', 'I', 'R']]
tp = ((S, '#02146b', 'S'), (I, '#ffb428', 'E'), (R, '#b4e830', 'I'))
fun.plotEpiDynamicsPop(tp, t, t[-1], 1000)
