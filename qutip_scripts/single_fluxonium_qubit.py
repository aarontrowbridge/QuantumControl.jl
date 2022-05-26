from qutip import sigmax, sigmaz, qsave

fq = 1 / (72 * 1e-9) 
# fq = 1

H_drift = fq * sigmaz()
H_drive = sigmax()

qsave(H_drift, 'qutip_saved_objects/single_fluxonium_qubit/H_drift.qu')
qsave(H_drive, 'qutip_saved_objects/single_fluxonium_qubit/H_drive.qu')